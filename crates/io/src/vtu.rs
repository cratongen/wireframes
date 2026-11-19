//! VTK XML: `.vtu` unstructured grids and `.vtp` poly data.

use std::io::{Read, Write};

use flate2::{Compression, read::ZlibDecoder, write::ZlibEncoder};
use wireframes_core::{Error, Mesh, Result};

use crate::vtk::{
    Array, Cells, PolyData, Ty, Values, concat_tables, data_to_poly, flatten, grid_problem,
    grid_to_data, mesh_to_grid, poly_to_data, to_points, unflatten,
};
use crate::{Data, UnstructuredGrid, WriteOptions};

struct Tag {
    name: String,
    attrs: Vec<(String, String)>,
    close: bool,
    empty: bool,
    start: usize,
    after: usize,
}

impl Tag {
    fn attr(&self, k: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(n, _)| n == k)
            .map(|(_, v)| v.as_str())
    }

    fn count(&self, k: &str) -> Result<usize> {
        match self.attr(k) {
            None => Ok(0),
            Some(v) => v
                .trim()
                .parse()
                .map_err(|_| Error::format(format!("{} {k}={v:?} is not a count", self.name))),
        }
    }
}

fn find(h: &[u8], n: &[u8]) -> Option<usize> {
    h.windows(n.len()).position(|w| w == n)
}

fn unescape(s: &str) -> String {
    s.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&")
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// One tag starting at `b[start] == b'<'`.
fn scan_tag(b: &[u8], start: usize) -> Result<Tag> {
    let bad = |m: &str| Error::format(format!("malformed xml tag at byte {start}: {m}"));
    let close = b.get(start + 1) == Some(&b'/');
    let mut j = start + 1 + close as usize;
    let ns = j;
    while j < b.len() && !b[j].is_ascii_whitespace() && b[j] != b'>' && b[j] != b'/' {
        j += 1;
    }
    let name = String::from_utf8_lossy(&b[ns..j]).into_owned();
    let mut attrs = Vec::new();
    let mut empty = false;
    loop {
        while j < b.len() && b[j].is_ascii_whitespace() {
            j += 1;
        }
        match b.get(j) {
            None => return Err(bad("unterminated")),
            Some(b'>') => {
                j += 1;
                break;
            }
            Some(b'/') => {
                if b.get(j + 1) != Some(&b'>') {
                    return Err(bad("stray '/'"));
                }
                empty = true;
                j += 2;
                break;
            }
            _ => {}
        }
        let ks = j;
        while j < b.len() && b[j] != b'=' && !b[j].is_ascii_whitespace() && b[j] != b'>' {
            j += 1;
        }
        let key = String::from_utf8_lossy(&b[ks..j]).into_owned();
        while j < b.len() && b[j].is_ascii_whitespace() {
            j += 1;
        }
        if b.get(j) != Some(&b'=') {
            return Err(bad("attribute without value"));
        }
        j += 1;
        while j < b.len() && b[j].is_ascii_whitespace() {
            j += 1;
        }
        let q = *b.get(j).ok_or_else(|| bad("unterminated"))?;
        if q != b'"' && q != b'\'' {
            return Err(bad("unquoted attribute"));
        }
        let vs = j + 1;
        let ve = vs
            + b[vs..]
                .iter()
                .position(|&c| c == q)
                .ok_or_else(|| bad("unterminated attribute"))?;
        attrs.push((key, unescape(&String::from_utf8_lossy(&b[vs..ve]))));
        j = ve + 1;
    }
    Ok(Tag {
        name,
        attrs,
        close,
        empty,
        start,
        after: j,
    })
}

fn scan(b: &[u8]) -> Result<Vec<Tag>> {
    let mut tags = Vec::new();
    let mut i = 0;
    while let Some(off) = b[i..].iter().position(|&c| c == b'<') {
        let start = i + off;
        let rest = &b[start..];
        let skip = |end: &[u8]| {
            find(rest, end)
                .map(|e| start + e + end.len())
                .ok_or_else(|| Error::format("unterminated xml declaration or comment"))
        };
        if rest.starts_with(b"<?") {
            i = skip(b"?>")?;
        } else if rest.starts_with(b"<!--") {
            i = skip(b"-->")?;
        } else if rest.starts_with(b"<!") {
            i = skip(b">")?;
        } else {
            let t = scan_tag(b, start)?;
            i = t.after;
            tags.push(t);
        }
    }
    Ok(tags)
}

fn b64_value(c: u8) -> Option<u32> {
    Some(match c {
        b'A'..=b'Z' => c - b'A',
        b'a'..=b'z' => c - b'a' + 26,
        b'0'..=b'9' => c - b'0' + 52,
        b'+' => 62,
        b'/' => 63,
        _ => return None,
    } as u32)
}

/// Decodes four characters at a time, so separately encoded chunks
/// (header, then data, each with its own padding) decode to their concatenation.
fn b64_decode(text: &[u8]) -> Result<Vec<u8>> {
    let s: Vec<u8> = text
        .iter()
        .copied()
        .filter(|c| !c.is_ascii_whitespace())
        .collect();
    if !s.len().is_multiple_of(4) {
        return Err(Error::format("base64 data length is not a multiple of 4"));
    }
    let mut out = Vec::with_capacity(s.len() / 4 * 3);
    for g in s.chunks_exact(4) {
        let mut n = 0u32;
        let mut pad = 0;
        for (i, &c) in g.iter().enumerate() {
            let v = if c == b'=' && i >= 2 {
                pad += 1;
                0
            } else if pad > 0 {
                return Err(Error::format("bad base64 padding"));
            } else {
                b64_value(c).ok_or_else(|| Error::format(format!("bad base64 byte {c:#x}")))?
            };
            n = n << 6 | v;
        }
        out.extend_from_slice(&n.to_be_bytes()[1..4 - pad]);
    }
    Ok(out)
}

struct Ctx<'a> {
    big: bool,
    header: usize,
    zlib: bool,
    appended: Option<(&'a [u8], bool)>,
    /// File byte offset of the first appended data byte, after the `_`.
    app_start: usize,
    offsets: Vec<usize>,
}

impl Ctx<'_> {
    fn header(&self, b: &[u8], i: usize) -> Result<u64> {
        let s = i
            .checked_mul(self.header)
            .and_then(|a| b.get(a..a + self.header))
            .ok_or_else(|| Error::format("binary array header ends early"))?;
        Ok(match (self.header, self.big) {
            (4, true) => u32::from_be_bytes(s.try_into().expect("4 bytes")) as u64,
            (4, false) => u32::from_le_bytes(s.try_into().expect("4 bytes")) as u64,
            (_, true) => u64::from_be_bytes(s.try_into().expect("8 bytes")),
            (_, false) => u64::from_le_bytes(s.try_into().expect("8 bytes")),
        })
    }

    /// The array bytes after its size header, decompressed.
    fn payload(&self, b: &[u8]) -> Result<Vec<u8>> {
        let short = || Error::format("binary array data ends early");
        let size = |x: u64| usize::try_from(x).map_err(|_| short());
        if !self.zlib {
            let n = size(self.header(b, 0)?)?;
            let end = self.header.checked_add(n).ok_or_else(short)?;
            return b
                .get(self.header..end)
                .map(<[u8]>::to_vec)
                .ok_or_else(short);
        }
        let nb = size(self.header(b, 0)?)?;
        if nb == 0 {
            return Ok(Vec::new());
        }
        let bs = self.header(b, 1)?;
        let lbs = self.header(b, 2)?;
        let mut pos = nb
            .checked_add(3)
            .and_then(|n| n.checked_mul(self.header))
            .filter(|&n| n <= b.len())
            .ok_or_else(short)?;
        let mut out = Vec::new();
        for k in 0..nb {
            let cs = size(self.header(b, 3 + k)?)?;
            let end = pos.checked_add(cs).ok_or_else(short)?;
            let block = b.get(pos..end).ok_or_else(short)?;
            pos = end;
            let expect = if k == nb - 1 && lbs != 0 { lbs } else { bs };
            let mut d = Vec::new();
            ZlibDecoder::new(block)
                .take(expect.saturating_add(1))
                .read_to_end(&mut d)
                .map_err(|e| Error::format(format!("bad zlib block: {e}")))?;
            if d.len() as u64 != expect {
                return Err(Error::format(format!(
                    "zlib block {k} holds {} bytes, header says {expect}",
                    d.len()
                )));
            }
            out.extend(d);
        }
        Ok(out)
    }
}

fn xml_ty(name: &str) -> Result<Ty> {
    Ok(match name {
        "Int8" | "Char" => Ty::I8,
        "UInt8" => Ty::U8,
        "Int16" => Ty::I16,
        "UInt16" => Ty::U16,
        "Int32" => Ty::I32,
        "UInt32" => Ty::U32,
        "Int64" => Ty::I64,
        "UInt64" => Ty::U64,
        "Float32" => Ty::F32,
        "Float64" => Ty::F64,
        other => return Err(Error::format(format!("unknown DataArray type {other:?}"))),
    })
}

fn decode_array(t: &Tag, content: &[u8], ctx: &Ctx, coords: bool) -> Result<Array> {
    let name = t.attr("Name").unwrap_or("").to_string();
    let ncomp = match t.attr("NumberOfComponents") {
        Some(_) => t.count("NumberOfComponents")?,
        None => 1,
    };
    let ty = xml_ty(t.attr("type").unwrap_or(""))?;
    let format = t.attr("format").unwrap_or("ascii");
    let at = match format {
        "appended" => ctx.app_start.saturating_add(t.count("offset")?),
        _ => t.after,
    };
    let what = |e: Error| match e {
        Error::Format { message, .. } => {
            Error::format(format!("array {name:?} with data at byte {at}: {message}"))
        }
        e => e,
    };
    let values = match format {
        "ascii" => {
            let s =
                std::str::from_utf8(content).map_err(|_| Error::format("non-text ascii array"))?;
            let toks: Vec<&str> = s.split_whitespace().collect();
            ty.parse_text(&toks).map_err(what)?
        }
        "binary" => {
            let bytes = b64_decode(content).map_err(what)?;
            ty.decode(&ctx.payload(&bytes).map_err(what)?, ctx.big)
                .map_err(what)?
        }
        "appended" => {
            let off = t.count("offset")?;
            let (app, raw) = ctx
                .appended
                .ok_or_else(|| what(Error::format("appended array without AppendedData")))?;
            let rest = app.get(off..).ok_or_else(|| {
                what(Error::format(format!(
                    "appended offset {off} is past the data, which ends at byte {}",
                    ctx.app_start + app.len()
                )))
            })?;
            let bytes = if raw {
                ctx.payload(rest).map_err(what)?
            } else {
                let end = ctx
                    .offsets
                    .iter()
                    .find(|&&o| o > off)
                    .map_or(app.len(), |&o| o);
                ctx.payload(&b64_decode(&app[off..end]).map_err(what)?)
                    .map_err(what)?
            };
            ty.decode(&bytes, ctx.big).map_err(what)?
        }
        other => return Err(Error::format(format!("unknown DataArray format {other:?}"))),
    };
    if coords && let Some(k) = values.to_f64().iter().position(|x| !x.is_finite()) {
        return Err(what(Error::format(format!(
            "point {} coordinate {} is not finite",
            k / 3,
            k % 3
        ))));
    }
    Ok(Array {
        name,
        ncomp,
        values,
    })
}

struct Piece {
    tag: Tag,
    arrays: Vec<(String, Array)>,
}

impl Piece {
    fn section(&mut self, sec: &str) -> Vec<Array> {
        let (mine, rest): (Vec<_>, Vec<_>) = std::mem::take(&mut self.arrays)
            .into_iter()
            .partition(|(s, _)| s == sec);
        self.arrays = rest;
        mine.into_iter().map(|(_, a)| a).collect()
    }

    fn take(&mut self, sec: &str, name: &str) -> Option<Values> {
        let i = self
            .arrays
            .iter()
            .position(|(s, a)| s == sec && a.name.eq_ignore_ascii_case(name))?;
        Some(self.arrays.remove(i).1.values)
    }

    fn points(&mut self) -> Result<Vec<[f64; 3]>> {
        let n = self.tag.count("NumberOfPoints")?;
        let arrays = self.section("Points");
        let Some(a) = arrays.into_iter().next() else {
            return if n == 0 {
                Ok(Vec::new())
            } else {
                Err(Error::format("piece has points but no Points array"))
            };
        };
        if a.ncomp != 3 || n.checked_mul(3) != Some(a.values.len()) {
            return Err(Error::format(format!(
                "Points array has {} values for {n} points",
                a.values.len()
            )));
        }
        Ok(to_points(a.values.to_f64()))
    }

    fn cells(&mut self, sec: &str, n: usize) -> Result<Cells> {
        let conn = self.take(sec, "connectivity");
        let ends = self.take(sec, "offsets");
        let (conn, ends) = match (conn, ends) {
            (Some(c), Some(e)) => (c.to_i64("connectivity")?, e.to_i64("offsets")?),
            (None, None) if n == 0 => return Ok(Cells::new()),
            _ => {
                return Err(Error::format(format!(
                    "{sec} needs connectivity and offsets"
                )));
            }
        };
        Ok(full_offsets(ends, Some(n), conn))
    }

    fn data(
        &mut self,
        rows_points: usize,
        rows_cells: usize,
    ) -> Result<(wireframes_core::Table, wireframes_core::Table)> {
        let pd = unflatten(self.section("PointData"), rows_points, "point")?;
        let cd = unflatten(self.section("CellData"), rows_cells, "cell")?;
        Ok((pd, cd))
    }
}

/// VTK XML offsets are end offsets; files from VTK 9.4 may carry a leading 0.
fn full_offsets(o: Vec<i64>, count: Option<usize>, conn: Vec<i64>) -> Cells {
    let full = match count {
        Some(n) => o.len().checked_sub(1) == Some(n),
        None => o.first() == Some(&0),
    };
    if full {
        Cells { offsets: o, conn }
    } else {
        Cells::from_ends(o, conn)
    }
}

fn polyhedra(p: &mut Piece, n: usize) -> Result<Option<(Vec<i64>, Vec<i64>)>> {
    let faces = p.take("Cells", "faces");
    let faceoffsets = p.take("Cells", "faceoffsets");
    if let (Some(f), Some(o)) = (faces, faceoffsets) {
        return Ok(Some((o.to_i64("faceoffsets")?, f.to_i64("faces")?)));
    }
    let fc = p.take("Cells", "face_connectivity");
    let fo = p.take("Cells", "face_offsets");
    let pf = p.take("Cells", "polyhedron_to_faces");
    let po = p.take("Cells", "polyhedron_offsets");
    let (Some(fc), Some(fo), Some(pf), Some(po)) = (fc, fo, pf, po) else {
        return Ok(None);
    };
    let faces = full_offsets(
        fo.to_i64("face_offsets")?,
        None,
        fc.to_i64("face_connectivity")?,
    );
    let cells = full_offsets(
        po.to_i64("polyhedron_offsets")?,
        Some(n),
        pf.to_i64("polyhedron_to_faces")?,
    );
    fn range(c: &Cells, i: i64) -> Result<&[i64]> {
        let bad = || Error::format("polyhedron face arrays are inconsistent");
        let i = usize::try_from(i).map_err(|_| bad())?;
        let (a, b) = (
            *c.offsets.get(i).ok_or_else(bad)?,
            *c.offsets.get(i + 1).ok_or_else(bad)?,
        );
        let (a, b) = (
            usize::try_from(a).map_err(|_| bad())?,
            usize::try_from(b).map_err(|_| bad())?,
        );
        c.conn.get(a..b).ok_or_else(bad)
    }
    let mut stream = Vec::new();
    let mut offsets = Vec::with_capacity(n);
    for c in 0..n {
        let fs = range(&cells, c as i64)?;
        if fs.is_empty() {
            offsets.push(-1);
            continue;
        }
        stream.push(fs.len() as i64);
        for &f in fs {
            let ids = range(&faces, f)?;
            stream.push(ids.len() as i64);
            stream.extend_from_slice(ids);
        }
        offsets.push(stream.len() as i64);
    }
    Ok(Some((offsets, stream)))
}

fn grid_piece(mut p: Piece) -> Result<UnstructuredGrid> {
    let nc = p.tag.count("NumberOfCells")?;
    let points = p.points()?;
    let cells = p.cells("Cells", nc)?;
    let types = match p.take("Cells", "types") {
        Some(t) => t.to_i64("types")?,
        None if nc == 0 => Vec::new(),
        None => return Err(Error::format("Cells has no types array")),
    };
    if types.len() != nc || cells.len() != nc {
        return Err(Error::format(format!(
            "piece declares {nc} cells but has {} types and {} offsets",
            types.len(),
            cells.len()
        )));
    }
    let cell_types = types
        .iter()
        .map(|&t| u8::try_from(t).map_err(|_| Error::format(format!("bad cell type {t}"))))
        .collect::<Result<Vec<u8>>>()?;
    let faces = polyhedra(&mut p, nc)?;
    let (point_data, cell_data) = p.data(points.len(), nc)?;
    Ok(UnstructuredGrid {
        points,
        cell_types,
        offsets: cells.offsets,
        connectivity: cells.conn,
        faces,
        point_data,
        cell_data,
    })
}

/// Shifts the point ids inside a legacy polyhedron face stream.
fn shift_stream(s: &[i64], by: i64) -> Result<Vec<i64>> {
    let bad = || Error::format("polyhedron face stream is malformed");
    let mut out = s.to_vec();
    let mut i = 0;
    while i < s.len() {
        let nf = usize::try_from(s[i]).map_err(|_| bad())?;
        i += 1;
        for _ in 0..nf {
            let np = usize::try_from(*s.get(i).ok_or_else(bad)?).map_err(|_| bad())?;
            let end = (i + 1)
                .checked_add(np)
                .filter(|&e| e <= s.len())
                .ok_or_else(bad)?;
            for v in &mut out[i + 1..end] {
                *v += by;
            }
            i = end;
        }
    }
    Ok(out)
}

fn merge_grids(mut gs: Vec<UnstructuredGrid>) -> Result<UnstructuredGrid> {
    if gs.len() == 1 {
        return Ok(gs.remove(0));
    }
    for g in &gs {
        if let Some(m) = grid_problem(g) {
            return Err(Error::format(m));
        }
    }
    let any_faces = gs.iter().any(|g| g.faces.is_some());
    let mut out = UnstructuredGrid {
        points: Vec::new(),
        cell_types: Vec::new(),
        offsets: vec![0],
        connectivity: Vec::new(),
        faces: any_faces.then(|| (Vec::new(), Vec::new())),
        point_data: concat_tables(&gs.iter().map(|g| g.point_data.clone()).collect::<Vec<_>>())?,
        cell_data: concat_tables(&gs.iter().map(|g| g.cell_data.clone()).collect::<Vec<_>>())?,
    };
    for g in gs {
        let base = out.points.len() as i64;
        let cbase = out.connectivity.len() as i64;
        out.offsets.extend(g.offsets[1..].iter().map(|o| o + cbase));
        out.connectivity
            .extend(g.connectivity.iter().map(|i| i + base));
        if let Some((fo, fs)) = &mut out.faces {
            let sbase = fs.len() as i64;
            match &g.faces {
                Some((go, gstream)) => {
                    fo.extend(go.iter().map(|&o| if o < 0 { o } else { o + sbase }));
                    fs.extend(shift_stream(gstream, base)?);
                }
                None => fo.extend(std::iter::repeat_n(-1, g.cell_types.len())),
            }
        }
        out.cell_types.extend(g.cell_types);
        out.points.extend(g.points);
    }
    Ok(out)
}

const POLY_SECTIONS: [(&str, &str); 4] = [
    ("Verts", "NumberOfVerts"),
    ("Lines", "NumberOfLines"),
    ("Polys", "NumberOfPolys"),
    ("Strips", "NumberOfStrips"),
];

fn poly_piece(mut p: Piece) -> Result<PolyData> {
    let points = p.points()?;
    let mut secs = Vec::with_capacity(4);
    for (sec, key) in POLY_SECTIONS {
        let n = p.tag.count(key)?;
        let c = p.cells(sec, n)?;
        if c.len() != n {
            return Err(Error::format(format!(
                "{sec} declares {n} cells but has {}",
                c.len()
            )));
        }
        secs.push(c);
    }
    let total = secs.iter().map(Cells::len).sum();
    let (point_data, cell_data) = p.data(points.len(), total)?;
    let [verts, lines, polys, strips]: [Cells; 4] = secs.try_into().ok().expect("four sections");
    Ok(PolyData {
        points,
        point_data,
        verts,
        lines,
        polys,
        strips,
        cell_data,
    })
}

fn merge_polys(mut ps: Vec<PolyData>) -> Result<PolyData> {
    if ps.len() == 1 {
        return Ok(ps.remove(0));
    }
    let mut points = Vec::new();
    let mut sections: Vec<Cells> = (0..4).map(|_| Cells::new()).collect();
    let mut tables: [Vec<wireframes_core::Table>; 4] = Default::default();
    for p in &ps {
        let base = points.len() as i64;
        points.extend_from_slice(&p.points);
        let mut row = 0u32;
        for (s, c) in p.sections().into_iter().enumerate() {
            let n = c.len() as u32;
            if row + n > p.cell_data.rows() as u32 {
                return Err(Error::format("piece cell data is shorter than its cells"));
            }
            tables[s].push(p.cell_data.take(&(row..row + n).collect::<Vec<_>>()));
            row += n;
            let cbase = sections[s].conn.len() as i64;
            sections[s]
                .offsets
                .extend(c.offsets[1..].iter().map(|o| o + cbase));
            sections[s].conn.extend(c.conn.iter().map(|i| i + base));
        }
    }
    let point_data = concat_tables(&ps.iter().map(|p| p.point_data.clone()).collect::<Vec<_>>())?;
    let per_section = tables
        .iter()
        .map(|t| concat_tables(t))
        .collect::<Result<Vec<_>>>()?;
    let cell_data = concat_tables(&per_section)?;
    let [verts, lines, polys, strips]: [Cells; 4] =
        sections.try_into().ok().expect("four sections");
    Ok(PolyData {
        points,
        point_data,
        verts,
        lines,
        polys,
        strips,
        cell_data,
    })
}

pub fn read(r: &[u8]) -> Result<Data> {
    let app_at = find(r, b"<AppendedData");
    let head = &r[..app_at.unwrap_or(r.len())];
    let tags = scan(head)?;
    let file = tags
        .iter()
        .find(|t| t.name == "VTKFile" && !t.close)
        .ok_or_else(|| Error::format("not a vtk xml file: no VTKFile element"))?;
    let kind = file.attr("type").unwrap_or("").to_string();
    let big = match file.attr("byte_order") {
        None | Some("LittleEndian") => false,
        Some("BigEndian") => true,
        Some(o) => return Err(Error::format(format!("unknown byte_order {o:?}"))),
    };
    let header = match file.attr("header_type") {
        None | Some("UInt32") => 4,
        Some("UInt64") => 8,
        Some(o) => return Err(Error::format(format!("unknown header_type {o:?}"))),
    };
    let zlib = match file.attr("compressor") {
        None | Some("") => false,
        Some("vtkZLibDataCompressor") => true,
        Some(o) => {
            return Err(Error::Unsupported(format!(
                "compressor {o} is not supported; only zlib is"
            )));
        }
    };
    let mut app_start = r.len();
    let appended = match app_at {
        None => None,
        Some(at) => {
            let t = scan_tag(r, at)?;
            let raw = match t.attr("encoding").unwrap_or("raw") {
                "raw" => true,
                "base64" => false,
                o => {
                    return Err(Error::format(format!(
                        "unknown AppendedData encoding {o:?}"
                    )));
                }
            };
            let mut j = t.after;
            while j < r.len() && r[j].is_ascii_whitespace() {
                j += 1;
            }
            if r.get(j) != Some(&b'_') {
                return Err(Error::format("AppendedData does not start with '_'"));
            }
            app_start = j + 1;
            let data = &r[j + 1..];
            let data = if raw {
                data
            } else {
                let end = data
                    .windows(b"</AppendedData".len())
                    .rposition(|w| w == b"</AppendedData")
                    .unwrap_or(data.len());
                &data[..end]
            };
            Some((data, raw))
        }
    };
    let mut offsets: Vec<usize> = tags
        .iter()
        .filter(|t| t.name == "DataArray" && t.attr("format") == Some("appended"))
        .filter_map(|t| t.attr("offset")?.trim().parse().ok())
        .collect();
    offsets.sort_unstable();
    offsets.dedup();
    let ctx = Ctx {
        big,
        header,
        zlib,
        appended,
        app_start,
        offsets,
    };
    let mut pieces: Vec<Piece> = Vec::new();
    let mut section: Option<String> = None;
    const SECTIONS: [&str; 9] = [
        "PointData",
        "CellData",
        "Points",
        "Cells",
        "Verts",
        "Lines",
        "Strips",
        "Polys",
        "FieldData",
    ];
    for (k, t) in tags.iter().enumerate() {
        if SECTIONS.contains(&t.name.as_str()) {
            section = (!t.close && !t.empty).then(|| t.name.clone());
            continue;
        }
        match (t.name.as_str(), t.close) {
            ("Piece", false) => pieces.push(Piece {
                tag: Tag {
                    name: t.name.clone(),
                    attrs: t.attrs.clone(),
                    close: false,
                    empty: t.empty,
                    start: t.start,
                    after: t.after,
                },
                arrays: Vec::new(),
            }),
            ("DataArray", false) => {
                let Some(sec) = section.as_deref().filter(|s| *s != "FieldData") else {
                    continue;
                };
                let content = if t.empty {
                    &head[0..0]
                } else {
                    let next = tags.get(k + 1).ok_or_else(|| {
                        Error::format(format!("DataArray at byte {} is not closed", t.start))
                    })?;
                    &head[t.after..next.start]
                };
                let a = decode_array(t, content, &ctx, sec == "Points")?;
                let piece = pieces
                    .last_mut()
                    .ok_or_else(|| Error::format("DataArray outside a Piece"))?;
                piece.arrays.push((sec.to_string(), a));
            }
            _ => {}
        }
    }
    if !r.trim_ascii_end().ends_with(b"</VTKFile>") {
        return Err(Error::format(format!(
            "file ends at byte {} without </VTKFile>; it is truncated",
            r.len()
        )));
    }
    match kind.as_str() {
        "UnstructuredGrid" => {
            let gs = pieces
                .into_iter()
                .map(grid_piece)
                .collect::<Result<Vec<_>>>()?;
            if gs.is_empty() {
                return Ok(Data::default());
            }
            grid_to_data(merge_grids(gs)?)
        }
        "PolyData" => {
            let ps = pieces
                .into_iter()
                .map(poly_piece)
                .collect::<Result<Vec<_>>>()?;
            if ps.is_empty() {
                return Ok(Data::default());
            }
            poly_to_data(merge_polys(ps)?)
        }
        other => Err(Error::Unsupported(format!(
            "vtk xml type {other:?} is not supported; use UnstructuredGrid or PolyData"
        ))),
    }
}

struct XmlOut {
    body: String,
    app: Vec<u8>,
    zlib: bool,
}

impl XmlOut {
    fn array(&mut self, name: &str, ncomp: usize, v: &Values) -> Result<()> {
        let (ty, raw): (&str, Vec<u8>) = match v {
            Values::F64(v) => ("Float64", v.iter().flat_map(|x| x.to_le_bytes()).collect()),
            Values::I64(v) => ("Int64", v.iter().flat_map(|x| x.to_le_bytes()).collect()),
            Values::U8(v) => ("UInt8", v.clone()),
        };
        let offset = self.app.len();
        let n = raw.len() as u64;
        if self.zlib {
            if raw.is_empty() {
                [0u64, 0, 0]
                    .iter()
                    .for_each(|h| self.app.extend(h.to_le_bytes()));
            } else {
                let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
                enc.write_all(&raw)?;
                let c = enc.finish()?;
                [1u64, n, n, c.len() as u64]
                    .iter()
                    .for_each(|h| self.app.extend(h.to_le_bytes()));
                self.app.extend(c);
            }
        } else {
            self.app.extend(n.to_le_bytes());
            self.app.extend(raw);
        }
        self.body.push_str(&format!(
            "        <DataArray type=\"{ty}\" Name=\"{}\" NumberOfComponents=\"{ncomp}\" format=\"appended\" offset=\"{offset}\"/>\n",
            escape(name)
        ));
        Ok(())
    }

    fn arrays(&mut self, tag: &str, arrays: &[Array]) -> Result<()> {
        if arrays.is_empty() {
            return Ok(());
        }
        self.body.push_str(&format!("      <{tag}>\n"));
        for a in arrays {
            self.array(&a.name, a.ncomp, &a.values)?;
        }
        self.body.push_str(&format!("      </{tag}>\n"));
        Ok(())
    }

    fn section(&mut self, tag: &str, f: impl FnOnce(&mut XmlOut) -> Result<()>) -> Result<()> {
        self.body.push_str(&format!("      <{tag}>\n"));
        f(self)?;
        self.body.push_str(&format!("      </{tag}>\n"));
        Ok(())
    }

    fn cells(&mut self, tag: &str, offsets: &[i64], conn: &[i64]) -> Result<()> {
        self.section(tag, |w| {
            w.array("connectivity", 1, &Values::I64(conn.to_vec()))?;
            w.array("offsets", 1, &Values::I64(offsets[1..].to_vec()))
        })
    }

    fn points(&mut self, p: &[[f64; 3]]) -> Result<()> {
        self.section("Points", |w| {
            w.array(
                "Points",
                3,
                &Values::F64(p.iter().flatten().copied().collect()),
            )
        })
    }

    fn finish(self, kind: &str, piece: &str) -> Vec<u8> {
        let compressor = if self.zlib {
            " compressor=\"vtkZLibDataCompressor\""
        } else {
            ""
        };
        let mut out = format!(
            "<?xml version=\"1.0\"?>\n<VTKFile type=\"{kind}\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"{compressor}>\n  <{kind}>\n    <Piece {piece}>\n{}    </Piece>\n  </{kind}>\n  <AppendedData encoding=\"raw\">\n   _",
            self.body
        )
        .into_bytes();
        out.extend(self.app);
        out.extend_from_slice(b"\n  </AppendedData>\n</VTKFile>\n");
        out
    }
}

/// Writes `.vtu`: the grid, or the mesh as triangle and polygon cells.
pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.polylines.is_some() || data.points.is_some() {
        return Err(Error::Unsupported(
            ".vtu holds a mesh or a grid; write polylines and points to .vtp".into(),
        ));
    }
    let owned;
    let g = match (&data.grid, &data.mesh) {
        (Some(_), Some(_)) => {
            return Err(Error::Unsupported(
                ".vtu holds a mesh or a grid, not both".into(),
            ));
        }
        (Some(g), None) => {
            if let Some(m) = grid_problem(g) {
                return Err(Error::invalid(m));
            }
            g
        }
        (None, Some(m)) => {
            owned = mesh_to_grid(m)?;
            &owned
        }
        (None, None) => {
            owned = mesh_to_grid(&Mesh::empty())?;
            &owned
        }
    };
    let pa = flatten(&g.point_data, opts, "vtu")?;
    let ca = flatten(&g.cell_data, opts, "vtu")?;
    let mut w = XmlOut {
        body: String::new(),
        app: Vec::new(),
        zlib: opts.zlib,
    };
    w.arrays("PointData", &pa)?;
    w.arrays("CellData", &ca)?;
    w.points(&g.points)?;
    w.section("Cells", |w| {
        w.array("connectivity", 1, &Values::I64(g.connectivity.clone()))?;
        w.array("offsets", 1, &Values::I64(g.offsets[1..].to_vec()))?;
        w.array("types", 1, &Values::U8(g.cell_types.clone()))?;
        if let Some((fo, fs)) = &g.faces {
            w.array("faces", 1, &Values::I64(fs.clone()))?;
            w.array("faceoffsets", 1, &Values::I64(fo.clone()))?;
        }
        Ok(())
    })?;
    let piece = format!(
        "NumberOfPoints=\"{}\" NumberOfCells=\"{}\"",
        g.points.len(),
        g.cell_types.len()
    );
    Ok(w.finish("UnstructuredGrid", &piece))
}

/// Writes `.vtp`: mesh as polys, polylines as lines, points as verts.
pub fn write_vtp(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.grid.is_some() {
        return Err(Error::Unsupported(
            ".vtp cannot hold an unstructured grid; write .vtu".into(),
        ));
    }
    let pd = data_to_poly(data)?;
    let pa = flatten(&pd.point_data, opts, "vtp")?;
    let ca = flatten(&pd.cell_data, opts, "vtp")?;
    let mut w = XmlOut {
        body: String::new(),
        app: Vec::new(),
        zlib: opts.zlib,
    };
    w.arrays("PointData", &pa)?;
    w.arrays("CellData", &ca)?;
    w.points(&pd.points)?;
    for (tag, c) in [
        ("Verts", &pd.verts),
        ("Lines", &pd.lines),
        ("Polys", &pd.polys),
    ] {
        if c.len() > 0 {
            w.cells(tag, &c.offsets, &c.conn)?;
        }
    }
    let piece = format!(
        "NumberOfPoints=\"{}\" NumberOfVerts=\"{}\" NumberOfLines=\"{}\" NumberOfStrips=\"0\" NumberOfPolys=\"{}\"",
        pd.points.len(),
        pd.verts.len(),
        pd.lines.len(),
        pd.polys.len()
    );
    Ok(w.finish("PolyData", &piece))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle(zlib: bool) -> Vec<u8> {
        let m = Mesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let data = Data {
            mesh: Some(m),
            ..Data::default()
        };
        let opts = WriteOptions {
            zlib,
            ..WriteOptions::default()
        };
        write(&data, &opts).unwrap()
    }

    #[test]
    fn a_truncated_appended_array_names_its_byte_offset() {
        for zlib in [false, true] {
            let full = triangle(zlib);
            let app = find(&full, b"<AppendedData").unwrap();
            let start = app + find(&full[app..], b"_").unwrap() + 1;
            for cut in [start, start + 3, full.len() - 30, full.len() - 1 - 30] {
                let e = read(&full[..cut]).unwrap_err().to_string();
                let at = e
                    .split("at byte ")
                    .nth(1)
                    .map(|r| r.split(|c: char| !c.is_ascii_digit()).next().unwrap_or(""))
                    .and_then(|n| n.parse::<usize>().ok());
                assert!(at.is_some_and(|a| a >= start && a <= full.len()), "{e}");
            }
        }
    }

    #[test]
    fn a_non_finite_point_names_its_array_byte_offset() {
        let mut b = triangle(false);
        let one = 1.0f64.to_le_bytes();
        let k = b.windows(8).position(|w| w == one).unwrap();
        b[k..k + 8].copy_from_slice(&f64::NAN.to_le_bytes());
        let e = read(&b).unwrap_err().to_string();
        assert!(
            e.contains("with data at byte") && e.ends_with("point 1 coordinate 0 is not finite"),
            "{e}"
        );
    }
}
