use std::collections::HashMap;
use std::io::Write;

use drusa_core::{Column, ColumnData, Error, Mesh, PointSet, Polylines, Result, Table};

use crate::common::{decode_name, encode_name, polygon_corners, polygon_groups};
use crate::{Data, WriteOptions};

pub(crate) fn is_ws(b: &u8) -> bool {
    b.is_ascii_whitespace()
}

pub(crate) fn tokens(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(is_ws).filter(|t| !t.is_empty())
}

pub(crate) fn trim(mut s: &[u8]) -> &[u8] {
    while let [first, rest @ ..] = s
        && first.is_ascii_whitespace()
    {
        s = rest;
    }
    while let [rest @ .., last] = s
        && last.is_ascii_whitespace()
    {
        s = rest;
    }
    s
}

pub(crate) fn strip_comment(line: &[u8]) -> &[u8] {
    match line.iter().position(|&b| b == b'#') {
        Some(i) => &line[..i],
        None => line,
    }
}

/// Calls `f(line_number, line)` for each line without its terminator. With
/// `continuation`, a line ending in `\` outside a `#` comment joins the next one;
/// the number is that of the first physical line.
pub(crate) fn for_each_line(
    data: &[u8],
    continuation: bool,
    mut f: impl FnMut(u64, &[u8]) -> Result<()>,
) -> Result<()> {
    let mut joined: Vec<u8> = Vec::new();
    let mut joined_line = 0u64;
    let mut line_no = 0u64;
    let mut rest = data;
    while !rest.is_empty() {
        line_no += 1;
        let (mut line, next) = match rest.iter().position(|&b| b == b'\n') {
            Some(i) => (&rest[..i], &rest[i + 1..]),
            None => (rest, &rest[rest.len()..]),
        };
        rest = next;
        if let [head @ .., b'\r'] = line {
            line = head;
        }
        if continuation {
            let t = trim(strip_comment(line));
            if let [head @ .., b'\\'] = t {
                if joined.is_empty() {
                    joined_line = line_no;
                }
                joined.extend_from_slice(head);
                joined.push(b' ');
                continue;
            }
            if !joined.is_empty() {
                joined.extend_from_slice(line);
                f(joined_line, &joined)?;
                joined.clear();
                continue;
            }
        }
        f(line_no, line)?;
    }
    if !joined.is_empty() {
        f(joined_line, &joined)?;
    }
    Ok(())
}

pub(crate) fn as_str(tok: &[u8]) -> &str {
    std::str::from_utf8(tok).unwrap_or("")
}

pub(crate) fn parse_f64(tok: &[u8], line: u64) -> Result<f64> {
    as_str(tok).parse::<f64>().map_err(|_| {
        Error::format_at(
            line,
            format!(
                "expected a number, found {:?}",
                String::from_utf8_lossy(tok)
            ),
        )
    })
}

pub(crate) fn parse_coord(tok: &[u8], line: u64) -> Result<f64> {
    let x = parse_f64(tok, line)?;
    if !x.is_finite() {
        return Err(Error::format_at(
            line,
            format!("coordinate {x} is not finite"),
        ));
    }
    Ok(x)
}

pub(crate) fn parse_i64(tok: &[u8], line: u64) -> Result<i64> {
    as_str(tok).parse::<i64>().map_err(|_| {
        Error::format_at(
            line,
            format!(
                "expected an integer, found {:?}",
                String::from_utf8_lossy(tok)
            ),
        )
    })
}

pub(crate) fn push_f64(out: &mut Vec<u8>, x: f64) -> Result<()> {
    write!(out, "{x:?}")?;
    Ok(())
}

/// Fans polygons from their first corner. Polygon ids are 0-based in push order;
/// the `polygon` column exists only when some polygon has more than 3 corners.
#[derive(Default)]
pub(crate) struct Polygons {
    pub tris: Vec<[u32; 3]>,
    ids: Vec<i64>,
    n: i64,
    big: bool,
    pub skipped: usize,
}

impl Polygons {
    /// Returns the number of triangles added; fewer than 3 corners adds none.
    pub fn push(&mut self, c: &[u32]) -> usize {
        if c.len() < 3 {
            self.skipped += 1;
            return 0;
        }
        self.big |= c.len() > 3;
        for k in 1..c.len() - 1 {
            self.tris.push([c[0], c[k], c[k + 1]]);
            self.ids.push(self.n);
        }
        self.n += 1;
        c.len() - 2
    }

    pub fn warn(&self, warnings: &mut Vec<String>) {
        if self.skipped > 0 {
            warnings.push(format!(
                "skipped {} faces with fewer than 3 corners",
                self.skipped
            ));
        }
    }

    pub fn finish(self) -> (Vec<[u32; 3]>, Option<Column>) {
        let col = self.big.then(|| Column::i64(self.ids));
        (self.tris, col)
    }
}

pub(crate) fn unsupported_column(
    format: &str,
    kind: &str,
    name: &str,
    opts: &WriteOptions,
) -> Result<()> {
    if opts.skip_unsupported_columns {
        Ok(())
    } else {
        Err(Error::Unsupported(format!(
            "{format} cannot store {kind} column {name:?}; set skip_unsupported_columns to drop it"
        )))
    }
}

pub(crate) fn check_no_columns(
    format: &str,
    kind: &str,
    t: &Table,
    opts: &WriteOptions,
) -> Result<()> {
    for name in t.names() {
        unsupported_column(format, kind, name, opts)?;
    }
    Ok(())
}

pub(crate) fn utf8_values<'a>(
    t: &'a Table,
    name: &str,
    format: &str,
) -> Result<Option<(&'a [String], &'a Column)>> {
    match t.column(name) {
        None => Ok(None),
        Some(c) => match &c.data {
            ColumnData::Utf8(v) => Ok(Some((v, c))),
            _ => Err(Error::Unsupported(format!(
                "{format} face column {name:?} must be utf8, found {}",
                c.type_name()
            ))),
        },
    }
}

fn index(tok: &[u8], n: usize, line: u64) -> Result<u32> {
    let field = tok.split(|&b| b == b'/').next().unwrap_or(tok);
    let i = parse_i64(field, line)?;
    let r = if i > 0 {
        i - 1
    } else if i < 0 {
        n as i64 + i
    } else {
        return Err(Error::format_at(
            line,
            "index 0 is not valid; obj indices start at 1",
        ));
    };
    if r < 0 || r > u32::MAX as i64 {
        return Err(Error::format_at(
            line,
            format!("index {i} is out of range with {n} vertices defined so far"),
        ));
    }
    Ok(r as u32)
}

const NULL_ID: u32 = u32::MAX;

struct Names {
    ids: HashMap<Vec<u8>, u32>,
    names: Vec<String>,
    buf: Vec<u8>,
}

impl Names {
    fn intern(&mut self, toks: &[u8]) -> u32 {
        self.buf.clear();
        for t in tokens(toks) {
            if !self.buf.is_empty() {
                self.buf.push(b' ');
            }
            self.buf.extend_from_slice(t);
        }
        if let Some(&id) = self.ids.get(&self.buf) {
            return id;
        }
        let name = match self.buf.as_slice() {
            b"" => Some("default".into()),
            raw => decode_name(raw),
        };
        let id = match name {
            None => NULL_ID,
            Some(n) => {
                self.names.push(n);
                self.names.len() as u32 - 1
            }
        };
        self.ids.insert(self.buf.clone(), id);
        id
    }

    fn column(&self, ids: &[u32]) -> Result<Column> {
        let values = ids
            .iter()
            .map(|&i| self.names.get(i as usize).cloned().unwrap_or_default())
            .collect();
        let c = Column::utf8(values);
        if ids.contains(&NULL_ID) {
            c.with_valid(ids.iter().map(|&i| i != NULL_ID).collect())
        } else {
            Ok(c)
        }
    }
}

pub fn read(r: &[u8]) -> Result<Data> {
    let mut verts: Vec<[f64; 3]> = Vec::new();
    let mut polys = Polygons::default();
    let mut corners: Vec<u32> = Vec::new();
    let mut names = Names {
        ids: HashMap::new(),
        names: Vec::new(),
        buf: Vec::new(),
    };
    let default = names.intern(b"");
    let (mut obj, mut grp) = (default, default);
    let (mut saw_o, mut saw_g) = (false, false);
    let (mut face_obj, mut face_grp): (Vec<u32>, Vec<u32>) = (Vec::new(), Vec::new());
    let mut line_idx: Vec<u32> = Vec::new();
    let mut line_off: Vec<usize> = vec![0];
    let mut point_idx: Vec<u32> = Vec::new();
    let mut short_lines = 0usize;
    let mut unknown = 0usize;

    for_each_line(r, true, |no, line| {
        let line = strip_comment(line);
        let mut it = tokens(line);
        let Some(kw) = it.next() else { return Ok(()) };
        match kw {
            b"v" => {
                let mut p = [0.0; 3];
                for c in &mut p {
                    let t = it
                        .next()
                        .ok_or_else(|| Error::format_at(no, "vertex needs 3 coordinates"))?;
                    *c = parse_coord(t, no)?;
                }
                verts.push(p);
            }
            b"f" | b"l" | b"p" => {
                corners.clear();
                for t in it {
                    corners.push(index(t, verts.len(), no)?);
                }
                match kw {
                    b"f" => {
                        let n = polys.push(&corners);
                        face_obj.extend(std::iter::repeat_n(obj, n));
                        face_grp.extend(std::iter::repeat_n(grp, n));
                    }
                    b"l" if corners.len() < 2 => short_lines += 1,
                    b"l" => {
                        line_idx.extend_from_slice(&corners);
                        line_off.push(line_idx.len());
                    }
                    _ => point_idx.extend_from_slice(&corners),
                }
            }
            b"o" | b"g" => {
                let rest = &line[line
                    .iter()
                    .position(|b| b.is_ascii_alphabetic())
                    .unwrap_or(0)
                    + 1..];
                let id = names.intern(rest);
                if kw == b"o" {
                    saw_o = true;
                    obj = id;
                } else {
                    saw_g = true;
                    grp = id;
                }
            }
            b"vt" | b"vn" | b"vp" | b"s" | b"mtllib" | b"usemtl" => {}
            _ => unknown += 1,
        }
        Ok(())
    })?;

    let n = verts.len();
    let check = |i: &u32| -> Result<()> {
        if *i as usize >= n {
            Err(Error::format(format!(
                "index {} is out of range; the file has {n} vertices",
                *i as u64 + 1
            )))
        } else {
            Ok(())
        }
    };
    polys.tris.iter().flatten().try_for_each(check)?;
    line_idx.iter().try_for_each(check)?;
    point_idx.iter().try_for_each(check)?;

    let mut warnings = Vec::new();
    polys.warn(&mut warnings);
    if short_lines > 0 {
        warnings.push(format!(
            "skipped {short_lines} lines with fewer than 2 vertices"
        ));
    }
    if unknown > 0 {
        warnings.push(format!("ignored {unknown} unknown statements"));
    }

    let mut data = Data::default();
    let (mut tris, polygon_col) = polys.finish();
    let n_faces = tris.len();

    let mut remap: Option<Vec<u32>> = None;
    if !line_idx.is_empty() || !point_idx.is_empty() {
        let mut in_mesh = vec![true; n];
        for &i in line_idx.iter().chain(&point_idx) {
            in_mesh[i as usize] = false;
        }
        for &i in tris.iter().flatten() {
            in_mesh[i as usize] = true;
        }
        let mut map = vec![u32::MAX; n];
        let mut k = 0u32;
        for i in 0..n {
            if in_mesh[i] {
                map[i] = k;
                k += 1;
            }
        }
        remap = Some(map);
    }
    let mesh_verts: Vec<[f64; 3]> = match &remap {
        None => verts.clone(),
        Some(map) => (0..n)
            .filter(|&i| map[i] != u32::MAX)
            .map(|i| verts[i])
            .collect(),
    };
    if let Some(map) = &remap {
        for t in &mut tris {
            for i in t.iter_mut() {
                *i = map[*i as usize];
            }
        }
    }
    if n_faces > 0 || !mesh_verts.is_empty() {
        let mut mesh = Mesh::new(mesh_verts, tris)?;
        if saw_o {
            mesh = mesh.with_face_column("object", names.column(&face_obj)?)?;
        }
        if saw_g {
            mesh = mesh.with_face_column("group", names.column(&face_grp)?)?;
        }
        if let Some(c) = polygon_col {
            mesh = mesh.with_face_column("polygon", c)?;
        }
        data.mesh = Some(mesh);
    }

    if line_off.len() > 1 {
        let mut parts = Vec::with_capacity(line_off.len() - 1);
        let mut closed = Vec::with_capacity(line_off.len() - 1);
        for w in line_off.windows(2) {
            let mut ids = &line_idx[w[0]..w[1]];
            let ring = ids.len() >= 4 && ids[0] == ids[ids.len() - 1];
            if ring {
                ids = &ids[..ids.len() - 1];
            }
            parts.push(ids.iter().map(|&i| verts[i as usize]).collect());
            closed.push(ring);
        }
        data.polylines = Some(Polylines::from_parts(parts, closed, None)?);
    }
    if !point_idx.is_empty() {
        data.points = Some(PointSet::new(
            point_idx.iter().map(|&i| verts[i as usize]).collect(),
        )?);
    }
    data.warnings = warnings;
    Ok(data)
}

fn push_vertex(out: &mut Vec<u8>, p: &[f64; 3]) -> Result<()> {
    out.extend_from_slice(b"v ");
    push_f64(out, p[0])?;
    out.push(b' ');
    push_f64(out, p[1])?;
    out.push(b' ');
    push_f64(out, p[2])?;
    out.push(b'\n');
    Ok(())
}

fn push_name(out: &mut Vec<u8>, kw: &[u8], name: Option<&str>) {
    out.extend_from_slice(kw);
    out.push(b' ');
    encode_name(name, true, out);
    out.push(b'\n');
}

pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.grid.is_some() {
        return Err(Error::Unsupported(
            "obj cannot store an unstructured grid".into(),
        ));
    }
    let mut out = Vec::new();
    let mut base = 0usize;
    if let Some(m) = &data.mesh {
        check_no_columns("obj", "vertex", m.vertex_attributes(), opts)?;
        let fa = m.face_attributes();
        let objects = utf8_values(fa, "object", "obj")?;
        let groups = utf8_values(fa, "group", "obj")?;
        for name in fa.names() {
            if !matches!(name.as_str(), "object" | "group" | "polygon") {
                unsupported_column("obj", "face", name, opts)?;
            }
        }
        out.reserve(m.n_vertices() * 40 + m.n_faces() * 24);
        for p in m.vertices() {
            push_vertex(&mut out, p)?;
        }
        let tris = m.triangles();
        let mut corners = Vec::new();
        let (mut last_o, mut last_g): (Option<Option<&str>>, Option<Option<&str>>) = (None, None);
        for g in polygon_groups(m)? {
            let f = g.start;
            for (vals, last, kw) in [(&objects, &mut last_o, b"o"), (&groups, &mut last_g, b"g")] {
                if let Some((v, c)) = vals {
                    let name = c.is_valid(f).then(|| v[f].as_str());
                    if *last != Some(name) {
                        push_name(&mut out, kw, name);
                        *last = Some(name);
                    }
                }
            }
            polygon_corners(tris, &g, &mut corners);
            out.push(b'f');
            for &i in &corners {
                write!(out, " {}", i as u64 + 1)?;
            }
            out.push(b'\n');
        }
        base = m.n_vertices();
    }
    if let Some(pl) = &data.polylines {
        check_no_columns("obj", "polyline vertex", pl.vertex_attributes(), opts)?;
        check_no_columns("obj", "polyline feature", pl.attributes(), opts)?;
        if !opts.skip_unsupported_columns && (0..pl.len()).any(|f| pl.feature_parts(f).len() != 1) {
            return Err(Error::Unsupported(
                "obj stores one part per polyline feature; set skip_unsupported_columns to write parts as features".into(),
            ));
        }
        for p in pl.vertices() {
            push_vertex(&mut out, p)?;
        }
        for p in 0..pl.n_parts() {
            let r = pl.part_range(p);
            out.push(b'l');
            for i in r.clone() {
                write!(out, " {}", base + i + 1)?;
            }
            if pl.closed()[p] {
                write!(out, " {}", base + r.start + 1)?;
            }
            out.push(b'\n');
        }
        base += pl.n_vertices();
    }
    if let Some(ps) = &data.points {
        check_no_columns("obj", "point", ps.attributes(), opts)?;
        for p in ps.coords() {
            push_vertex(&mut out, p)?;
        }
        for i in 0..ps.len() {
            writeln!(out, "p {}", base + i + 1)?;
        }
    }
    Ok(out)
}
