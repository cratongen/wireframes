use std::io::Write;

use wireframes_core::{Error, Mesh, Result};

use crate::common::{polygon_corners, polygon_groups};
use crate::obj::{
    Polygons, check_no_columns, parse_coord, parse_i64, push_f64, strip_comment, tokens, trim,
    unsupported_column,
};
use crate::{Data, WriteOptions};

struct Lines<'a> {
    rest: &'a [u8],
    no: u64,
}

impl<'a> Lines<'a> {
    fn next_content(&mut self) -> Option<(u64, &'a [u8])> {
        while !self.rest.is_empty() {
            self.no += 1;
            let (line, next) = match self.rest.iter().position(|&b| b == b'\n') {
                Some(i) => (&self.rest[..i], &self.rest[i + 1..]),
                None => (self.rest, &self.rest[self.rest.len()..]),
            };
            self.rest = next;
            let line = trim(strip_comment(line));
            if !line.is_empty() {
                return Some((self.no, line));
            }
        }
        None
    }
}

struct Keyword {
    homogeneous: bool,
    ndim: bool,
}

fn keyword(tok: &[u8], line: u64) -> Result<Keyword> {
    let Some(mut p) = tok.strip_suffix(b"OFF") else {
        return Err(Error::format_at(
            line,
            format!(
                "expected an OFF keyword, found {:?}",
                String::from_utf8_lossy(tok)
            ),
        ));
    };
    let mut kw = Keyword {
        homogeneous: false,
        ndim: false,
    };
    for prefix in [&b"ST"[..], b"C", b"N", b"4", b"n"] {
        if let Some(rest) = p.strip_prefix(prefix) {
            p = rest;
            match prefix {
                b"4" => kw.homogeneous = true,
                b"n" => kw.ndim = true,
                _ => {}
            }
        }
    }
    if !p.is_empty() {
        return Err(Error::format_at(
            line,
            format!("unknown OFF keyword {:?}", String::from_utf8_lossy(tok)),
        ));
    }
    Ok(kw)
}

fn count(tok: &[u8], line: u64, what: &str) -> Result<usize> {
    let n = parse_i64(tok, line)?;
    usize::try_from(n).map_err(|_| Error::format_at(line, format!("{what} count {n} is negative")))
}

pub fn read(r: &[u8]) -> Result<Data> {
    let mut lines = Lines { rest: r, no: 0 };
    let (no, first) = lines
        .next_content()
        .ok_or_else(|| Error::format("empty off file"))?;
    let mut toks = tokens(first);
    let kw = keyword(toks.next().unwrap_or_default(), no)?;
    let mut header: Vec<(u64, &[u8])> = toks.map(|t| (no, t)).collect();
    let need = if kw.ndim { 4 } else { 3 };
    if header
        .first()
        .is_some_and(|t| t.1.eq_ignore_ascii_case(b"BINARY"))
    {
        return Err(Error::Unsupported("binary off is not supported".into()));
    }
    while header.len() < need {
        let (no, line) = lines
            .next_content()
            .ok_or_else(|| Error::format("off header ends before the element counts"))?;
        if line.eq_ignore_ascii_case(b"BINARY") {
            return Err(Error::Unsupported("binary off is not supported".into()));
        }
        header.extend(tokens(line).map(|t| (no, t)));
    }
    let mut h = header.into_iter();
    let mut next = || h.next().expect("length checked");
    if kw.ndim {
        let (no, t) = next();
        let dim = parse_i64(t, no)?;
        if dim != 3 {
            return Err(Error::Unsupported(format!(
                "nOFF with dimension {dim}; only 3 is supported"
            )));
        }
    }
    let (no, t) = next();
    let nv = count(t, no, "vertex")?;
    let (no, t) = next();
    let nf = count(t, no, "face")?;

    let width = if kw.homogeneous { 4 } else { 3 };
    let mut verts = Vec::with_capacity(nv.min(r.len() / 6));
    for _ in 0..nv {
        let (no, line) = lines.next_content().ok_or_else(|| {
            Error::format(format!("file ends after {} of {nv} vertices", verts.len()))
        })?;
        let mut p = [0.0; 4];
        let mut it = tokens(line);
        for c in p.iter_mut().take(width) {
            let t = it
                .next()
                .ok_or_else(|| Error::format_at(no, format!("vertex needs {width} values")))?;
            *c = parse_coord(t, no)?;
        }
        if kw.homogeneous {
            if p[3] == 0.0 {
                return Err(Error::format_at(no, "homogeneous coordinate w is zero"));
            }
            let q = [p[0] / p[3], p[1] / p[3], p[2] / p[3]];
            if q.iter().any(|x| !x.is_finite()) {
                return Err(Error::format_at(no, "vertex divided by w is not finite"));
            }
            verts.push(q);
        } else {
            verts.push([p[0], p[1], p[2]]);
        }
    }

    let mut polys = Polygons::default();
    let mut corners = Vec::new();
    for f in 0..nf {
        let (no, line) = lines
            .next_content()
            .ok_or_else(|| Error::format(format!("file ends after {f} of {nf} faces")))?;
        let mut it = tokens(line);
        let n = count(it.next().unwrap_or_default(), no, "corner")?;
        corners.clear();
        for _ in 0..n {
            let t = it.next().ok_or_else(|| {
                Error::format_at(no, format!("face lists fewer than {n} corners"))
            })?;
            let i = parse_i64(t, no)?;
            if i < 0 || i as u64 >= nv as u64 {
                return Err(Error::format_at(
                    no,
                    format!("face index {i} is out of range for {nv} vertices"),
                ));
            }
            corners.push(i as u32);
        }
        polys.push(&corners);
    }

    let mut data = Data::default();
    polys.warn(&mut data.warnings);
    let (tris, col) = polys.finish();
    let mut mesh = Mesh::new(verts, tris)?;
    if let Some(c) = col {
        mesh = mesh.with_face_column("polygon", c)?;
    }
    data.mesh = Some(mesh);
    Ok(data)
}

pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.polylines.is_some() || data.points.is_some() || data.grid.is_some() {
        return Err(Error::Unsupported(
            "off stores only a mesh, not polylines, points or a grid".into(),
        ));
    }
    let empty = Mesh::empty();
    let m = data.mesh.as_ref().unwrap_or(&empty);
    check_no_columns("off", "vertex", m.vertex_attributes(), opts)?;
    for name in m.face_attributes().names() {
        if name != "polygon" {
            unsupported_column("off", "face", name, opts)?;
        }
    }
    let groups = polygon_groups(m)?;
    let mut out = Vec::with_capacity(m.n_vertices() * 40 + m.n_faces() * 20 + 32);
    writeln!(out, "OFF\n{} {} 0", m.n_vertices(), groups.len())?;
    for p in m.vertices() {
        push_f64(&mut out, p[0])?;
        out.push(b' ');
        push_f64(&mut out, p[1])?;
        out.push(b' ');
        push_f64(&mut out, p[2])?;
        out.push(b'\n');
    }
    let mut corners = Vec::new();
    for g in &groups {
        polygon_corners(m.triangles(), g, &mut corners);
        write!(out, "{}", corners.len())?;
        for &i in &corners {
            write!(out, " {i}")?;
        }
        out.push(b'\n');
    }
    Ok(out)
}
