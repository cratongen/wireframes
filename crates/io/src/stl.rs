use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::io::Write;

use wireframes_core::{Column, ColumnData, Error, Mesh, Result};

use crate::common::{decode_name, encode_name};
use crate::obj::{
    as_str, check_no_columns, parse_coord, push_f64, trim, unsupported_column, utf8_values,
};
use crate::{Data, WriteOptions};

/// Exact weld: corners with bitwise equal coordinates (after mapping `-0.0` to
/// `+0.0`) share one vertex, numbered and valued by first appearance.
#[derive(Default)]
struct Welder {
    ids: HashMap<[u64; 3], u32>,
    verts: Vec<[f64; 3]>,
    tris: Vec<[u32; 3]>,
}

impl Welder {
    fn corner(&mut self, p: [f64; 3]) -> u32 {
        let key = p.map(|x| if x == 0.0 { 0u64 } else { x.to_bits() });
        match self.ids.entry(key) {
            Entry::Occupied(e) => *e.get(),
            Entry::Vacant(e) => {
                let id = self.verts.len() as u32;
                self.verts.push(p);
                e.insert(id);
                id
            }
        }
    }

    fn polygon(&mut self, c: &[[f64; 3]]) {
        let a = self.corner(c[0]);
        let mut prev = self.corner(c[1]);
        for p in &c[2..] {
            let next = self.corner(*p);
            self.tris.push([a, prev, next]);
            prev = next;
        }
    }
}

fn binary_count(r: &[u8]) -> Option<u64> {
    (r.len() >= 84).then(|| u32::from_le_bytes([r[80], r[81], r[82], r[83]]) as u64)
}

fn starts_with_solid(r: &[u8]) -> bool {
    let t = trim(&r[..r.len().min(512)]);
    t.len() >= 5 && t[..5].eq_ignore_ascii_case(b"solid")
}

pub fn read(r: &[u8]) -> Result<Data> {
    let len = r.len() as u64;
    let padded = |n: u64| {
        let extra = len - 84 - 50 * n;
        read_binary(
            r,
            n as usize,
            vec![format!(
                "binary stl declares {n} triangles; ignored {extra} trailing bytes"
            )],
        )
    };
    match binary_count(r) {
        Some(n) if 84 + 50 * n == len => read_binary(r, n as usize, Vec::new()),
        count if starts_with_solid(r) => match (read_ascii(r), count) {
            (Ok(d), _) => Ok(d),
            (Err(_), Some(n)) if 84 + 50 * n < len => padded(n),
            (Err(e), _) => Err(e),
        },
        Some(n) if 84 + 50 * n < len => padded(n),
        Some(n) => Err(Error::format(format!(
            "binary stl declares {n} triangles ({} bytes) but the file has {len} bytes",
            84 + 50 * n
        ))),
        None => Err(Error::format(
            "not an stl file: shorter than a binary header and does not start with \"solid\"",
        )),
    }
}

fn f32_at(b: &[u8], i: usize) -> f64 {
    f32::from_le_bytes([b[i], b[i + 1], b[i + 2], b[i + 3]]) as f64
}

fn read_binary(r: &[u8], n: usize, warnings: Vec<String>) -> Result<Data> {
    let mut w = Welder {
        tris: Vec::with_capacity(n),
        ..Welder::default()
    };
    let mut attr = Vec::with_capacity(n);
    for (f, rec) in r[84..84 + 50 * n].chunks_exact(50).enumerate() {
        let mut c = [[0.0; 3]; 3];
        for (k, p) in c.iter_mut().enumerate() {
            for (a, x) in p.iter_mut().enumerate() {
                *x = f32_at(rec, 12 + 12 * k + 4 * a);
                if !x.is_finite() {
                    return Err(Error::format(format!(
                        "triangle {f} has a non-finite coordinate"
                    )));
                }
            }
        }
        w.polygon(&c);
        attr.push(u16::from_le_bytes([rec[48], rec[49]]) as i64);
    }
    let mut mesh = Mesh::new(w.verts, w.tris)?;
    if attr.iter().any(|&a| a != 0) {
        mesh = mesh.with_face_column("stl_attribute", Column::i64(attr))?;
    }
    Ok(Data {
        mesh: Some(mesh),
        warnings,
        ..Data::default()
    })
}

struct Cursor<'a> {
    d: &'a [u8],
    pos: usize,
    line: u64,
}

impl<'a> Cursor<'a> {
    fn skip_ws(&mut self) {
        while let Some(&b) = self.d.get(self.pos) {
            if !b.is_ascii_whitespace() {
                break;
            }
            if b == b'\n' {
                self.line += 1;
            }
            self.pos += 1;
        }
    }

    fn token(&mut self) -> Option<&'a [u8]> {
        self.skip_ws();
        let start = self.pos;
        while self.pos < self.d.len() && !self.d[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        (self.pos > start).then(|| &self.d[start..self.pos])
    }

    fn rest_of_line(&mut self) -> &'a [u8] {
        let start = self.pos;
        while self.pos < self.d.len() && self.d[self.pos] != b'\n' {
            self.pos += 1;
        }
        trim(&self.d[start..self.pos])
    }

    fn expect(&mut self, kw: &str) -> Result<()> {
        match self.token() {
            Some(t) if t.eq_ignore_ascii_case(kw.as_bytes()) => Ok(()),
            t => Err(Error::format_at(
                self.line,
                format!(
                    "expected {kw:?}, found {:?}",
                    String::from_utf8_lossy(t.unwrap_or(b"end of file"))
                ),
            )),
        }
    }

    fn number(&mut self) -> Result<f64> {
        let line = self.line;
        let t = self
            .token()
            .ok_or_else(|| Error::format_at(line, "file ends inside a facet"))?;
        parse_coord(t, self.line)
    }
}

fn read_ascii(r: &[u8]) -> Result<Data> {
    let mut c = Cursor {
        d: r,
        pos: 0,
        line: 1,
    };
    let mut w = Welder::default();
    let mut names: Vec<&[u8]> = Vec::new();
    let mut solid_of_face: Vec<u32> = Vec::new();
    let mut corners: Vec<[f64; 3]> = Vec::new();
    let mut current = None;
    let mut short = 0usize;
    let mut open = None;
    while let Some(t) = c.token() {
        if t.eq_ignore_ascii_case(b"solid") {
            if let Some(at) = open {
                return Err(Error::format_at(
                    c.line,
                    format!("solid opened on line {at} has no endsolid before this solid"),
                ));
            }
            open = Some(c.line);
            current = Some(names.len() as u32);
            names.push(c.rest_of_line());
        } else if t.eq_ignore_ascii_case(b"endsolid") {
            c.rest_of_line();
            current = None;
            open = None;
        } else if t.eq_ignore_ascii_case(b"facet") {
            c.expect("normal")?;
            for _ in 0..3 {
                c.number().or_else(|e| match e {
                    Error::Format { .. } => Ok(0.0),
                    e => Err(e),
                })?;
            }
            c.expect("outer")?;
            c.expect("loop")?;
            corners.clear();
            loop {
                let line = c.line;
                match c.token() {
                    Some(t) if t.eq_ignore_ascii_case(b"vertex") => {
                        corners.push([c.number()?, c.number()?, c.number()?]);
                    }
                    Some(t) if t.eq_ignore_ascii_case(b"endloop") => break,
                    t => {
                        return Err(Error::format_at(
                            line,
                            format!(
                                "expected \"vertex\" or \"endloop\", found {:?}",
                                as_str(t.unwrap_or(b"end of file"))
                            ),
                        ));
                    }
                }
            }
            c.expect("endfacet")?;
            if corners.len() < 3 {
                short += 1;
                continue;
            }
            let before = w.tris.len();
            w.polygon(&corners);
            let id = *current.get_or_insert_with(|| {
                names.push(b"");
                names.len() as u32 - 1
            });
            solid_of_face.extend(std::iter::repeat_n(id, w.tris.len() - before));
        } else {
            return Err(Error::format_at(
                c.line,
                format!(
                    "expected \"solid\", \"facet\" or \"endsolid\", found {:?}",
                    String::from_utf8_lossy(t)
                ),
            ));
        }
    }
    if let Some(at) = open {
        return Err(Error::format_at(
            c.line,
            format!("solid opened on line {at} ends without endsolid; the file is truncated"),
        ));
    }
    let mut data = Data::default();
    if short > 0 {
        data.warnings
            .push(format!("skipped {short} facets with fewer than 3 vertices"));
    }
    let mut mesh = Mesh::new(w.verts, w.tris)?;
    let first = solid_of_face.first().map(|&i| names[i as usize]);
    let one_default = matches!(first, Some(b"" | DEFAULT_SOLID))
        && solid_of_face
            .iter()
            .all(|&i| Some(names[i as usize]) == first);
    if !solid_of_face.is_empty() && !one_default {
        let decoded: Vec<Option<String>> = names.iter().map(|n| decode_name(n)).collect();
        let valid: Vec<bool> = solid_of_face
            .iter()
            .map(|&i| decoded[i as usize].is_some())
            .collect();
        let col = solid_of_face
            .iter()
            .map(|&i| decoded[i as usize].clone().unwrap_or_default())
            .collect();
        mesh = mesh.with_face_column("solid", Column::utf8(col).with_valid(valid)?)?;
    }
    data.mesh = Some(mesh);
    Ok(data)
}

fn normal(c: [[f64; 3]; 3]) -> [f64; 3] {
    let u = [c[1][0] - c[0][0], c[1][1] - c[0][1], c[1][2] - c[0][2]];
    let v = [c[2][0] - c[0][0], c[2][1] - c[0][1], c[2][2] - c[0][2]];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len > 0.0 && len.is_finite() {
        n.map(|x| x / len)
    } else {
        [0.0; 3]
    }
}

/// The name a mesh without a `solid` column is written under; a column value
/// equal to it is written escaped so the reader keeps the column.
const DEFAULT_SOLID: &[u8] = b"wireframes";

fn solid_name(name: Option<&str>) -> String {
    let mut out = Vec::new();
    match name {
        Some(n) if n.as_bytes() == DEFAULT_SOLID => out.extend_from_slice(b"%77ireframes"),
        _ => encode_name(name, false, &mut out),
    }
    String::from_utf8_lossy(&out).into_owned()
}

pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.polylines.is_some() || data.points.is_some() || data.grid.is_some() {
        return Err(Error::Unsupported(
            "stl stores only a mesh, not polylines, points or a grid".into(),
        ));
    }
    let empty = Mesh::empty();
    let m = data.mesh.as_ref().unwrap_or(&empty);
    check_no_columns("stl", "vertex", m.vertex_attributes(), opts)?;
    let fa = m.face_attributes();
    let keep = if opts.ascii { "solid" } else { "stl_attribute" };
    for name in fa.names() {
        if name != keep {
            unsupported_column("stl", "face", name, opts)?;
        }
    }
    if opts.ascii {
        write_ascii(m, utf8_values(fa, "solid", "stl")?)
    } else {
        write_binary(m, fa.column("stl_attribute"))
    }
}

fn write_ascii(m: &Mesh, solids: Option<(&[String], &Column)>) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(m.n_faces() * 260 + 32);
    let name_of = |f: usize| -> String {
        match solids {
            Some((v, c)) => solid_name(c.is_valid(f).then(|| v[f].as_str())),
            None => "wireframes".into(),
        }
    };
    let mut current: Option<String> = None;
    for f in 0..m.n_faces() {
        let name = name_of(f);
        if current.as_ref() != Some(&name) {
            if let Some(prev) = &current {
                writeln!(out, "endsolid {prev}")?;
            }
            writeln!(out, "solid {name}")?;
            current = Some(name);
        }
        let c = m.corners(f as u32);
        out.extend_from_slice(b"  facet normal");
        for x in normal(c) {
            out.push(b' ');
            push_f64(&mut out, x)?;
        }
        out.extend_from_slice(b"\n    outer loop\n");
        for p in c {
            out.extend_from_slice(b"      vertex");
            for x in p {
                out.push(b' ');
                push_f64(&mut out, x)?;
            }
            out.push(b'\n');
        }
        out.extend_from_slice(b"    endloop\n  endfacet\n");
    }
    match current {
        Some(prev) => writeln!(out, "endsolid {prev}")?,
        None => out.extend_from_slice(b"solid wireframes\nendsolid wireframes\n"),
    }
    Ok(out)
}

/// How many vertices binary STL moves by storing them as f32, and the largest
/// coordinate change.
pub fn f32_rounding(m: &Mesh) -> (usize, f64) {
    let mut moved = 0;
    let mut worst = 0.0f64;
    for p in m.vertices() {
        let d = p
            .iter()
            .map(|&x| (x - x as f32 as f64).abs())
            .fold(0.0, f64::max);
        if d != 0.0 {
            moved += 1;
            worst = worst.max(d);
        }
    }
    (moved, worst)
}

fn to_f32(x: f64) -> Result<f32> {
    let y = x as f32;
    if y.is_finite() {
        Ok(y)
    } else {
        Err(Error::invalid(format!(
            "coordinate {x} does not fit in the f32 range of binary stl"
        )))
    }
}

fn write_binary(m: &Mesh, attr: Option<&Column>) -> Result<Vec<u8>> {
    let attr = match attr {
        None => None,
        Some(c) => match &c.data {
            ColumnData::Int64(v) => Some((v, c)),
            _ => {
                return Err(Error::invalid(format!(
                    "face column \"stl_attribute\" must be int64, found {}",
                    c.type_name()
                )));
            }
        },
    };
    let n = u32::try_from(m.n_faces()).map_err(|_| {
        Error::invalid(format!(
            "binary stl holds at most {} triangles, mesh has {}",
            u32::MAX,
            m.n_faces()
        ))
    })?;
    let mut out = Vec::with_capacity(84 + 50 * m.n_faces());
    let mut header = [0u8; 80];
    header[..10].copy_from_slice(b"wireframes");
    out.extend_from_slice(&header);
    out.extend_from_slice(&n.to_le_bytes());
    for f in 0..m.n_faces() {
        let c = m.corners(f as u32);
        for x in normal(c) {
            out.extend_from_slice(&(x as f32).to_le_bytes());
        }
        for p in c {
            for x in p {
                out.extend_from_slice(&to_f32(x)?.to_le_bytes());
            }
        }
        let a = match attr {
            Some((v, c)) if c.is_valid(f) => u16::try_from(v[f]).map_err(|_| {
                Error::invalid(format!(
                    "stl_attribute {} on face {f} does not fit in u16",
                    v[f]
                ))
            })?,
            _ => 0,
        };
        out.extend_from_slice(&a.to_le_bytes());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f32_rounding_counts_the_vertices_binary_stl_moves() {
        let m = Mesh::new(
            vec![
                [0.5, 1.0, 2.0],
                [700_000.1, 7_500_000.0, 0.0],
                [0.0, 7_500_000.3, 0.0],
            ],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let (moved, worst) = f32_rounding(&m);
        assert_eq!(moved, 2);
        let expect = (7_500_000.3f64 - 7_500_000.3f32 as f64).abs();
        assert_eq!(
            worst,
            expect.max((700_000.1f64 - 700_000.1f32 as f64).abs())
        );
        let exact = Mesh::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.25, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        assert_eq!(f32_rounding(&exact), (0, 0.0));
    }

    const FACET: &str = "facet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\n";

    #[test]
    fn an_ascii_solid_without_endsolid_is_truncated() {
        let whole = format!("solid a\n{FACET}{FACET}endsolid a\n");
        assert!(read(whole.as_bytes()).is_ok());
        let cut = format!("solid a\n{FACET}");
        let e = read(cut.as_bytes()).unwrap_err().to_string();
        assert!(e.contains("line 1 ends without endsolid"), "{e}");
        let unclosed = format!("solid a\n{FACET}solid b\n{FACET}endsolid b\n");
        assert!(read(unclosed.as_bytes()).is_err());
    }

    #[test]
    fn a_binary_file_with_a_solid_header_and_padding_reads_as_binary() {
        let mut b = b"solid exported by package".to_vec();
        b.resize(80, b' ');
        b.extend(1u32.to_le_bytes());
        for x in [0f32, 0., 1., 0., 0., 0., 1., 0., 0., 0., 1., 0.] {
            b.extend(x.to_le_bytes());
        }
        b.extend([0u8; 2]);
        b.extend([0u8; 7]);
        let d = read(&b).unwrap();
        assert_eq!(d.mesh.unwrap().n_faces(), 1);
        assert_eq!(
            d.warnings,
            ["binary stl declares 1 triangles; ignored 7 trailing bytes"]
        );
    }
}
