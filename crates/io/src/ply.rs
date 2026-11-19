use std::io::Write;

use wireframes_core::{Column, ColumnData, Error, Mesh, Result, Table};

use crate::common::{polygon_corners, polygon_groups};
use crate::obj::{Polygons, as_str, push_f64, tokens, unsupported_column};
use crate::{Data, WriteOptions};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    F32,
    F64,
}

impl Ty {
    fn parse(t: &[u8], line: u64) -> Result<Ty> {
        Ok(match t {
            b"char" | b"int8" => Ty::I8,
            b"uchar" | b"uint8" => Ty::U8,
            b"short" | b"int16" => Ty::I16,
            b"ushort" | b"uint16" => Ty::U16,
            b"int" | b"int32" => Ty::I32,
            b"uint" | b"uint32" => Ty::U32,
            b"float" | b"float32" => Ty::F32,
            b"double" | b"float64" => Ty::F64,
            b"int64" | b"uint64" | b"long" | b"ulong" => {
                return Err(Error::Unsupported(format!(
                    "ply property type {} is not supported",
                    as_str(t)
                )));
            }
            _ => {
                return Err(Error::format_at(
                    line,
                    format!("unknown ply type {:?}", String::from_utf8_lossy(t)),
                ));
            }
        })
    }

    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::F64 => 8,
        }
    }

    fn is_float(self) -> bool {
        matches!(self, Ty::F32 | Ty::F64)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Scalar(Ty),
    List(Ty, Ty),
}

struct Prop {
    name: String,
    kind: Kind,
}

struct Element {
    name: String,
    count: usize,
    props: Vec<Prop>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Enc {
    Ascii,
    Le,
    Be,
}

#[derive(Clone, Copy)]
enum Num {
    I(i64),
    F(f64),
}

impl Num {
    fn f64(self) -> f64 {
        match self {
            Num::I(i) => i as f64,
            Num::F(x) => x,
        }
    }
}

struct Body<'a> {
    d: &'a [u8],
    pos: usize,
    line: u64,
    enc: Enc,
}

impl<'a> Body<'a> {
    fn err(&self, msg: impl Into<String>) -> Error {
        match self.enc {
            Enc::Ascii => Error::format_at(self.line, msg),
            _ => Error::format(format!("{} (at byte {})", msg.into(), self.pos)),
        }
    }

    fn token(&mut self) -> Result<&'a [u8]> {
        while let Some(&b) = self.d.get(self.pos) {
            if !b.is_ascii_whitespace() {
                break;
            }
            if b == b'\n' {
                self.line += 1;
            }
            self.pos += 1;
        }
        let start = self.pos;
        while self.pos < self.d.len() && !self.d[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        if self.pos == start {
            return Err(self.err("unexpected end of data"));
        }
        Ok(&self.d[start..self.pos])
    }

    fn bytes<const N: usize>(&mut self) -> Result<[u8; N]> {
        let b = self
            .d
            .get(self.pos..self.pos + N)
            .ok_or_else(|| self.err("unexpected end of data"))?;
        self.pos += N;
        Ok(b.try_into().expect("sized slice"))
    }

    fn num(&mut self, t: Ty) -> Result<Num> {
        macro_rules! bin {
            ($ty:ty, $n:literal) => {{
                let b = self.bytes::<$n>()?;
                if self.enc == Enc::Le {
                    <$ty>::from_le_bytes(b)
                } else {
                    <$ty>::from_be_bytes(b)
                }
            }};
        }
        if self.enc == Enc::Ascii {
            let tok = self.token()?;
            let s = as_str(tok);
            if t.is_float() {
                return s.parse::<f64>().map(Num::F).map_err(|_| {
                    self.err(format!(
                        "expected a number, found {:?}",
                        String::from_utf8_lossy(tok)
                    ))
                });
            }
            if let Ok(i) = s.parse::<i64>() {
                return Ok(Num::I(i));
            }
            return match s.parse::<f64>() {
                Ok(x) if x.fract() == 0.0 && x.abs() < 9.0e15 => Ok(Num::I(x as i64)),
                _ => Err(self.err(format!(
                    "expected an integer, found {:?}",
                    String::from_utf8_lossy(tok)
                ))),
            };
        }
        Ok(match t {
            Ty::I8 => Num::I(self.bytes::<1>()?[0] as i8 as i64),
            Ty::U8 => Num::I(self.bytes::<1>()?[0] as i64),
            Ty::I16 => Num::I(bin!(i16, 2) as i64),
            Ty::U16 => Num::I(bin!(u16, 2) as i64),
            Ty::I32 => Num::I(bin!(i32, 4) as i64),
            Ty::U32 => Num::I(bin!(u32, 4) as i64),
            Ty::F32 => Num::F(bin!(f32, 4) as f64),
            Ty::F64 => Num::F(bin!(f64, 8)),
        })
    }

    fn int(&mut self, t: Ty) -> Result<i64> {
        match self.num(t)? {
            Num::I(i) => Ok(i),
            Num::F(_) => unreachable!("header admits only integer types here"),
        }
    }

    fn count(&mut self, t: Ty, item: Ty) -> Result<usize> {
        let n = self.int(t)?;
        let n = usize::try_from(n).map_err(|_| self.err(format!("list length {n} is negative")))?;
        if self.enc != Enc::Ascii && n.saturating_mul(item.size()) > self.d.len() - self.pos {
            return Err(self.err(format!("list of {n} items runs past the end of data")));
        }
        Ok(n)
    }

    fn skip(&mut self, kind: Kind) -> Result<()> {
        match kind {
            Kind::Scalar(t) => {
                self.num(t)?;
            }
            Kind::List(c, t) => {
                let n = self.count(c, t)?;
                if self.enc == Enc::Ascii {
                    for _ in 0..n {
                        self.token()?;
                    }
                } else {
                    self.pos += n * t.size();
                }
            }
        }
        Ok(())
    }
}

struct Header {
    enc: Enc,
    elements: Vec<Element>,
    body: usize,
    lines: u64,
}

fn header(r: &[u8]) -> Result<Header> {
    let mut pos = 0;
    let mut no = 0u64;
    let mut enc = None;
    let mut elements: Vec<Element> = Vec::new();
    loop {
        let Some(len) = r[pos..].iter().position(|&b| b == b'\n') else {
            return Err(Error::format("ply header has no end_header line"));
        };
        let line = &r[pos..pos + len];
        pos += len + 1;
        no += 1;
        let mut t = tokens(line);
        let kw = t.next().unwrap_or_default();
        if no == 1 {
            if kw != b"ply" {
                return Err(Error::format_at(
                    1,
                    "not a ply file: first line must be \"ply\"",
                ));
            }
            continue;
        }
        match kw {
            b"format" => {
                enc = Some(match t.next().unwrap_or_default() {
                    b"ascii" => Enc::Ascii,
                    b"binary_little_endian" => Enc::Le,
                    b"binary_big_endian" => Enc::Be,
                    f => {
                        return Err(Error::format_at(
                            no,
                            format!("unknown ply format {:?}", String::from_utf8_lossy(f)),
                        ));
                    }
                });
            }
            b"element" => {
                let name = t
                    .next()
                    .ok_or_else(|| Error::format_at(no, "element needs a name"))?;
                let count = t
                    .next()
                    .and_then(|c| as_str(c).parse::<usize>().ok())
                    .ok_or_else(|| Error::format_at(no, "element needs a non-negative count"))?;
                elements.push(Element {
                    name: String::from_utf8_lossy(name).into_owned(),
                    count,
                    props: Vec::new(),
                });
            }
            b"property" => {
                let el = elements
                    .last_mut()
                    .ok_or_else(|| Error::format_at(no, "property before any element"))?;
                let first = t
                    .next()
                    .ok_or_else(|| Error::format_at(no, "property needs a type"))?;
                let kind = if first == b"list" {
                    let c = Ty::parse(t.next().unwrap_or_default(), no)?;
                    let i = Ty::parse(t.next().unwrap_or_default(), no)?;
                    if c.is_float() {
                        return Err(Error::format_at(no, "list count type must be an integer"));
                    }
                    Kind::List(c, i)
                } else {
                    Kind::Scalar(Ty::parse(first, no)?)
                };
                let name = t
                    .next()
                    .ok_or_else(|| Error::format_at(no, "property needs a name"))?;
                el.props.push(Prop {
                    name: String::from_utf8_lossy(name).into_owned(),
                    kind,
                });
            }
            b"end_header" => break,
            b"comment" | b"obj_info" | b"" => {}
            _ => {
                return Err(Error::format_at(
                    no,
                    format!("unknown ply header line {:?}", String::from_utf8_lossy(kw)),
                ));
            }
        }
    }
    let enc = enc.ok_or_else(|| Error::format("ply header has no format line"))?;
    Ok(Header {
        enc,
        elements,
        body: pos,
        lines: no,
    })
}

enum Buf {
    F(Vec<f64>),
    I(Vec<i64>),
}

impl Buf {
    fn new(t: Ty) -> Buf {
        if t.is_float() {
            Buf::F(Vec::new())
        } else {
            Buf::I(Vec::new())
        }
    }

    fn push(&mut self, v: Num, times: usize) {
        match (self, v) {
            (Buf::F(b), v) => b.extend(std::iter::repeat_n(v.f64(), times)),
            (Buf::I(b), Num::I(i)) => b.extend(std::iter::repeat_n(i, times)),
            (Buf::I(_), Num::F(_)) => unreachable!("integer property decoded as float"),
        }
    }

    fn column(self) -> Result<Column> {
        match self {
            Buf::I(v) => Ok(Column::i64(v)),
            Buf::F(v) => {
                let valid = v.iter().map(|x| !x.is_nan()).collect();
                Column::f64(v).with_valid(valid)
            }
        }
    }
}

#[derive(Clone, Copy)]
enum Target {
    Coord(usize),
    Attr(usize),
    Indices,
    Skip,
}

fn table(rows: usize, bufs: Vec<(String, Buf)>) -> Result<Table> {
    let mut t = Table::new(rows);
    for (name, b) in bufs {
        t = t.with_column(&name, b.column()?)?;
    }
    Ok(t)
}

fn index_ty(p: &Prop) -> Option<Ty> {
    match p.kind {
        Kind::List(_, t) if !t.is_float() => Some(t),
        _ => None,
    }
}

pub fn read(r: &[u8]) -> Result<Data> {
    let h = header(r)?;
    let mut b = Body {
        d: r,
        pos: h.body,
        line: h.lines + 1,
        enc: h.enc,
    };
    let mut data = Data::default();
    let mut verts: Option<Vec<[f64; 3]>> = None;
    let mut vbufs: Vec<(String, Buf)> = Vec::new();
    let mut fbufs: Vec<(String, Buf)> = Vec::new();
    let mut polys = Polygons::default();
    let mut face_rows: Vec<Option<u32>> = Vec::new();
    let mut n_face_tris = 0usize;
    let mut saw_strips = false;
    let mut saw_face = false;
    let mut corners: Vec<u32> = Vec::new();
    let mut rec: Vec<Num> = Vec::new();
    let mut skipped_lists: Vec<String> = Vec::new();

    for el in &h.elements {
        let is = |n: &str| el.name == n;
        let role = if is("vertex") {
            0
        } else if is("face") {
            1
        } else if is("tristrips") {
            2
        } else {
            3
        };
        if (role == 0 && verts.is_some()) || (role == 1 && saw_face) {
            return Err(Error::format(format!(
                "duplicate ply element {:?}",
                el.name
            )));
        }
        let mut targets = Vec::with_capacity(el.props.len());
        let mut have_indices = false;
        let mut coords = [false; 3];
        for p in &el.props {
            let t = match (role, p.kind) {
                (0, Kind::Scalar(_)) if matches!(p.name.as_str(), "x" | "y" | "z") => {
                    let a = (p.name.as_bytes()[0] - b'x') as usize;
                    coords[a] = true;
                    Target::Coord(a)
                }
                (0, Kind::Scalar(t)) => {
                    vbufs.push((p.name.clone(), Buf::new(t)));
                    Target::Attr(vbufs.len() - 1)
                }
                (1, Kind::Scalar(t)) => {
                    fbufs.push((p.name.clone(), Buf::new(t)));
                    Target::Attr(fbufs.len() - 1)
                }
                (1 | 2, Kind::List(..))
                    if !have_indices
                        && matches!(p.name.as_str(), "vertex_indices" | "vertex_index") =>
                {
                    if index_ty(p).is_none() {
                        return Err(Error::format(format!(
                            "ply {} {} must hold integers",
                            el.name, p.name
                        )));
                    }
                    have_indices = true;
                    Target::Indices
                }
                (0 | 1, Kind::List(..)) => {
                    skipped_lists.push(format!("{}.{}", el.name, p.name));
                    Target::Skip
                }
                _ => Target::Skip,
            };
            targets.push(t);
        }
        if role == 0
            && let Some(a) = coords.iter().position(|c| !c)
        {
            return Err(Error::format(format!(
                "ply vertex element lacks property {}",
                ["x", "y", "z"][a]
            )));
        }
        if (role == 1 || role == 2) && !have_indices {
            return Err(Error::format(format!(
                "ply {} element lacks a vertex_indices list",
                el.name
            )));
        }
        saw_face |= role == 1;
        if el.props.is_empty() {
            continue;
        }
        let cap = el.count.min(b.d.len());
        let mut vs = Vec::with_capacity(if role == 0 { cap } else { 0 });
        for _ in 0..el.count {
            let mut p = [0.0; 3];
            corners.clear();
            rec.clear();
            for (prop, t) in el.props.iter().zip(&targets) {
                match (*t, prop.kind) {
                    (Target::Coord(a), Kind::Scalar(ty)) => {
                        let at = b.pos;
                        let x = b.num(ty)?.f64();
                        if !x.is_finite() {
                            if b.enc != Enc::Ascii {
                                b.pos = at;
                            }
                            return Err(b.err(format!("vertex coordinate {x} is not finite")));
                        }
                        p[a] = x;
                    }
                    (Target::Attr(i), Kind::Scalar(ty)) => {
                        let v = b.num(ty)?;
                        if role == 0 {
                            vbufs[i].1.push(v, 1);
                        } else {
                            rec.push(v);
                        }
                    }
                    (Target::Indices, Kind::List(c, ty)) => {
                        let n = b.count(c, ty)?;
                        for _ in 0..n {
                            let i = b.int(ty)?;
                            if role == 2 && i == -1 {
                                corners.push(u32::MAX);
                            } else if i < 0 || i >= u32::MAX as i64 {
                                return Err(b.err(format!("vertex index {i} is out of range")));
                            } else {
                                corners.push(i as u32);
                            }
                        }
                    }
                    (_, kind) => b.skip(kind)?,
                }
            }
            match role {
                0 => vs.push(p),
                1 => {
                    let n = polys.push(&corners);
                    let start = n_face_tris as u32;
                    face_rows.extend((start..start + n as u32).map(Some));
                    n_face_tris += n;
                    let mut vals = rec.iter();
                    for (_, buf) in fbufs.iter_mut() {
                        buf.push(*vals.next().expect("one value per attribute"), n);
                    }
                }
                2 => {
                    saw_strips = true;
                    for strip in corners.split(|&i| i == u32::MAX) {
                        for k in 0..strip.len().saturating_sub(2) {
                            let t = if k % 2 == 0 {
                                [strip[k], strip[k + 1], strip[k + 2]]
                            } else {
                                [strip[k + 1], strip[k], strip[k + 2]]
                            };
                            if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
                                polys.push(&t);
                                face_rows.push(None);
                            }
                        }
                    }
                }
                _ => {}
            }
        }
        if role == 0 {
            verts = Some(vs);
        }
    }

    if !skipped_lists.is_empty() {
        data.warnings.push(format!(
            "skipped list properties {}",
            skipped_lists.join(", ")
        ));
    }
    polys.warn(&mut data.warnings);
    let Some(verts) = verts else {
        if !polys.tris.is_empty() {
            return Err(Error::format("ply has faces but no vertex element"));
        }
        return Ok(data);
    };
    let n = verts.len();
    if let Some(&i) = polys.tris.iter().flatten().find(|&&i| i as usize >= n) {
        return Err(Error::format(format!(
            "face index {i} is out of range for {n} vertices"
        )));
    }
    let (tris, polygon) = polys.finish();
    let vt = table(n, vbufs)?;
    let mut ft = table(n_face_tris, fbufs)?;
    if saw_strips {
        ft = ft.take_or_null(&face_rows);
    }
    let mut mesh = Mesh::new(verts, tris)?.with_attributes(vt, ft)?;
    match mesh.face_attributes().column("polygon") {
        Some(c) => {
            if let Some(ids) = integral(c) {
                mesh = mesh.with_face_column("polygon", ids)?;
            }
        }
        None => {
            if let Some(c) = polygon {
                mesh = mesh.with_face_column("polygon", c)?;
            }
        }
    }
    data.mesh = Some(mesh);
    Ok(data)
}

pub(crate) const F64_EXACT_INT: i64 = 1 << 53;

/// A float column whose valid values are all integers that `f64` holds exactly,
/// as Int64; the writer stores an Int64 `polygon` column with nulls as double.
fn integral(c: &Column) -> Option<Column> {
    let ColumnData::Float64(v) = &c.data else {
        return None;
    };
    let ok = (0..v.len())
        .all(|i| !c.is_valid(i) || (v[i].fract() == 0.0 && v[i].abs() <= F64_EXACT_INT as f64));
    if !ok {
        return None;
    }
    let ids = v
        .iter()
        .map(|&x| if x.is_nan() { 0 } else { x as i64 })
        .collect();
    let valid = (0..v.len()).map(|i| c.is_valid(i)).collect();
    Column::i64(ids).with_valid(valid).ok()
}

enum Out<'a> {
    Double(&'a Column),
    Int(&'a [i64]),
    UChar(&'a Column),
}

impl Out<'_> {
    fn ty(&self) -> &'static str {
        match self {
            Out::Double(_) => "double",
            Out::Int(_) => "int",
            Out::UChar(_) => "uchar",
        }
    }
}

fn plan<'a>(
    t: &'a Table,
    reserved: &[&str],
    kind: &str,
    opts: &WriteOptions,
) -> Result<Vec<(&'a str, Out<'a>)>> {
    let mut out = Vec::new();
    for (name, c) in t.names().iter().zip(t.columns()) {
        let bad_name = name.is_empty()
            || name
                .bytes()
                .any(|b| b.is_ascii_whitespace() || !b.is_ascii())
            || reserved.contains(&name.as_str());
        if bad_name {
            unsupported_column("ply", kind, name, opts)?;
            continue;
        }
        let o = match &c.data {
            ColumnData::Float64(_) => Out::Double(c),
            ColumnData::Int64(v) => {
                if c.valid.is_none() && v.iter().all(|&x| i32::try_from(x).is_ok()) {
                    Out::Int(v)
                } else if v.iter().all(|&x| x.unsigned_abs() <= F64_EXACT_INT as u64) {
                    Out::Double(c)
                } else {
                    if !opts.skip_unsupported_columns {
                        return Err(Error::Unsupported(format!(
                            "ply cannot store {kind} column {name:?}: it has nulls or values outside int32, and values beyond 2^53 do not fit a double exactly; set skip_unsupported_columns to drop it"
                        )));
                    }
                    continue;
                }
            }
            ColumnData::Bool(_) => Out::UChar(c),
            ColumnData::Utf8(_) => {
                unsupported_column("ply", kind, name, opts)?;
                continue;
            }
        };
        out.push((name.as_str(), o));
    }
    Ok(out)
}

fn put_value(out: &mut Vec<u8>, o: &Out, row: usize, ascii: bool) -> Result<()> {
    match o {
        Out::Double(c) => {
            let x = if !c.is_valid(row) {
                f64::NAN
            } else {
                match &c.data {
                    ColumnData::Float64(v) => v[row],
                    ColumnData::Int64(v) => v[row] as f64,
                    _ => unreachable!("planned as double"),
                }
            };
            if ascii {
                push_f64(out, x)?;
            } else {
                out.extend_from_slice(&x.to_le_bytes());
            }
        }
        Out::Int(v) => {
            let x = v[row] as i32;
            if ascii {
                write!(out, "{x}")?;
            } else {
                out.extend_from_slice(&x.to_le_bytes());
            }
        }
        Out::UChar(c) => {
            let x = match &c.data {
                ColumnData::Bool(v) => c.is_valid(row) && v[row],
                _ => unreachable!("planned as uchar"),
            } as u8;
            if ascii {
                write!(out, "{x}")?;
            } else {
                out.push(x);
            }
        }
    }
    Ok(())
}

/// The reader numbers polygons 0, 1, ... and records them only when some polygon
/// has more than 3 corners; a `polygon` column of exactly that form is left out
/// of the file, any other is written as a face property the reader keeps.
pub(crate) fn regenerated_by_reader(m: &Mesh, groups: &[std::ops::Range<usize>]) -> bool {
    let Some(c) = m.face_attributes().column("polygon") else {
        return false;
    };
    let ColumnData::Int64(v) = &c.data else {
        return false;
    };
    c.valid.is_none()
        && m.face_attributes()
            .names()
            .last()
            .is_some_and(|n| n == "polygon")
        && groups.iter().any(|g| g.len() > 1)
        && groups
            .iter()
            .enumerate()
            .all(|(k, g)| v[g.clone()].iter().all(|&id| id == k as i64))
}

pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.polylines.is_some() || data.points.is_some() || data.grid.is_some() {
        return Err(Error::Unsupported(
            "ply stores only a mesh, not polylines, points or a grid".into(),
        ));
    }
    let empty = Mesh::empty();
    let m = data.mesh.as_ref().unwrap_or(&empty);
    if m.n_vertices() > i32::MAX as usize {
        return Err(Error::Unsupported(format!(
            "ply indices are int; {} vertices exceed {}",
            m.n_vertices(),
            i32::MAX
        )));
    }
    let vcols = plan(m.vertex_attributes(), &["x", "y", "z"], "vertex", opts)?;
    let groups = polygon_groups(m)?;
    let mut fcols = plan(
        m.face_attributes(),
        &["vertex_indices", "vertex_index"],
        "face",
        opts,
    )?;
    if regenerated_by_reader(m, &groups) {
        fcols.retain(|(name, _)| *name != "polygon");
    }
    let tris = m.triangles();
    let wide = groups.iter().any(|g| g.len() + 2 > 255);
    let ascii = opts.ascii;

    let mut out =
        Vec::with_capacity(256 + m.n_vertices() * (24 + 8 * vcols.len()) + groups.len() * 16);
    writeln!(
        out,
        "ply\nformat {} 1.0\ncomment wireframes\nelement vertex {}",
        if ascii {
            "ascii"
        } else {
            "binary_little_endian"
        },
        m.n_vertices()
    )?;
    out.extend_from_slice(b"property double x\nproperty double y\nproperty double z\n");
    for (name, o) in &vcols {
        writeln!(out, "property {} {name}", o.ty())?;
    }
    writeln!(
        out,
        "element face {}\nproperty list {} int vertex_indices",
        groups.len(),
        if wide { "int" } else { "uchar" }
    )?;
    for (name, o) in &fcols {
        writeln!(out, "property {} {name}", o.ty())?;
    }
    out.extend_from_slice(b"end_header\n");

    for (i, p) in m.vertices().iter().enumerate() {
        for (k, x) in p.iter().enumerate() {
            if ascii {
                if k > 0 {
                    out.push(b' ');
                }
                push_f64(&mut out, *x)?;
            } else {
                out.extend_from_slice(&x.to_le_bytes());
            }
        }
        for (_, o) in &vcols {
            if ascii {
                out.push(b' ');
            }
            put_value(&mut out, o, i, ascii)?;
        }
        if ascii {
            out.push(b'\n');
        }
    }
    let mut corners = Vec::new();
    for g in &groups {
        polygon_corners(tris, g, &mut corners);
        let n = corners.len();
        if ascii {
            write!(out, "{n}")?;
            for &i in &corners {
                write!(out, " {i}")?;
            }
        } else {
            if wide {
                out.extend_from_slice(&(n as i32).to_le_bytes());
            } else {
                out.push(n as u8);
            }
            for &i in &corners {
                out.extend_from_slice(&(i as i32).to_le_bytes());
            }
        }
        for (_, o) in &fcols {
            if ascii {
                out.push(b' ');
            }
            put_value(&mut out, o, g.start, ascii)?;
        }
        if ascii {
            out.push(b'\n');
        }
    }
    Ok(out)
}
