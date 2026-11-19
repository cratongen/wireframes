//! AutoCAD DXF: group-code reader for ASCII and binary files, hand-written
//! ASCII R2000 writer.

use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;

use wireframes_core::{Column, ColumnData, Error, PointSet, Polylines, Result, Table};

use crate::common::polygons;
use crate::vtk::{Fan, fmt_f64};
use crate::{Data, DxfEntity, WriteOptions};

const CP1252_HIGH: [char; 32] = [
    '\u{20AC}', '\u{81}', '\u{201A}', '\u{192}', '\u{201E}', '\u{2026}', '\u{2020}', '\u{2021}',
    '\u{2C6}', '\u{2030}', '\u{160}', '\u{2039}', '\u{152}', '\u{8D}', '\u{17D}', '\u{8F}',
    '\u{90}', '\u{2018}', '\u{2019}', '\u{201C}', '\u{201D}', '\u{2022}', '\u{2013}', '\u{2014}',
    '\u{2DC}', '\u{2122}', '\u{161}', '\u{203A}', '\u{153}', '\u{9D}', '\u{17E}', '\u{178}',
];

/// UTF-8 when valid, else Windows-1252.
pub(crate) fn decode_text(b: &[u8]) -> String {
    match std::str::from_utf8(b) {
        Ok(s) => s.to_string(),
        Err(_) => b
            .iter()
            .map(|&c| match c {
                0x80..=0x9F => CP1252_HIGH[(c - 0x80) as usize],
                _ => c as char,
            })
            .collect(),
    }
}

/// AC1021 (AutoCAD 2007) and later store text as UTF-8; earlier versions use a code page.
fn utf8_only(b: &[u8]) -> bool {
    let Some(i) = b.windows(8).position(|w| w == b"$ACADVER") else {
        return false;
    };
    let rest = &b[i + 8..b.len().min(i + 48)];
    rest.windows(6)
        .find(|w| w.starts_with(b"AC") && w[2..].iter().all(u8::is_ascii_digit))
        .is_some_and(|v| v >= &b"AC1021"[..])
}

fn dxf_text(b: &[u8], utf8_only: bool) -> Result<String> {
    if !utf8_only {
        return Ok(decode_text(b));
    }
    String::from_utf8(b.to_vec())
        .map_err(|_| Error::format("dxf from AC1021 on must be UTF-8, but the text is not"))
}

/// Replaces AutoCAD `\U+XXXX` escapes in one left-to-right pass, so a decoded
/// backslash never starts another escape.
fn unescape_unicode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('\\') {
        out.push_str(&rest[..i]);
        let tail = &rest[i..];
        let c = tail
            .get(1..7)
            .filter(|t| {
                (t.starts_with("U+") || t.starts_with("u+"))
                    && t[2..].bytes().all(|b| b.is_ascii_hexdigit())
            })
            .and_then(|t| u32::from_str_radix(&t[2..], 16).ok())
            .and_then(char::from_u32);
        match c {
            Some(c) => {
                out.push(c);
                rest = &tail[7..];
            }
            None => {
                out.push('\\');
                rest = &tail[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// The inverse of `unescape_unicode` for a layer name: backslashes, edge spaces
/// and anything outside printable ASCII become `\U+XXXX`, so the name survives
/// the reader's trimming and an R2000 code page.
fn escape_layer(l: &str, what: &str) -> Result<String> {
    if l.is_empty() || l.contains(['\n', '\r']) {
        return Err(Error::invalid(format!(
            "{what} layer {l:?} is not a valid dxf layer name"
        )));
    }
    let n = l.chars().count();
    let mut out = String::with_capacity(l.len());
    for (i, c) in l.chars().enumerate() {
        let edge = c == ' ' && (i == 0 || i + 1 == n);
        if c == '\\' || edge || !(c.is_ascii_graphic() || c == ' ') {
            if c as u32 > 0xFFFF {
                return Err(Error::invalid(format!(
                    "{what} layer {l:?} has {c:?}, outside what dxf \\U+XXXX escapes can hold"
                )));
            }
            let _ = write!(out, "\\U+{:04X}", c as u32);
        } else {
            out.push(c);
        }
    }
    Ok(out)
}

struct Group {
    code: i32,
    value: String,
    line: Option<u64>,
}

fn group_err(g: &Group, m: String) -> Error {
    match g.line {
        Some(l) => Error::format_at(l, m),
        None => Error::format(m),
    }
}

impl Group {
    fn f64(&self) -> Result<f64> {
        match self.value.trim().parse::<f64>() {
            Ok(x) if x.is_finite() => Ok(x),
            Ok(_) => Err(group_err(
                self,
                format!("group {} value {:?} is not finite", self.code, self.value),
            )),
            Err(_) => Err(group_err(
                self,
                format!("group {} value {:?} is not a number", self.code, self.value),
            )),
        }
    }

    fn i64(&self) -> Result<i64> {
        let v = self.value.trim();
        v.parse::<i64>()
            .or_else(|_| match v.parse::<f64>() {
                Ok(x) if x.fract() == 0.0 && x.abs() < 9.0e15 => Ok(x as i64),
                _ => Err(()),
            })
            .map_err(|_| {
                group_err(
                    self,
                    format!(
                        "group {} value {:?} is not an integer",
                        self.code, self.value
                    ),
                )
            })
    }
}

fn ascii_groups(b: &[u8]) -> Result<Vec<Group>> {
    let text = dxf_text(b, utf8_only(b))?;
    let lines: Vec<&str> = text.lines().collect();
    let mut out = Vec::with_capacity(lines.len() / 2);
    let mut i = 0;
    while i < lines.len() {
        let c = lines[i].trim();
        if c.is_empty() && lines[i..].iter().all(|l| l.trim().is_empty()) {
            break;
        }
        let code = c.parse::<i32>().map_err(|_| {
            Error::format_at(i as u64 + 1, format!("expected a group code, found {c:?}"))
        })?;
        let value = lines
            .get(i + 1)
            .ok_or_else(|| Error::format_at(i as u64 + 1, "group code without a value"))?;
        out.push(Group {
            code,
            value: unescape_unicode(value.trim()),
            line: Some(i as u64 + 2),
        });
        if code == 0 && value.trim() == "EOF" {
            break;
        }
        i += 2;
    }
    Ok(out)
}

const BINARY_SENTINEL: &[u8] = b"AutoCAD Binary DXF\r\n\x1a\0";

fn binary_groups(b: &[u8]) -> Result<Vec<Group>> {
    let short = || Error::format("binary dxf ends inside a group");
    let two_byte = b.get(1) == Some(&0);
    let utf8 = utf8_only(b);
    let mut pos = 0;
    let mut out = Vec::new();
    let take = |pos: &mut usize, n: usize| -> Result<&[u8]> {
        let s = b
            .get(*pos..pos.checked_add(n).ok_or_else(short)?)
            .ok_or_else(short)?;
        *pos += n;
        Ok(s)
    };
    while pos < b.len() {
        let code = if two_byte {
            i16::from_le_bytes(take(&mut pos, 2)?.try_into().expect("2 bytes")) as i32
        } else {
            match take(&mut pos, 1)?[0] {
                255 => i16::from_le_bytes(take(&mut pos, 2)?.try_into().expect("2 bytes")) as i32,
                c => c as i32,
            }
        };
        let value = match code {
            10..=59 | 110..=149 | 210..=239 | 460..=469 | 1010..=1059 => fmt_f64(
                f64::from_le_bytes(take(&mut pos, 8)?.try_into().expect("8 bytes")),
            ),
            60..=79 | 170..=179 | 270..=289 | 370..=389 | 400..=409 | 1060..=1070 => {
                i16::from_le_bytes(take(&mut pos, 2)?.try_into().expect("2 bytes")).to_string()
            }
            90..=99 | 420..=429 | 440..=459 | 1071 => {
                i32::from_le_bytes(take(&mut pos, 4)?.try_into().expect("4 bytes")).to_string()
            }
            160..=169 => {
                i64::from_le_bytes(take(&mut pos, 8)?.try_into().expect("8 bytes")).to_string()
            }
            290..=299 => take(&mut pos, 1)?[0].to_string(),
            310..=319 | 1004 => {
                let n = take(&mut pos, 1)?[0] as usize;
                take(&mut pos, n)?
                    .iter()
                    .map(|x| format!("{x:02X}"))
                    .collect()
            }
            0..=9
            | 100..=105
            | 300..=309
            | 320..=369
            | 390..=399
            | 410..=419
            | 430..=439
            | 470..=481
            | 999
            | 1000..=1003
            | 1005..=1009 => {
                let n = b[pos..].iter().position(|&c| c == 0).ok_or_else(short)?;
                let s = unescape_unicode(&dxf_text(&b[pos..pos + n], utf8)?);
                pos += n + 1;
                s
            }
            _ => {
                return Err(Error::format(format!(
                    "unknown binary dxf group code {code}"
                )));
            }
        };
        let eof = code == 0 && value == "EOF";
        out.push(Group {
            code,
            value,
            line: None,
        });
        if eof {
            break;
        }
    }
    Ok(out)
}

fn first(body: &[Group], code: i32) -> Option<&Group> {
    body.iter().find(|g| g.code == code)
}

fn f_or(body: &[Group], code: i32, default: f64) -> Result<f64> {
    first(body, code).map_or(Ok(default), Group::f64)
}

fn i_or(body: &[Group], code: i32, default: i64) -> Result<i64> {
    first(body, code).map_or(Ok(default), Group::i64)
}

fn point(body: &[Group], code: i32) -> Result<[f64; 3]> {
    Ok([
        f_or(body, code, 0.0)?,
        f_or(body, code + 10, 0.0)?,
        f_or(body, code + 20, 0.0)?,
    ])
}

/// A point whose x, y and z groups must all be present.
fn required_point(kind: &str, head: &Group, body: &[Group], code: i32) -> Result<[f64; 3]> {
    let mut p = [0.0; 3];
    for (k, c) in [code, code + 10, code + 20].into_iter().enumerate() {
        p[k] = first(body, c)
            .ok_or_else(|| group_err(head, format!("{kind} has no group {c}")))?
            .f64()?;
    }
    Ok(p)
}

fn layer(body: &[Group]) -> String {
    first(body, 8).map_or_else(|| "0".to_string(), |g| g.value.clone())
}

const BYLAYER: i64 = 256;

fn color(body: &[Group]) -> Result<i64> {
    i_or(body, 62, BYLAYER)
}

/// OCS to WCS by the arbitrary axis algorithm; the default extrusion is exact identity.
fn ocs(n: [f64; 3]) -> Option<[[f64; 3]; 3]> {
    if n == [0.0, 0.0, 1.0] {
        return None;
    }
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len == 0.0 || !len.is_finite() {
        return None;
    }
    let n = n.map(|x| x / len);
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let unit = |v: [f64; 3]| {
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        v.map(|x| x / l)
    };
    let ax = if n[0].abs() < 1.0 / 64.0 && n[1].abs() < 1.0 / 64.0 {
        unit(cross([0.0, 1.0, 0.0], n))
    } else {
        unit(cross([0.0, 0.0, 1.0], n))
    };
    let ay = unit(cross(n, ax));
    Some([ax, ay, n])
}

fn to_wcs(m: Option<[[f64; 3]; 3]>, p: [f64; 3]) -> [f64; 3] {
    match m {
        None => p,
        Some([ax, ay, n]) => std::array::from_fn(|k| p[0] * ax[k] + p[1] * ay[k] + p[2] * n[k]),
    }
}

#[derive(Default)]
struct Out {
    weld: HashMap<[u64; 3], u32>,
    vertices: Vec<[f64; 3]>,
    fan: Fan,
    face_layer: Vec<String>,
    face_color: Vec<i64>,
    parts: Vec<Vec<[f64; 3]>>,
    closed: Vec<bool>,
    part_layer: Vec<String>,
    part_color: Vec<i64>,
    part_feature: Vec<Option<i64>>,
    points: Vec<[f64; 3]>,
    point_layer: Vec<String>,
    point_color: Vec<i64>,
    skipped: BTreeMap<String, usize>,
    bulges: usize,
    paper: usize,
    smoothed: usize,
    short_parts: usize,
    short_faces: usize,
}

impl Out {
    fn vid(&mut self, p: [f64; 3]) -> u32 {
        let key = p.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
        let n = self.vertices.len() as u32;
        *self.weld.entry(key).or_insert_with(|| {
            self.vertices.push(p);
            n
        })
    }

    fn face(&mut self, corners: &[u32], layer: &str, color: i64) {
        let before = self.fan.tris.len();
        self.fan.polygon(corners, 0);
        for _ in before..self.fan.tris.len() {
            self.face_layer.push(layer.to_string());
            self.face_color.push(color);
        }
    }

    fn part(&mut self, v: Vec<[f64; 3]>, mut closed: bool, body: &[Group]) -> Result<()> {
        let (layer, color, feature) = (layer(body), color(body)?, xdata_feature(body)?);
        if closed && v.len() < 3 {
            closed = false;
        }
        if v.is_empty() {
            self.short_parts += 1;
            return Ok(());
        }
        self.parts.push(v);
        self.closed.push(closed);
        self.part_layer.push(layer);
        self.part_color.push(color);
        self.part_feature.push(feature);
        Ok(())
    }
}

const WIREFRAMES_APP: &str = "WIREFRAMES";

/// The feature index in a `WIREFRAMES` extended-data block, written on every part of
/// a feature with several parts.
fn xdata_feature(body: &[Group]) -> Result<Option<i64>> {
    let Some(k) = body
        .iter()
        .position(|g| g.code == 1001 && g.value == WIREFRAMES_APP)
    else {
        return Ok(None);
    };
    body[k + 1..]
        .iter()
        .take_while(|g| g.code != 1001)
        .find(|g| g.code == 1071)
        .map(Group::i64)
        .transpose()
}

fn polyline(out: &mut Out, start: &Group, head: &[Group], verts: &[&[Group]]) -> Result<()> {
    let flags = i_or(head, 70, 0)?;
    let lay = layer(head);
    let col = color(head)?;
    let vflag = |v: &[Group]| i_or(v, 70, 0);
    if flags & 64 != 0 {
        let mut local = Vec::new();
        for v in verts {
            let f = vflag(v)?;
            if f & 128 != 0 && f & 64 != 0 {
                let p = point(v, 10)?;
                local.push(out.vid(p));
            }
        }
        for v in verts {
            let f = vflag(v)?;
            if f & 128 == 0 || f & 64 != 0 {
                continue;
            }
            let mut ids = Vec::with_capacity(4);
            for code in 71..=74 {
                let Some(g) = first(v, code) else {
                    continue;
                };
                let k = g.i64()?.unsigned_abs() as usize;
                if k == 0 {
                    continue;
                }
                let id = *local.get(k - 1).ok_or_else(|| {
                    group_err(
                        g,
                        format!("polyface face uses vertex {k} of {}", local.len()),
                    )
                })?;
                ids.push(id);
            }
            if ids.len() == 4 && ids[3] == ids[2] {
                ids.pop();
            }
            if ids.len() < 3 {
                out.short_faces += 1;
                continue;
            }
            let c = match first(v, 62) {
                Some(g) => g.i64()?,
                None => col,
            };
            out.face(&ids, &lay, c);
        }
        return Ok(());
    }
    if flags & 16 != 0 {
        let fitted = i_or(head, 75, 0)? != 0 && flags & 4 != 0;
        let (m, n) = if fitted {
            (i_or(head, 73, 0)?, i_or(head, 74, 0)?)
        } else {
            (i_or(head, 71, 0)?, i_or(head, 72, 0)?)
        };
        let mut grid = Vec::new();
        for v in verts {
            if fitted && vflag(v)? & 16 != 0 {
                continue;
            }
            let p = point(v, 10)?;
            grid.push(out.vid(p));
        }
        let (m, n) = (m.max(0) as usize, n.max(0) as usize);
        if m.checked_mul(n) != Some(grid.len()) {
            return Err(group_err(
                start,
                format!(
                    "polygon mesh declares {m} x {n} vertices but has {}",
                    grid.len()
                ),
            ));
        }
        let (wrap_m, wrap_n) = (flags & 1 != 0, flags & 32 != 0);
        let rows = if wrap_m { m } else { m.saturating_sub(1) };
        let cols = if wrap_n { n } else { n.saturating_sub(1) };
        for i in 0..rows {
            for j in 0..cols {
                let (i1, j1) = ((i + 1) % m, (j + 1) % n);
                let q = [
                    grid[i * n + j],
                    grid[i * n + j1],
                    grid[i1 * n + j1],
                    grid[i1 * n + j],
                ];
                out.face(&q, &lay, col);
            }
        }
        return Ok(());
    }
    let mut pts = Vec::with_capacity(verts.len());
    if flags & 8 != 0 {
        for v in verts {
            if vflag(v)? & 16 == 0 {
                pts.push(point(v, 10)?);
            }
        }
    } else {
        let elevation = f_or(head, 30, 0.0)?;
        let m = ocs(extrusion(head)?);
        for v in verts {
            if vflag(v)? & 16 != 0 {
                continue;
            }
            if f_or(v, 42, 0.0)? != 0.0 {
                out.bulges += 1;
            }
            pts.push(to_wcs(m, [f_or(v, 10, 0.0)?, f_or(v, 20, 0.0)?, elevation]));
        }
    }
    out.part(pts, flags & 1 != 0, head)
}

fn extrusion(body: &[Group]) -> Result<[f64; 3]> {
    Ok([
        f_or(body, 210, 0.0)?,
        f_or(body, 220, 0.0)?,
        f_or(body, 230, 1.0)?,
    ])
}

fn lwpolyline(out: &mut Out, body: &[Group]) -> Result<()> {
    let elevation = f_or(body, 38, 0.0)?;
    let m = ocs(extrusion(body)?);
    let mut xy: Vec<[f64; 2]> = Vec::new();
    for g in body {
        match g.code {
            10 => xy.push([g.f64()?, 0.0]),
            20 => {
                if let Some(p) = xy.last_mut() {
                    p[1] = g.f64()?;
                }
            }
            42 if g.f64()? != 0.0 => out.bulges += 1,
            _ => {}
        }
    }
    let pts = xy
        .iter()
        .map(|p| to_wcs(m, [p[0], p[1], elevation]))
        .collect();
    out.part(pts, i_or(body, 70, 0)? & 1 != 0, body)
}

fn mesh_entity(out: &mut Out, head: &Group, body: &[Group]) -> Result<()> {
    let bad = |g: &Group, m: &str| group_err(g, format!("MESH entity: {m}"));
    let mut verts: Option<Vec<u32>> = None;
    let mut list: Option<&[Group]> = None;
    let mut smoothed = false;
    let mut j = 0;
    while j < body.len() {
        let at = &body[j];
        match at.code {
            91 if verts.is_none() => {
                smoothed |= at.i64()? > 0;
                j += 1;
            }
            92 if verts.is_none() => {
                let n = usize::try_from(at.i64()?).map_err(|_| bad(at, "negative vertex count"))?;
                j += 1;
                if n.checked_mul(3).is_none_or(|k| j + k > body.len()) {
                    return Err(bad(at, "vertex list ends early"));
                }
                let mut v = Vec::with_capacity(n);
                for _ in 0..n {
                    if body[j].code != 10 || body[j + 1].code != 20 || body[j + 2].code != 30 {
                        return Err(bad(&body[j], "expected vertex groups 10, 20, 30"));
                    }
                    let p = [body[j].f64()?, body[j + 1].f64()?, body[j + 2].f64()?];
                    v.push(out.vid(p));
                    j += 3;
                }
                verts = Some(v);
            }
            93 if list.is_none() => {
                let n =
                    usize::try_from(at.i64()?).map_err(|_| bad(at, "negative face list size"))?;
                j += 1;
                if j + n > body.len() {
                    return Err(bad(at, "face list ends early"));
                }
                for g in &body[j..j + n] {
                    if g.code != 90 {
                        return Err(bad(g, "expected face list group 90"));
                    }
                }
                list = Some(&body[j..j + n]);
                j += n;
            }
            _ => j += 1,
        }
    }
    out.smoothed += usize::from(smoothed);
    let verts = verts.unwrap_or_default();
    let list = list.unwrap_or_default();
    let (lay, col) = (layer(body), color(body)?);
    let mut i = 0;
    while i < list.len() {
        let n = usize::try_from(list[i].i64()?).map_err(|_| bad(&list[i], "negative face size"))?;
        let ids = list
            .get(i + 1..i + 1 + n)
            .ok_or_else(|| bad(list.last().unwrap_or(head), "face list ends inside a face"))?
            .iter()
            .map(|g| {
                let k = g.i64()?;
                usize::try_from(k)
                    .ok()
                    .and_then(|k| verts.get(k).copied())
                    .ok_or_else(|| bad(g, &format!("face uses vertex {k} of {}", verts.len())))
            })
            .collect::<Result<Vec<u32>>>()?;
        if ids.len() < 3 {
            out.short_faces += 1;
        } else {
            out.face(&ids, &lay, col);
        }
        i += 1 + n;
    }
    Ok(())
}

fn unit_name(code: i64) -> Option<&'static str> {
    Some(match code {
        1 => "in",
        2 => "ft",
        4 => "mm",
        5 => "cm",
        6 => "m",
        7 => "km",
        _ => return None,
    })
}

fn unit_code(name: &str) -> i64 {
    match name.trim().to_ascii_lowercase().as_str() {
        "in" => 1,
        "ft" | "foot" | "feet" => 2,
        "mm" => 4,
        "cm" => 5,
        "m" | "metre" | "meter" | "metres" | "meters" => 6,
        "km" => 7,
        _ => 0,
    }
}

fn attr_table(layers: Vec<String>, colors: Vec<i64>) -> Result<Table> {
    let mut t = Table::new(layers.len()).with_column("layer", Column::utf8(layers))?;
    if colors.iter().any(|&c| c != BYLAYER) {
        t = t.with_column("color", Column::i64(colors))?;
    }
    Ok(t)
}

pub fn read(r: &[u8]) -> Result<Data> {
    let binary = r.starts_with(BINARY_SENTINEL);
    let groups = match r.strip_prefix(BINARY_SENTINEL) {
        Some(rest) => binary_groups(rest)?,
        None => ascii_groups(r.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(r))?,
    };
    let ends_with_eof = groups
        .last()
        .is_some_and(|g| g.code == 0 && g.value == "EOF");
    if !ends_with_eof {
        let m = "file ends without the 0/EOF group; it is truncated";
        return Err(match groups.last() {
            Some(g) => group_err(g, m.to_string()),
            None => Error::format(m),
        });
    }
    let mut entities: Vec<(&str, &Group, &[Group])> = Vec::new();
    let mut section = String::new();
    let mut length_unit = None;
    let mut i = 0;
    while i < groups.len() {
        if groups[i].code != 0 {
            i += 1;
            continue;
        }
        let kind = groups[i].value.as_str();
        let start = i + 1;
        let mut end = start;
        while end < groups.len() && groups[end].code != 0 {
            end += 1;
        }
        let body = &groups[start..end];
        match kind {
            "SECTION" => {
                section = first(body, 2).map(|g| g.value.clone()).unwrap_or_default();
                if section == "HEADER" {
                    for (k, g) in body.iter().enumerate() {
                        if g.code == 9
                            && g.value == "$INSUNITS"
                            && let Some(v) = body.get(k + 1).filter(|v| v.code == 70)
                        {
                            length_unit = unit_name(v.i64()?).map(String::from);
                        }
                    }
                }
            }
            "ENDSEC" => section.clear(),
            "EOF" => break,
            _ if section == "ENTITIES" => entities.push((kind, &groups[i], body)),
            _ => {}
        }
        i = end;
    }
    let mut out = Out::default();
    let mut k = 0;
    while k < entities.len() {
        let (kind, head, body) = entities[k];
        k += 1;
        if i_or(body, 67, 0)? == 1 && !matches!(kind, "VERTEX" | "SEQEND") {
            out.paper += 1;
            if kind == "POLYLINE" {
                while k < entities.len() && entities[k].0 == "VERTEX" {
                    k += 1;
                }
                if k < entities.len() && entities[k].0 == "SEQEND" {
                    k += 1;
                }
            }
            continue;
        }
        match kind {
            "3DFACE" => {
                let c = [
                    required_point(kind, head, body, 10)?,
                    required_point(kind, head, body, 11)?,
                    required_point(kind, head, body, 12)?,
                ];
                let mut ids = vec![out.vid(c[0]), out.vid(c[1]), out.vid(c[2])];
                if first(body, 13).is_some()
                    || first(body, 23).is_some()
                    || first(body, 33).is_some()
                {
                    let d = out.vid(required_point(kind, head, body, 13)?);
                    if d != ids[2] {
                        ids.push(d);
                    }
                }
                out.face(&ids, &layer(body), color(body)?);
            }
            "POLYLINE" => {
                let mut verts = Vec::new();
                while k < entities.len() && entities[k].0 == "VERTEX" {
                    verts.push(entities[k].2);
                    k += 1;
                }
                if k < entities.len() && entities[k].0 == "SEQEND" {
                    k += 1;
                }
                polyline(&mut out, head, body, &verts)?;
            }
            "LWPOLYLINE" => lwpolyline(&mut out, body)?,
            "LINE" => {
                let v = vec![point(body, 10)?, point(body, 11)?];
                out.part(v, false, body)?;
            }
            "POINT" => {
                out.points.push(point(body, 10)?);
                out.point_layer.push(layer(body));
                out.point_color.push(color(body)?);
            }
            "MESH" if !binary => mesh_entity(&mut out, head, body)?,
            "VERTEX" | "SEQEND" => {}
            other => *out.skipped.entry(other.to_string()).or_default() += 1,
        }
    }
    let mut warnings: Vec<String> = out
        .skipped
        .iter()
        .map(|(kind, n)| format!("skipped {n} {kind} entities"))
        .collect();
    if out.bulges > 0 {
        warnings.push(format!(
            "ignored nonzero bulge on {} polyline vertices; arcs read as straight segments",
            out.bulges
        ));
    }
    if out.paper > 0 {
        warnings.push(format!(
            "skipped {} paper space entities (group 67 = 1); only model space is read",
            out.paper
        ));
    }
    if out.smoothed > 0 {
        warnings.push(format!(
            "read the level-0 cage of {} MESH entities with a smoothing level above 0; \
             CAD shows the subdivided surface, which is not read",
            out.smoothed
        ));
    }
    if out.short_parts > 0 {
        warnings.push(format!(
            "skipped {} polylines with no vertices",
            out.short_parts
        ));
    }
    if out.short_faces > 0 {
        warnings.push(format!(
            "skipped {} faces with fewer than 3 corners",
            out.short_faces
        ));
    }
    let mut data = Data {
        length_unit: length_unit.clone(),
        warnings,
        ..Data::default()
    };
    if !out.fan.tris.is_empty() || !out.vertices.is_empty() {
        let n = out.vertices.len();
        let faces = attr_table(out.face_layer, out.face_color)?;
        let mesh = out.fan.mesh(out.vertices, Table::new(n), faces)?;
        data.mesh = Some(mesh.with_length_unit(length_unit.clone()));
    }
    if !out.parts.is_empty() {
        let tags = &out.part_feature;
        let mut feature_of_part = Vec::with_capacity(tags.len());
        let mut firsts = Vec::new();
        for (p, tag) in tags.iter().enumerate() {
            if p == 0 || tag.is_none() || *tag != tags[p - 1] {
                firsts.push(p as u32);
            }
            feature_of_part.push(firsts.len() as u32 - 1);
        }
        let parts = attr_table(out.part_layer, out.part_color)?;
        let pl = Polylines::from_parts(out.parts, out.closed, Some(feature_of_part))?
            .with_attributes(parts.take(&firsts))?;
        data.polylines = Some(pl.with_length_unit(length_unit.clone()));
    }
    if !out.points.is_empty() {
        let ps = PointSet::new(out.points)?
            .with_attributes(attr_table(out.point_layer, out.point_color)?)?;
        data.points = Some(ps.with_length_unit(length_unit));
    }
    Ok(data)
}

/// Per-row layer and color (`None` means BYLAYER) from a table holding at most
/// `layer`, `color` and the names in `extra`.
fn layer_color(
    table: &Table,
    extra: &[&str],
    what: &str,
    opts: &WriteOptions,
) -> Result<(Vec<String>, Vec<Option<i64>>)> {
    let rows = table.rows();
    let mut layers = vec!["0".to_string(); rows];
    let mut colors = vec![None; rows];
    for (name, col) in table.names().iter().zip(table.columns()) {
        match (name.as_str(), &col.data) {
            ("layer", ColumnData::Utf8(v)) => {
                for i in 0..rows {
                    if col.is_valid(i) {
                        layers[i] = escape_layer(&v[i], what)?;
                    }
                }
            }
            ("color", ColumnData::Int64(v)) => {
                for i in 0..rows {
                    if col.is_valid(i) && v[i] != BYLAYER {
                        colors[i] = Some(v[i]);
                    }
                }
            }
            (n, _) if extra.contains(&n) => {}
            _ if opts.skip_unsupported_columns => {}
            _ => {
                return Err(Error::Unsupported(format!(
                    "dxf cannot store {what} column {name:?} ({}); it stores layer (utf8) and color (int64)",
                    col.type_name()
                )));
            }
        }
    }
    Ok((layers, colors))
}

fn no_columns(table: &Table, what: &str, opts: &WriteOptions) -> Result<()> {
    match table.names().first() {
        Some(n) if !opts.skip_unsupported_columns => Err(Error::Unsupported(format!(
            "dxf cannot store {what} column {n:?}"
        ))),
        _ => Ok(()),
    }
}

struct W(String);

impl W {
    fn g(&mut self, code: i32, value: impl std::fmt::Display) {
        let _ = write!(self.0, "{code:>3}\n{value}\n");
    }

    fn f(&mut self, code: i32, x: f64) {
        self.g(code, fmt_f64(x));
    }

    fn p(&mut self, code: i32, p: [f64; 3]) {
        self.f(code, p[0]);
        self.f(code + 10, p[1]);
        self.f(code + 20, p[2]);
    }

    fn entity(&mut self, kind: &str, layer: &str, color: Option<i64>) {
        self.g(0, kind);
        self.g(8, layer);
        if let Some(c) = color {
            self.g(62, c);
        }
    }
}

const POLYFACE_MAX: usize = 32767;

pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    if data.grid.is_some() {
        return Err(Error::Unsupported(
            "dxf cannot store an unstructured grid; write .vtu".into(),
        ));
    }
    let mut ents = W(String::new());
    let mut used: Vec<String> = vec!["0".to_string()];
    let mut note = |l: &str| {
        if !used.iter().any(|u| u == l) {
            used.push(l.to_string());
        }
    };
    if let Some(m) = &data.mesh {
        no_columns(m.vertex_attributes(), "vertex", opts)?;
        let (polys, table) = polygons(m)?;
        let (layers, colors) = layer_color(&table, &[], "face", opts)?;
        let mut faces: Vec<(Vec<u32>, usize)> = Vec::with_capacity(polys.len());
        for (k, p) in polys.into_iter().enumerate() {
            if p.len() <= 4 {
                faces.push((p, k));
            } else {
                for i in 1..p.len() - 1 {
                    faces.push((vec![p[0], p[i], p[i + 1]], k));
                }
            }
        }
        let v = m.vertices();
        match opts.dxf_entity {
            DxfEntity::Face3d => {
                for (f, k) in &faces {
                    note(&layers[*k]);
                    ents.entity("3DFACE", &layers[*k], colors[*k]);
                    for (c, &id) in f.iter().enumerate() {
                        ents.p(10 + c as i32, v[id as usize]);
                    }
                    if f.len() == 3 {
                        ents.p(13, v[f[2] as usize]);
                    }
                }
            }
            DxfEntity::Polyface => {
                let mut s = 0;
                while s < faces.len() {
                    let lay = &layers[faces[s].1];
                    note(lay);
                    let mut local: HashMap<u32, usize> = HashMap::new();
                    let mut order: Vec<u32> = Vec::new();
                    let mut e = s;
                    while e < faces.len() && layers[faces[e].1] == *lay && e - s < POLYFACE_MAX {
                        let new = faces[e]
                            .0
                            .iter()
                            .filter(|id| !local.contains_key(id))
                            .count();
                        if order.len() + new > POLYFACE_MAX {
                            break;
                        }
                        for &id in &faces[e].0 {
                            local.entry(id).or_insert_with(|| {
                                order.push(id);
                                order.len()
                            });
                        }
                        e += 1;
                    }
                    ents.entity("POLYLINE", lay, None);
                    ents.g(66, 1);
                    ents.p(10, [0.0; 3]);
                    ents.g(70, 64);
                    ents.g(71, order.len());
                    ents.g(72, e - s);
                    for &id in &order {
                        ents.entity("VERTEX", lay, None);
                        ents.p(10, v[id as usize]);
                        ents.g(70, 192);
                    }
                    for (f, k) in &faces[s..e] {
                        ents.entity("VERTEX", lay, colors[*k]);
                        ents.p(10, [0.0; 3]);
                        ents.g(70, 128);
                        for (c, id) in f.iter().enumerate() {
                            ents.g(71 + c as i32, local[id]);
                        }
                    }
                    ents.entity("SEQEND", lay, None);
                    s = e;
                }
            }
        }
    }
    if let Some(pl) = &data.polylines {
        no_columns(pl.vertex_attributes(), "polyline vertex", opts)?;
        let (layers, colors) = layer_color(pl.attributes(), &[], "polyline", opts)?;
        for f in 0..pl.len() {
            note(&layers[f]);
            let parts = pl.feature_parts(f);
            if parts.is_empty() {
                return Err(Error::Unsupported(format!(
                    "dxf cannot store polyline feature {f}: it has no parts"
                )));
            }
            let multipart = parts.len() > 1;
            for p in parts {
                ents.entity("POLYLINE", &layers[f], colors[f]);
                ents.g(66, 1);
                ents.p(10, [0.0; 3]);
                ents.g(70, 8 | pl.closed()[p] as i32);
                if multipart {
                    ents.g(1001, WIREFRAMES_APP);
                    ents.g(1071, f);
                }
                for &q in pl.part(p) {
                    ents.entity("VERTEX", &layers[f], None);
                    ents.p(10, q);
                    ents.g(70, 32);
                }
                ents.entity("SEQEND", &layers[f], None);
            }
        }
    }
    if let Some(ps) = &data.points {
        let (layers, colors) = layer_color(ps.attributes(), &[], "point", opts)?;
        for (i, &q) in ps.coords().iter().enumerate() {
            note(&layers[i]);
            ents.entity("POINT", &layers[i], colors[i]);
            ents.p(10, q);
        }
    }
    let unit = data
        .length_unit
        .as_deref()
        .or_else(|| data.mesh.as_ref().and_then(|m| m.length_unit.as_deref()))
        .or_else(|| {
            data.polylines
                .as_ref()
                .and_then(|m| m.length_unit.as_deref())
        })
        .or_else(|| data.points.as_ref().and_then(|m| m.length_unit.as_deref()))
        .map_or(0, unit_code);
    let mut w = W(String::new());
    w.g(0, "SECTION");
    w.g(2, "HEADER");
    w.g(9, "$ACADVER");
    w.g(1, "AC1015");
    w.g(9, "$INSUNITS");
    w.g(70, unit);
    w.g(0, "ENDSEC");
    w.g(0, "SECTION");
    w.g(2, "TABLES");
    w.g(0, "TABLE");
    w.g(2, "LAYER");
    w.g(70, used.len());
    for l in &used {
        w.g(0, "LAYER");
        w.g(2, l);
        w.g(70, 0);
        w.g(62, 7);
        w.g(6, "CONTINUOUS");
    }
    w.g(0, "ENDTAB");
    w.g(0, "TABLE");
    w.g(2, "APPID");
    w.g(70, 1);
    w.g(0, "APPID");
    w.g(2, WIREFRAMES_APP);
    w.g(70, 0);
    w.g(0, "ENDTAB");
    w.g(0, "ENDSEC");
    w.g(0, "SECTION");
    w.g(2, "ENTITIES");
    w.0.push_str(&ents.0);
    w.g(0, "ENDSEC");
    w.g(0, "EOF");
    Ok(w.0.into_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn face(corners: &str) -> String {
        format!("0\nSECTION\n2\nENTITIES\n0\n3DFACE\n8\nA\n{corners}0\nENDSEC\n0\nEOF\n")
    }

    const CORNERS: &str = "10\n0\n20\n0\n30\n0\n11\n1\n21\n0\n31\n0\n12\n0\n22\n1\n32\n0\n";

    #[test]
    fn a_file_without_eof_is_an_error_naming_its_last_line() {
        let full = face(CORNERS);
        assert!(read(full.as_bytes()).is_ok());
        let cut = &full[..full.find("0\nENDSEC").unwrap()];
        let e = read(cut.as_bytes()).unwrap_err().to_string();
        assert!(e.contains("line 26") && e.contains("EOF"), "{e}");
        assert!(read(b"").is_err());
    }

    #[test]
    fn paper_space_entities_are_skipped_with_a_warning() {
        let two = face(CORNERS).replace(
            "0\nENDSEC",
            &format!("0\n3DFACE\n8\nTITLE\n67\n1\n{CORNERS}0\nPOLYLINE\n67\n1\n70\n8\n0\nVERTEX\n10\n0\n20\n0\n30\n0\n0\nSEQEND\n0\nENDSEC"),
        );
        let d = read(two.as_bytes()).unwrap();
        assert_eq!(d.mesh.unwrap().n_faces(), 1);
        assert!(d.polylines.is_none());
        assert_eq!(
            d.warnings,
            ["skipped 2 paper space entities (group 67 = 1); only model space is read"]
        );
    }

    #[test]
    fn a_utf8_byte_order_mark_is_skipped() {
        let mut b = b"\xEF\xBB\xBF".to_vec();
        b.extend(face(CORNERS).bytes());
        assert_eq!(read(&b).unwrap().mesh.unwrap().n_faces(), 1);
    }

    fn mesh(level: i64, count: &str, face: &str) -> String {
        let v = "10\n0\n20\n0\n30\n0\n10\n1\n20\n0\n30\n0\n10\n0\n20\n1\n30\n0\n";
        format!(
            "0\nSECTION\n2\nENTITIES\n0\nMESH\n8\nA\n71\n2\n72\n0\n91\n{level}\n92\n{count}\n{v}93\n4\n90\n3\n{face}0\nENDSEC\n0\nEOF\n"
        )
    }

    #[test]
    fn a_smoothed_mesh_warns_that_only_its_cage_is_read() {
        let flat = read(mesh(0, "3", "90\n0\n90\n1\n90\n2\n").as_bytes()).unwrap();
        assert!(flat.warnings.is_empty());
        let smooth = read(mesh(2, "3", "90\n0\n90\n1\n90\n2\n").as_bytes()).unwrap();
        assert_eq!(smooth.mesh.unwrap().n_faces(), 1);
        assert!(smooth.warnings[0].contains("level-0 cage of 1 MESH"));
    }

    #[test]
    fn entity_errors_name_the_line_of_the_group() {
        let e = |s: String| read(s.as_bytes()).unwrap_err().to_string();
        let count = e(mesh(0, "2000000000", "90\n0\n90\n1\n90\n2\n"));
        assert!(
            count.starts_with("line 16: MESH entity: vertex list ends early"),
            "{count}"
        );
        let index = e(mesh(0, "3", "90\n0\n90\n9\n90\n2\n"));
        assert!(
            index.starts_with("line 42: MESH entity: face uses vertex 9"),
            "{index}"
        );
        let nan = e(face(&CORNERS.replace("21\n0\n", "21\nnan\n")));
        assert!(
            nan.starts_with("line 18: group 21 value \"nan\" is not finite"),
            "{nan}"
        );
    }

    #[test]
    fn a_3dface_missing_a_corner_group_is_an_error() {
        let e = read(face(&CORNERS.replace("22\n1\n", "")).as_bytes())
            .unwrap_err()
            .to_string();
        assert!(e.contains("3DFACE has no group 22"), "{e}");
    }

    #[test]
    fn cp1252_fallback_decodes_high_bytes() {
        assert_eq!(decode_text(b"caf\xe9 \x80"), "caf\u{e9} \u{20AC}");
        assert_eq!(decode_text("caf\u{e9}".as_bytes()), "caf\u{e9}");
    }

    #[test]
    fn unicode_escapes_are_replaced() {
        assert_eq!(unescape_unicode("A\\U+00E9B"), "A\u{e9}B");
        assert_eq!(unescape_unicode("bad\\U+zz"), "bad\\U+zz");
        assert_eq!(unescape_unicode("\\u+0041\\U+0042"), "AB");
        assert_eq!(unescape_unicode("\\U++041"), "\\U++041");
        assert_eq!(unescape_unicode("\\U+005CU+0041"), "\\U+0041");
    }

    #[test]
    fn layer_escape_inverts_unescape() {
        for l in [" A ", "X\\U+0041", "caf\u{e9}", "a\tb", "\\", "0"] {
            let e = escape_layer(l, "t").unwrap();
            assert!(e.bytes().all(|b| b.is_ascii_graphic()) || e.contains(' '));
            assert_eq!(e.trim(), e);
            assert_eq!(unescape_unicode(&e), l);
        }
        assert!(escape_layer("", "t").is_err());
        assert!(escape_layer("a\nb", "t").is_err());
        assert!(escape_layer("\u{1F600}", "t").is_err());
    }

    #[test]
    fn version_gates_the_cp1252_fallback() {
        let head = |v: &str| format!("  0\nSECTION\n  2\nHEADER\n  9\n$ACADVER\n  1\n{v}\n");
        assert!(!utf8_only(head("AC1015").as_bytes()));
        assert!(utf8_only(head("AC1021").as_bytes()));
        assert!(utf8_only(head("AC1032").as_bytes()));
        assert!(!utf8_only(b"no header"));
        assert!(dxf_text(b"caf\xe9", true).is_err());
        assert_eq!(dxf_text(b"caf\xe9", false).unwrap(), "caf\u{e9}");
    }

    #[test]
    fn arbitrary_axis_matches_identity_and_flip() {
        assert!(ocs([0.0, 0.0, 1.0]).is_none());
        let m = ocs([0.0, 0.0, -1.0]);
        let p = to_wcs(m, [1.0, 2.0, 3.0]);
        assert_eq!(p, [-1.0, 2.0, -3.0]);
    }
}
