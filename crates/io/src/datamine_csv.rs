//! Datamine wireframe CSV pairs: `_pt.csv` points and `_tr.csv` triangles.

use std::collections::HashMap;

use drusa_core::{Column, ColumnData, Error, Mesh, Result, Table};

use crate::dxf::decode_text;
use crate::vtk::fmt_f64;
use crate::{Data, WriteOptions};

const NULL: f64 = -1.0e30;
const TRIANGLE_IDS: [&str; 4] = ["TRIANGLE", "TRIANG", "TRI", "TRISOL"];

fn unquote(f: &str) -> String {
    let f = f.trim();
    match f.strip_prefix('"').and_then(|s| s.strip_suffix('"')) {
        Some(inner) => inner.replace("\"\"", "\""),
        None => f.to_string(),
    }
}

/// `;` when the header record has a `;` outside quotes, else `,`.
fn sniff(text: &str) -> char {
    let mut quoted = false;
    for c in text.chars() {
        match c {
            '"' => quoted = !quoted,
            ';' if !quoted => return ';',
            '\n' if !quoted => break,
            _ => {}
        }
    }
    ','
}

/// Records with their first line number; delimiters and line breaks inside
/// quotes are field data. Blank records are skipped.
fn records(text: &str, d: char, file: &str) -> Result<Vec<(u64, Vec<String>)>> {
    let mut out = Vec::new();
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    let mut line = 1u64;
    let mut start = 1u64;
    let mut end_record = |fields: &mut Vec<String>, cur: &mut String, start: u64| {
        if cur.ends_with('\r') {
            cur.pop();
        }
        fields.push(unquote(&std::mem::take(cur)));
        let fields = std::mem::take(fields);
        if !(fields.len() == 1 && fields[0].is_empty()) {
            out.push((start, fields));
        }
    };
    for c in text.chars() {
        match c {
            '"' => {
                quoted = !quoted;
                cur.push(c);
            }
            '\n' if !quoted => {
                end_record(&mut fields, &mut cur, start);
                line += 1;
                start = line;
            }
            '\n' => {
                line += 1;
                cur.push(c);
            }
            c if c == d && !quoted => fields.push(unquote(&std::mem::take(&mut cur))),
            c => cur.push(c),
        }
    }
    if quoted {
        return Err(Error::format_at(
            start,
            format!("{file} has a quote that is never closed"),
        ));
    }
    end_record(&mut fields, &mut cur, start);
    Ok(out)
}

struct Csv {
    file: &'static str,
    delim: char,
    header: Vec<String>,
    rows: Vec<(u64, Vec<String>)>,
}

impl Csv {
    fn parse(b: &[u8], file: &'static str) -> Result<Csv> {
        let text = decode_text(b);
        let text = text.strip_prefix('\u{feff}').unwrap_or(&text);
        let delim = sniff(text);
        let mut lines = records(text, delim, file)?.into_iter();
        let (_, header) = lines
            .next()
            .ok_or_else(|| Error::format(format!("{file} is empty")))?;
        let mut rows = Vec::new();
        for (n, mut f) in lines {
            while f.len() > header.len() && f.last().is_some_and(String::is_empty) {
                f.pop();
            }
            if f.len() != header.len() {
                return Err(Error::format_at(
                    n,
                    format!(
                        "{file} row has {} fields, header has {}",
                        f.len(),
                        header.len()
                    ),
                ));
            }
            rows.push((n, f));
        }
        Ok(Csv {
            file,
            delim,
            header,
            rows,
        })
    }

    fn col(&self, names: &[&str]) -> Option<usize> {
        self.header
            .iter()
            .position(|h| names.iter().any(|n| h.eq_ignore_ascii_case(n)))
    }

    fn need(&self, name: &str) -> Result<usize> {
        self.col(&[name]).ok_or_else(|| {
            Error::format(format!(
                "{} has no {name} column; columns are {:?}",
                self.file, self.header
            ))
        })
    }

    fn num(&self, s: &str) -> Option<f64> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        if self.delim == ';' {
            s.replace(',', ".").parse().ok()
        } else {
            s.parse().ok()
        }
    }

    /// Drops a type or unit row right after the header.
    fn skip_type_row(&mut self, probe: usize) {
        if self
            .rows
            .first()
            .is_some_and(|(_, r)| self.num(&r[probe]).is_none())
        {
            self.rows.remove(0);
        }
    }

    fn value(&self, line: u64, row: &[String], c: usize) -> Result<f64> {
        match self.num(&row[c]) {
            Some(v) if v != NULL => Ok(v),
            _ => Err(Error::format_at(
                line,
                format!(
                    "{} {} {:?} is not a number",
                    self.file, self.header[c], row[c]
                ),
            )),
        }
    }

    fn id(&self, line: u64, row: &[String], c: usize) -> Result<i64> {
        let v = self.value(line, row, c)?;
        if v.fract() != 0.0 || v.abs() >= 9.0e15 {
            return Err(Error::format_at(
                line,
                format!("{} {} {v} is not an integer", self.file, self.header[c]),
            ));
        }
        Ok(v as i64)
    }

    /// Extra columns: Float64 when every non-empty value is numeric, else Utf8.
    /// Empty fields and `-1.0E30` are null.
    fn extras(&self, skip: &[usize]) -> Result<Table> {
        let mut t = Table::new(self.rows.len());
        for (c, name) in self.header.iter().enumerate() {
            if skip.contains(&c) {
                continue;
            }
            let raw: Vec<&str> = self.rows.iter().map(|(_, r)| r[c].as_str()).collect();
            let nums: Option<Vec<Option<f64>>> = raw
                .iter()
                .map(|s| {
                    if s.is_empty() {
                        Some(None)
                    } else {
                        self.num(s).map(|v| (v != NULL).then_some(v))
                    }
                })
                .collect();
            let valid: Vec<bool>;
            let col = match nums {
                Some(v) => {
                    valid = v.iter().map(Option::is_some).collect();
                    Column::f64(v.iter().map(|x| x.unwrap_or(f64::NAN)).collect())
                }
                None => {
                    valid = raw.iter().map(|s| !s.is_empty()).collect();
                    Column::utf8(raw.iter().map(|s| s.to_string()).collect())
                }
            };
            let col = if valid.iter().all(|&v| v) {
                col
            } else {
                col.with_valid(valid)?
            };
            t = t.with_column(name, col)?;
        }
        Ok(t)
    }
}

/// Reads a `_pt.csv` and `_tr.csv` pair.
pub fn read(pt: &[u8], tr: &[u8]) -> Result<Data> {
    let mut p = Csv::parse(pt, "_pt.csv")?;
    let cols = [p.need("XP")?, p.need("YP")?, p.need("ZP")?, p.need("PID")?];
    p.skip_type_row(cols[0]);
    let mut vertices = Vec::with_capacity(p.rows.len());
    let mut index: HashMap<i64, u32> = HashMap::with_capacity(p.rows.len());
    for (k, (line, row)) in p.rows.iter().enumerate() {
        let xyz = [
            p.value(*line, row, cols[0])?,
            p.value(*line, row, cols[1])?,
            p.value(*line, row, cols[2])?,
        ];
        let pid = p.id(*line, row, cols[3])?;
        if index.insert(pid, k as u32).is_some() {
            return Err(Error::format_at(
                *line,
                format!("duplicate PID {pid} in _pt.csv"),
            ));
        }
        vertices.push(xyz);
    }
    let vertex_table = p.extras(&cols)?;

    let mut t = Csv::parse(tr, "_tr.csv")?;
    let tcols = [t.need("PID1")?, t.need("PID2")?, t.need("PID3")?];
    t.skip_type_row(tcols[0]);
    let mut skip = tcols.to_vec();
    skip.extend(t.col(&TRIANGLE_IDS));
    let mut triangles = Vec::with_capacity(t.rows.len());
    for (line, row) in &t.rows {
        let mut tri = [0u32; 3];
        for (k, &c) in tcols.iter().enumerate() {
            let pid = t.id(*line, row, c)?;
            tri[k] = *index.get(&pid).ok_or_else(|| {
                Error::format_at(
                    *line,
                    format!("_tr.csv uses PID {pid}, missing from _pt.csv"),
                )
            })?;
        }
        triangles.push(tri);
    }
    let face_table = t.extras(&skip)?;
    let mesh = Mesh::new(vertices, triangles)?.with_attributes(vertex_table, face_table)?;
    Ok(Data {
        mesh: Some(mesh),
        ..Data::default()
    })
}

fn field(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r', ';']) || s.trim() != s {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

fn cell(col: &Column, i: usize) -> String {
    let valid = col.is_valid(i);
    match &col.data {
        ColumnData::Float64(v) => fmt_f64(if valid { v[i] } else { NULL }),
        ColumnData::Int64(v) if valid => v[i].to_string(),
        ColumnData::Bool(v) if valid => (v[i] as u8).to_string(),
        ColumnData::Utf8(v) if valid => field(&v[i]),
        ColumnData::Int64(_) | ColumnData::Bool(_) => fmt_f64(NULL),
        ColumnData::Utf8(_) => String::new(),
    }
}

fn header(fixed: &[&str], t: &Table, reserved: &[&str], what: &str) -> Result<String> {
    let mut h: Vec<String> = fixed.iter().map(|s| s.to_string()).collect();
    for name in t.names() {
        if name.is_empty() || reserved.iter().any(|r| r.eq_ignore_ascii_case(name)) {
            return Err(Error::invalid(format!(
                "{what} column {name:?} clashes with a datamine field"
            )));
        }
        h.push(field(name));
    }
    Ok(h.join(","))
}

/// Writes the `(pt, tr)` pair. PIDs are 1-based in vertex order.
pub fn write(data: &Data, _opts: &WriteOptions) -> Result<(Vec<u8>, Vec<u8>)> {
    if data.polylines.is_some() || data.points.is_some() || data.grid.is_some() {
        return Err(Error::Unsupported(
            "a datamine wireframe pair holds only a mesh".into(),
        ));
    }
    let mesh = data
        .mesh
        .as_ref()
        .ok_or_else(|| Error::Unsupported("a datamine wireframe pair needs a mesh".into()))?;
    let vt = mesh.vertex_attributes();
    let ft = mesh.face_attributes();
    let mut reserved = vec!["XP", "YP", "ZP", "PID"];
    let mut pt = header(&["XP", "YP", "ZP", "PID"], vt, &reserved, "vertex")?;
    pt.push('\n');
    for (i, v) in mesh.vertices().iter().enumerate() {
        let mut row = vec![
            fmt_f64(v[0]),
            fmt_f64(v[1]),
            fmt_f64(v[2]),
            (i + 1).to_string(),
        ];
        row.extend(vt.columns().iter().map(|c| cell(c, i)));
        pt.push_str(&row.join(","));
        pt.push('\n');
    }
    reserved = vec!["PID1", "PID2", "PID3"];
    reserved.extend(TRIANGLE_IDS);
    let mut tr = header(&["TRIANGLE", "PID1", "PID2", "PID3"], ft, &reserved, "face")?;
    tr.push('\n');
    for (i, t) in mesh.triangles().iter().enumerate() {
        let mut row = vec![(i + 1).to_string()];
        row.extend(t.iter().map(|&k| (k + 1).to_string()));
        row.extend(ft.columns().iter().map(|c| cell(c, i)));
        tr.push_str(&row.join(","));
        tr.push('\n');
    }
    Ok((pt.into_bytes(), tr.into_bytes()))
}
