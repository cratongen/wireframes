//! Datamine wireframe CSV pairs: `_pt.csv` points and `_tr.csv` triangles.

use std::collections::HashMap;
use std::path::Path;

use wireframes_core::{Column, ColumnData, Error, Mesh, Result, Table};

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

/// A field's text and whether it was quoted.
type Field = (String, bool);

fn field_of(raw: &str) -> Field {
    (unquote(raw), raw.trim().starts_with('"'))
}

/// Records with their first line number; delimiters and line breaks inside
/// quotes are field data. Blank records are skipped.
fn records(text: &str, d: char, file: &str) -> Result<Vec<(u64, Vec<Field>)>> {
    let mut out = Vec::new();
    let mut fields: Vec<Field> = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    let mut line = 1u64;
    let mut start = 1u64;
    let mut end_record = |fields: &mut Vec<Field>, cur: &mut String, start: u64| {
        if cur.ends_with('\r') {
            cur.pop();
        }
        fields.push(field_of(&std::mem::take(cur)));
        let fields = std::mem::take(fields);
        if !(fields.len() == 1 && fields[0].0.is_empty()) {
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
            c if c == d && !quoted => fields.push(field_of(&std::mem::take(&mut cur))),
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
    /// Per column: some value was quoted, so the column is text.
    quoted: Vec<bool>,
}

impl Csv {
    fn parse(b: &[u8], file: &'static str) -> Result<Csv> {
        let text = decode_text(b);
        let text = text.strip_prefix('\u{feff}').unwrap_or(&text);
        let delim = sniff(text);
        let mut lines = records(text, delim, file)?.into_iter();
        let (_, mut header) = lines
            .next()
            .ok_or_else(|| Error::format(format!("{file} is empty")))?;
        while header.last().is_some_and(|(h, q)| h.is_empty() && !q) {
            header.pop();
        }
        let header: Vec<String> = header.into_iter().map(|(h, _)| h).collect();
        let mut quoted = vec![false; header.len()];
        let mut rows = Vec::new();
        for (n, mut f) in lines {
            while f.len() > header.len() && f.last().is_some_and(|(s, q)| s.is_empty() && !q) {
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
            for (q, (_, fq)) in quoted.iter_mut().zip(&f) {
                *q |= fq;
            }
            rows.push((n, f.into_iter().map(|(s, _)| s).collect()));
        }
        Ok(Csv {
            file,
            delim,
            header,
            rows,
            quoted,
        })
    }

    /// The column of the first name in `names` that the header holds.
    fn col(&self, names: &[&str]) -> Option<usize> {
        names
            .iter()
            .find_map(|n| self.header.iter().position(|h| h.eq_ignore_ascii_case(n)))
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
            Some(v) if v != NULL && v.is_finite() => Ok(v),
            Some(v) if v != NULL => Err(Error::format_at(
                line,
                format!("{} {} {v} is not finite", self.file, self.header[c]),
            )),
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

    /// Extra columns: Float64 when no value is quoted and every non-empty value
    /// is numeric, else Utf8. Empty fields and `-1.0E30` are null.
    fn extras(&self, skip: &[usize]) -> Result<Table> {
        let mut t = Table::new(self.rows.len());
        for (c, name) in self.header.iter().enumerate() {
            if skip.contains(&c) {
                continue;
            }
            let raw: Vec<&str> = self.rows.iter().map(|(_, r)| r[c].as_str()).collect();
            let nums: Option<Vec<Option<f64>>> = if self.quoted[c] {
                None
            } else {
                raw.iter()
                    .map(|s| {
                        if s.is_empty() {
                            Some(None)
                        } else {
                            self.num(s).map(|v| (v != NULL).then_some(v))
                        }
                    })
                    .collect()
            };
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
    read_pair(pt, tr, None)
}

/// As [`read`]; with `paths = Some((pt, tr))` a format error names the file of
/// the pair that holds it.
pub(crate) fn read_pair(pt: &[u8], tr: &[u8], paths: Option<(&Path, &Path)>) -> Result<Data> {
    let in_pt = |e: Error| match paths {
        Some((p, _)) => e.with_path(p),
        None => e,
    };
    let in_tr = |e: Error| match paths {
        Some((_, t)) => e.with_path(t),
        None => e,
    };
    let (vertices, index, vertex_table) = read_points(pt).map_err(in_pt)?;
    let (triangles, face_table) = read_triangles(tr, &index).map_err(in_tr)?;
    let mesh = Mesh::new(vertices, triangles)?.with_attributes(vertex_table, face_table)?;
    Ok(Data {
        mesh: Some(mesh),
        ..Data::default()
    })
}

type Points = (Vec<[f64; 3]>, HashMap<i64, u32>, Table);

fn read_points(pt: &[u8]) -> Result<Points> {
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
    Ok((vertices, index, p.extras(&cols)?))
}

fn read_triangles(tr: &[u8], index: &HashMap<i64, u32>) -> Result<(Vec<[u32; 3]>, Table)> {
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
    Ok((triangles, t.extras(&skip)?))
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
        ColumnData::Utf8(v) if valid => format!("\"{}\"", v[i].replace('"', "\"\"")),
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

/// `t` without the Int64 columns that a double cannot hold exactly, or an error
/// naming the first one unless `skip_unsupported_columns` is set.
fn storable(t: &Table, what: &str, opts: &WriteOptions) -> Result<Table> {
    let mut out = t.clone();
    for (name, col) in t.names().iter().zip(t.columns()) {
        let ColumnData::Int64(v) = &col.data else {
            continue;
        };
        let exact = (0..v.len()).all(|i| !col.is_valid(i) || v[i].unsigned_abs() <= 1 << 53);
        if exact {
            continue;
        }
        if !opts.skip_unsupported_columns {
            return Err(Error::Unsupported(format!(
                "datamine csv cannot store {what} column {name:?}: datamine numbers are doubles, and values beyond 2^53 do not fit a double exactly; set skip_unsupported_columns to drop it"
            )));
        }
        out = out.without_column(name);
    }
    Ok(out)
}

/// Writes the `(pt, tr)` pair. PIDs are 1-based in vertex order.
pub fn write(data: &Data, opts: &WriteOptions) -> Result<(Vec<u8>, Vec<u8>)> {
    if data.polylines.is_some() || data.points.is_some() || data.grid.is_some() {
        return Err(Error::Unsupported(
            "a datamine wireframe pair holds only a mesh".into(),
        ));
    }
    let mesh = data
        .mesh
        .as_ref()
        .ok_or_else(|| Error::Unsupported("a datamine wireframe pair needs a mesh".into()))?;
    let vt = &storable(mesh.vertex_attributes(), "vertex", opts)?;
    let ft = &storable(mesh.face_attributes(), "face", opts)?;
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
    reserved = vec!["TRIANGLE", "PID1", "PID2", "PID3"];
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

#[cfg(test)]
mod tests {
    use super::*;

    const PT: &str = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0,2\n0,1,0,3\n";

    #[test]
    fn triangle_id_alias_is_chosen_by_priority_not_column_order() {
        let t = Csv::parse(b"TRISOL,TRI,TRIANGLE,PID1\n1,2,3,4\n", "_tr.csv").unwrap();
        assert_eq!(t.col(&TRIANGLE_IDS), Some(2));
        let t = Csv::parse(b"TRISOL,PID1\n1,2\n", "_tr.csv").unwrap();
        assert_eq!(t.col(&TRIANGLE_IDS), Some(0));
    }

    #[test]
    fn a_trisol_column_beside_triangle_round_trips() {
        let tr = "TRISOL,TRIANGLE,PID1,PID2,PID3\n7,1,1,2,3\n";
        let data = read(PT.as_bytes(), tr.as_bytes()).unwrap();
        let (pt2, tr2) = write(&data, &WriteOptions::default()).unwrap();
        let again = read(&pt2, &tr2).unwrap().mesh.unwrap();
        let col = again.face_attributes().column("TRISOL").unwrap();
        assert_eq!(col.data, ColumnData::Float64(vec![7.0]));
    }

    const TR: &str = "TRIANGLE,PID1,PID2,PID3\n1,1,2,3\n";

    #[test]
    fn a_trailing_delimiter_adds_no_columns() {
        let pt = PT.replace('\n', ",\n");
        let tr = TR.replace('\n', ",\n");
        let m = read(pt.as_bytes(), tr.as_bytes()).unwrap().mesh.unwrap();
        assert_eq!(m.vertex_attributes().num_columns(), 0);
        assert_eq!(m.face_attributes().num_columns(), 0);
        let filled = "TRIANGLE,PID1,PID2,PID3,\n1,1,2,3,X\n";
        assert!(read(PT.as_bytes(), filled.as_bytes()).is_err());
    }

    #[test]
    fn quoted_text_that_looks_numeric_stays_text() {
        let m = Mesh::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap()
        .with_face_column("CODE", Column::utf8(vec!["007".into()]))
        .unwrap();
        let data = Data {
            mesh: Some(m),
            ..Data::default()
        };
        let (pt, tr) = write(&data, &WriteOptions::default()).unwrap();
        let back = read(&pt, &tr).unwrap().mesh.unwrap();
        let col = back.face_attributes().column("CODE").unwrap();
        assert_eq!(col.data, ColumnData::Utf8(vec!["007".into()]));
        let plain = "TRIANGLE,PID1,PID2,PID3,CODE\n1,1,2,3,007\n";
        let m = read(PT.as_bytes(), plain.as_bytes()).unwrap().mesh.unwrap();
        let col = m.face_attributes().column("CODE").unwrap();
        assert_eq!(col.data, ColumnData::Float64(vec![7.0]));
    }

    #[test]
    fn int64_beyond_2_pow_53_is_refused_or_skipped() {
        let tri = |v: i64| {
            let m = Mesh::new(
                vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                vec![[0, 1, 2]],
            )
            .unwrap()
            .with_face_column("ID", Column::i64(vec![v]))
            .unwrap();
            Data {
                mesh: Some(m),
                ..Data::default()
            }
        };
        assert!(write(&tri(1 << 53), &WriteOptions::default()).is_ok());
        let e = write(&tri((1 << 53) + 1), &WriteOptions::default()).unwrap_err();
        assert!(matches!(e, Error::Unsupported(_)), "{e}");
        let skip = WriteOptions {
            skip_unsupported_columns: true,
            ..WriteOptions::default()
        };
        let (pt, tr) = write(&tri((1 << 53) + 1), &skip).unwrap();
        let back = read(&pt, &tr).unwrap().mesh.unwrap();
        assert_eq!(back.face_attributes().num_columns(), 0);
    }

    #[test]
    fn a_non_finite_coordinate_names_its_line() {
        let pt = PT.replace("1,0,0,2", "1,nan,0,2");
        let e = read_pair(
            pt.as_bytes(),
            TR.as_bytes(),
            Some((Path::new("a_pt.csv"), Path::new("a_tr.csv"))),
        )
        .unwrap_err();
        assert_eq!(e.to_string(), "a_pt.csv:3: _pt.csv YP NaN is not finite");
    }

    #[test]
    fn each_error_names_the_file_of_the_pair_that_holds_it() {
        let (pt, tr) = (Path::new("a_pt.csv"), Path::new("a_tr.csv"));
        let bad_pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0\n";
        let good_tr = "TRIANGLE,PID1,PID2,PID3\n1,1,2,3\n";
        let e = read_pair(bad_pt.as_bytes(), good_tr.as_bytes(), Some((pt, tr))).unwrap_err();
        assert!(e.to_string().starts_with("a_pt.csv:3: "), "{e}");
        let bad_tr = "TRIANGLE,PID1,PID2,PID3\n1,1,2,9\n";
        let e = read_pair(PT.as_bytes(), bad_tr.as_bytes(), Some((pt, tr))).unwrap_err();
        assert!(e.to_string().starts_with("a_tr.csv:2: "), "{e}");
    }
}
