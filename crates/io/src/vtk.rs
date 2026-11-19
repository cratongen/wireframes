//! VTK legacy files, plus the cell, array and polygon model that `vtu` and `dxf` share.

use wireframes_core::{Column, ColumnData, Error, Mesh, PointSet, Polylines, Result, Table};

use crate::common::polygons;
use crate::{Data, UnstructuredGrid, WriteOptions};

pub(crate) struct Cells {
    /// `cells + 1` long, starting at 0.
    pub offsets: Vec<i64>,
    pub conn: Vec<i64>,
}

impl Cells {
    pub(crate) fn new() -> Cells {
        Cells {
            offsets: vec![0],
            conn: Vec::new(),
        }
    }

    pub(crate) fn from_ends(ends: Vec<i64>, conn: Vec<i64>) -> Cells {
        let mut offsets = Vec::with_capacity(ends.len() + 1);
        offsets.push(0);
        offsets.extend(ends);
        Cells { offsets, conn }
    }

    pub(crate) fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub(crate) fn push(&mut self, ids: impl IntoIterator<Item = i64>) {
        self.conn.extend(ids);
        self.offsets.push(self.conn.len() as i64);
    }

    /// Valid only after `check`.
    pub(crate) fn cell(&self, i: usize) -> &[i64] {
        &self.conn[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }

    pub(crate) fn check(&self, what: &str, n_points: usize) -> Result<()> {
        match offsets_problem(&self.offsets, &self.conn, n_points) {
            Some(m) => Err(Error::format(format!("{what}: {m}"))),
            None => Ok(()),
        }
    }
}

pub(crate) fn offsets_problem(offsets: &[i64], conn: &[i64], n_points: usize) -> Option<String> {
    let Some(&first) = offsets.first() else {
        return Some("offsets are empty".into());
    };
    if first != 0 {
        return Some("offsets must start at 0".into());
    }
    if offsets.windows(2).any(|w| w[1] < w[0]) {
        return Some("offsets must not decrease".into());
    }
    let last = offsets[offsets.len() - 1];
    if last != conn.len() as i64 {
        return Some(format!(
            "offsets end at {last} but connectivity has {} ids",
            conn.len()
        ));
    }
    if let Some(id) = conn.iter().find(|&&id| id < 0 || id >= n_points as i64) {
        return Some(format!("point id {id} out of range for {n_points} points"));
    }
    None
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Values {
    F64(Vec<f64>),
    I64(Vec<i64>),
    U8(Vec<u8>),
}

impl Values {
    pub(crate) fn len(&self) -> usize {
        match self {
            Values::F64(v) => v.len(),
            Values::I64(v) => v.len(),
            Values::U8(v) => v.len(),
        }
    }

    pub(crate) fn to_f64(&self) -> Vec<f64> {
        match self {
            Values::F64(v) => v.clone(),
            Values::I64(v) => v.iter().map(|&x| x as f64).collect(),
            Values::U8(v) => v.iter().map(|&x| x as f64).collect(),
        }
    }

    pub(crate) fn to_i64(&self, what: &str) -> Result<Vec<i64>> {
        match self {
            Values::I64(v) => Ok(v.clone()),
            Values::U8(v) => Ok(v.iter().map(|&x| x as i64).collect()),
            Values::F64(v) => v
                .iter()
                .map(|&x| {
                    if x.fract() == 0.0 && x.abs() < 9.0e15 {
                        Ok(x as i64)
                    } else {
                        Err(Error::format(format!("{what} holds a non-integer {x}")))
                    }
                })
                .collect(),
        }
    }

    fn component(&self, c: usize, n: usize) -> Values {
        match self {
            Values::F64(v) => Values::F64(v.iter().skip(c).step_by(n).copied().collect()),
            Values::I64(v) => Values::I64(v.iter().skip(c).step_by(n).copied().collect()),
            Values::U8(v) => Values::U8(v.iter().skip(c).step_by(n).copied().collect()),
        }
    }

    fn into_column(self) -> Column {
        match self {
            Values::F64(v) => Column::f64(v),
            Values::I64(v) => Column::i64(v),
            Values::U8(v) => Column::i64(v.into_iter().map(i64::from).collect()),
        }
    }
}

pub(crate) struct Array {
    pub name: String,
    pub ncomp: usize,
    pub values: Values,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    I64,
    U64,
    F32,
    F64,
}

macro_rules! decode_as {
    ($t:ty, $bytes:expr, $big:expr) => {
        $bytes.chunks_exact(std::mem::size_of::<$t>()).map(|c| {
            let a = c.try_into().expect("chunk has the type's size");
            if $big {
                <$t>::from_be_bytes(a)
            } else {
                <$t>::from_le_bytes(a)
            }
        })
    };
}

impl Ty {
    pub(crate) fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::I64 | Ty::U64 | Ty::F64 => 8,
        }
    }

    pub(crate) fn decode(self, b: &[u8], big: bool) -> Result<Values> {
        if !b.len().is_multiple_of(self.size()) {
            return Err(Error::format(format!(
                "binary array of {} bytes is not a whole number of {}-byte values",
                b.len(),
                self.size()
            )));
        }
        Ok(match self {
            Ty::I8 => Values::I64(decode_as!(i8, b, big).map(i64::from).collect()),
            Ty::U8 => Values::U8(b.to_vec()),
            Ty::I16 => Values::I64(decode_as!(i16, b, big).map(i64::from).collect()),
            Ty::U16 => Values::I64(decode_as!(u16, b, big).map(i64::from).collect()),
            Ty::I32 => Values::I64(decode_as!(i32, b, big).map(i64::from).collect()),
            Ty::U32 => Values::I64(decode_as!(u32, b, big).map(i64::from).collect()),
            Ty::I64 => Values::I64(decode_as!(i64, b, big).collect()),
            Ty::U64 => Values::I64(
                decode_as!(u64, b, big)
                    .map(|x| {
                        i64::try_from(x).map_err(|_| {
                            Error::format(format!("unsigned 64-bit value {x} exceeds i64"))
                        })
                    })
                    .collect::<Result<_>>()?,
            ),
            Ty::F32 => Values::F64(decode_as!(f32, b, big).map(f64::from).collect()),
            Ty::F64 => Values::F64(decode_as!(f64, b, big).collect()),
        })
    }

    pub(crate) fn parse_text(self, tokens: &[&str]) -> Result<Values> {
        let bad = |t: &str| Error::format(format!("cannot parse {t:?} as a number"));
        Ok(match self {
            Ty::F32 | Ty::F64 => Values::F64(
                tokens
                    .iter()
                    .map(|t| t.parse::<f64>().map_err(|_| bad(t)))
                    .collect::<Result<_>>()?,
            ),
            Ty::U8 => Values::U8(
                tokens
                    .iter()
                    .map(|t| t.parse::<u8>().map_err(|_| bad(t)))
                    .collect::<Result<_>>()?,
            ),
            Ty::U64 => Values::I64(
                tokens
                    .iter()
                    .map(|t| {
                        t.parse::<u64>()
                            .ok()
                            .and_then(|x| i64::try_from(x).ok())
                            .ok_or_else(|| bad(t))
                    })
                    .collect::<Result<_>>()?,
            ),
            _ => Values::I64(
                tokens
                    .iter()
                    .map(|t| t.parse::<i64>().map_err(|_| bad(t)))
                    .collect::<Result<_>>()?,
            ),
        })
    }
}

pub(crate) fn fmt_f64(x: f64) -> String {
    if x.is_nan() {
        "nan".into()
    } else if x.is_infinite() {
        if x > 0.0 { "inf".into() } else { "-inf".into() }
    } else {
        format!("{x:?}")
    }
}

/// Columns as arrays; a column with nulls gets a `<name>_valid` UInt8 companion.
pub(crate) fn flatten(table: &Table, opts: &WriteOptions, format: &str) -> Result<Vec<Array>> {
    let mut out = Vec::new();
    for (name, col) in table.names().iter().zip(table.columns()) {
        let companion = format!("{name}_valid");
        let is_companion = matches!(col.data, ColumnData::Bool(_))
            && name
                .strip_suffix("_valid")
                .is_some_and(|base| table.column(base).is_some());
        if is_companion || (col.null_count() > 0 && table.column(&companion).is_some()) {
            return Err(Error::Unsupported(format!(
                "{format} stores the nulls of a column x in a UInt8 column x_valid, so column {name:?} would be ambiguous; rename it"
            )));
        }
        let values = match &col.data {
            ColumnData::Float64(v) => Values::F64(v.clone()),
            ColumnData::Int64(v) => Values::I64(v.clone()),
            ColumnData::Bool(v) => Values::U8(v.iter().map(|&b| b as u8).collect()),
            ColumnData::Utf8(_) if opts.skip_unsupported_columns => continue,
            ColumnData::Utf8(_) => {
                return Err(Error::Unsupported(format!(
                    "{format} cannot store utf8 column {name:?}; drop it or set skip_unsupported_columns"
                )));
            }
        };
        out.push(Array {
            name: name.clone(),
            ncomp: 1,
            values,
        });
        if col.null_count() > 0 {
            let valid = (0..col.len()).map(|i| col.is_valid(i) as u8).collect();
            out.push(Array {
                name: companion,
                ncomp: 1,
                values: Values::U8(valid),
            });
        }
    }
    Ok(out)
}

/// Arrays as columns: components split into `<name>_0..`, UInt8 `<name>_valid`
/// companions applied as validity. A multi-component array over zero rows adds
/// no columns, so the column count stays bounded by the input size.
pub(crate) fn unflatten(arrays: Vec<Array>, rows: usize, what: &str) -> Result<Table> {
    let mut t = Table::new(rows);
    let mut companions = Vec::new();
    for a in arrays {
        if a.ncomp == 0 || rows.checked_mul(a.ncomp) != Some(a.values.len()) {
            return Err(Error::format(format!(
                "{what} array {:?} has {} values for {rows} rows of {} components",
                a.name,
                a.values.len(),
                a.ncomp
            )));
        }
        if a.ncomp == 1 {
            if let Values::U8(v) = &a.values
                && v.iter().all(|&f| f <= 1)
                && let Some(base) = a.name.strip_suffix("_valid")
            {
                companions.push((base.to_string(), v.iter().map(|&f| f == 1).collect()));
            }
            t = t.with_column(&a.name, a.values.into_column())?;
        } else if rows > 0 {
            for c in 0..a.ncomp {
                let col = a.values.component(c, a.ncomp).into_column();
                t = t.with_column(&format!("{}_{c}", a.name), col)?;
            }
        }
    }
    for (base, valid) in companions {
        let Some(col) = t.column(&base) else {
            continue;
        };
        let col = col.clone().with_valid(valid)?;
        t = t
            .with_column(&base, col)?
            .without_column(&format!("{base}_valid"));
    }
    Ok(t)
}

/// Drops the validity mask of columns without nulls and, when `drop_all_null`,
/// columns that are null on every row.
pub(crate) fn tidy(t: Table, drop_all_null: bool) -> Table {
    let rows = t.rows();
    let mut out = Table::new(rows);
    for (name, col) in t.names().iter().zip(t.columns()) {
        let nulls = col.null_count();
        if drop_all_null && nulls == rows {
            continue;
        }
        let mut col = col.clone();
        if nulls == 0 {
            col.valid = None;
        }
        out = out.with_column(name, col).expect("same row count");
    }
    out
}

pub(crate) fn concat_tables(tables: &[Table]) -> Result<Table> {
    let mut it = tables.iter();
    let Some(first) = it.next() else {
        return Ok(Table::new(0));
    };
    it.try_fold(first.clone(), |a, b| Table::concat(&a, b))
}

/// Triangles fanned from polygons, with the source cell of each triangle and a
/// `polygon` id column when some polygon had more than three corners.
#[derive(Default)]
pub(crate) struct Fan {
    pub tris: Vec<[u32; 3]>,
    pub source: Vec<u32>,
    polygon: Vec<i64>,
    next: i64,
    wide: bool,
}

impl Fan {
    pub(crate) fn polygon(&mut self, ids: &[u32], source: u32) {
        for i in 1..ids.len().saturating_sub(1) {
            self.tris.push([ids[0], ids[i], ids[i + 1]]);
            self.source.push(source);
            self.polygon.push(self.next);
        }
        self.wide |= ids.len() > 3;
        self.next += 1;
    }

    pub(crate) fn strip(&mut self, ids: &[u32], source: u32) {
        for i in 0..ids.len().saturating_sub(2) {
            let t = if i % 2 == 0 {
                [ids[i], ids[i + 1], ids[i + 2]]
            } else {
                [ids[i + 1], ids[i], ids[i + 2]]
            };
            self.tris.push(t);
            self.source.push(source);
            self.polygon.push(self.next);
            self.next += 1;
        }
    }

    /// `face_table` has one row per triangle.
    pub(crate) fn mesh(
        self,
        vertices: Vec<[f64; 3]>,
        vertex_table: Table,
        face_table: Table,
    ) -> Result<Mesh> {
        let mut face = face_table;
        if self.wide {
            face = face.with_column("polygon", Column::i64(self.polygon))?;
        }
        Mesh::new(vertices, self.tris)?.with_attributes(vertex_table, face)
    }
}

pub(crate) struct PolyData {
    pub points: Vec<[f64; 3]>,
    pub point_data: Table,
    pub verts: Cells,
    pub lines: Cells,
    pub polys: Cells,
    pub strips: Cells,
    /// Rows in cell order: verts, lines, polys, strips.
    pub cell_data: Table,
}

impl PolyData {
    pub(crate) fn sections(&self) -> [&Cells; 4] {
        [&self.verts, &self.lines, &self.polys, &self.strips]
    }
}

const FEATURE: &str = "wireframes_feature";

/// Splits poly data into points (verts), polylines (lines) and a mesh (polys
/// and strips). The mesh keeps the points that faces use or that no cell uses.
/// Columns null on every row of a section are dropped when the file mixes kinds.
pub(crate) fn poly_to_data(pd: PolyData) -> Result<Data> {
    let n = pd.points.len();
    for (c, what) in pd
        .sections()
        .iter()
        .zip(["verts", "lines", "polys", "strips"])
    {
        c.check(what, n)?;
    }
    let [nv, nl, np, ns] = pd.sections().map(|c| c.len());
    if pd.cell_data.rows() != nv + nl + np + ns {
        return Err(Error::format(format!(
            "cell data has {} rows for {} cells",
            pd.cell_data.rows(),
            nv + nl + np + ns
        )));
    }
    if pd.point_data.rows() != n {
        return Err(Error::format(format!(
            "point data has {} rows for {n} points",
            pd.point_data.rows()
        )));
    }
    let multi = [nv > 0, nl > 0, np + ns > 0].iter().filter(|b| **b).count() > 1;
    let mut data = Data::default();
    if nv > 0 {
        let mut ids = Vec::new();
        let mut cell = Vec::new();
        for c in 0..nv {
            for &id in pd.verts.cell(c) {
                ids.push(id as u32);
                cell.push(c as u32);
            }
        }
        let coords = ids.iter().map(|&i| pd.points[i as usize]).collect();
        let mut attrs = tidy(pd.point_data.take(&ids), multi);
        let cells = tidy(pd.cell_data.take(&cell), multi);
        for (name, col) in cells.names().iter().zip(cells.columns()) {
            let name = match attrs.column(name) {
                Some(_) => format!("{name}_cell"),
                None => name.clone(),
            };
            attrs = attrs.with_column(&name, col.clone())?;
        }
        data.points = Some(PointSet::new(coords)?.with_attributes(attrs)?);
    }
    if nl > 0 {
        let rows: Vec<u32> = (nv as u32..(nv + nl) as u32).collect();
        let mut table = tidy(pd.cell_data.take(&rows), multi);
        let feature_of_cell: Vec<i64> = match table.column(FEATURE).map(|c| &c.data) {
            Some(ColumnData::Int64(v)) => {
                let f = v.clone();
                table = table.without_column(FEATURE);
                f
            }
            _ => (0..nl as i64).collect(),
        };
        let mut vertices = Vec::new();
        let mut vids: Vec<u32> = Vec::new();
        let mut parts = vec![0u32];
        let mut features = vec![0u32];
        let mut closed = Vec::new();
        let mut firsts = Vec::new();
        for c in 0..nl {
            let f = feature_of_cell[c];
            let prev = c.checked_sub(1).map(|p| feature_of_cell[p]);
            match prev {
                None if f == 0 => firsts.push(c as u32),
                Some(q) if f == q => {}
                Some(q) if Some(f) == q.checked_add(1) => {
                    features.push(parts.len() as u32 - 1);
                    firsts.push(c as u32);
                }
                _ => {
                    return Err(Error::format(format!(
                        "{FEATURE} must start at 0 and increase by at most 1 per line cell"
                    )));
                }
            }
            let ids = pd.lines.cell(c);
            if ids.len() == 1 {
                return Err(Error::format(format!(
                    "line cell {c} has 1 point, needs 0 (an empty feature) or at least 2"
                )));
            }
            if ids.is_empty() {
                continue;
            }
            let is_closed = ids.len() >= 4 && ids[0] == ids[ids.len() - 1];
            let used = if is_closed {
                &ids[..ids.len() - 1]
            } else {
                ids
            };
            vertices.extend(used.iter().map(|&i| pd.points[i as usize]));
            vids.extend(used.iter().map(|&i| i as u32));
            parts.push(vertices.len() as u32);
            closed.push(is_closed);
        }
        features.push(parts.len() as u32 - 1);
        let pl = Polylines::new(vertices, parts, features, closed, table.take(&firsts))
            .map_err(|e| Error::format(e.to_string()))?
            .with_vertex_attributes(tidy(pd.point_data.take(&vids), multi))?;
        data.polylines = Some(pl);
    }
    if np + ns > 0 || nv + nl == 0 {
        let keep: Vec<u32> = if nv + nl == 0 {
            (0..n as u32).collect()
        } else {
            let mut by_line = vec![false; n];
            let mut by_face = vec![false; n];
            for &id in pd.verts.conn.iter().chain(&pd.lines.conn) {
                by_line[id as usize] = true;
            }
            for &id in pd.polys.conn.iter().chain(&pd.strips.conn) {
                by_face[id as usize] = true;
            }
            (0..n)
                .filter(|&i| by_face[i] || !by_line[i])
                .map(|i| i as u32)
                .collect()
        };
        let mut remap = vec![u32::MAX; n];
        for (new, &old) in keep.iter().enumerate() {
            remap[old as usize] = new as u32;
        }
        let mut fan = Fan::default();
        let base = (nv + nl) as u32;
        for (sec, cells, start) in [
            ("polygon", &pd.polys, base),
            ("strip", &pd.strips, base + np as u32),
        ] {
            for c in 0..cells.len() {
                let ids: Vec<u32> = cells.cell(c).iter().map(|&i| remap[i as usize]).collect();
                if ids.len() < 3 {
                    return Err(Error::format(format!(
                        "{sec} cell {c} has {} points, needs at least 3",
                        ids.len()
                    )));
                }
                if sec == "polygon" {
                    fan.polygon(&ids, start + c as u32);
                } else {
                    fan.strip(&ids, start + c as u32);
                }
            }
        }
        let vertices = keep.iter().map(|&i| pd.points[i as usize]).collect();
        let vt = tidy(pd.point_data.take(&keep), multi);
        let ft = tidy(pd.cell_data.take(&fan.source), multi);
        data.mesh = Some(fan.mesh(vertices, vt, ft)?);
    }
    Ok(data)
}

/// Mesh vertices, then polyline vertices, then points; cells as verts (one per
/// point), lines (closed parts repeat their first id) and polygons. Features
/// with several parts carry a `wireframes_feature` line column.
pub(crate) fn data_to_poly(data: &Data) -> Result<PolyData> {
    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut vtables = Vec::new();
    let mut verts = Cells::new();
    let mut lines = Cells::new();
    let mut polys = Cells::new();
    if let Some(m) = &data.mesh {
        points.extend_from_slice(m.vertices());
        vtables.push(m.vertex_attributes().clone());
    }
    let pl_base = points.len() as i64;
    if let Some(pl) = &data.polylines {
        points.extend_from_slice(pl.vertices());
        vtables.push(pl.vertex_attributes().clone());
    }
    if let Some(ps) = &data.points {
        let base = points.len() as i64;
        points.extend_from_slice(ps.coords());
        vtables.push(ps.attributes().clone());
        for i in 0..ps.len() {
            verts.push([base + i as i64]);
        }
    }
    let mut ctables = vec![Table::new(verts.len())];
    if let Some(pl) = &data.polylines {
        let mesh_faces = data.mesh.as_ref().map(Mesh::face_attributes);
        if [Some(pl.attributes()), mesh_faces]
            .iter()
            .flatten()
            .any(|t| t.column(FEATURE).is_some())
        {
            return Err(Error::Unsupported(format!(
                "column {FEATURE:?} is reserved for the feature of each polyline part; rename it"
            )));
        }
        let mut feature = Vec::new();
        let mut rows = Vec::new();
        let mut multipart = false;
        for f in 0..pl.len() {
            let parts = pl.feature_parts(f);
            multipart |= parts.len() > 1;
            if parts.is_empty() {
                lines.push(std::iter::empty());
                feature.push(f as i64);
                rows.push(f as u32);
            }
            for p in parts {
                let mut ids: Vec<i64> = pl.part_range(p).map(|i| pl_base + i as i64).collect();
                if pl.closed()[p] {
                    ids.push(ids[0]);
                }
                lines.push(ids);
                feature.push(f as i64);
                rows.push(f as u32);
            }
        }
        let mut t = pl.attributes().take(&rows);
        if multipart {
            t = t.with_column(FEATURE, Column::i64(feature))?;
        }
        ctables.push(t);
    }
    if let Some(m) = &data.mesh {
        let (ps, t) = polygons(m)?;
        for p in ps {
            polys.push(p.into_iter().map(i64::from));
        }
        ctables.push(t);
    }
    Ok(PolyData {
        point_data: concat_tables(&vtables)?,
        cell_data: concat_tables(&ctables)?,
        points,
        verts,
        lines,
        polys,
        strips: Cells::new(),
    })
}

pub(crate) fn grid_problem(g: &UnstructuredGrid) -> Option<String> {
    let n = g.cell_types.len();
    if g.offsets.len() != n + 1 {
        return Some(format!(
            "grid offsets have {} values for {n} cells; expected {}",
            g.offsets.len(),
            n + 1
        ));
    }
    if let Some(m) = offsets_problem(&g.offsets, &g.connectivity, g.points.len()) {
        return Some(format!("grid cells: {m}"));
    }
    if g.point_data.rows() != g.points.len() {
        return Some(format!(
            "grid point data has {} rows for {} points",
            g.point_data.rows(),
            g.points.len()
        ));
    }
    if g.cell_data.rows() != n {
        return Some(format!(
            "grid cell data has {} rows for {n} cells",
            g.cell_data.rows()
        ));
    }
    if let Some((fo, _)) = &g.faces
        && fo.len() != n
    {
        return Some(format!(
            "grid face offsets have {} values for {n} cells",
            fo.len()
        ));
    }
    None
}

/// A grid of only triangle, polygon and quad cells becomes a mesh; anything
/// else stays a grid.
pub(crate) fn grid_to_data(g: UnstructuredGrid) -> Result<Data> {
    if let Some(m) = grid_problem(&g) {
        return Err(Error::format(m));
    }
    let surface = g.faces.is_none() && g.cell_types.iter().all(|t| matches!(t, 5 | 7 | 9));
    if !surface {
        let grid = UnstructuredGrid {
            point_data: tidy(g.point_data, false),
            cell_data: tidy(g.cell_data, false),
            ..g
        };
        return Ok(Data {
            grid: Some(grid),
            ..Data::default()
        });
    }
    let mut fan = Fan::default();
    for c in 0..g.cell_types.len() {
        let ids: Vec<u32> = g.connectivity[g.offsets[c] as usize..g.offsets[c + 1] as usize]
            .iter()
            .map(|&i| i as u32)
            .collect();
        if ids.len() < 3 {
            return Err(Error::format(format!(
                "cell {c} has {} points, needs at least 3",
                ids.len()
            )));
        }
        fan.polygon(&ids, c as u32);
    }
    let ft = tidy(g.cell_data.take(&fan.source), false);
    let mesh = fan.mesh(g.points, tidy(g.point_data, false), ft)?;
    Ok(Data {
        mesh: Some(mesh),
        ..Data::default()
    })
}

/// Triangles as VTK_TRIANGLE (5), rebuilt polygons as VTK_POLYGON (7).
pub(crate) fn mesh_to_grid(m: &Mesh) -> Result<UnstructuredGrid> {
    let (ps, cell_data) = polygons(m)?;
    let mut cells = Cells::new();
    let mut types = Vec::with_capacity(ps.len());
    for p in ps {
        types.push(if p.len() == 3 { 5 } else { 7 });
        cells.push(p.into_iter().map(i64::from));
    }
    Ok(UnstructuredGrid {
        points: m.vertices().to_vec(),
        cell_types: types,
        offsets: cells.offsets,
        connectivity: cells.conn,
        faces: None,
        point_data: m.vertex_attributes().clone(),
        cell_data,
    })
}

pub(crate) fn to_points(v: Vec<f64>) -> Vec<[f64; 3]> {
    v.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

struct Legacy<'a> {
    b: &'a [u8],
    pos: usize,
    line: u64,
    binary: bool,
}

impl<'a> Legacy<'a> {
    fn err(&self, m: impl Into<String>) -> Error {
        Error::format_at(self.line, m)
    }

    fn raw_line(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.b.len() {
            return None;
        }
        let start = self.pos;
        while self.pos < self.b.len() && self.b[self.pos] != b'\n' {
            self.pos += 1;
        }
        let mut l = &self.b[start..self.pos];
        if self.pos < self.b.len() {
            self.pos += 1;
        }
        self.line += 1;
        if l.last() == Some(&b'\r') {
            l = &l[..l.len() - 1];
        }
        Some(l)
    }

    fn line_tokens(&mut self) -> Option<Vec<String>> {
        loop {
            let l = self.raw_line()?;
            let toks: Vec<String> = String::from_utf8_lossy(l)
                .split_whitespace()
                .map(String::from)
                .collect();
            if !toks.is_empty() {
                return Some(toks);
            }
        }
    }

    fn peek_tokens(&mut self) -> Option<Vec<String>> {
        let (pos, line) = (self.pos, self.line);
        let t = self.line_tokens();
        self.pos = pos;
        self.line = line;
        t
    }

    fn need_line(&mut self, what: &str) -> Result<Vec<String>> {
        self.line_tokens()
            .ok_or_else(|| self.err(format!("file ends before {what}")))
    }

    fn token(&mut self) -> Option<&'a [u8]> {
        while self.pos < self.b.len() && self.b[self.pos].is_ascii_whitespace() {
            if self.b[self.pos] == b'\n' {
                self.line += 1;
            }
            self.pos += 1;
        }
        let start = self.pos;
        while self.pos < self.b.len() && !self.b[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        (self.pos > start).then(|| &self.b[start..self.pos])
    }

    fn values(&mut self, n: usize, ty: Ty) -> Result<Values> {
        if self.binary {
            let size = n
                .checked_mul(ty.size())
                .filter(|&s| s <= self.b.len() - self.pos)
                .ok_or_else(|| self.err(format!("binary data for {n} values ends early")))?;
            let v = ty.decode(&self.b[self.pos..self.pos + size], true)?;
            self.pos += size;
            return Ok(v);
        }
        let mut toks = Vec::with_capacity(n.min(self.b.len() / 2 + 1));
        for _ in 0..n {
            let t = self
                .token()
                .ok_or_else(|| self.err(format!("file ends before {n} values")))?;
            toks.push(std::str::from_utf8(t).map_err(|_| self.err("non-text value"))?);
        }
        ty.parse_text(&toks).map_err(|e| self.err(e.to_string()))
    }

    fn cells(&mut self, a: usize, b: usize) -> Result<Cells> {
        if self
            .peek_tokens()
            .is_some_and(|t| t[0].eq_ignore_ascii_case("OFFSETS"))
        {
            let t = self.need_line("OFFSETS")?;
            let ty = legacy_ty(&t, 1, self.line)?;
            let offsets = self.values(a, ty)?.to_i64("offsets")?;
            let t = self.need_line("CONNECTIVITY")?;
            if !t[0].eq_ignore_ascii_case("CONNECTIVITY") {
                return Err(self.err(format!("expected CONNECTIVITY, found {:?}", t[0])));
            }
            let ty = legacy_ty(&t, 1, self.line)?;
            let conn = self.values(b, ty)?.to_i64("connectivity")?;
            let offsets = if offsets.is_empty() { vec![0] } else { offsets };
            return Ok(Cells { offsets, conn });
        }
        let v = self.values(b, Ty::I32)?.to_i64("cells")?;
        let mut cells = Cells::new();
        let mut i = 0usize;
        for k in 0..a {
            let short = || Error::format(format!("cell list ends inside cell {k}"));
            let n = usize::try_from(*v.get(i).ok_or_else(short)?).map_err(|_| short())?;
            let end = i
                .checked_add(1)
                .and_then(|x| x.checked_add(n))
                .filter(|&e| e <= v.len())
                .ok_or_else(short)?;
            cells.push(v[i + 1..end].iter().copied());
            i = end;
        }
        if i != v.len() {
            return Err(self.err(format!(
                "cell list size {} does not match its {a} cells ({i} values)",
                v.len()
            )));
        }
        Ok(cells)
    }

    fn skip_metadata(&mut self) {
        while let Some(l) = self.raw_line() {
            if l.iter().all(|c| c.is_ascii_whitespace()) {
                break;
            }
        }
    }
}

fn legacy_ty(t: &[String], i: usize, line: u64) -> Result<Ty> {
    let name = t
        .get(i)
        .ok_or_else(|| Error::format_at(line, "missing data type"))?;
    Ok(match name.to_ascii_lowercase().as_str() {
        "char" | "vtktypeint8" => Ty::I8,
        "unsigned_char" | "vtktypeuint8" => Ty::U8,
        "short" | "vtktypeint16" => Ty::I16,
        "unsigned_short" | "vtktypeuint16" => Ty::U16,
        "int" | "vtktypeint32" => Ty::I32,
        "unsigned_int" | "vtktypeuint32" => Ty::U32,
        "long" | "vtktypeint64" | "vtkidtype" => Ty::I64,
        "unsigned_long" | "vtktypeuint64" => Ty::U64,
        "float" | "vtktypefloat32" => Ty::F32,
        "double" | "vtktypefloat64" => Ty::F64,
        "bit" => {
            return Err(Error::Unsupported(
                "vtk bit arrays are not supported".into(),
            ));
        }
        other => {
            return Err(Error::format_at(
                line,
                format!("unknown data type {other:?}"),
            ));
        }
    })
}

fn count(t: &[String], i: usize, line: u64) -> Result<usize> {
    t.get(i)
        .and_then(|s| s.parse::<usize>().ok())
        .ok_or_else(|| Error::format_at(line, format!("{} needs a count", t[0])))
}

fn decode_name(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%'
            && i + 3 <= b.len()
            && let Some(v) = std::str::from_utf8(&b[i + 1..i + 3])
                .ok()
                .and_then(|h| u8::from_str_radix(h, 16).ok())
        {
            out.push(v);
            i += 3;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn encode_name(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        if c.is_whitespace() || c == '%' || c.is_control() {
            let mut buf = [0u8; 4];
            for b in c.encode_utf8(&mut buf).bytes() {
                out.push_str(&format!("%{b:02X}"));
            }
        } else {
            out.push(c);
        }
    }
    out
}

pub fn read(r: &[u8]) -> Result<Data> {
    let mut c = Legacy {
        b: r,
        pos: 0,
        line: 0,
        binary: false,
    };
    let head = c.raw_line().unwrap_or_default();
    if !String::from_utf8_lossy(head)
        .to_ascii_lowercase()
        .starts_with("# vtk datafile version")
    {
        return Err(Error::format_at(1, "not a vtk legacy file"));
    }
    c.raw_line();
    let mode = c.need_line("ASCII or BINARY")?;
    c.binary = match mode[0].to_ascii_uppercase().as_str() {
        "ASCII" => false,
        "BINARY" => true,
        other => return Err(c.err(format!("expected ASCII or BINARY, found {other:?}"))),
    };
    let mut kind: Option<String> = None;
    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut sections: [Option<Cells>; 4] = [None, None, None, None];
    let mut cells: Option<Cells> = None;
    let mut types: Option<Vec<i64>> = None;
    let mut point_arrays: Vec<Array> = Vec::new();
    let mut cell_arrays: Vec<Array> = Vec::new();
    let mut point_rows: Option<usize> = None;
    let mut cell_rows: Option<usize> = None;
    let mut target: Option<bool> = None;
    while let Some(t) = c.line_tokens() {
        let line = c.line;
        let key = t[0].to_ascii_uppercase();
        let rows = match target {
            Some(true) => point_rows.unwrap_or(0),
            Some(false) => cell_rows.unwrap_or(0),
            None => 0,
        };
        let mut push = |a: Array, c: &Legacy| -> Result<()> {
            match target {
                Some(true) => point_arrays.push(a),
                Some(false) => cell_arrays.push(a),
                None => return Err(c.err("attribute data before POINT_DATA or CELL_DATA")),
            }
            Ok(())
        };
        match key.as_str() {
            "DATASET" => {
                let k = t.get(1).ok_or_else(|| c.err("DATASET needs a type"))?;
                kind = Some(k.to_ascii_uppercase());
            }
            "POINTS" => {
                let n = count(&t, 1, line)?;
                let ty = legacy_ty(&t, 2, line)?;
                let n3 = n.checked_mul(3).ok_or_else(|| c.err("too many points"))?;
                points = to_points(c.values(n3, ty)?.to_f64());
            }
            "VERTICES" | "LINES" | "POLYGONS" | "TRIANGLE_STRIPS" | "CELLS" => {
                let a = count(&t, 1, line)?;
                let b = count(&t, 2, line)?;
                let cs = c.cells(a, b)?;
                match key.as_str() {
                    "VERTICES" => sections[0] = Some(cs),
                    "LINES" => sections[1] = Some(cs),
                    "POLYGONS" => sections[2] = Some(cs),
                    "TRIANGLE_STRIPS" => sections[3] = Some(cs),
                    _ => cells = Some(cs),
                }
            }
            "CELL_TYPES" => {
                let n = count(&t, 1, line)?;
                types = Some(c.values(n, Ty::I32)?.to_i64("cell types")?);
            }
            "POINT_DATA" | "CELL_DATA" => {
                let n = count(&t, 1, line)?;
                let (slot, is_point) = if key == "POINT_DATA" {
                    (&mut point_rows, true)
                } else {
                    (&mut cell_rows, false)
                };
                if slot.is_some_and(|r| r != n) {
                    return Err(c.err(format!("{key} {n} repeats with a different count")));
                }
                *slot = Some(n);
                target = Some(is_point);
            }
            "SCALARS" => {
                let name = decode_name(t.get(1).ok_or_else(|| c.err("SCALARS needs a name"))?);
                let ty = legacy_ty(&t, 2, line)?;
                let ncomp = match t.get(3) {
                    Some(s) => s
                        .parse::<usize>()
                        .map_err(|_| c.err("bad component count"))?,
                    None => 1,
                };
                if c.peek_tokens()
                    .is_some_and(|p| p[0].eq_ignore_ascii_case("LOOKUP_TABLE"))
                {
                    c.line_tokens();
                }
                let n = rows
                    .checked_mul(ncomp)
                    .ok_or_else(|| c.err("array too large"))?;
                let values = c.values(n, ty)?;
                push(
                    Array {
                        name,
                        ncomp,
                        values,
                    },
                    &c,
                )?;
            }
            "COLOR_SCALARS" => {
                let name = decode_name(
                    t.get(1)
                        .ok_or_else(|| c.err("COLOR_SCALARS needs a name"))?,
                );
                let ncomp = count(&t, 2, line)?;
                let n = rows
                    .checked_mul(ncomp)
                    .ok_or_else(|| c.err("array too large"))?;
                let values = if c.binary {
                    Values::F64(
                        c.values(n, Ty::U8)?
                            .to_f64()
                            .iter()
                            .map(|v| v / 255.0)
                            .collect(),
                    )
                } else {
                    c.values(n, Ty::F64)?
                };
                push(
                    Array {
                        name,
                        ncomp,
                        values,
                    },
                    &c,
                )?;
            }
            "LOOKUP_TABLE" => {
                let size = count(&t, 2, line)?;
                let n = size
                    .checked_mul(4)
                    .ok_or_else(|| c.err("table too large"))?;
                c.values(n, if c.binary { Ty::U8 } else { Ty::F64 })?;
            }
            "VECTORS" | "NORMALS" | "TENSORS" | "TENSORS6" | "TEXTURE_COORDINATES" => {
                let name = decode_name(
                    t.get(1)
                        .ok_or_else(|| c.err(format!("{key} needs a name")))?,
                );
                let (ncomp, ty) = match key.as_str() {
                    "VECTORS" | "NORMALS" => (3, legacy_ty(&t, 2, line)?),
                    "TENSORS" => (9, legacy_ty(&t, 2, line)?),
                    "TENSORS6" => (6, legacy_ty(&t, 2, line)?),
                    _ => (count(&t, 2, line)?, legacy_ty(&t, 3, line)?),
                };
                let n = rows
                    .checked_mul(ncomp)
                    .ok_or_else(|| c.err("array too large"))?;
                let values = c.values(n, ty)?;
                push(
                    Array {
                        name,
                        ncomp,
                        values,
                    },
                    &c,
                )?;
            }
            "FIELD" => {
                let k = count(&t, 2, line)?;
                for _ in 0..k {
                    let mut a = c.need_line("a field array")?;
                    if a[0].eq_ignore_ascii_case("METADATA") {
                        c.skip_metadata();
                        a = c.need_line("a field array")?;
                    }
                    let line = c.line;
                    let ncomp = count(&a, 1, line)?;
                    let ntuples = count(&a, 2, line)?;
                    let ty = legacy_ty(&a, 3, line)?;
                    let n = ncomp
                        .checked_mul(ntuples)
                        .ok_or_else(|| c.err("array too large"))?;
                    let values = c.values(n, ty)?;
                    if target.is_some() {
                        if ntuples != rows {
                            return Err(c.err(format!(
                                "field array {:?} has {ntuples} tuples for {rows} rows",
                                a[0]
                            )));
                        }
                        push(
                            Array {
                                name: decode_name(&a[0]),
                                ncomp,
                                values,
                            },
                            &c,
                        )?;
                    }
                }
            }
            "METADATA" => c.skip_metadata(),
            other => return Err(c.err(format!("unknown keyword {other:?}"))),
        }
    }
    let n_points = points.len();
    if point_rows.is_some_and(|r| r != n_points) {
        return Err(Error::format(format!(
            "POINT_DATA has {} rows for {n_points} points",
            point_rows.unwrap_or(0)
        )));
    }
    let point_data = unflatten(point_arrays, n_points, "point")?;
    match kind.as_deref() {
        Some("POLYDATA") => {
            let [v, l, p, s] = sections.map(|c| c.unwrap_or_else(Cells::new));
            let n_cells = v.len() + l.len() + p.len() + s.len();
            if cell_rows.is_some_and(|r| r != n_cells) {
                return Err(Error::format(format!(
                    "CELL_DATA has {} rows for {n_cells} cells",
                    cell_rows.unwrap_or(0)
                )));
            }
            poly_to_data(PolyData {
                points,
                point_data,
                verts: v,
                lines: l,
                polys: p,
                strips: s,
                cell_data: unflatten(cell_arrays, n_cells, "cell")?,
            })
        }
        Some("UNSTRUCTURED_GRID") => {
            let cells = cells.unwrap_or_else(Cells::new);
            let types = types.unwrap_or_default();
            if types.len() != cells.len() {
                return Err(Error::format(format!(
                    "CELL_TYPES has {} values for {} cells",
                    types.len(),
                    cells.len()
                )));
            }
            let cell_types = types
                .iter()
                .map(|&t| u8::try_from(t).map_err(|_| Error::format(format!("bad cell type {t}"))))
                .collect::<Result<Vec<u8>>>()?;
            if cell_rows.is_some_and(|r| r != cells.len()) {
                return Err(Error::format(format!(
                    "CELL_DATA has {} rows for {} cells",
                    cell_rows.unwrap_or(0),
                    cells.len()
                )));
            }
            let cell_data = unflatten(cell_arrays, cells.len(), "cell")?;
            grid_to_data(UnstructuredGrid {
                points,
                cell_types,
                offsets: cells.offsets,
                connectivity: cells.conn,
                faces: None,
                point_data,
                cell_data,
            })
        }
        Some(other) => Err(Error::Unsupported(format!(
            "vtk dataset {other} is not supported; use POLYDATA or UNSTRUCTURED_GRID"
        ))),
        None => Err(Error::format("vtk file has no DATASET")),
    }
}

struct LegacyWriter {
    out: Vec<u8>,
    binary: bool,
}

impl LegacyWriter {
    fn text(&mut self, s: &str) {
        self.out.extend_from_slice(s.as_bytes());
    }

    fn values(&mut self, v: &Values, int32: bool) {
        if self.binary {
            match v {
                Values::F64(v) => v.iter().for_each(|x| self.out.extend(x.to_be_bytes())),
                Values::I64(v) if int32 => v
                    .iter()
                    .for_each(|&x| self.out.extend((x as i32).to_be_bytes())),
                Values::I64(v) => v.iter().for_each(|x| self.out.extend(x.to_be_bytes())),
                Values::U8(v) => self.out.extend_from_slice(v),
            }
            self.out.push(b'\n');
            return;
        }
        let toks: Vec<String> = match v {
            Values::F64(v) => v.iter().map(|&x| fmt_f64(x)).collect(),
            Values::I64(v) => v.iter().map(|x| x.to_string()).collect(),
            Values::U8(v) => v.iter().map(|x| x.to_string()).collect(),
        };
        for chunk in toks.chunks(9) {
            self.text(&chunk.join(" "));
            self.out.push(b'\n');
        }
    }

    fn cells(&mut self, key: &str, offsets: &[i64], conn: &[i64]) {
        self.text(&format!(
            "{key} {} {}\nOFFSETS vtktypeint64\n",
            offsets.len(),
            conn.len()
        ));
        self.values(&Values::I64(offsets.to_vec()), false);
        self.text("CONNECTIVITY vtktypeint64\n");
        self.values(&Values::I64(conn.to_vec()), false);
    }

    fn arrays(&mut self, key: &str, rows: usize, arrays: &[Array]) {
        if arrays.is_empty() {
            return;
        }
        self.text(&format!("{key} {rows}\nFIELD FieldData {}\n", arrays.len()));
        for a in arrays {
            let ty = match a.values {
                Values::F64(_) => "double",
                Values::I64(_) => "vtktypeint64",
                Values::U8(_) => "unsigned_char",
            };
            self.text(&format!(
                "{} {} {rows} {ty}\n",
                encode_name(&a.name),
                a.ncomp
            ));
            self.values(&a.values, false);
        }
    }

    fn points(&mut self, p: &[[f64; 3]]) {
        self.text(&format!("POINTS {} double\n", p.len()));
        self.values(&Values::F64(p.iter().flatten().copied().collect()), false);
    }
}

pub fn write(data: &Data, opts: &WriteOptions) -> Result<Vec<u8>> {
    let mut w = LegacyWriter {
        out: Vec::new(),
        binary: !opts.ascii,
    };
    w.text("# vtk DataFile Version 5.1\nwireframes\n");
    w.text(if w.binary { "BINARY\n" } else { "ASCII\n" });
    if let Some(g) = &data.grid {
        if data.mesh.is_some() || data.polylines.is_some() || data.points.is_some() {
            return Err(Error::Unsupported(
                "a vtk file holds either a grid or surface data, not both".into(),
            ));
        }
        if let Some(m) = grid_problem(g) {
            return Err(Error::invalid(m));
        }
        if g.faces.is_some() {
            return Err(Error::Unsupported(
                "legacy vtk cannot store polyhedron faces; write .vtu".into(),
            ));
        }
        let pa = flatten(&g.point_data, opts, "vtk")?;
        let ca = flatten(&g.cell_data, opts, "vtk")?;
        w.text("DATASET UNSTRUCTURED_GRID\n");
        w.points(&g.points);
        w.cells("CELLS", &g.offsets, &g.connectivity);
        w.text(&format!("CELL_TYPES {}\n", g.cell_types.len()));
        w.values(
            &Values::I64(g.cell_types.iter().map(|&t| t as i64).collect()),
            true,
        );
        w.arrays("POINT_DATA", g.points.len(), &pa);
        w.arrays("CELL_DATA", g.cell_types.len(), &ca);
        return Ok(w.out);
    }
    let pd = data_to_poly(data)?;
    let pa = flatten(&pd.point_data, opts, "vtk")?;
    let ca = flatten(&pd.cell_data, opts, "vtk")?;
    w.text("DATASET POLYDATA\n");
    w.points(&pd.points);
    for (key, c) in [
        ("VERTICES", &pd.verts),
        ("LINES", &pd.lines),
        ("POLYGONS", &pd.polys),
    ] {
        if c.len() > 0 {
            w.cells(key, &c.offsets, &c.conn);
        }
    }
    w.arrays("POINT_DATA", pd.points.len(), &pa);
    w.arrays("CELL_DATA", pd.cell_data.rows(), &ca);
    Ok(w.out)
}
