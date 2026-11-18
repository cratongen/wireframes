use crate::{Error, Result};

/// Equality compares `f64` by bits, except that every NaN equals every NaN,
/// so tables with null floats compare equal after a round trip.
#[derive(Clone, Debug)]
pub enum ColumnData {
    Float64(Vec<f64>),
    Int64(Vec<i64>),
    Bool(Vec<bool>),
    Utf8(Vec<String>),
}

/// A typed column with an optional validity mask and unit label. Null rows hold
/// `NaN`, `0`, `false` or `""`.
#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    pub data: ColumnData,
    /// None means every row is valid.
    pub valid: Option<Vec<bool>>,
    pub unit: Option<String>,
}

macro_rules! map_data {
    ($data:expr, $v:ident => $body:expr) => {
        match $data {
            ColumnData::Float64($v) => ColumnData::Float64($body),
            ColumnData::Int64($v) => ColumnData::Int64($body),
            ColumnData::Bool($v) => ColumnData::Bool($body),
            ColumnData::Utf8($v) => ColumnData::Utf8($body),
        }
    };
}

impl PartialEq for ColumnData {
    fn eq(&self, o: &ColumnData) -> bool {
        match (self, o) {
            (ColumnData::Float64(a), ColumnData::Float64(b)) => {
                a.len() == b.len()
                    && a.iter()
                        .zip(b)
                        .all(|(x, y)| x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()))
            }
            (ColumnData::Int64(a), ColumnData::Int64(b)) => a == b,
            (ColumnData::Bool(a), ColumnData::Bool(b)) => a == b,
            (ColumnData::Utf8(a), ColumnData::Utf8(b)) => a == b,
            _ => false,
        }
    }
}

impl ColumnData {
    pub fn len(&self) -> usize {
        match self {
            ColumnData::Float64(v) => v.len(),
            ColumnData::Int64(v) => v.len(),
            ColumnData::Bool(v) => v.len(),
            ColumnData::Utf8(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn same_type(&self, o: &ColumnData) -> bool {
        std::mem::discriminant(self) == std::mem::discriminant(o)
    }
}

impl Column {
    fn new(data: ColumnData) -> Column {
        Column {
            data,
            valid: None,
            unit: None,
        }
    }

    pub fn f64(v: Vec<f64>) -> Column {
        Column::new(ColumnData::Float64(v))
    }

    pub fn i64(v: Vec<i64>) -> Column {
        Column::new(ColumnData::Int64(v))
    }

    pub fn bool(v: Vec<bool>) -> Column {
        Column::new(ColumnData::Bool(v))
    }

    pub fn utf8(v: Vec<String>) -> Column {
        Column::new(ColumnData::Utf8(v))
    }

    /// An all-null column of `n` rows with the type of `kind`.
    pub fn nulls(kind: &ColumnData, n: usize) -> Column {
        let data = match kind {
            ColumnData::Float64(_) => ColumnData::Float64(vec![f64::NAN; n]),
            ColumnData::Int64(_) => ColumnData::Int64(vec![0; n]),
            ColumnData::Bool(_) => ColumnData::Bool(vec![false; n]),
            ColumnData::Utf8(_) => ColumnData::Utf8(vec![String::new(); n]),
        };
        Column {
            data,
            valid: Some(vec![false; n]),
            unit: None,
        }
    }

    /// The same column with validity `valid`; null rows get the null payload.
    pub fn with_valid(mut self, valid: Vec<bool>) -> Result<Column> {
        if valid.len() != self.len() {
            return Err(Error::invalid(format!(
                "validity mask has {} rows, column has {}",
                valid.len(),
                self.len()
            )));
        }
        match &mut self.data {
            ColumnData::Float64(v) => null_out(v, &valid, f64::NAN),
            ColumnData::Int64(v) => null_out(v, &valid, 0),
            ColumnData::Bool(v) => null_out(v, &valid, false),
            ColumnData::Utf8(v) => null_out(v, &valid, String::new()),
        }
        self.valid = if valid.iter().all(|&b| b) {
            None
        } else {
            Some(valid)
        };
        Ok(self)
    }

    pub fn with_unit(mut self, unit: Option<String>) -> Column {
        self.unit = unit;
        self
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn is_valid(&self, i: usize) -> bool {
        self.valid.as_ref().is_none_or(|v| v[i])
    }

    pub fn null_count(&self) -> usize {
        self.valid
            .as_ref()
            .map_or(0, |v| v.iter().filter(|&&b| !b).count())
    }

    pub fn type_name(&self) -> &'static str {
        match self.data {
            ColumnData::Float64(_) => "float64",
            ColumnData::Int64(_) => "int64",
            ColumnData::Bool(_) => "bool",
            ColumnData::Utf8(_) => "utf8",
        }
    }

    fn gather(&self, rows: &[Option<u32>]) -> Column {
        let data = map_data!(&self.data, v => gather_vec(v, rows));
        let valid: Vec<bool> = rows
            .iter()
            .map(|r| r.is_some_and(|i| self.is_valid(i as usize)))
            .collect();
        Column {
            data,
            valid: (!valid.iter().all(|&b| b)).then_some(valid),
            unit: self.unit.clone(),
        }
    }
}

trait NullValue {
    fn null() -> Self;
}
impl NullValue for f64 {
    fn null() -> f64 {
        f64::NAN
    }
}
impl NullValue for i64 {
    fn null() -> i64 {
        0
    }
}
impl NullValue for bool {
    fn null() -> bool {
        false
    }
}
impl NullValue for String {
    fn null() -> String {
        String::new()
    }
}

fn gather_vec<T: Clone + NullValue>(v: &[T], rows: &[Option<u32>]) -> Vec<T> {
    rows.iter()
        .map(|r| r.map_or_else(T::null, |i| v[i as usize].clone()))
        .collect()
}

fn null_out<T: Clone>(v: &mut [T], valid: &[bool], null: T) {
    for (x, &ok) in v.iter_mut().zip(valid) {
        if !ok {
            *x = null.clone();
        }
    }
}

/// Named columns of equal length.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Table {
    names: Vec<String>,
    columns: Vec<Column>,
    rows: usize,
}

impl Table {
    pub fn new(rows: usize) -> Table {
        Table {
            names: vec![],
            columns: vec![],
            rows,
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn num_columns(&self) -> usize {
        self.columns.len()
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn columns(&self) -> &[Column] {
        &self.columns
    }

    pub fn column(&self, name: &str) -> Option<&Column> {
        self.names
            .iter()
            .position(|n| n == name)
            .map(|i| &self.columns[i])
    }

    pub fn require(&self, name: &str) -> Result<&Column> {
        self.column(name).ok_or_else(|| Error::MissingColumn {
            name: name.to_string(),
            available: self.names.clone(),
        })
    }

    /// Adds `column` as `name`, replacing an existing column in its position.
    pub fn with_column(mut self, name: &str, column: Column) -> Result<Table> {
        if column.len() != self.rows {
            return Err(Error::invalid(format!(
                "column {name:?} has {} rows, table has {}",
                column.len(),
                self.rows
            )));
        }
        if let Some(v) = &column.valid
            && v.len() != column.len()
        {
            return Err(Error::invalid(format!(
                "column {name:?} validity mask has the wrong length"
            )));
        }
        match self.names.iter().position(|n| n == name) {
            Some(i) => self.columns[i] = column,
            None => {
                self.names.push(name.to_string());
                self.columns.push(column);
            }
        }
        Ok(self)
    }

    pub fn without_column(mut self, name: &str) -> Table {
        if let Some(i) = self.names.iter().position(|n| n == name) {
            self.names.remove(i);
            self.columns.remove(i);
        }
        self
    }

    /// Rows `rows` in that order. Every index must be `< self.rows()`.
    pub fn take(&self, rows: &[u32]) -> Table {
        let rows: Vec<Option<u32>> = rows.iter().map(|&r| Some(r)).collect();
        self.take_or_null(&rows)
    }

    /// Like `take`, with `None` giving a null row.
    pub fn take_or_null(&self, rows: &[Option<u32>]) -> Table {
        Table {
            names: self.names.clone(),
            columns: self.columns.iter().map(|c| c.gather(rows)).collect(),
            rows: rows.len(),
        }
    }

    pub fn filter(&self, keep: &[bool]) -> Result<Table> {
        if keep.len() != self.rows {
            return Err(Error::invalid(format!(
                "mask has {} rows, table has {}",
                keep.len(),
                self.rows
            )));
        }
        let rows: Vec<u32> = (0..self.rows as u32)
            .filter(|&i| keep[i as usize])
            .collect();
        Ok(self.take(&rows))
    }

    /// Rows of `a` then `b`. Columns are the union in first-seen order; a column
    /// missing on one side is null there.
    pub fn concat(a: &Table, b: &Table) -> Result<Table> {
        let mut names = a.names.clone();
        for n in &b.names {
            if !names.contains(n) {
                names.push(n.clone());
            }
        }
        let mut out = Table::new(a.rows + b.rows);
        for name in &names {
            let (ca, cb) = (a.column(name), b.column(name));
            let kind = &ca.or(cb).expect("name comes from a or b").data;
            let ca = ca.cloned().unwrap_or_else(|| Column::nulls(kind, a.rows));
            let cb = cb.cloned().unwrap_or_else(|| Column::nulls(kind, b.rows));
            if !ca.data.same_type(&cb.data) {
                return Err(Error::invalid(format!(
                    "column {name:?} is {} in one table and {} in the other",
                    ca.type_name(),
                    cb.type_name()
                )));
            }
            let data = match (ca.data, cb.data) {
                (ColumnData::Float64(mut x), ColumnData::Float64(y)) => {
                    x.extend(y);
                    ColumnData::Float64(x)
                }
                (ColumnData::Int64(mut x), ColumnData::Int64(y)) => {
                    x.extend(y);
                    ColumnData::Int64(x)
                }
                (ColumnData::Bool(mut x), ColumnData::Bool(y)) => {
                    x.extend(y);
                    ColumnData::Bool(x)
                }
                (ColumnData::Utf8(mut x), ColumnData::Utf8(y)) => {
                    x.extend(y);
                    ColumnData::Utf8(x)
                }
                _ => unreachable!("types checked above"),
            };
            let valid = match (ca.valid, cb.valid) {
                (None, None) => None,
                (va, vb) => {
                    let mut v = va.unwrap_or_else(|| vec![true; a.rows]);
                    v.extend(vb.unwrap_or_else(|| vec![true; b.rows]));
                    Some(v)
                }
            };
            let unit = ca.unit.or(cb.unit);
            out = out.with_column(name, Column { data, valid, unit })?;
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn sample() -> Table {
        Table::new(3)
            .with_column("x", Column::f64(vec![1.0, 2.0, 3.0]))
            .unwrap()
            .with_column(
                "n",
                Column::i64(vec![4, 5, 6])
                    .with_valid(vec![true, false, true])
                    .unwrap(),
            )
            .unwrap()
            .with_column("s", Column::utf8(vec!["a".into(), "b".into(), "c".into()]))
            .unwrap()
    }

    #[test]
    fn columns_are_length_checked_and_replaced_in_place() {
        let t = sample();
        assert!(t.clone().with_column("x", Column::f64(vec![1.0])).is_err());
        let t = t.with_column("x", Column::bool(vec![true; 3])).unwrap();
        assert_eq!(t.names(), ["x", "n", "s"]);
        assert_eq!(t.column("x").unwrap().type_name(), "bool");
        assert_eq!(t.clone().without_column("n").names(), ["x", "s"]);
        assert!(matches!(t.require("zz"), Err(Error::MissingColumn { .. })));
    }

    #[test]
    fn nulls_hold_null_payloads() {
        let t = sample();
        let n = t.column("n").unwrap();
        assert_eq!(n.null_count(), 1);
        assert!(!n.is_valid(1));
        assert_eq!(n.data, ColumnData::Int64(vec![4, 0, 6]));
        let c = Column::nulls(&ColumnData::Float64(vec![]), 2);
        assert_eq!(c.null_count(), 2);
        assert!(matches!(&c.data, ColumnData::Float64(v) if v.iter().all(|x| x.is_nan())));
    }

    #[test]
    fn take_or_null_and_filter() {
        let t = sample();
        let g = t.take_or_null(&[Some(2), None, Some(0)]);
        assert_eq!(g.rows(), 3);
        assert_eq!(
            g.column("s").unwrap().data,
            ColumnData::Utf8(vec!["c".into(), "".into(), "a".into()])
        );
        assert_eq!(g.column("s").unwrap().valid, Some(vec![true, false, true]));
        let f = t.filter(&[false, true, true]).unwrap();
        assert_eq!(
            f.column("x").unwrap().data,
            ColumnData::Float64(vec![2.0, 3.0])
        );
        assert!(t.filter(&[true]).is_err());
    }

    #[test]
    fn concat_unions_columns_and_checks_types() {
        let a = sample();
        let b = Table::new(1)
            .with_column("y", Column::f64(vec![9.0]))
            .unwrap();
        let c = Table::concat(&a, &b).unwrap();
        assert_eq!(c.names(), ["x", "n", "s", "y"]);
        assert_eq!(c.rows(), 4);
        assert_eq!(
            c.column("x").unwrap().valid,
            Some(vec![true, true, true, false])
        );
        assert_eq!(
            c.column("y").unwrap().valid,
            Some(vec![false, false, false, true])
        );
        let bad = Table::new(1)
            .with_column("x", Column::i64(vec![1]))
            .unwrap();
        assert!(Table::concat(&a, &bad).is_err());
    }

    proptest! {
        #[test]
        fn take_then_concat_of_parts_is_identity(values in proptest::collection::vec(-1e6..1e6f64, 0..40), cut in 0usize..40) {
            let n = values.len();
            let cut = cut.min(n);
            let t = Table::new(n).with_column("v", Column::f64(values)).unwrap();
            let head: Vec<u32> = (0..cut as u32).collect();
            let tail: Vec<u32> = (cut as u32..n as u32).collect();
            let joined = Table::concat(&t.take(&head), &t.take(&tail)).unwrap();
            prop_assert_eq!(joined, t);
        }
    }
}
