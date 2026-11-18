use crate::mesh::check_finite;
use crate::{Aabb, Column, Error, Result, Table};

/// Points with one attribute row each.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct PointSet {
    coords: Vec<[f64; 3]>,
    attributes: Table,
    pub crs: Option<String>,
    pub length_unit: Option<String>,
}

impl PointSet {
    pub fn new(coords: Vec<[f64; 3]>) -> Result<PointSet> {
        check_finite("coords", &coords)?;
        let n = coords.len();
        Ok(PointSet {
            coords,
            attributes: Table::new(n),
            crs: None,
            length_unit: None,
        })
    }

    pub fn with_attributes(mut self, table: Table) -> Result<PointSet> {
        if table.rows() != self.coords.len() {
            return Err(Error::invalid(format!(
                "attribute table has {} rows, there are {} points",
                table.rows(),
                self.coords.len()
            )));
        }
        self.attributes = table;
        Ok(self)
    }

    pub fn with_column(mut self, name: &str, column: Column) -> Result<PointSet> {
        self.attributes = self.attributes.with_column(name, column)?;
        Ok(self)
    }

    pub fn coords(&self) -> &[[f64; 3]] {
        &self.coords
    }

    pub fn attributes(&self) -> &Table {
        &self.attributes
    }

    pub fn len(&self) -> usize {
        self.coords.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }

    pub fn bounds(&self) -> Option<Aabb> {
        (!self.coords.is_empty()).then(|| Aabb::from_points(&self.coords))
    }

    pub fn filter(&self, keep: &[bool]) -> Result<PointSet> {
        let attributes = self.attributes.filter(keep)?;
        let coords = self
            .coords
            .iter()
            .zip(keep)
            .filter(|(_, k)| **k)
            .map(|(p, _)| *p)
            .collect();
        Ok(PointSet {
            coords,
            attributes,
            crs: self.crs.clone(),
            length_unit: self.length_unit.clone(),
        })
    }

    pub fn with_crs(mut self, crs: Option<String>) -> PointSet {
        self.crs = crs;
        self
    }

    pub fn with_length_unit(mut self, unit: Option<String>) -> PointSet {
        self.length_unit = unit;
        self
    }

    pub fn with_vertices(mut self, coords: Vec<[f64; 3]>) -> Result<PointSet> {
        if coords.len() != self.coords.len() {
            return Err(Error::invalid(format!(
                "coords has {} rows, point set has {}",
                coords.len(),
                self.coords.len()
            )));
        }
        check_finite("coords", &coords)?;
        self.coords = coords;
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn points_filter_with_attributes() {
        let p = PointSet::new(vec![[0.0; 3], [1.0, 2.0, 3.0], [4.0; 3]])
            .unwrap()
            .with_column("g", Column::f64(vec![0.1, 0.2, 0.3]))
            .unwrap()
            .with_crs(Some("EPSG:31983".into()));
        let f = p.filter(&[false, true, true]).unwrap();
        assert_eq!(f.coords(), &[[1.0, 2.0, 3.0], [4.0; 3]]);
        assert_eq!(f.attributes().rows(), 2);
        assert_eq!(f.crs, p.crs);
        assert_eq!(p.bounds().unwrap().max, [4.0; 3]);
        assert!(PointSet::new(vec![[f64::INFINITY, 0.0, 0.0]]).is_err());
        assert!(p.clone().with_vertices(vec![]).is_err());
        assert!(p.with_attributes(Table::new(1)).is_err());
    }
}
