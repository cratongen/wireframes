use std::ops::Range;

use crate::mesh::check_finite;
use crate::predicates::{Sign, orient2d};
use crate::{Aabb, Column, Error, Result, Table, units};

/// 3D strings grouped into features. A closed part does not repeat its first
/// vertex; in plan, a feature's inside is even-odd over its closed parts.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Polylines {
    vertices: Vec<[f64; 3]>,
    /// Offsets into `vertices`, `n_parts + 1` long, starting at 0.
    parts: Vec<u32>,
    /// Offsets into parts, `n_features + 1` long, starting at 0.
    features: Vec<u32>,
    closed: Vec<bool>,
    vertex_attributes: Table,
    /// One row per feature.
    attributes: Table,
    pub crs: Option<String>,
    pub length_unit: Option<String>,
}

fn check_offsets(what: &str, offsets: &[u32], total: usize) -> Result<()> {
    if offsets.first() != Some(&0) || *offsets.last().expect("non-empty") as usize != total {
        return Err(Error::invalid(format!(
            "{what} offsets must start at 0 and end at {total}"
        )));
    }
    if offsets.windows(2).any(|w| w[0] > w[1]) {
        return Err(Error::invalid(format!(
            "{what} offsets must be non-decreasing"
        )));
    }
    Ok(())
}

fn in_ring(xy: [f64; 2], ring: &[[f64; 3]]) -> bool {
    let mut inside = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        let (a, b) = ([a[0], a[1]], [b[0], b[1]]);
        let o = Sign::of(orient2d(a, b, xy));
        if o == Sign::Zero
            && xy[0] >= a[0].min(b[0])
            && xy[0] <= a[0].max(b[0])
            && xy[1] >= a[1].min(b[1])
            && xy[1] <= a[1].max(b[1])
        {
            return true;
        }
        if (a[1] > xy[1]) != (b[1] > xy[1]) {
            let up = b[1] > a[1];
            if (o == Sign::Positive) == up {
                inside = !inside;
            }
        }
    }
    inside
}

fn ring_area(ring: &[[f64; 3]]) -> f64 {
    let o = ring[0];
    let s: f64 = (1..ring.len().saturating_sub(1))
        .map(|i| {
            let (a, b) = (ring[i], ring[i + 1]);
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        })
        .sum();
    0.5 * s.abs()
}

fn segment_distance(p: [f64; 2], a: [f64; 3], b: [f64; 3]) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len2 = dx * dx + dy * dy;
    let t = if len2 > 0.0 {
        (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let (x, y) = (a[0] + t * dx - p[0], a[1] + t * dy - p[1]);
    (x * x + y * y).sqrt()
}

impl Polylines {
    pub fn new(
        vertices: Vec<[f64; 3]>,
        parts: Vec<u32>,
        features: Vec<u32>,
        closed: Vec<bool>,
        attributes: Table,
    ) -> Result<Polylines> {
        check_finite("vertices", &vertices)?;
        check_offsets("part", &parts, vertices.len())?;
        let n_parts = parts.len() - 1;
        check_offsets("feature", &features, n_parts)?;
        if closed.len() != n_parts {
            return Err(Error::invalid(format!(
                "closed has {} values, there are {n_parts} parts",
                closed.len()
            )));
        }
        for p in 0..n_parts {
            let n = (parts[p + 1] - parts[p]) as usize;
            let need = if closed[p] { 3 } else { 1 };
            if n < need {
                return Err(Error::invalid(format!(
                    "part {p} has {n} vertices; {} parts need at least {need}",
                    if closed[p] { "closed" } else { "open" }
                )));
            }
        }
        let n_features = features.len() - 1;
        if attributes.rows() != n_features {
            return Err(Error::invalid(format!(
                "attribute table has {} rows, there are {n_features} features",
                attributes.rows()
            )));
        }
        let nv = vertices.len();
        Ok(Polylines {
            vertices,
            parts,
            features,
            closed,
            vertex_attributes: Table::new(nv),
            attributes,
            crs: None,
            length_unit: None,
        })
    }

    /// `feature_of_part` is non-decreasing from 0; `None` makes each part its own feature.
    pub fn from_parts(
        parts: Vec<Vec<[f64; 3]>>,
        closed: Vec<bool>,
        feature_of_part: Option<Vec<u32>>,
    ) -> Result<Polylines> {
        let n = parts.len();
        let feature_of_part = feature_of_part.unwrap_or_else(|| (0..n as u32).collect());
        if feature_of_part.len() != n {
            return Err(Error::invalid(format!(
                "feature_of_part has {} values, there are {n} parts",
                feature_of_part.len()
            )));
        }
        let mut features = vec![0u32];
        let mut prev: Option<u32> = None;
        for (p, &f) in feature_of_part.iter().enumerate() {
            let ok = match prev {
                None => f == 0,
                Some(q) => f == q || f == q + 1,
            };
            if !ok {
                return Err(Error::invalid(
                    "feature_of_part must start at 0 and increase by at most 1",
                ));
            }
            if prev.is_some_and(|q| f != q) {
                features.push(p as u32);
            }
            prev = Some(f);
        }
        if n > 0 {
            features.push(n as u32);
        }
        let mut offsets = vec![0u32];
        let mut vertices = vec![];
        for part in parts {
            vertices.extend(part);
            offsets.push(vertices.len() as u32);
        }
        let nf = features.len() - 1;
        Polylines::new(vertices, offsets, features, closed, Table::new(nf))
    }

    pub fn with_vertex_attributes(mut self, table: Table) -> Result<Polylines> {
        if table.rows() != self.vertices.len() {
            return Err(Error::invalid(format!(
                "vertex attribute table has {} rows, there are {} vertices",
                table.rows(),
                self.vertices.len()
            )));
        }
        self.vertex_attributes = table;
        Ok(self)
    }

    /// Number of features.
    pub fn len(&self) -> usize {
        self.features.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn n_parts(&self) -> usize {
        self.parts.len() - 1
    }

    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }

    pub fn part_offsets(&self) -> &[u32] {
        &self.parts
    }

    pub fn feature_offsets(&self) -> &[u32] {
        &self.features
    }

    pub fn part(&self, p: usize) -> &[[f64; 3]] {
        &self.vertices[self.part_range(p)]
    }

    pub fn part_range(&self, p: usize) -> Range<usize> {
        self.parts[p] as usize..self.parts[p + 1] as usize
    }

    pub fn feature_parts(&self, f: usize) -> Range<usize> {
        self.features[f] as usize..self.features[f + 1] as usize
    }

    pub fn closed(&self) -> &[bool] {
        &self.closed
    }

    /// Closed parts of feature `f`.
    pub fn rings(&self, f: usize) -> Vec<usize> {
        self.feature_parts(f).filter(|&p| self.closed[p]).collect()
    }

    pub fn attributes(&self) -> &Table {
        &self.attributes
    }

    pub fn vertex_attributes(&self) -> &Table {
        &self.vertex_attributes
    }

    pub fn with_column(mut self, name: &str, column: Column) -> Result<Polylines> {
        self.attributes = self.attributes.with_column(name, column)?;
        Ok(self)
    }

    pub fn with_vertex_column(mut self, name: &str, column: Column) -> Result<Polylines> {
        self.vertex_attributes = self.vertex_attributes.with_column(name, column)?;
        Ok(self)
    }

    pub fn with_attributes(mut self, table: Table) -> Result<Polylines> {
        if table.rows() != self.len() {
            return Err(Error::invalid(format!(
                "attribute table has {} rows, there are {} features",
                table.rows(),
                self.len()
            )));
        }
        self.attributes = table;
        Ok(self)
    }

    /// Whole features where `keep` is true.
    pub fn filter(&self, keep: &[bool]) -> Result<Polylines> {
        let attributes = self.attributes.filter(keep)?;
        let (mut vertices, mut parts, mut features, mut closed, mut rows) =
            (vec![], vec![0u32], vec![0u32], vec![], vec![]);
        for f in (0..self.len()).filter(|&f| keep[f]) {
            for p in self.feature_parts(f) {
                vertices.extend_from_slice(self.part(p));
                rows.extend(self.part_range(p).map(|i| i as u32));
                parts.push(vertices.len() as u32);
                closed.push(self.closed[p]);
            }
            features.push(parts.len() as u32 - 1);
        }
        let mut out = Polylines::new(vertices, parts, features, closed, attributes)?
            .with_vertex_attributes(self.vertex_attributes.take(&rows))?;
        out.crs.clone_from(&self.crs);
        out.length_unit.clone_from(&self.length_unit);
        Ok(out)
    }

    pub fn bounds(&self) -> Option<Aabb> {
        (!self.vertices.is_empty()).then(|| Aabb::from_points(&self.vertices))
    }

    pub fn concat(items: &[&Polylines]) -> Result<Polylines> {
        let (mut vertices, mut parts, mut features, mut closed) =
            (vec![], vec![0u32], vec![0u32], vec![]);
        let (mut attributes, mut vattr) = (Table::new(0), Table::new(0));
        let mut labels = (None, None);
        for (k, it) in items.iter().enumerate() {
            labels = units::check_same(
                &format!("polylines {k}"),
                (&it.crs, &it.length_unit),
                (&labels.0, &labels.1),
            )?;
            let (v0, p0) = (vertices.len() as u32, parts.len() as u32 - 1);
            vertices.extend_from_slice(&it.vertices);
            parts.extend(it.parts[1..].iter().map(|&o| o + v0));
            features.extend(it.features[1..].iter().map(|&o| o + p0));
            closed.extend_from_slice(&it.closed);
            attributes = Table::concat(&attributes, &it.attributes)?;
            vattr = Table::concat(&vattr, &it.vertex_attributes)?;
        }
        let mut out = Polylines::new(vertices, parts, features, closed, attributes)?
            .with_vertex_attributes(vattr)?;
        out.crs = labels.0;
        out.length_unit = labels.1;
        Ok(out)
    }

    /// 3D length per feature, closed parts including their closing segment.
    pub fn length(&self) -> Vec<f64> {
        let dist = |a: [f64; 3], b: [f64; 3]| {
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        };
        (0..self.len())
            .map(|f| {
                self.feature_parts(f)
                    .map(|p| {
                        let v = self.part(p);
                        let open: f64 = v.windows(2).map(|w| dist(w[0], w[1])).sum();
                        open + if self.closed[p] {
                            dist(v[v.len() - 1], v[0])
                        } else {
                            0.0
                        }
                    })
                    .sum()
            })
            .collect()
    }

    /// Plan area per feature, even-odd: a ring inside an odd number of the
    /// feature's other rings is a hole. Rings of a feature must not cross.
    pub fn area(&self) -> Vec<f64> {
        (0..self.len())
            .map(|f| {
                let rings = self.rings(f);
                rings
                    .iter()
                    .map(|&r| {
                        let p = self.part(r)[0];
                        let depth = rings
                            .iter()
                            .filter(|&&o| o != r && in_ring([p[0], p[1]], self.part(o)))
                            .count();
                        let a = ring_area(self.part(r));
                        if depth % 2 == 0 { a } else { -a }
                    })
                    .sum()
            })
            .collect()
    }

    /// Plan even-odd containment; a point on a ring edge counts as inside.
    pub fn feature_contains(&self, f: usize, xy: [f64; 2]) -> bool {
        let rings = self.rings(f);
        if rings.iter().any(|&r| on_ring(xy, self.part(r))) {
            return true;
        }
        rings.iter().filter(|&&r| in_ring(xy, self.part(r))).count() % 2 == 1
    }

    /// First feature containing `xy` in plan.
    pub fn locate(&self, xy: [f64; 2]) -> Option<usize> {
        (0..self.len()).find(|&f| self.feature_contains(f, xy))
    }

    /// Plan distance to the nearest ring edge of `feature`, or of any feature;
    /// negative inside when `signed`. `None` without rings.
    pub fn distance(&self, xy: [f64; 2], feature: Option<usize>, signed: bool) -> Option<f64> {
        let fs = feature.map_or(0..self.len(), |f| f..f + 1);
        let d = fs
            .clone()
            .flat_map(|f| self.rings(f))
            .flat_map(|r| {
                let v = self.part(r);
                (0..v.len()).map(move |i| segment_distance(xy, v[i], v[(i + 1) % v.len()]))
            })
            .reduce(f64::min)?;
        let inside = signed && fs.into_iter().any(|f| self.feature_contains(f, xy));
        Some(if inside { -d } else { d })
    }

    pub fn with_crs(mut self, crs: Option<String>) -> Polylines {
        self.crs = crs;
        self
    }

    pub fn with_length_unit(mut self, unit: Option<String>) -> Polylines {
        self.length_unit = unit;
        self
    }

    pub fn with_vertices(mut self, vertices: Vec<[f64; 3]>) -> Result<Polylines> {
        if vertices.len() != self.vertices.len() {
            return Err(Error::invalid(format!(
                "vertices has {} rows, polylines have {}",
                vertices.len(),
                self.vertices.len()
            )));
        }
        check_finite("vertices", &vertices)?;
        self.vertices = vertices;
        Ok(self)
    }
}

fn on_ring(xy: [f64; 2], ring: &[[f64; 3]]) -> bool {
    (0..ring.len()).any(|i| {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        orient2d([a[0], a[1]], [b[0], b[1]], xy) == 0.0
            && xy[0] >= a[0].min(b[0])
            && xy[0] <= a[0].max(b[0])
            && xy[1] >= a[1].min(b[1])
            && xy[1] <= a[1].max(b[1])
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn square(x0: f64, s: f64) -> Vec<[f64; 3]> {
        vec![
            [x0, x0, 1.0],
            [x0 + s, x0, 1.0],
            [x0 + s, x0 + s, 1.0],
            [x0, x0 + s, 1.0],
        ]
    }

    fn donut() -> Polylines {
        Polylines::from_parts(
            vec![
                square(0.0, 4.0),
                square(1.0, 2.0),
                vec![[9.0, 9.0, 0.0], [9.0, 10.0, 0.0]],
            ],
            vec![true, true, false],
            Some(vec![0, 0, 1]),
        )
        .unwrap()
    }

    #[test]
    fn structure_and_validation() {
        let p = donut();
        assert_eq!((p.len(), p.n_parts(), p.n_vertices()), (2, 3, 10));
        assert_eq!(p.rings(0), vec![0, 1]);
        assert_eq!(p.feature_parts(1), 2..3);
        assert!(Polylines::from_parts(vec![vec![[0.0; 3]; 2]], vec![true], None).is_err());
        assert!(
            Polylines::from_parts(vec![vec![[0.0; 3]; 2]], vec![false], Some(vec![1])).is_err()
        );
        assert!(
            Polylines::new(
                vec![[0.0; 3]; 2],
                vec![0, 2],
                vec![0, 1],
                vec![false],
                Table::new(2)
            )
            .is_err()
        );
        assert!(p.clone().with_column("id", Column::i64(vec![1, 2])).is_ok());
        assert!(p.clone().with_column("id", Column::i64(vec![1])).is_err());
    }

    #[test]
    fn an_open_part_may_be_a_single_point() {
        let p = Polylines::from_parts(vec![vec![[1.0, 2.0, 3.0]]], vec![false], None).unwrap();
        assert_eq!((p.n_parts(), p.part(0)), (1, &[[1.0, 2.0, 3.0]][..]));
        assert_eq!(p.length(), vec![0.0]);
        assert!(Polylines::from_parts(vec![vec![[0.0; 3]]], vec![true], None).is_err());
        assert!(Polylines::from_parts(vec![vec![]], vec![false], None).is_err());
    }

    #[test]
    fn plan_measures_of_a_donut() {
        let p = donut();
        assert_eq!(p.area(), vec![12.0, 0.0]);
        assert_eq!(p.length(), vec![24.0, 1.0]);
        assert!(p.feature_contains(0, [0.5, 0.5]));
        assert!(!p.feature_contains(0, [2.0, 2.0]));
        assert!(p.feature_contains(0, [1.0, 2.0]));
        assert!(p.feature_contains(0, [0.0, 0.0]));
        assert!(!p.feature_contains(0, [5.0, 2.0]));
        assert_eq!(p.locate([0.5, 3.5]), Some(0));
        assert_eq!(p.locate([2.0, 2.0]), None);
        assert_eq!(p.distance([2.0, 2.0], None, true), Some(1.0));
        assert_eq!(p.distance([0.5, 2.0], Some(0), true), Some(-0.5));
        assert_eq!(p.distance([0.5, 2.0], Some(1), true), None);
    }

    #[test]
    fn filter_and_concat_keep_vertex_rows() {
        let p = donut()
            .with_vertex_column("i", Column::i64((0..10).collect()))
            .unwrap();
        let f = p.filter(&[false, true]).unwrap();
        assert_eq!(f.n_vertices(), 2);
        assert_eq!(
            f.vertex_attributes().column("i").unwrap().data,
            crate::ColumnData::Int64(vec![8, 9])
        );
        let c = Polylines::concat(&[&f, &p]).unwrap();
        assert_eq!((c.len(), c.n_parts(), c.n_vertices()), (3, 4, 12));
        assert_eq!(c.feature_parts(2), 3..4);
        assert_eq!(c.part(0), f.part(0));
    }

    proptest! {
        #[test]
        fn convex_polygon_contains_its_centroid(n in 3usize..20, r in 0.1..100.0f64, cx in -1e5..1e5f64) {
            let ring: Vec<[f64; 3]> = (0..n)
                .map(|i| {
                    let a = i as f64 / n as f64 * std::f64::consts::TAU;
                    [cx + r * a.cos(), r * a.sin(), 0.0]
                })
                .collect();
            let p = Polylines::from_parts(vec![ring], vec![true], None).unwrap();
            prop_assert!(p.feature_contains(0, [cx, 0.0]));
            prop_assert!(!p.feature_contains(0, [cx + 2.0 * r, 0.0]));
            let want = 0.5 * n as f64 * r * r * (std::f64::consts::TAU / n as f64).sin();
            prop_assert!((p.area()[0] - want).abs() <= 1e-9 * want.max(1.0) * (1.0 + cx.abs() / r));
        }
    }
}
