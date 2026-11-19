use std::collections::HashMap;

use wireframes_core::{Error, Mesh, Polylines, Result, tolerance};

use crate::geom::{Cdt, check_finite, ekey, even_odd};

#[derive(Clone, Debug, Default)]
pub struct DtmOptions<'a> {
    pub breaklines: Option<&'a Polylines>,
    /// Closed rings, even-odd; triangles outside are dropped.
    pub boundary: Option<&'a Polylines>,
    /// Plan duplicate merge distance; 0 merges exact duplicates only.
    pub min_spacing: f64,
    /// Drops boundary triangles with a plan edge longer than this, repeatedly.
    pub max_edge: Option<f64>,
}

/// An input vertex: `what`, then the part for strings, then the index.
#[derive(Clone, Copy)]
struct Src(&'static str, Option<usize>, usize);

impl std::fmt::Display for Src {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.1 {
            None => write!(f, "{}[{}]", self.0, self.2),
            Some(part) => write!(f, "{} part {part} vertex {}", self.0, self.2),
        }
    }
}

/// Merges plan positions closer than `spacing` into the earliest. At
/// `spacing = 0` only exact 3D duplicates merge: two heights at one plan
/// position are an error, never a silent choice.
struct Merger {
    spacing: f64,
    origin: [f64; 2],
    kept: Vec<[f64; 3]>,
    /// Where each kept point came from, for error messages.
    source: Vec<Src>,
    grid: HashMap<(u64, u64), Vec<u32>>,
}

fn cell_key(x: f64, y: f64) -> (u64, u64) {
    let c = tolerance::canonical;
    (c(x).to_bits(), c(y).to_bits())
}

impl Merger {
    /// Cell coordinates stay `f64` so a spacing tiny against the coordinates
    /// cannot overflow; neighbour cells that round onto the same value are
    /// visited once.
    fn cell(&self, p: [f64; 3]) -> [f64; 2] {
        if self.spacing > 0.0 {
            [0, 1].map(|k| ((p[k] - self.origin[k]) / self.spacing).floor())
        } else {
            [p[0], p[1]]
        }
    }

    fn add(&mut self, p: [f64; 3], source: Src) -> Result<u32> {
        let [cx, cy] = self.cell(p);
        let steps: &[f64] = if self.spacing > 0.0 {
            &[-1.0, 0.0, 1.0]
        } else {
            &[0.0]
        };
        let mut cells: Vec<(u64, u64)> = steps
            .iter()
            .flat_map(|&dx| steps.iter().map(move |&dy| cell_key(cx + dx, cy + dy)))
            .collect();
        cells.sort_unstable();
        cells.dedup();
        let s2 = self.spacing * self.spacing;
        let mut found: Option<u32> = None;
        for &i in cells.iter().flat_map(|c| self.grid.get(c)).flatten() {
            let q = self.kept[i as usize];
            let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2);
            let hit = if self.spacing > 0.0 {
                d2 < s2
            } else {
                p[0] == q[0] && p[1] == q[1]
            };
            if hit && found.is_none_or(|f| i < f) {
                found = Some(i);
            }
        }
        if let Some(i) = found {
            let q = self.kept[i as usize];
            if self.spacing == 0.0 && p[2] != q[2] {
                return Err(Error::invalid(format!(
                    "{} and {} share the plan position ({}, {}) with different z ({} and {}): \
                     a dtm holds one height per position; remove one, or pass min_spacing > 0 \
                     to keep the earlier",
                    self.source[i as usize], source, p[0], p[1], q[2], p[2]
                )));
            }
            return Ok(i);
        }
        let i = self.kept.len() as u32;
        self.kept.push(p);
        self.source.push(source);
        self.grid.entry(cell_key(cx, cy)).or_default().push(i);
        Ok(i)
    }
}

pub fn dtm(points: &[[f64; 3]], opts: &DtmOptions) -> Result<Mesh> {
    check_finite("points", points)?;
    let spacing = tolerance::check("min_spacing", opts.min_spacing)?;
    if let Some(m) = opts.max_edge {
        tolerance::check("max_edge", m)?;
    }
    for (what, pl) in [("breaklines", opts.breaklines), ("boundary", opts.boundary)] {
        if let Some(pl) = pl {
            check_finite(what, pl.vertices())?;
        }
    }
    let mut merger = Merger {
        spacing,
        origin: points.first().map_or([0.0; 2], |p| [p[0], p[1]]),
        kept: vec![],
        source: vec![],
        grid: HashMap::new(),
    };
    for (i, &p) in points.iter().enumerate() {
        merger.add(p, Src("points", None, i))?;
    }
    let n_points = merger.kept.len();
    let mut strings: Vec<(Vec<u32>, bool)> = vec![];
    if let Some(b) = opts.breaklines {
        for part in 0..b.n_parts() {
            let ids = (b.part(part).iter().enumerate())
                .map(|(k, &p)| merger.add(p, Src("breaklines", Some(part), k)))
                .collect::<Result<_>>()?;
            strings.push((ids, b.closed()[part]));
        }
    }
    let mut rings2: Vec<Vec<[f64; 2]>> = vec![];
    if let Some(b) = opts.boundary {
        for part in 0..b.n_parts() {
            if !b.closed()[part] {
                return Err(Error::invalid(format!(
                    "boundary must hold closed rings: part {part} is open"
                )));
            }
            let ring = b.part(part);
            rings2.push(ring.iter().map(|p| [p[0], p[1]]).collect());
            let ids = (ring.iter().enumerate())
                .map(|(k, &p)| merger.add(p, Src("boundary", Some(part), k)))
                .collect::<Result<_>>()?;
            strings.push((ids, true));
        }
    }
    let mut cdt = Cdt::new(merger.kept)?;
    for id in 0..n_points as u32 {
        cdt.pin(id);
    }
    for (ids, closed) in &strings {
        let n = ids.len();
        let segs = if *closed { n } else { n.saturating_sub(1) };
        for k in 0..segs {
            cdt.constrain(ids[k], ids[(k + 1) % n]);
        }
    }
    let tris = cdt.triangles();
    if tris.is_empty() {
        return Err(Error::invalid(
            "points must hold at least 3 points not all collinear in plan",
        ));
    }
    let verts = cdt.verts;
    let mut keep: Vec<bool> = if opts.boundary.is_some() {
        tris.iter()
            .map(|t| {
                let c = t.map(|v| verts[v as usize]);
                let m = [
                    (c[0][0] + c[1][0] + c[2][0]) / 3.0,
                    (c[0][1] + c[1][1] + c[2][1]) / 3.0,
                ];
                even_odd(&rings2, m)
            })
            .collect()
    } else {
        vec![true; tris.len()]
    };
    if let Some(max_edge) = opts.max_edge {
        peel_long_edges(&verts, &tris, &mut keep, max_edge);
    }
    let tris = tris
        .into_iter()
        .zip(keep)
        .filter(|(_, k)| *k)
        .map(|(t, _)| t);
    Mesh::new(verts, tris.collect())
}

fn plan_len(verts: &[[f64; 3]], a: u32, b: u32) -> f64 {
    let (p, q) = (verts[a as usize], verts[b as usize]);
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt()
}

fn peel_long_edges(verts: &[[f64; 3]], tris: &[[u32; 3]], keep: &mut [bool], max_edge: f64) {
    let mut edge_faces: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    for (f, t) in tris.iter().enumerate() {
        if keep[f] {
            for e in 0..3 {
                edge_faces
                    .entry(ekey(t[e], t[(e + 1) % 3]))
                    .or_default()
                    .push(f as u32);
            }
        }
    }
    let exposed = |f: u32, keep: &[bool], edge_faces: &HashMap<(u32, u32), Vec<u32>>| {
        let t = tris[f as usize];
        (0..3).any(|e| {
            let (a, b) = (t[e], t[(e + 1) % 3]);
            let alive = edge_faces[&ekey(a, b)]
                .iter()
                .filter(|&&g| keep[g as usize])
                .count();
            alive == 1 && plan_len(verts, a, b) > max_edge
        })
    };
    let mut stack: Vec<u32> = (0..tris.len() as u32)
        .rev()
        .filter(|&f| keep[f as usize])
        .collect();
    while let Some(f) = stack.pop() {
        if !keep[f as usize] || !exposed(f, keep, &edge_faces) {
            continue;
        }
        keep[f as usize] = false;
        let t = tris[f as usize];
        for e in 0..3 {
            for &g in &edge_faces[&ekey(t[e], t[(e + 1) % 3])] {
                if keep[g as usize] {
                    stack.push(g);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::collections::HashSet;
    use wireframes_core::fixtures::random_points;
    use wireframes_core::predicates::{incircle, orient2d};

    fn xy(p: [f64; 3]) -> [f64; 2] {
        [p[0], p[1]]
    }

    fn edges(m: &Mesh) -> HashMap<(u32, u32), Vec<u32>> {
        let mut e: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (f, t) in m.triangles().iter().enumerate() {
            for k in 0..3 {
                e.entry(ekey(t[k], t[(k + 1) % 3]))
                    .or_default()
                    .push(f as u32);
            }
        }
        e
    }

    fn on_segment(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
        orient2d(a, b, p) == 0.0
            && p[0] >= a[0].min(b[0])
            && p[0] <= a[0].max(b[0])
            && p[1] >= a[1].min(b[1])
            && p[1] <= a[1].max(b[1])
    }

    /// Mesh edges lying on each breakline segment, and a check that they chain from end to end.
    fn constrained_edges(m: &Mesh, lines: &Polylines) -> HashSet<(u32, u32)> {
        let v = m.vertices();
        let es = edges(m);
        let mut out = HashSet::new();
        for part in 0..lines.n_parts() {
            let pts = lines.part(part);
            for k in 0..pts.len() - 1 {
                let (a, b) = (xy(pts[k]), xy(pts[k + 1]));
                let d = [b[0] - a[0], b[1] - a[1]];
                let mut on: Vec<(f64, u32)> = (0..v.len() as u32)
                    .filter(|&i| on_segment(a, b, xy(v[i as usize])))
                    .map(|i| {
                        let p = v[i as usize];
                        ((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1], i)
                    })
                    .collect();
                on.sort_by(|x, y| x.0.total_cmp(&y.0));
                assert_eq!(xy(v[on[0].1 as usize]), a);
                assert_eq!(xy(v[on.last().unwrap().1 as usize]), b);
                for w in on.windows(2) {
                    let e = ekey(w[0].1, w[1].1);
                    assert!(es.contains_key(&e), "breakline piece {e:?} is not an edge");
                    out.insert(e);
                }
            }
        }
        out
    }

    fn assert_delaunay(m: &Mesh, constrained: &HashSet<(u32, u32)>) {
        let v = m.vertices();
        let t = m.triangles();
        for (e, fs) in edges(m) {
            if fs.len() != 2 || constrained.contains(&e) {
                continue;
            }
            let a = t[fs[0] as usize];
            let opp = t[fs[1] as usize]
                .into_iter()
                .find(|&x| x != e.0 && x != e.1)
                .unwrap();
            let c = a.map(|i| xy(v[i as usize]));
            assert!(
                incircle(c[0], c[1], c[2], xy(v[opp as usize])) <= 0.0,
                "edge {e:?} fails the empty circle test"
            );
        }
    }

    fn assert_up(m: &Mesh) {
        for f in 0..m.n_faces() as u32 {
            let c = m.corners(f);
            assert!(orient2d(xy(c[0]), xy(c[1]), xy(c[2])) > 0.0);
        }
    }

    fn strip_breaklines(seed: u64, k: usize) -> Polylines {
        let mut rng = wireframes_core::fixtures::Rng::new(seed);
        let parts = (0..k)
            .map(|s| {
                let n = 2 + rng.below(4);
                let mut xs: Vec<f64> = (0..n).map(|_| rng.uniform(0.0, 10.0)).collect();
                xs.sort_by(f64::total_cmp);
                xs.dedup();
                if xs.len() < 2 {
                    xs = vec![1.0, 9.0];
                }
                xs.iter()
                    .map(|&x| [x, s as f64 + rng.uniform(0.1, 0.9), rng.uniform(0.0, 5.0)])
                    .collect()
            })
            .collect::<Vec<Vec<[f64; 3]>>>();
        let closed = vec![false; parts.len()];
        Polylines::from_parts(parts, closed, None).unwrap()
    }

    proptest! {
        #[test]
        fn delaunay_and_breaklines_are_edges(seed in 0u64..10_000, n in 3usize..120, k in 0usize..5) {
            let pts = random_points(seed, n, [0.0, 0.0, 0.0], [10.0, k.max(1) as f64, 5.0]);
            let lines = strip_breaklines(seed + 1, k);
            let opts = DtmOptions { breaklines: (k > 0).then_some(&lines), ..Default::default() };
            let m = dtm(&pts, &opts).unwrap();
            let constrained = if k > 0 { constrained_edges(&m, &lines) } else { HashSet::new() };
            assert_delaunay(&m, &constrained);
            assert_up(&m);
            if k == 0 {
                prop_assert_eq!(m.vertices(), &pts[..]);
            } else {
                prop_assert_eq!(&m.vertices()[..n], &pts[..]);
            }
        }

        #[test]
        fn duplicates_merge_into_the_earlier_point(seed in 0u64..10_000, n in 3usize..60) {
            let mut pts = random_points(seed, n, [0.0; 3], [1.0; 3]);
            pts.extend(pts.clone());
            let m = dtm(&pts, &DtmOptions::default()).unwrap();
            prop_assert_eq!(m.vertices(), &pts[..n]);
            let k = (seed as usize) % n;
            pts[n + k][2] += 1.0;
            let e = dtm(&pts, &DtmOptions::default()).unwrap_err().to_string();
            let want = format!("points[{k}] and points[{}] share the plan position", n + k);
            prop_assert!(e.contains(&want), "{}", e);
            let m = dtm(&pts, &DtmOptions { min_spacing: 1e-9, ..Default::default() }).unwrap();
            prop_assert_eq!(m.vertices(), &pts[..n]);
        }
    }

    #[test]
    fn crossing_breaklines_split_with_z_of_the_first() {
        let pts = [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
        ];
        let lines = Polylines::from_parts(
            vec![
                vec![[1.0, 2.0, 10.0], [3.0, 2.0, 20.0]],
                vec![[2.0, 1.0, -5.0], [2.0, 3.0, -5.0]],
            ],
            vec![false, false],
            None,
        )
        .unwrap();
        let m = dtm(
            &pts,
            &DtmOptions {
                breaklines: Some(&lines),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(m.n_vertices(), 9);
        assert_eq!(m.vertices()[8], [2.0, 2.0, 15.0]);
        let es = edges(&m);
        for (a, b) in [(4, 8), (8, 5), (6, 8), (8, 7)] {
            assert!(es.contains_key(&ekey(a, b)));
        }
        assert_up(&m);
    }

    #[test]
    fn boundary_drops_outside_and_holes() {
        let pts = random_points(7, 400, [0.0; 3], [10.0, 10.0, 1.0]);
        let sq = |a: f64, b: f64| vec![[a, a, 0.0], [b, a, 0.0], [b, b, 0.0], [a, b, 0.0]];
        let boundary = Polylines::from_parts(
            vec![sq(1.0, 9.0), sq(4.0, 6.0)],
            vec![true, true],
            Some(vec![0, 0]),
        )
        .unwrap();
        let m = dtm(
            &pts,
            &DtmOptions {
                boundary: Some(&boundary),
                ..Default::default()
            },
        )
        .unwrap();
        let area: f64 = (0..m.n_faces() as u32)
            .map(|f| {
                let c = m.corners(f);
                orient2d(xy(c[0]), xy(c[1]), xy(c[2])) / 2.0
            })
            .sum();
        assert!((area - 60.0).abs() < 1e-9, "area {area}");
    }

    #[test]
    fn a_breakline_vertex_on_a_spot_height_of_another_z_is_refused() {
        let mut pts = vec![[0.0, 0.0, 300.0], [10.0, 0.0, 300.0], [0.0, 10.0, 300.0]];
        let crest = Polylines::from_parts(
            vec![vec![[2.0, 2.0, 305.0], [8.0, 1.0, 305.0]]],
            vec![false],
            None,
        )
        .unwrap();
        let opts = DtmOptions {
            breaklines: Some(&crest),
            ..Default::default()
        };
        pts.push([2.0, 2.0, 305.0]);
        assert_eq!(dtm(&pts, &opts).unwrap().n_vertices(), 5);
        pts[3][2] = 300.0;
        let e = dtm(&pts, &opts).unwrap_err().to_string();
        assert!(
            e.contains("points[3] and breaklines part 0 vertex 0 share the plan position (2, 2)"),
            "{e}"
        );
    }

    #[test]
    fn min_spacing_keeps_the_earlier_point() {
        let pts = [
            [0.0, 0.0, 1.0],
            [0.05, 0.0, 2.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let m = dtm(
            &pts,
            &DtmOptions {
                min_spacing: 0.1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(m.vertices(), &[pts[0], pts[2], pts[3]]);
    }

    #[test]
    fn max_edge_peels_the_long_hull() {
        let mut pts = vec![];
        for j in 0..5 {
            for i in 0..5 {
                pts.push([i as f64, j as f64, 0.0]);
            }
        }
        pts.push([20.0, 2.0, 0.0]);
        let m = dtm(
            &pts,
            &DtmOptions {
                max_edge: Some(2.0),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(m.n_faces(), 32);
    }

    #[test]
    fn collinear_and_non_finite_points_error() {
        let pts = [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 2.0, 0.0]];
        assert!(dtm(&pts, &DtmOptions::default()).is_err());
        assert!(dtm(&[[f64::NAN, 0.0, 0.0]], &DtmOptions::default()).is_err());
    }
}
