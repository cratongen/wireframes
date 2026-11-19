//! Constrained triangulation driven by exact orientation and lexicographic
//! order callbacks, so it runs on implicit points. Validity is decided
//! exactly. Quality comes from Lawson flips that are chosen by exact
//! `incircle` on a fixed `f64` approximation of the points: a flip is applied
//! only when its quad is strictly convex both exactly and in the
//! approximation, so every flip is a legal Delaunay flip of the approximate
//! point set (which guarantees termination) and keeps the exact
//! triangulation valid.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use wireframes_core::predicates::{Sign, incircle, orient2d};
use wireframes_core::{Error, Result};

pub(crate) trait Geom {
    /// Counter-clockwise orientation of points `a, b, c`.
    fn orient(&self, a: u32, b: u32, c: u32) -> Sign;
    /// A total order that is monotone along every line (lexicographic).
    fn cmp(&self, a: u32, b: u32) -> Ordering;
    /// A fixed `f64` approximation of point `a`, used only to choose flips.
    fn approx(&self, a: u32) -> [f64; 2];
}

pub(crate) struct Cdt {
    pub tris: Vec<[u32; 3]>,
    pub constrained: HashSet<(u32, u32)>,
}

fn key(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

struct Tm {
    t: Vec<[u32; 3]>,
    alive: Vec<bool>,
    e: HashMap<(u32, u32), u32>,
    vt: Vec<u32>,
}

impl Tm {
    fn add(&mut self, t: [u32; 3]) {
        let id = self.t.len() as u32;
        for k in 0..3 {
            self.e.insert((t[k], t[(k + 1) % 3]), id);
            self.vt[t[k] as usize] = id;
        }
        self.t.push(t);
        self.alive.push(true);
    }

    fn kill(&mut self, id: u32) {
        let t = self.t[id as usize];
        self.alive[id as usize] = false;
        for k in 0..3 {
            self.e.remove(&(t[k], t[(k + 1) % 3]));
        }
    }

    fn with_edge(&self, a: u32, b: u32) -> Option<u32> {
        self.e.get(&(a, b)).copied()
    }

    fn rot(&self, id: u32, a: u32) -> [u32; 3] {
        let t = self.t[id as usize];
        let k = t.iter().position(|&v| v == a).expect("vertex of triangle");
        [t[k], t[(k + 1) % 3], t[(k + 2) % 3]]
    }

    fn star(&self, a: u32) -> Vec<u32> {
        let start = self.vt[a as usize];
        if start == u32::MAX || !self.alive[start as usize] {
            return vec![];
        }
        let mut out = vec![start];
        let mut cur = start;
        loop {
            let [_, _, y] = self.rot(cur, a);
            match self.with_edge(a, y) {
                Some(n) if n == start => return out,
                Some(n) => {
                    out.push(n);
                    cur = n;
                }
                None => break,
            }
        }
        cur = start;
        loop {
            let [_, x, _] = self.rot(cur, a);
            match self.with_edge(x, a) {
                Some(n) if n != start => {
                    out.push(n);
                    cur = n;
                }
                _ => return out,
            }
        }
    }
}

fn ear_clip<G: Geom>(g: &G, poly: &[u32]) -> Result<Vec<[u32; 3]>> {
    let mut v = poly.to_vec();
    let mut out = Vec::with_capacity(v.len().saturating_sub(2));
    while v.len() > 3 {
        let n = v.len();
        let ear = (0..n).find(|&i| {
            let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
            g.orient(a, b, c) == Sign::Positive
                && v.iter().all(|&p| {
                    p == a
                        || p == b
                        || p == c
                        || g.orient(a, b, p) == Sign::Negative
                        || g.orient(b, c, p) == Sign::Negative
                        || g.orient(c, a, p) == Sign::Negative
                })
        });
        let i = ear.ok_or_else(|| Error::kernel("no ear in a constraint cavity"))?;
        out.push([v[(i + n - 1) % n], v[i], v[(i + 1) % n]]);
        v.remove(i);
    }
    if v.len() == 3 {
        if g.orient(v[0], v[1], v[2]) != Sign::Positive {
            return Err(Error::kernel("flat constraint cavity"));
        }
        out.push([v[0], v[1], v[2]]);
    }
    Ok(out)
}

/// Triangulates the convex hull of points `0..n` (pairwise distinct) and
/// inserts `segments` as chains of edges, splitting them at collinear points.
/// Two constraints whose relative interiors cross give an error.
pub(crate) fn cdt<G: Geom>(g: &G, n: u32, segments: &[[u32; 2]]) -> Result<Cdt> {
    let mut tm = Tm {
        t: vec![],
        alive: vec![],
        e: HashMap::new(),
        vt: vec![u32::MAX; n as usize],
    };
    let mut order: Vec<u32> = (0..n).collect();
    order.sort_by(|&a, &b| g.cmp(a, b).then(a.cmp(&b)));
    if order
        .windows(2)
        .any(|w| g.cmp(w[0], w[1]) == Ordering::Equal)
    {
        return Err(Error::kernel("duplicate point in triangulation"));
    }
    let mut constrained = HashSet::new();
    if n < 3 {
        return Ok(Cdt {
            tris: vec![],
            constrained,
        });
    }
    let mut k = 2;
    while k < order.len() && g.orient(order[0], order[k - 1], order[k]) == Sign::Zero {
        k += 1;
    }
    if k == order.len() {
        return Ok(Cdt {
            tris: vec![],
            constrained,
        });
    }
    let line = &order[..k];
    let p = order[k];
    let mut hull: Vec<u32>;
    if g.orient(line[0], line[k - 1], p) == Sign::Positive {
        for w in line.windows(2) {
            tm.add([w[0], w[1], p]);
        }
        hull = line.to_vec();
    } else {
        for w in line.windows(2) {
            tm.add([w[1], w[0], p]);
        }
        hull = line.iter().rev().copied().collect();
    }
    hull.push(p);
    for &p in &order[k + 1..] {
        let m = hull.len();
        let vis: Vec<bool> = (0..m)
            .map(|i| g.orient(hull[i], hull[(i + 1) % m], p) == Sign::Negative)
            .collect();
        let s = (0..m)
            .find(|&i| vis[i] && !vis[(i + m - 1) % m])
            .ok_or_else(|| Error::kernel("point not outside hull"))?;
        let mut c = 0;
        while vis[(s + c) % m] {
            let i = (s + c) % m;
            tm.add([hull[(i + 1) % m], hull[i], p]);
            c += 1;
        }
        let e = (s + c) % m;
        let mut next = Vec::with_capacity(m - c + 2);
        let mut i = e;
        loop {
            next.push(hull[i]);
            if i == s {
                break;
            }
            i = (i + 1) % m;
        }
        next.push(p);
        hull = next;
    }

    let mut stack: Vec<[u32; 2]> = segments.iter().rev().copied().collect();
    while let Some([a, b]) = stack.pop() {
        if a == b {
            continue;
        }
        if tm.with_edge(a, b).is_some() || tm.with_edge(b, a).is_some() {
            constrained.insert(key(a, b));
            continue;
        }
        let toward = |x: u32| g.cmp(a, x) == g.cmp(a, b);
        let mut start = None;
        for id in tm.star(a) {
            let [_, x, y] = tm.rot(id, a);
            let (ox, oy) = (g.orient(a, x, b), g.orient(a, y, b));
            if ox == Sign::Zero && toward(x) {
                start = Some(Err(x));
                break;
            }
            if oy == Sign::Zero && toward(y) {
                start = Some(Err(y));
                break;
            }
            if ox == Sign::Positive && oy == Sign::Negative {
                start = Some(Ok((id, x, y)));
                break;
            }
        }
        let (t0, mut r, mut l) = match start {
            None => return Err(Error::kernel("constraint leaves the hull")),
            Some(Err(mid)) => {
                constrained.insert(key(a, mid));
                stack.push([mid, b]);
                continue;
            }
            Some(Ok(s)) => s,
        };
        let mut removed = vec![t0];
        let (mut left, mut right) = (vec![l], vec![r]);
        let end = loop {
            if constrained.contains(&key(r, l)) {
                return Err(Error::invalid("constraints cross"));
            }
            let t = tm
                .with_edge(l, r)
                .ok_or_else(|| Error::kernel("constraint leaves the hull"))?;
            removed.push(t);
            let [_, _, z] = tm.rot(t, l);
            if z == b {
                break b;
            }
            match g.orient(a, b, z) {
                Sign::Zero => break z,
                Sign::Positive => {
                    left.push(z);
                    l = z;
                }
                Sign::Negative => {
                    right.push(z);
                    r = z;
                }
            }
        };
        for &t in &removed {
            tm.kill(t);
        }
        let mut lp = vec![a, end];
        lp.extend(left.iter().rev());
        let mut rp = vec![end, a];
        rp.extend(right.iter());
        for t in ear_clip(g, &lp)?.into_iter().chain(ear_clip(g, &rp)?) {
            tm.add(t);
        }
        constrained.insert(key(a, end));
        if end != b {
            stack.push([end, b]);
        }
    }
    let mut stack: Vec<(u32, u32)> = tm.e.keys().filter(|(a, b)| a < b).copied().collect();
    stack.sort_unstable();
    while let Some((a, b)) = stack.pop() {
        if constrained.contains(&key(a, b)) {
            continue;
        }
        let (Some(t1), Some(t2)) = (tm.with_edge(a, b), tm.with_edge(b, a)) else {
            continue;
        };
        let [_, _, c] = tm.rot(t1, a);
        let [_, _, d] = tm.rot(t2, b);
        let ap = [a, b, c, d].map(|i| g.approx(i));
        let pos = |x: usize, y: usize, z: usize| orient2d(ap[x], ap[y], ap[z]) > 0.0;
        if !(pos(0, 1, 2)
            && incircle(ap[0], ap[1], ap[2], ap[3]) > 0.0
            && pos(0, 3, 2)
            && pos(3, 1, 2))
        {
            continue;
        }
        if g.orient(a, d, c) != Sign::Positive || g.orient(d, b, c) != Sign::Positive {
            continue;
        }
        tm.kill(t1);
        tm.kill(t2);
        tm.add([a, d, c]);
        tm.add([d, b, c]);
        stack.extend([key(a, d), key(d, b), key(b, c), key(c, a)]);
    }
    let tris =
        tm.t.iter()
            .zip(&tm.alive)
            .filter(|(_, a)| **a)
            .map(|(t, _)| *t)
            .collect();
    Ok(Cdt { tris, constrained })
}

/// Triangles inside the constraints by even-odd parity, counted from the hull.
pub(crate) fn even_odd(c: &Cdt) -> Vec<[u32; 3]> {
    let mut e: HashMap<(u32, u32), usize> = HashMap::new();
    for (i, t) in c.tris.iter().enumerate() {
        for k in 0..3 {
            e.insert((t[k], t[(k + 1) % 3]), i);
        }
    }
    let mut depth = vec![u8::MAX; c.tris.len()];
    let mut queue = std::collections::VecDeque::new();
    for (i, t) in c.tris.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if !e.contains_key(&(b, a)) && depth[i] == u8::MAX {
                depth[i] = c.constrained.contains(&key(a, b)) as u8;
                queue.push_back(i);
            }
        }
    }
    while let Some(i) = queue.pop_front() {
        let t = c.tris[i];
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if let Some(&j) = e.get(&(b, a))
                && depth[j] == u8::MAX
            {
                depth[j] = depth[i] ^ c.constrained.contains(&key(a, b)) as u8;
                queue.push_back(j);
            }
        }
    }
    c.tris
        .iter()
        .zip(&depth)
        .filter(|(_, d)| **d == 1)
        .map(|(t, _)| *t)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    struct P(Vec<[f64; 2]>);

    impl Geom for P {
        fn orient(&self, a: u32, b: u32, c: u32) -> Sign {
            Sign::of(orient2d(
                self.0[a as usize],
                self.0[b as usize],
                self.0[c as usize],
            ))
        }
        fn cmp(&self, a: u32, b: u32) -> Ordering {
            let (p, q) = (self.0[a as usize], self.0[b as usize]);
            p[0].partial_cmp(&q[0])
                .unwrap()
                .then(p[1].partial_cmp(&q[1]).unwrap())
        }
        fn approx(&self, a: u32) -> [f64; 2] {
            self.0[a as usize]
        }
    }

    fn area(p: &P, t: &[[u32; 3]]) -> f64 {
        t.iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| p.0[i as usize]);
                0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
            })
            .sum()
    }

    fn has_edge(t: &[[u32; 3]], a: u32, b: u32) -> bool {
        t.iter()
            .any(|t| (0..3).any(|k| key(t[k], t[(k + 1) % 3]) == key(a, b)))
    }

    #[test]
    fn square_with_hole_even_odd() {
        let p = P(vec![
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [1.0, 1.0],
            [3.0, 1.0],
            [3.0, 3.0],
            [1.0, 3.0],
        ]);
        let segs = [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 4],
        ];
        let c = cdt(&p, 8, &segs).unwrap();
        assert!((area(&p, &c.tris) - 16.0).abs() < 1e-12);
        let inside = even_odd(&c);
        assert!((area(&p, &inside) - 12.0).abs() < 1e-12);
    }

    #[test]
    fn collinear_points_split_constraints() {
        let p = P(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [1.5, 2.0],
            [1.5, -2.0],
        ]);
        let c = cdt(&p, 6, &[[0, 3]]).unwrap();
        for k in 0..3 {
            assert!(has_edge(&c.tris, k, k + 1));
        }
        assert!(cdt(&p, 6, &[[4, 5], [0, 3]]).is_err());
    }

    proptest! {
        #[test]
        fn triangulates_hull_with_random_constraints(
            pts in proptest::collection::vec((0i32..12, 0i32..12), 3..40),
            segs in proptest::collection::vec((0usize..40, 0usize..40), 0..6),
        ) {
            let mut v: Vec<[f64; 2]> = pts.iter().map(|&(x, y)| [x as f64, y as f64]).collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v.dedup();
            prop_assume!(v.len() >= 3);
            let p = P(v);
            let n = p.0.len() as u32;
            let mut ok_segs: Vec<[u32; 2]> = vec![];
            let crosses = |s: [u32; 2], t: [u32; 2]| {
                let o = |a, b, c| p.orient(a, b, c);
                o(s[0], s[1], t[0]).times(o(s[0], s[1], t[1])) == Sign::Negative
                    && o(t[0], t[1], s[0]).times(o(t[0], t[1], s[1])) == Sign::Negative
            };
            for (a, b) in segs {
                let s = [(a as u32) % n, (b as u32) % n];
                if s[0] != s[1] && ok_segs.iter().all(|&t| !crosses(s, t)) {
                    ok_segs.push(s);
                }
            }
            let c = cdt(&p, n, &ok_segs).unwrap();
            let mut hull_area = 0.0;
            {
                let mut idx: Vec<u32> = (0..n).collect();
                idx.sort_by(|&a, &b| p.cmp(a, b));
                let mut h: Vec<u32> = vec![];
                for pass in 0..2 {
                    let start = h.len();
                    let it: Vec<u32> = if pass == 0 { idx.clone() } else { idx.iter().rev().copied().collect() };
                    for q in it {
                        while h.len() >= start + 2 && p.orient(h[h.len() - 2], h[h.len() - 1], q) != Sign::Positive { h.pop(); }
                        h.push(q);
                    }
                    h.pop();
                }
                for i in 1..h.len().saturating_sub(1) {
                    hull_area += area(&p, &[[h[0], h[i], h[i + 1]]]);
                }
            }
            prop_assert!((area(&p, &c.tris) - hull_area).abs() < 1e-9);
            if hull_area == 0.0 {
                prop_assert!(c.tris.is_empty());
                return Ok(());
            }
            for t in &c.tris {
                prop_assert_eq!(p.orient(t[0], t[1], t[2]), Sign::Positive);
            }
            let mut by_edge: HashMap<(u32, u32), u32> = HashMap::new();
            for t in &c.tris {
                for k in 0..3 {
                    by_edge.insert((t[k], t[(k + 1) % 3]), t[(k + 2) % 3]);
                }
            }
            for (&(a, b), &x) in &by_edge {
                if let Some(&y) = by_edge.get(&(b, a))
                    && !c.constrained.contains(&key(a, b))
                {
                    let q = [a, b, x, y].map(|i| p.0[i as usize]);
                    let convex = orient2d(q[0], q[3], q[2]) > 0.0 && orient2d(q[3], q[1], q[2]) > 0.0;
                    prop_assert!(!convex || incircle(q[0], q[1], q[2], q[3]) <= 0.0, "edge {:?} is not locally Delaunay", (a, b));
                }
            }
            let used: HashSet<u32> = c.tris.iter().flatten().copied().collect();
            prop_assert_eq!(used.len(), n as usize);
            for s in ok_segs {
                let mut chain: Vec<u32> = (0..n).filter(|&q| q == s[0] || q == s[1] || (p.orient(s[0], s[1], q) == Sign::Zero && p.cmp(s[0], q) == p.cmp(q, s[1]))).collect();
                chain.sort_by(|&a, &b| p.cmp(a, b));
                for w in chain.windows(2) {
                    prop_assert!(has_edge(&c.tris, w[0], w[1]));
                }
            }
        }
    }
}
