//! Bounding volume hierarchy over item boxes. Box tests use floating comparisons
//! and only prune; callers make every final decision with exact predicates.

use crate::Aabb;

const LEAF: usize = 4;

#[derive(Clone, Debug)]
struct Node {
    bounds: Aabb,
    /// Leaf: `start..end` into `items`. Inner: `start` is the left child, `end`
    /// the right child.
    start: u32,
    end: u32,
    leaf: bool,
}

#[derive(Clone, Debug, Default)]
pub struct Bvh {
    nodes: Vec<Node>,
    items: Vec<u32>,
    boxes: Vec<Aabb>,
}

impl Bvh {
    /// Median split on the longest axis of the item centroids, leaves of up to
    /// four items, ties ordered by item index. Empty boxes are left out.
    pub fn new(boxes: &[Aabb]) -> Bvh {
        let mut items: Vec<u32> = (0..boxes.len() as u32)
            .filter(|&i| !boxes[i as usize].is_empty())
            .collect();
        let centers: Vec<[f64; 3]> = boxes.iter().map(|b| b.center()).collect();
        let mut nodes = vec![];
        if !items.is_empty() {
            let n = items.len();
            build(&mut nodes, &mut items, 0, n, boxes, &centers);
        }
        Bvh {
            nodes,
            items,
            boxes: boxes.to_vec(),
        }
    }

    /// Number of indexed items.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn item_box(&self, i: u32) -> &Aabb {
        &self.boxes[i as usize]
    }

    fn walk(&self, mut enter: impl FnMut(&Aabb) -> bool, mut f: impl FnMut(u32)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = vec![0u32];
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n as usize];
            if !enter(&node.bounds) {
                continue;
            }
            if node.leaf {
                for &i in &self.items[node.start as usize..node.end as usize] {
                    if enter(&self.boxes[i as usize]) {
                        f(i);
                    }
                }
            } else {
                stack.push(node.end);
                stack.push(node.start);
            }
        }
    }

    /// Items whose box meets `b`.
    pub fn query(&self, b: &Aabb, f: impl FnMut(u32)) {
        self.walk(|x| x.intersects(b), f);
    }

    /// Pairs `(i, j)` of an item of `self` and an item of `other` whose boxes meet.
    pub fn pairs(&self, other: &Bvh, mut f: impl FnMut(u32, u32)) {
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return;
        }
        let mut stack = vec![(0u32, 0u32)];
        while let Some((a, b)) = stack.pop() {
            let (na, nb) = (&self.nodes[a as usize], &other.nodes[b as usize]);
            if !na.bounds.intersects(&nb.bounds) {
                continue;
            }
            match (na.leaf, nb.leaf) {
                (true, true) => {
                    for &i in &self.items[na.start as usize..na.end as usize] {
                        for &j in &other.items[nb.start as usize..nb.end as usize] {
                            if self.boxes[i as usize].intersects(&other.boxes[j as usize]) {
                                f(i, j);
                            }
                        }
                    }
                }
                (false, true) => {
                    stack.push((na.end, b));
                    stack.push((na.start, b));
                }
                (true, false) => {
                    stack.push((a, nb.end));
                    stack.push((a, nb.start));
                }
                (false, false) => {
                    for (x, y) in [
                        (na.end, nb.end),
                        (na.end, nb.start),
                        (na.start, nb.end),
                        (na.start, nb.start),
                    ] {
                        stack.push((x, y));
                    }
                }
            }
        }
    }

    /// Pairs `(i, j)` with `i < j` of items whose boxes meet, each pair once.
    pub fn self_pairs(&self, mut f: impl FnMut(u32, u32)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = vec![(0u32, 0u32)];
        while let Some((a, b)) = stack.pop() {
            let (na, nb) = (&self.nodes[a as usize], &self.nodes[b as usize]);
            if !na.bounds.intersects(&nb.bounds) {
                continue;
            }
            if a == b {
                if na.leaf {
                    let items = &self.items[na.start as usize..na.end as usize];
                    for (k, &i) in items.iter().enumerate() {
                        for &j in &items[k + 1..] {
                            if self.boxes[i as usize].intersects(&self.boxes[j as usize]) {
                                f(i.min(j), i.max(j));
                            }
                        }
                    }
                } else {
                    stack.push((na.end, na.end));
                    stack.push((na.start, na.end));
                    stack.push((na.start, na.start));
                }
                continue;
            }
            match (na.leaf, nb.leaf) {
                (true, true) => {
                    for &i in &self.items[na.start as usize..na.end as usize] {
                        for &j in &self.items[nb.start as usize..nb.end as usize] {
                            if self.boxes[i as usize].intersects(&self.boxes[j as usize]) {
                                f(i.min(j), i.max(j));
                            }
                        }
                    }
                }
                (false, _) => {
                    stack.push((na.end, b));
                    stack.push((na.start, b));
                }
                (true, false) => {
                    stack.push((a, nb.end));
                    stack.push((a, nb.start));
                }
            }
        }
    }

    /// Items whose box the ray `origin + t * dir`, `0 <= t <= t_max`, may meet.
    /// Boxes are padded by a few ulps so rounding never drops a hit.
    pub fn ray(&self, origin: [f64; 3], dir: [f64; 3], t_max: f64, f: impl FnMut(u32)) {
        let hit = |b: &Aabb| {
            let (mut lo, mut hi) = (0.0f64, t_max);
            for k in 0..3 {
                let pad = 8.0
                    * f64::EPSILON
                    * (b.min[k].abs().max(b.max[k].abs()).max(origin[k].abs()) + 1e-300);
                let (bmin, bmax) = (b.min[k] - pad, b.max[k] + pad);
                if dir[k] == 0.0 {
                    if origin[k] < bmin || origin[k] > bmax {
                        return false;
                    }
                    continue;
                }
                let (t0, t1) = ((bmin - origin[k]) / dir[k], (bmax - origin[k]) / dir[k]);
                let (t0, t1) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
                lo = lo.max(t0 * (1.0 - 4.0 * f64::EPSILON) - f64::MIN_POSITIVE);
                hi = hi.min(t1 * (1.0 + 4.0 * f64::EPSILON) + f64::MIN_POSITIVE);
                if lo > hi {
                    return false;
                }
            }
            true
        };
        self.walk(hit, f);
    }

    /// The item with the smallest `dist2(item)`, smaller index on ties, visiting
    /// only items whose box could hold a closer point.
    pub fn nearest(&self, p: [f64; 3], mut dist2: impl FnMut(u32) -> f64) -> Option<(u32, f64)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<(u32, f64)> = None;
        let mut stack = vec![(0u32, self.nodes[0].bounds.distance2(p))];
        while let Some((n, d)) = stack.pop() {
            if best.is_some_and(|(_, b)| d > b) {
                continue;
            }
            let node = &self.nodes[n as usize];
            if node.leaf {
                for &i in &self.items[node.start as usize..node.end as usize] {
                    if best.is_some_and(|(_, b)| self.boxes[i as usize].distance2(p) > b) {
                        continue;
                    }
                    let di = dist2(i);
                    best = match best {
                        Some((j, b)) if b < di || (b == di && j < i) => Some((j, b)),
                        _ => Some((i, di)),
                    };
                }
            } else {
                let (l, r) = (node.start, node.end);
                let (dl, dr) = (
                    self.nodes[l as usize].bounds.distance2(p),
                    self.nodes[r as usize].bounds.distance2(p),
                );
                if dl <= dr {
                    stack.push((r, dr));
                    stack.push((l, dl));
                } else {
                    stack.push((l, dl));
                    stack.push((r, dr));
                }
            }
        }
        best
    }
}

fn build(
    nodes: &mut Vec<Node>,
    items: &mut [u32],
    start: usize,
    end: usize,
    boxes: &[Aabb],
    centers: &[[f64; 3]],
) -> u32 {
    let slice = &mut items[start..end];
    let bounds = slice
        .iter()
        .fold(Aabb::empty(), |b, &i| b.union(&boxes[i as usize]));
    let id = nodes.len() as u32;
    nodes.push(Node {
        bounds,
        start: start as u32,
        end: end as u32,
        leaf: true,
    });
    if slice.len() <= LEAF {
        return id;
    }
    let cb = slice.iter().fold(Aabb::empty(), |b, &i| {
        let c = centers[i as usize];
        b.union(&Aabb { min: c, max: c })
    });
    let axis = (0..3)
        .max_by(|&a, &b| {
            (cb.max[a] - cb.min[a])
                .total_cmp(&(cb.max[b] - cb.min[b]))
                .then(b.cmp(&a))
        })
        .expect("three axes");
    let mid = slice.len() / 2;
    slice.select_nth_unstable_by(mid, |&a, &b| {
        centers[a as usize][axis]
            .total_cmp(&centers[b as usize][axis])
            .then(a.cmp(&b))
    });
    let left = build(nodes, items, start, start + mid, boxes, centers);
    let right = build(nodes, items, start + mid, end, boxes, centers);
    nodes[id as usize] = Node {
        bounds,
        start: left,
        end: right,
        leaf: false,
    };
    id
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn boxes() -> impl Strategy<Value = Vec<Aabb>> {
        proptest::collection::vec(
            (
                proptest::array::uniform3(-10.0..10.0f64),
                proptest::array::uniform3(0.0..3.0f64),
            ),
            0..60,
        )
        .prop_map(|v| {
            v.into_iter()
                .map(|(c, s)| Aabb {
                    min: c,
                    max: [c[0] + s[0], c[1] + s[1], c[2] + s[2]],
                })
                .collect()
        })
    }

    #[test]
    fn empty_tree_answers_nothing() {
        let t = Bvh::new(&[]);
        assert!(t.is_empty());
        t.query(
            &Aabb {
                min: [0.0; 3],
                max: [1.0; 3],
            },
            |_| panic!(),
        );
        assert_eq!(t.nearest([0.0; 3], |_| 0.0), None);
    }

    proptest! {
        #[test]
        fn query_matches_brute_force(bs in boxes(), q in boxes()) {
            let t = Bvh::new(&bs);
            for b in &q {
                let mut got = vec![];
                t.query(b, |i| got.push(i));
                got.sort_unstable();
                let want: Vec<u32> = (0..bs.len() as u32).filter(|&i| bs[i as usize].intersects(b)).collect();
                prop_assert_eq!(got, want);
            }
        }

        #[test]
        fn pairs_match_brute_force(a in boxes(), b in boxes()) {
            let (ta, tb) = (Bvh::new(&a), Bvh::new(&b));
            let mut got = vec![];
            ta.pairs(&tb, |i, j| got.push((i, j)));
            got.sort_unstable();
            let mut want = vec![];
            for i in 0..a.len() {
                for j in 0..b.len() {
                    if a[i].intersects(&b[j]) {
                        want.push((i as u32, j as u32));
                    }
                }
            }
            prop_assert_eq!(got, want);
            let mut got = vec![];
            ta.self_pairs(|i, j| got.push((i, j)));
            got.sort_unstable();
            let mut want = vec![];
            for i in 0..a.len() {
                for j in i + 1..a.len() {
                    if a[i].intersects(&a[j]) {
                        want.push((i as u32, j as u32));
                    }
                }
            }
            prop_assert_eq!(got, want);
        }

        #[test]
        fn nearest_matches_brute_force(bs in boxes(), p in proptest::array::uniform3(-12.0..12.0f64)) {
            let t = Bvh::new(&bs);
            let got = t.nearest(p, |i| bs[i as usize].distance2(p));
            let want = (0..bs.len() as u32)
                .map(|i| (i, bs[i as usize].distance2(p)))
                .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            prop_assert_eq!(got, want);
        }

        #[test]
        fn ray_never_misses_a_box_it_crosses(bs in boxes(), o in proptest::array::uniform3(-12.0..12.0f64), d in proptest::array::uniform3(-1.0..1.0f64)) {
            let t = Bvh::new(&bs);
            let mut got = vec![];
            t.ray(o, d, f64::INFINITY, |i| got.push(i));
            for (i, b) in bs.iter().enumerate() {
                let inside_at = |s: f64| b.contains([o[0] + s * d[0], o[1] + s * d[1], o[2] + s * d[2]]);
                if (0..200).any(|k| inside_at(k as f64 * 0.25)) {
                    prop_assert!(got.contains(&(i as u32)));
                }
            }
        }
    }
}
