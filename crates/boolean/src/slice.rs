use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use rayon::prelude::*;
use wireframes_core::predicates::{Sign, orient2d_axis};
use wireframes_core::topology::Shells;
use wireframes_core::validate::degenerate_mask;
use wireframes_core::{Column, EdgeMap, Error, Mesh, Polylines, Result, Table};

use crate::hp::{self, Pt};
use crate::ops::{Assembled, collapse_rounded, is_closed_result};
use crate::triangulate::{Geom, cdt, even_odd};

/// A section point: an input vertex on the plane, or the exact crossing of an
/// edge whose endpoints lie strictly on opposite sides.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Key {
    V(u32),
    E(u32, u32),
}

/// The plane `(p - origin) . normal = 0`, perturbed toward its negative side
/// so that a vertex on the plane counts as positive.
struct Section {
    negative: Vec<bool>,
    keys: Vec<Key>,
    index: HashMap<Key, u32>,
    /// Per face crossed, the directed segment (inside on its left seen from
    /// the normal side when the face is outward) and the face.
    segs: Vec<(u32, u32, u32)>,
}

fn check_plane(origin: [f64; 3], normal: [f64; 3]) -> Result<()> {
    if origin.iter().chain(&normal).any(|x| !x.is_finite()) {
        return Err(Error::invalid("origin and normal must be finite"));
    }
    if normal.iter().all(|x| *x == 0.0) {
        return Err(Error::invalid("normal must be nonzero"));
    }
    Ok(())
}

fn crossing_key(negative: &[bool], on: &[bool], u: u32, v: u32) -> Key {
    let p = if negative[u as usize] { v } else { u };
    if on[p as usize] {
        Key::V(p)
    } else {
        Key::E(u.min(v), u.max(v))
    }
}

fn section(mesh: &Mesh, origin: [f64; 3], normal: [f64; 3]) -> Section {
    let sides: Vec<Sign> = mesh
        .vertices()
        .par_iter()
        .map(|&p| hp::plane_side(p, origin, normal))
        .collect();
    let negative: Vec<bool> = sides.iter().map(|s| *s == Sign::Negative).collect();
    let on: Vec<bool> = sides.iter().map(|s| *s == Sign::Zero).collect();
    let mut sec = Section {
        negative,
        keys: vec![],
        index: HashMap::new(),
        segs: vec![],
    };
    for (f, t) in mesh.triangles().iter().enumerate() {
        let neg = t.map(|v| sec.negative[v as usize]);
        if neg.iter().all(|&x| x) || neg.iter().all(|&x| !x) {
            continue;
        }
        let mut into = None;
        let mut out = None;
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            match (neg[k], neg[(k + 1) % 3]) {
                (false, true) => into = Some(crossing_key(&sec.negative, &on, a, b)),
                (true, false) => out = Some(crossing_key(&sec.negative, &on, a, b)),
                _ => {}
            }
        }
        let (i, o) = (into.expect("mixed face"), out.expect("mixed face"));
        if i == o {
            continue;
        }
        let mut id = |k: Key| {
            let n = sec.keys.len() as u32;
            *sec.index.entry(k).or_insert_with(|| {
                sec.keys.push(k);
                n
            })
        };
        let (a, b) = (id(i), id(o));
        sec.segs.push((a, b, f as u32));
    }
    sec
}

fn key_point(mesh: &Mesh, k: Key, origin: [f64; 3], normal: [f64; 3]) -> Pt {
    match k {
        Key::V(v) => Pt::E(mesh.vertices()[v as usize]),
        Key::E(u, v) => Pt::plane_cross(
            mesh.vertices()[u as usize],
            mesh.vertices()[v as usize],
            origin,
            normal,
        ),
    }
}

struct CapGeom<'a> {
    pts: &'a [Pt],
    drop: usize,
}

impl Geom for CapGeom<'_> {
    fn orient(&self, a: u32, b: u32, c: u32) -> Sign {
        hp::orient2d(
            &self.pts[a as usize],
            &self.pts[b as usize],
            &self.pts[c as usize],
            self.drop,
        )
    }

    fn cmp(&self, a: u32, b: u32) -> Ordering {
        hp::cmp_proj(&self.pts[a as usize], &self.pts[b as usize], self.drop)
    }

    fn approx(&self, a: u32) -> [f64; 2] {
        let (i, j) = wireframes_core::predicates::kept_axes(self.drop);
        let p = self.pts[a as usize].approx();
        [p[i], p[j]]
    }
}

fn dominant(normal: [f64; 3]) -> usize {
    (0..3)
        .max_by(|&a, &b| normal[a].abs().total_cmp(&normal[b].abs()).then(b.cmp(&a)))
        .expect("three axes")
}

/// For each section point, the representative of the points whose rounded
/// coordinates project onto the same point along `drop`: one that is an
/// input vertex (`snap`) when there is one, else the smallest index.
///
/// Rounded section points lie within a few ulps of the plane, and `drop` is
/// the dominant axis of its normal, so points with one projection are a few
/// ulps apart in 3D. They must be one cap point because the cap is
/// triangulated in that projection.
fn coincident(pts: &[Pt], snap: &[Option<u32>], drop: usize) -> Vec<u32> {
    let edge = |a: u32| snap[a as usize].is_none();
    let mut order: Vec<u32> = (0..pts.len() as u32).collect();
    order.sort_by(|&a, &b| {
        hp::cmp_proj(&pts[a as usize], &pts[b as usize], drop)
            .then(edge(a).cmp(&edge(b)))
            .then(a.cmp(&b))
    });
    let mut same: Vec<u32> = (0..pts.len() as u32).collect();
    for w in order.windows(2) {
        if hp::cmp_proj(&pts[w[0] as usize], &pts[w[1] as usize], drop) == Ordering::Equal {
            same[w[1] as usize] = same[w[0] as usize];
        }
    }
    same
}

/// `list` after the cleanup booleans apply to rounded output
/// ([`collapse_rounded`]): vertices with equal coordinates become one, faces
/// that repeat a vertex or cancel go, and collinear corners and folds are
/// split away. Also returns `merged` with the vertices that absorbed others,
/// which may now be shared by touching shells.
fn clean(
    list: &[([u32; 3], Option<u32>)],
    coords: &dyn Fn(u32) -> [f64; 3],
    nv: u32,
    merged: &HashSet<u32>,
) -> (Vec<([u32; 3], Option<u32>)>, HashSet<u32>) {
    let mut ids: Vec<u32> = list.iter().flat_map(|(t, _)| *t).collect();
    ids.sort_unstable();
    ids.dedup();
    let at = |g: u32| ids.binary_search(&g).expect("listed") as u32;
    let asm = Assembled {
        vertices: ids.iter().map(|&g| coords(g)).collect(),
        triangles: list.iter().map(|(t, _)| t.map(at)).collect(),
        vertex_rows: ids.iter().map(|&g| (g < nv).then_some(g)).collect(),
        face_source: list
            .iter()
            .map(|(_, f)| f.map_or((1, 0), |f| (0, f)))
            .collect(),
        normal: vec![],
    };
    let (mut out, mut point) = collapse_rounded(asm, &ids);
    loop {
        let n = out.triangles.len();
        (out, point) = collapse_rounded(out, &point);
        if out.triangles.len() == n {
            break;
        }
    }
    let bits = |g: u32| coords(g).map(|x| (x + 0.0).to_bits());
    let kept: HashMap<[u64; 3], u32> = point.iter().map(|&r| (bits(r), r)).collect();
    let mut merged = merged.clone();
    merged.extend(
        ids.iter()
            .filter_map(|&g| kept.get(&bits(g)).filter(|&&r| r != g))
            .copied(),
    );
    let list = out
        .triangles
        .iter()
        .zip(&out.face_source)
        .map(|(t, &(m, f))| (t.map(|i| point[i as usize]), (m == 0).then_some(f)))
        .collect();
    (list, merged)
}

/// Shell label of each face of `list`, the faces joined as [`Shells`] joins
/// them.
fn shells_of(
    list: &[([u32; 3], Option<u32>)],
    coords: &dyn Fn(u32) -> [f64; 3],
) -> Result<Vec<u32>> {
    let mut ids: Vec<u32> = list.iter().flat_map(|(t, _)| *t).collect();
    ids.sort_unstable();
    ids.dedup();
    let at = |g: u32| ids.binary_search(&g).expect("listed") as u32;
    let m = Mesh::new(
        ids.iter().map(|&g| coords(g)).collect(),
        list.iter().map(|(t, _)| t.map(at)).collect(),
    )?;
    let skip = degenerate_mask(&m);
    Ok(Shells::new(&m, &EdgeMap::new(&m, Some(&skip))).of)
}

/// `(negative, positive)`; points on the plane go to positive.
/// Closed meshes are capped with a triangulation of the section; cap faces
/// have null attributes. Faces lying in the plane go to the positive side,
/// except faces of a closed mesh whose normal points along `normal`: those
/// bound material on the negative side only, so they go to the negative
/// side and are cut out of the cap.
pub fn split_by_plane(mesh: &Mesh, origin: [f64; 3], normal: [f64; 3]) -> Result<(Mesh, Mesh)> {
    check_plane(origin, normal)?;
    let normal = hp::unit_scale(normal);
    let mut sec = section(mesh, origin, normal);
    let drop = dominant(normal);
    let closed = mesh.is_closed();
    let on: Vec<bool> = mesh
        .vertices()
        .par_iter()
        .map(|&p| hp::plane_side(p, origin, normal) == Sign::Zero)
        .collect();
    let roof: Vec<bool> = mesh
        .triangles()
        .iter()
        .map(|t| {
            closed
                && t.iter().all(|&v| on[v as usize])
                && Sign::of(orient2d_axis(
                    mesh.vertices()[t[0] as usize],
                    mesh.vertices()[t[1] as usize],
                    mesh.vertices()[t[2] as usize],
                    drop,
                )) == Sign::of(normal[drop])
        })
        .collect();
    let mut cap_segs: Vec<[u32; 2]> = sec.segs.iter().map(|s| [s.0, s.1]).collect();
    for (t, _) in mesh.triangles().iter().zip(&roof).filter(|(_, r)| **r) {
        for k in 0..3 {
            let mut id = |v: u32| {
                let n = sec.keys.len() as u32;
                *sec.index.entry(Key::V(v)).or_insert_with(|| {
                    sec.keys.push(Key::V(v));
                    n
                })
            };
            cap_segs.push([id(t[k]), id(t[(k + 1) % 3])]);
        }
    }
    let mut touch_segs: Vec<[u32; 2]> = vec![];
    if closed {
        let mut touch = on.clone();
        for t in mesh.triangles() {
            if t.iter().any(|&v| sec.negative[v as usize]) || t.iter().all(|&v| on[v as usize]) {
                for &v in t {
                    touch[v as usize] = false;
                }
            }
        }
        let mut id = |v: u32| {
            let n = sec.keys.len() as u32;
            *sec.index.entry(Key::V(v)).or_insert_with(|| {
                sec.keys.push(Key::V(v));
                n
            })
        };
        for v in 0..mesh.n_vertices() as u32 {
            if touch[v as usize] {
                id(v);
            }
        }
        for t in mesh.triangles() {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if a < b && touch[a as usize] && touch[b as usize] {
                    touch_segs.push([id(a), id(b)]);
                }
            }
        }
    }
    let pts: Vec<Pt> = sec
        .keys
        .par_iter()
        .map(|&k| Pt::E(key_point(mesh, k, origin, normal).to_f64()))
        .collect();
    let snap: Vec<Option<u32>> = sec
        .keys
        .iter()
        .zip(&pts)
        .map(|(&k, p)| match k {
            Key::V(v) => Some(v),
            Key::E(u, v) => {
                let p = p.to_f64();
                [u, v]
                    .into_iter()
                    .find(|&w| mesh.vertices()[w as usize] == p)
            }
        })
        .collect();
    let same = if closed {
        coincident(&pts, &snap, drop)
    } else {
        (0..pts.len() as u32).collect()
    };
    for s in cap_segs.iter_mut() {
        *s = s.map(|i| same[i as usize]);
    }
    cap_segs.retain(|s| s[0] != s[1]);
    let mut count: HashMap<(u32, u32), usize> = HashMap::new();
    for s in &cap_segs {
        *count.entry((s[0].min(s[1]), s[0].max(s[1]))).or_default() += 1;
    }
    cap_segs.retain(|s| count[&(s[0].min(s[1]), s[0].max(s[1]))] % 2 == 1);
    cap_segs.sort_unstable_by_key(|s| (s[0].min(s[1]), s[0].max(s[1])));
    cap_segs.dedup_by_key(|s| (s[0].min(s[1]), s[0].max(s[1])));
    let nv = mesh.n_vertices() as u32;
    let own = |i: u32| snap[i as usize].unwrap_or(nv + i);
    let gid = |k: Key| own(same[sec.index[&k] as usize]);
    let mut vmap: Vec<u32> = (0..nv).collect();
    for (i, k) in sec.keys.iter().enumerate() {
        if let Key::V(v) = *k {
            vmap[v as usize] = own(same[i]);
        }
    }
    let merged: HashSet<u32> = (0..pts.len() as u32)
        .filter(|&i| same[i as usize] != i)
        .flat_map(|i| [own(i), own(same[i as usize])])
        .map(|g| if g < nv { vmap[g as usize] } else { g })
        .collect();
    let mut faces: [Vec<([u32; 3], Option<u32>)>; 2] = [vec![], vec![]];
    for (f, input) in mesh.triangles().iter().enumerate() {
        let neg = input.map(|v| sec.negative[v as usize]);
        let t = &input.map(|v| vmap[v as usize]);
        if roof[f] {
            faces[0].push((*t, Some(f as u32)));
            continue;
        }
        if neg.iter().all(|&x| x) || neg.iter().all(|&x| !x) {
            faces[usize::from(!neg[0])].push((*t, Some(f as u32)));
            continue;
        }
        let mut poly: [Vec<u32>; 2] = [vec![], vec![]];
        for k in 0..3 {
            poly[usize::from(!neg[k])].push(t[k]);
            if neg[k] != neg[(k + 1) % 3] {
                let x = match crossing_key(&sec.negative, &on, input[k], input[(k + 1) % 3]) {
                    Key::V(v) => vmap[v as usize],
                    e => gid(e),
                };
                poly[0].push(x);
                poly[1].push(x);
            }
        }
        for (side, mut p) in poly.into_iter().enumerate() {
            p.dedup();
            while p.len() > 1 && p[0] == p[p.len() - 1] {
                p.pop();
            }
            for i in 1..p.len().saturating_sub(1) {
                faces[side].push(([p[0], p[i], p[i + 1]], Some(f as u32)));
            }
        }
    }
    if closed && !cap_segs.is_empty() {
        let reps: Vec<u32> = (0..pts.len() as u32)
            .filter(|&i| same[i as usize] == i)
            .collect();
        let local: HashMap<u32, u32> = reps
            .iter()
            .enumerate()
            .map(|(k, &i)| (i, k as u32))
            .collect();
        let upts: Vec<Pt> = reps.iter().map(|&i| pts[i as usize].clone()).collect();
        let edge = |s: &[u32; 2]| {
            let s = s.map(|i| local[&same[i as usize]]);
            (s[0].min(s[1]), s[0].max(s[1]))
        };
        let parity: HashSet<(u32, u32)> = cap_segs.iter().map(edge).collect();
        let mut touching: Vec<(u32, u32)> = touch_segs
            .iter()
            .map(edge)
            .filter(|e| e.0 != e.1 && !parity.contains(e))
            .collect();
        touching.sort_unstable();
        touching.dedup();
        let usegs: Vec<[u32; 2]> = cap_segs
            .iter()
            .map(|s| s.map(|i| local[&i]))
            .chain(touching.iter().map(|e| [e.0, e.1]))
            .collect();
        let mut c = cdt(&CapGeom { pts: &upts, drop }, upts.len() as u32, &usegs)
            .map_err(|e| Error::invalid(format!("cannot cap the section: {e}")))?;
        for e in &touching {
            c.constrained.remove(e);
        }
        let up = Sign::of(normal[drop]) == Sign::Positive;
        for t in even_odd(&c) {
            let g = t.map(|i| own(reps[i as usize]));
            let g = g.map(|x| if x < nv { vmap[x as usize] } else { x });
            let (n, p) = if up {
                (g, [g[0], g[2], g[1]])
            } else {
                ([g[0], g[2], g[1]], g)
            };
            faces[0].push((n, None));
            faces[1].push((p, None));
        }
    }
    let coords = |id: u32| {
        if id < nv {
            mesh.vertices()[id as usize]
        } else {
            pts[(id - nv) as usize].to_f64()
        }
    };
    let build = |list: &[([u32; 3], Option<u32>)]| -> Result<Mesh> {
        let (list, merged) = clean(list, &coords, nv, &merged);
        let list = &list[..];
        let shell = if merged.is_empty() {
            vec![0; list.len()]
        } else {
            shells_of(list, &coords)?
        };
        let mut map: HashMap<(u32, u32), u32> = HashMap::new();
        let (mut v, mut rows, mut tris) = (vec![], vec![], vec![]);
        for (f, (t, _)) in list.iter().enumerate() {
            tris.push(t.map(|id| {
                let part = if merged.contains(&id) { shell[f] } else { 0 };
                *map.entry((id, part)).or_insert_with(|| {
                    v.push(coords(id));
                    rows.push((id < nv).then_some(id));
                    v.len() as u32 - 1
                })
            }));
        }
        let frows: Vec<Option<u32>> = list.iter().map(|x| x.1).collect();
        let m = if tris.is_empty() {
            Mesh::empty().with_attributes(
                mesh.vertex_attributes().take(&[]),
                mesh.face_attributes().take(&[]),
            )?
        } else {
            Mesh::new(v, tris)?.with_attributes(
                mesh.vertex_attributes().take_or_null(&rows),
                mesh.face_attributes().take_or_null(&frows),
            )?
        };
        Ok(m.with_crs(mesh.crs.clone())
            .with_length_unit(mesh.length_unit.clone()))
    };
    let parts = (build(&faces[0])?, build(&faces[1])?);
    if closed && !(is_closed_result(&parts.0)? && is_closed_result(&parts.1)?) {
        return Err(Error::invalid(
            "internal: rounding the split to f64 left it open or non-manifold",
        ));
    }
    Ok(parts)
}

/// Section lines at each offset along `normal`; empty `offsets` means `[0.0]`.
/// The plane for offset `t` passes through `origin + t * normal / |normal|`
/// (rounded to f64 once); every side decision on that plane is exact.
pub fn slice_by_plane(
    mesh: &Mesh,
    origin: [f64; 3],
    normal: [f64; 3],
    offsets: &[f64],
) -> Result<Polylines> {
    check_plane(origin, normal)?;
    if mesh.face_attributes().column("offset").is_some() {
        return Err(Error::invalid(
            "mesh has a face column named offset, which slice_by_plane writes",
        ));
    }
    if offsets.iter().any(|x| !x.is_finite()) {
        return Err(Error::invalid("offsets must be finite"));
    }
    let offsets = if offsets.is_empty() {
        &[0.0][..]
    } else {
        offsets
    };
    let unit = hp::unit_scale(normal);
    let len = unit.iter().map(|x| x * x).sum::<f64>().sqrt();
    let ft = mesh.face_attributes();
    let same_attr =
        |f: u32, g: u32| f == g || ft.num_columns() == 0 || ft.take(&[f]) == ft.take(&[g]);
    let mut parts: Vec<Vec<[f64; 3]>> = vec![];
    let mut closed = vec![];
    let mut feat_offset = vec![];
    let mut feat_face = vec![];
    for &t in offsets {
        let o = [0, 1, 2].map(|k| origin[k] + t * unit[k] / len);
        let sec = section(mesh, o, unit);
        let pts: Vec<[f64; 3]> = sec
            .keys
            .par_iter()
            .map(|&k| key_point(mesh, k, o, unit).to_f64())
            .collect();
        let mut next: HashMap<u32, Vec<usize>> = HashMap::new();
        let mut has_in = vec![false; sec.keys.len()];
        for (i, s) in sec.segs.iter().enumerate() {
            next.entry(s.0).or_default().push(i);
            has_in[s.1 as usize] = true;
        }
        let mut used = vec![false; sec.segs.len()];
        let mut chains: Vec<(Vec<usize>, bool)> = vec![];
        let starts: Vec<usize> = (0..sec.segs.len())
            .filter(|&i| !has_in[sec.segs[i].0 as usize])
            .chain(0..sec.segs.len())
            .collect();
        for s in starts {
            if used[s] {
                continue;
            }
            let mut chain = vec![s];
            used[s] = true;
            let first = sec.segs[s].0;
            let mut cur = sec.segs[s].1;
            let mut ring = false;
            loop {
                if cur == first {
                    ring = true;
                    break;
                }
                let Some(n) = next
                    .get(&cur)
                    .and_then(|l| l.iter().copied().find(|&i| !used[i]))
                else {
                    break;
                };
                used[n] = true;
                chain.push(n);
                cur = sec.segs[n].1;
            }
            chains.push((chain, ring));
        }
        for (mut chain, ring) in chains {
            let face = |i: usize| sec.segs[chain[i]].2;
            let breaks: Vec<usize> = (1..chain.len())
                .filter(|&i| !same_attr(face(i - 1), face(i)))
                .collect();
            let wrap_break = ring && !same_attr(face(chain.len() - 1), face(0));
            if ring && breaks.is_empty() && !wrap_break {
                let ring_pts: Vec<[f64; 3]> =
                    chain.iter().map(|&i| pts[sec.segs[i].0 as usize]).collect();
                if ring_pts.len() >= 3 {
                    parts.push(ring_pts);
                    closed.push(true);
                    feat_offset.push(t);
                    feat_face.push(sec.segs[chain[0]].2);
                }
                continue;
            }
            if ring && !wrap_break {
                let b = breaks[0];
                chain.rotate_left(b);
            }
            let mut start = 0;
            for i in 1..=chain.len() {
                let cut =
                    i == chain.len() || !same_attr(sec.segs[chain[i - 1]].2, sec.segs[chain[i]].2);
                if cut {
                    let run = &chain[start..i];
                    let mut line: Vec<[f64; 3]> =
                        run.iter().map(|&s| pts[sec.segs[s].0 as usize]).collect();
                    line.push(pts[sec.segs[run[run.len() - 1]].1 as usize]);
                    parts.push(line);
                    closed.push(false);
                    feat_offset.push(t);
                    feat_face.push(sec.segs[run[0]].2);
                    start = i;
                }
            }
        }
    }
    let n = parts.len();
    let pl = Polylines::from_parts(parts, closed, None)?;
    let taken = ft.take(&feat_face);
    let attrs = taken.names().iter().zip(taken.columns()).try_fold(
        Table::new(n).with_column("offset", Column::f64(feat_offset))?,
        |acc, (name, c)| acc.with_column(name, c.clone()),
    )?;
    Ok(pl
        .with_attributes(attrs)?
        .with_crs(mesh.crs.clone())
        .with_length_unit(mesh.length_unit.clone()))
}
