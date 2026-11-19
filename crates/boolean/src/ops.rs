use std::collections::HashMap;

use wireframes_core::predicates::Sign;
use wireframes_core::topology::{Radial, radial_pairs};
use wireframes_core::{Error, Mesh, Result, Table, units};
use wireframes_ops::{RepairOptions, RepairReport, repair};

use crate::arrangement::{self, Arrangement};
use crate::classify::{self, Status};
use crate::conform::conform;
use crate::hp;
use crate::round::{self, Assembled, Draft, RoundingReport};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Union,
    Intersection,
    Difference,
    SymmetricDifference,
}

/// `repair`: repair each input with `RepairOptions::default()` first.
#[derive(Clone, Debug, Default)]
pub struct BooleanOptions {
    pub repair: bool,
}

#[derive(Clone, Debug)]
pub struct BooleanResult {
    pub mesh: Mesh,
    /// Per output face: (0 for a, 1 for b, input face index after repair).
    pub face_source: Vec<(u8, u32)>,
    pub repairs: [Option<RepairReport>; 2],
    /// What rounding the exact result to `f64` did; see [`RoundingReport`].
    pub rounding: RoundingReport,
}

/// Repairs `mesh` when asked, then refuses it unless it is a valid solid or
/// has no faces (the empty solid a boolean returns).
pub(crate) fn prepare(
    mesh: &Mesh,
    context: &str,
    fix: bool,
) -> Result<(Mesh, Option<RepairReport>)> {
    let (m, r) = if fix {
        let (m, r) = repair(mesh, &RepairOptions::default())?;
        (m, Some(r))
    } else {
        (mesh.clone(), None)
    };
    let report = m.validate(0.0)?;
    if m.n_faces() > 0 && !report.summary.is_solid() {
        return Err(Error::InvalidMesh {
            context: context.to_string(),
            report: Box::new(report),
        });
    }
    Ok((m, r))
}

/// Sub-triangle ids with a flip flag.
pub(crate) type Piece = Vec<(u32, bool)>;

/// One output shell: sub-triangle ids with a flip flag and, per corner, the
/// umbrella of faces about that corner's point it belongs to. Shells, and
/// umbrellas within a shell, get separate vertex copies, so shells may touch
/// each other and a shell may touch itself without making the output
/// non-manifold.
pub(crate) type Shell = Vec<(u32, bool, [u32; 3])>;

/// The exact result of `pieces` rounded to `f64` by [`round::materialize`],
/// with its report and whether the rounding certified it.
pub(crate) fn assemble(
    arr: &Arrangement,
    pieces: &[Shell],
) -> Result<(Assembled, RoundingReport, bool)> {
    let r = round::materialize(draft(arr, pieces))?;
    Ok((r.mesh, r.report, r.certified))
}

/// The exact result of `pieces`: a vertex copy per exact point and umbrella
/// of each shell.
pub(crate) fn draft(arr: &Arrangement, pieces: &[Shell]) -> Draft {
    let mut d = Draft {
        point: vec![],
        home: vec![],
        err: vec![],
        rows: vec![],
        triangles: vec![],
        face_source: vec![],
        normal: vec![],
        shells: vec![],
        closed: true,
    };
    for piece in pieces {
        let mut map: HashMap<(u32, u32), u32> = HashMap::new();
        let mut faces = vec![];
        for &(s, flip, umbrella) in piece {
            let sub = &arr.subs[s as usize];
            let v = [0, 1, 2].map(|c| {
                let p = sub.v[c];
                *map.entry((p, umbrella[c])).or_insert_with(|| {
                    let pt = &arr.pts[p as usize];
                    let home = pt.to_f64();
                    d.err.push(if pt.explicit().is_some() {
                        [0.0; 3]
                    } else {
                        home.map(round::half_ulp)
                    });
                    d.home.push(home);
                    d.rows.push(arr.source_row(p));
                    d.point.push(p);
                    d.home.len() as u32 - 1
                })
            });
            faces.push(d.triangles.len() as u32);
            d.triangles.push(if flip { [v[0], v[2], v[1]] } else { v });
            let f = &arr.faces[sub.face as usize];
            d.face_source.push((f.mesh, f.index));
            let n = arr.planes[sub.face as usize]
                .normal_exp()
                .clone()
                .map(|x| x.estimate());
            d.normal.push(if flip { n.map(|x| -x) } else { n });
        }
        d.shells.push(faces);
    }
    d
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// The shells of a closed piece, each ordered as in `piece` and the shells
/// ordered by first sub-triangle. Sub-triangles join across edges with two of
/// them, and across edges with more in the pairs that bound wedges of solid,
/// found by exact radial order about the edge. Corners of joined
/// sub-triangles at the ends of the edge share an umbrella.
pub(crate) fn shells(arr: &Arrangement, piece: &Piece) -> Result<Vec<Shell>> {
    let tri = |k: usize| {
        let (s, flip) = piece[k];
        let v = arr.subs[s as usize].v;
        if flip { [v[0], v[2], v[1]] } else { v }
    };
    let mut at: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    for k in 0..piece.len() {
        let t = tri(k);
        for e in 0..3 {
            at.entry(ekey(t[e], t[(e + 1) % 3]))
                .or_default()
                .push(k as u32);
        }
    }
    let mut parent: Vec<u32> = (0..piece.len() as u32).collect();
    fn root(parent: &mut [u32], mut i: u32) -> u32 {
        while parent[i as usize] != i {
            parent[i as usize] = parent[parent[i as usize] as usize];
            i = parent[i as usize];
        }
        i
    }
    let join = |parent: &mut [u32], a: u32, b: u32| {
        let (a, b) = (root(parent, a), root(parent, b));
        parent[a.max(b) as usize] = a.min(b);
    };
    let mut corners: Vec<u32> = (0..3 * piece.len() as u32).collect();
    let corner = |k: u32, p: u32| {
        let v = arr.subs[piece[k as usize].0 as usize].v;
        3 * k + v.iter().position(|&x| x == p).expect("corner on the edge") as u32
    };
    let link = |parent: &mut [u32], corners: &mut [u32], f: u32, g: u32, p: u32, q: u32| {
        join(parent, f, g);
        join(corners, corner(f, p), corner(g, p));
        join(corners, corner(f, q), corner(g, q));
    };
    let mut keys: Vec<(u32, u32)> = at.keys().copied().collect();
    keys.sort_unstable();
    for (p, q) in keys {
        let ks = &at[&(p, q)];
        if ks.len() == 2 {
            link(&mut parent, &mut corners, ks[0], ks[1], p, q);
            continue;
        }
        if ks.len() < 2 {
            continue;
        }
        let d = hp::direction(&arr.pts[p as usize], &arr.pts[q as usize]);
        let forward: Vec<bool> = ks
            .iter()
            .map(|&k| {
                let t = tri(k as usize);
                (0..3).any(|e| t[e] == p && t[(e + 1) % 3] == q)
            })
            .collect();
        let m: Vec<[wireframes_core::exact::Expansion; 3]> = ks
            .iter()
            .zip(&forward)
            .map(|(&k, &fwd)| {
                let (s, flip) = piece[k as usize];
                let n = arr.planes[arr.subs[s as usize].face as usize].normal_exp();
                if flip == fwd {
                    n.clone().map(|x| x.neg())
                } else {
                    n.clone()
                }
            })
            .collect();
        let turn = |a: usize, b: usize| hp::det3(&d, &m[a], &m[b]);
        let same = |a: usize, b: usize| hp::dot3(&m[a], &m[b]) == Sign::Positive;
        match radial_pairs(&forward, turn, same) {
            Radial::Paired(pairs) => {
                for [a, b] in pairs {
                    link(&mut parent, &mut corners, ks[a], ks[b], p, q);
                }
            }
            _ => {
                return Err(Error::kernel(
                    "the faces around an edge of the result do not bound wedges of solid",
                ));
            }
        }
    }
    let mut index: HashMap<u32, usize> = HashMap::new();
    let mut out: Vec<Shell> = vec![];
    for (k, &(s, flip)) in piece.iter().enumerate() {
        let r = root(&mut parent, k as u32);
        let next = out.len();
        let i = *index.entry(r).or_insert(next);
        if i == out.len() {
            out.push(vec![]);
        }
        let umbrella = [0, 1, 2].map(|c| root(&mut corners, 3 * k as u32 + c));
        out[i].push((s, flip, umbrella));
    }
    Ok(out)
}

/// Every edge used once in each direction and every vertex with a single fan.
pub(crate) fn is_closed_manifold(triangles: &[[u32; 3]]) -> bool {
    let mut dir: HashMap<(u32, u32), u32> = HashMap::new();
    for t in triangles {
        for k in 0..3 {
            *dir.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
        }
    }
    if dir
        .iter()
        .any(|(&(a, b), &n)| n != 1 || dir.get(&(b, a)) != Some(&1))
    {
        return false;
    }
    let mut link: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
    for t in triangles {
        for k in 0..3 {
            link.entry(t[k])
                .or_default()
                .push((t[(k + 1) % 3], t[(k + 2) % 3]));
        }
    }
    link.values().all(|arcs| {
        let next: HashMap<u32, u32> = arcs.iter().copied().collect();
        if next.len() != arcs.len() {
            return false;
        }
        let start = arcs[0].0;
        let mut c = arcs[0].1;
        let mut n = 1;
        while c != start {
            match next.get(&c) {
                Some(&d) => c = d,
                None => return false,
            }
            n += 1;
            if n > arcs.len() {
                return false;
            }
        }
        n == arcs.len()
    })
}

pub(crate) fn build_mesh(
    asm: Assembled,
    vertex_tables: (&Table, &Table),
    face_tables: (&Table, &Table),
    labels: (Option<String>, Option<String>),
) -> Result<Mesh> {
    let vt = Table::concat(vertex_tables.0, vertex_tables.1)?;
    let ft = Table::concat(face_tables.0, face_tables.1)?;
    let n0 = face_tables.0.rows() as u32;
    let rows: Vec<u32> = asm
        .face_source
        .iter()
        .map(|&(m, f)| if m == 0 { f } else { n0 + f })
        .collect();
    let mesh = if asm.triangles.is_empty() {
        Mesh::empty().with_attributes(vt.take(&[]), ft.take(&[]))?
    } else {
        Mesh::new(asm.vertices, asm.triangles)?
            .with_attributes(vt.take_or_null(&asm.vertex_rows), ft.take(&rows))?
    };
    Ok(mesh.with_crs(labels.0).with_length_unit(labels.1))
}

/// Whether a split of a closed mesh that is not a solid is still closed and
/// manifold: its edges pair into wedges of solid and its vertices are manifold.
pub(crate) fn is_closed_result(mesh: &Mesh) -> Result<bool> {
    if mesh.is_empty() || is_closed_manifold(mesh.triangles()) {
        return Ok(true);
    }
    let s = mesh.validate(0.0)?.summary;
    Ok(s.is_closed
        && s.boundary_edges == 0
        && s.non_manifold_edges == 0
        && s.non_manifold_vertices == 0
        && s.inconsistent_edges == 0)
}

/// `mesh` when the rounding certified it as a solid. Exact results of solids
/// are solids; a rounded result that the rounding operations could not
/// certify is an error that carries the report and the mesh, never a return.
pub(crate) fn solid_result(mesh: Mesh, certified: bool, what: &str) -> Result<Mesh> {
    round::certify(mesh, certified, what)
}

fn keep_rule(op: Op, mesh: usize, st: Status) -> Option<bool> {
    use Status::*;
    match (op, mesh, st) {
        (Op::Union, _, Out) | (Op::Union, 0, OnSame) => Some(false),
        (Op::Intersection, _, In) | (Op::Intersection, 0, OnSame) => Some(false),
        (Op::Difference, 0, Out) | (Op::Difference, 0, OnOpp) => Some(false),
        (Op::Difference, 1, In) => Some(true),
        (Op::SymmetricDifference, _, Out) => Some(false),
        (Op::SymmetricDifference, _, In) => Some(true),
        _ => None,
    }
}

/// The sub-triangles `op` keeps, each with a flip flag.
pub(crate) fn kept(arr: &Arrangement, op: Op) -> Result<Piece> {
    let p = classify::patches(arr);
    let st = classify::solid_status(arr, &p)?;
    Ok(arr
        .subs
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            keep_rule(op, classify::mesh_of(arr, s), st[i]).map(|flip| (i as u32, flip))
        })
        .collect())
}

pub fn boolean(a: &Mesh, b: &Mesh, op: Op, opts: &BooleanOptions) -> Result<BooleanResult> {
    let labels = units::check_same("b", (&a.crs, &a.length_unit), (&b.crs, &b.length_unit))?;
    let (a, ra) = prepare(a, "a", opts.repair)?;
    let (b, rb) = prepare(b, "b", opts.repair)?;
    let run = |a: &Mesh, b: &Mesh| -> Result<(Mesh, Vec<(u8, u32)>, RoundingReport, bool)> {
        let arr = arrangement::build(a, b)?;
        let kept = kept(&arr, op)?;
        let (asm, rounding, certified) = assemble(&arr, &shells(&arr, &kept)?)?;
        let face_source = asm.face_source.clone();
        let mesh = build_mesh(
            asm,
            (a.vertex_attributes(), b.vertex_attributes()),
            (a.face_attributes(), b.face_attributes()),
            labels.clone(),
        )?;
        Ok((mesh, face_source, rounding, certified))
    };
    let (mut mesh, mut face_source, mut rounding, mut certified) = run(&a, &b)?;
    if !certified && let Some(c) = conform(&a, &b)? {
        let (m, fs, mut r, ok) = run(&c.meshes[0], &c.meshes[1])?;
        if ok {
            r.inserted_points += c.inserted;
            r.flipped_edges += c.flips;
            r.steps += c.steps;
            r.step_bound = r.step_bound.saturating_add(c.step_bound);
            r.volume_bound = round::up(r.volume_bound + c.volume);
            mesh = m;
            face_source = fs
                .into_iter()
                .map(|(k, f)| (k, c.src[k as usize][f as usize]))
                .collect();
            rounding = r;
            certified = true;
        }
    }
    let mesh = solid_result(mesh, certified, "the exact result")?;
    Ok(BooleanResult {
        mesh,
        face_source,
        repairs: [ra, rb],
        rounding,
    })
}

pub fn union(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Union, opts)
}

pub fn intersection(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Intersection, opts)
}

pub fn difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Difference, opts)
}

/// The parts of `a` outside `b` and of `b` outside `a`. Where they meet along
/// the intersection curve they are separate shells in contact.
pub fn symmetric_difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::SymmetricDifference, opts)
}
