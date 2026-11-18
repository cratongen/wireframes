use drusa_core::predicates::Sign;
use drusa_core::{Error, Mesh, Result, units};
use drusa_ops::{RepairOptions, RepairReport, repair};

use crate::BooleanOptions;
use crate::arrangement::{self, Arrangement};
use crate::classify::{self, Patches, Sample, Target};
use crate::hp;
use crate::ops::{Piece, assemble, build_mesh, is_closed_manifold, prepare};

pub struct ClipResult {
    pub below: Mesh,
    pub above: Mesh,
    pub repairs: [Option<RepairReport>; 2],
}

fn prepare_surface(surface: &Mesh, fix: bool) -> Result<(Mesh, Option<RepairReport>)> {
    let (m, r) = if fix {
        let (m, r) = repair(surface, &RepairOptions::default())?;
        (m, Some(r))
    } else {
        (surface.clone(), None)
    };
    let report = m.validate(0.0)?;
    let s = &report.summary;
    let bad = s.degenerate_faces
        + s.duplicate_faces
        + s.duplicate_vertices
        + s.non_manifold_edges
        + s.non_manifold_vertices
        + s.inconsistent_edges
        + s.self_intersections
        + s.coplanar_overlaps;
    if bad > 0 || m.is_empty() {
        return Err(Error::InvalidMesh {
            context: "surface".to_string(),
            report: Box::new(report),
        });
    }
    Ok((m, r))
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// Side of the surface holding the solid sub-triangle `s` next to its edge
/// `e`, which lies on the surface: `true` for the side the normals point to.
fn wedge_above(arr: &Arrangement, p: &Patches, s: u32, e: (u32, u32)) -> Result<bool> {
    let sub = &arr.subs[s as usize];
    let w = *sub
        .v
        .iter()
        .find(|&&v| v != e.0 && v != e.1)
        .expect("edge of the sub-triangle");
    let w = &arr.pts[w as usize];
    let mut faces: Vec<u32> = p.edges[1][&e]
        .iter()
        .map(|&t| arr.subs[t as usize].face)
        .collect();
    faces.sort_unstable();
    faces.dedup();
    let side = |f: u32| hp::side(w, &arr.planes[f as usize]);
    let flat = || {
        Error::invalid("surface does not cut through the solid: it touches the solid's boundary")
    };
    match faces[..] {
        [g] => match side(g) {
            Sign::Zero => Err(flat()),
            s => Ok(s == Sign::Positive),
        },
        [g1, g2] => {
            let c1 = arr.faces[g1 as usize].v;
            let far = arr.faces[g2 as usize]
                .v
                .iter()
                .position(|v| !c1.contains(v))
                .expect("faces sharing an edge");
            let far = hp::Pt::E(arr.faces[g2 as usize].c[far]);
            let bend = hp::side(&far, &arr.planes[g1 as usize]);
            let (s1, s2) = (side(g1), side(g2));
            match bend {
                Sign::Positive => Ok(s1 == Sign::Positive && s2 == Sign::Positive),
                Sign::Negative => Ok(!(s1 == Sign::Negative && s2 == Sign::Negative)),
                Sign::Zero => match s1 {
                    Sign::Zero => Err(flat()),
                    s => Ok(s == Sign::Positive),
                },
            }
        }
        _ => Err(Error::invalid(
            "internal: curve edge on more than two surface faces",
        )),
    }
}

/// Splits a solid by an open surface; `above` is on the side the surface normals point to.
/// Both pieces are closed and outward; the cap faces of `above` are the
/// flipped surface pieces inside the solid, those of `below` keep the
/// surface's orientation.
pub fn clip_by_surface(solid: &Mesh, surface: &Mesh, opts: &BooleanOptions) -> Result<ClipResult> {
    let labels = units::check_same(
        "surface",
        (&solid.crs, &solid.length_unit),
        (&surface.crs, &surface.length_unit),
    )?;
    let (a, ra) = prepare(solid, "solid", opts.repair)?;
    let (b, rb) = prepare_surface(surface, opts.repair)?;
    let arr = arrangement::build(&a, &b)?;
    let p = classify::patches(&arr);
    let solid_target = Target::new(
        arr.faces
            .iter()
            .filter(|f| f.mesh == 0)
            .map(|f| f.c)
            .collect(),
    );
    let surface_target = Target::new(
        arr.faces
            .iter()
            .filter(|f| f.mesh == 1)
            .map(|f| f.c)
            .collect(),
    );

    let mut above_patch = vec![false; p.members.len()];
    let mut cap_patch = vec![false; p.members.len()];
    for (k, members) in p.members.iter().enumerate() {
        let m = classify::mesh_of(&arr, &arr.subs[members[0] as usize]);
        if m == 1 {
            cap_patch[k] = solid_target.inside(&classify::sample(&arr, members)?)?;
            continue;
        }
        let curve = members.iter().find_map(|&s| {
            let v = arr.subs[s as usize].v;
            (0..3)
                .map(|e| ekey(v[e], v[(e + 1) % 3]))
                .find(|e| p.edges[1].contains_key(e))
                .map(|e| (s, e))
        });
        above_patch[k] = match curve {
            Some((s, e)) => wedge_above(&arr, &p, s, e)?,
            None => {
                let sm: Sample = classify::sample(&arr, members)?;
                let (up, down) = surface_target.crossings(&sm)?;
                if let Some(h) = surface_target.nearest(&sm, &up, true) {
                    up[h].facing != Sign::Positive
                } else if let Some(h) = surface_target.nearest(&sm, &down, false) {
                    down[h].facing == Sign::Positive
                } else {
                    false
                }
            }
        };
    }

    let (mut below, mut above): (Piece, Piece) = (vec![], vec![]);
    for (i, s) in arr.subs.iter().enumerate() {
        let i = i as u32;
        let m = classify::mesh_of(&arr, s);
        match (m, p.on[i as usize]) {
            (0, Some(true)) => below.push((i, false)),
            (0, Some(false)) => above.push((i, false)),
            (0, None) => {
                if above_patch[p.of[i as usize] as usize] {
                    above.push((i, false));
                } else {
                    below.push((i, false));
                }
            }
            (_, None) if cap_patch[p.of[i as usize] as usize] => {
                for e in 0..3 {
                    let k = ekey(s.v[e], s.v[(e + 1) % 3]);
                    if p.edges[1][&k].len() == 1 && !p.edges[0].contains_key(&k) {
                        return Err(Error::invalid("surface does not cut through the solid"));
                    }
                }
                below.push((i, false));
                above.push((i, true));
            }
            _ => {}
        }
    }
    let mut out = vec![];
    for piece in [below, above] {
        let asm = assemble(&arr, &[piece]);
        if !asm.triangles.is_empty() && !is_closed_manifold(&asm.triangles) {
            return Err(Error::invalid(
                "clip result is not manifold: the surface touches the solid along an edge or vertex",
            ));
        }
        out.push(build_mesh(
            asm,
            (a.vertex_attributes(), b.vertex_attributes()),
            (a.face_attributes(), b.face_attributes()),
            labels.clone(),
        )?);
    }
    let above = out.pop().expect("two pieces");
    let below = out.pop().expect("two pieces");
    Ok(ClipResult {
        below,
        above,
        repairs: [ra, rb],
    })
}
