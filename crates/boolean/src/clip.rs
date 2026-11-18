use drusa_core::predicates::Sign;
use drusa_core::{Error, Mesh, Result, units};
use drusa_ops::{RepairOptions, RepairReport, repair};

use crate::BooleanOptions;
use crate::arrangement;
use crate::classify::{self, Sample, Target};
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
            cap_patch[k] = classify::patch_inside(&arr, &p, members, &solid_target)?;
            continue;
        }
        above_patch[k] = match classify::curve_side(&arr, &p, members) {
            Some(above) => above,
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
