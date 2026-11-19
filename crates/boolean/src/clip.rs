use wireframes_core::predicates::{Sign, kept_axes};
use wireframes_core::{Error, Mesh, Result, units};
use wireframes_ops::{RepairOptions, RepairReport, repair};

use crate::BooleanOptions;
use crate::arrangement::{self, Arrangement};
use crate::classify::{self, Sample, Target};
use crate::conform::conform;
use crate::hp::Pt;
use crate::ops::{Piece, assemble, build_mesh, prepare, shells, solid_result};
use crate::round::{self, RoundingReport};

pub struct ClipResult {
    pub below: Mesh,
    pub above: Mesh,
    pub repairs: [Option<RepairReport>; 2],
    /// What rounding `below` and `above` to `f64` did.
    pub rounding: [RoundingReport; 2],
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

/// Side of the surface holding a sample, from the nearest crossing on the
/// sample's axis line, then on the other axis lines through an explicit sample.
/// A sample whose lines all miss the surface lies beyond its extent and has no
/// side.
fn ray_side(surface: &Target, sm: Sample) -> Result<Option<bool>> {
    let mut lines = vec![sm];
    if let Some(q) = lines[0].p.explicit() {
        lines.extend([0, 1].map(|axis| {
            let (i, j) = kept_axes(axis);
            Sample {
                p: Pt::E(q),
                axis,
                uv: [q[i], q[j]],
            }
        }));
    }
    for sm in &lines {
        let (up, down) = surface.crossings(sm)?;
        if let Some(h) = surface.nearest(sm, &up, true) {
            return Ok(Some(up[h].facing != Sign::Positive));
        }
        if let Some(h) = surface.nearest(sm, &down, false) {
            return Ok(Some(down[h].facing == Sign::Positive));
        }
    }
    Ok(None)
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// The sub-triangles of the part below and the part above the surface.
pub(crate) fn pieces(arr: &Arrangement) -> Result<[Piece; 2]> {
    let p = classify::patches(arr);
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

    let mut above_patch = vec![None; p.members.len()];
    let mut cap_patch = vec![false; p.members.len()];
    for (k, members) in p.members.iter().enumerate() {
        let m = classify::mesh_of(arr, &arr.subs[members[0] as usize]);
        if m == 1 {
            cap_patch[k] = classify::patch_inside(arr, &p, members, &solid_target)?;
            continue;
        }
        above_patch[k] = match classify::curve_side(arr, &p, members) {
            Some(above) => Some(above),
            None => ray_side(&surface_target, classify::sample(arr, members)?)?,
        };
    }

    let mut beyond = false;
    let (mut below, mut above): (Piece, Piece) = (vec![], vec![]);
    for (i, s) in arr.subs.iter().enumerate() {
        let i = i as u32;
        let m = classify::mesh_of(arr, s);
        match (m, p.on[i as usize]) {
            (0, Some(true)) => below.push((i, false)),
            (0, Some(false)) => above.push((i, false)),
            (0, None) => match above_patch[p.of[i as usize] as usize] {
                Some(true) => above.push((i, false)),
                Some(false) => below.push((i, false)),
                None => beyond = true,
            },
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
    if beyond {
        return Err(Error::invalid(
            "surface does not reach part of the solid: no line along x, y or z through it meets \
             the surface, so it is neither above nor below; extend the surface past the solid",
        ));
    }
    Ok([below, above])
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
    type Piece2 = (Mesh, RoundingReport, bool);
    let run = |a: &Mesh, b: &Mesh| -> Result<[Piece2; 2]> {
        let arr = arrangement::build(a, b)?;
        let mut out = vec![];
        for piece in pieces(&arr)? {
            let (asm, r, certified) = assemble(&arr, &shells(&arr, &piece)?)?;
            let mesh = build_mesh(
                asm,
                (a.vertex_attributes(), b.vertex_attributes()),
                (a.face_attributes(), b.face_attributes()),
                labels.clone(),
            )?;
            out.push((mesh, r, certified));
        }
        Ok(out.try_into().map_err(|_| ()).expect("two pieces"))
    };
    let mut parts = run(&a, &b)?;
    if !(parts[0].2 && parts[1].2)
        && let Some(c) = conform(&a, &b)?
    {
        let mut again = run(&c.meshes[0], &c.meshes[1])?;
        if again[0].2 && again[1].2 {
            for (k, part) in again.iter_mut().enumerate() {
                let r = &mut part.1;
                if k == 0 {
                    r.inserted_points += c.inserted;
                    r.flipped_edges += c.flips;
                    r.steps += c.steps;
                }
                r.step_bound = r.step_bound.saturating_add(c.step_bound);
                r.volume_bound = round::up(r.volume_bound + c.volume);
            }
            parts = again;
        }
    }
    let [(lo, rlo, clo), (hi, rhi, chi)] = parts;
    let below = solid_result(lo, clo, "the exact clip")?;
    let above = solid_result(hi, chi, "the exact clip")?;
    let rounding = [rlo, rhi];
    Ok(ClipResult {
        below,
        above,
        repairs: [ra, rb],
        rounding,
    })
}
