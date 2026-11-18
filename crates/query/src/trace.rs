use drusa_core::{Error, Result};

use crate::distance::check_finite;
use crate::ray::{HitKind, Line, Param};
use crate::{Location, SolidIndex};

/// Desurveyed trace vertices and their measured depths, non-decreasing.
#[derive(Clone, Debug, PartialEq)]
pub struct Trace {
    pub xyz: Vec<[f64; 3]>,
    pub depth: Vec<f64>,
}

/// `entry_face` and `exit_face` are the faces the trace meets at `from` and
/// at `to` (the smallest face index there), `None` where it meets none.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    pub from: f64,
    pub to: f64,
    pub inside: bool,
    pub entry_face: Option<u32>,
    pub exit_face: Option<u32>,
}

fn check(trace: &Trace) -> Result<()> {
    let n = trace.xyz.len();
    if n != trace.depth.len() {
        return Err(Error::invalid(format!(
            "trace xyz and depth must have the same length: {n} vs {}",
            trace.depth.len()
        )));
    }
    if n < 2 {
        return Err(Error::invalid(format!(
            "trace needs at least 2 vertices: it has {n}"
        )));
    }
    check_finite("xyz", &trace.xyz)?;
    if let Some(i) = trace.depth.iter().position(|d| !d.is_finite()) {
        return Err(Error::invalid(format!(
            "depth must be finite: depth {i} is {}",
            trace.depth[i]
        )));
    }
    if let Some(i) = trace.depth.windows(2).position(|w| w[1] < w[0]) {
        return Err(Error::invalid(format!(
            "depth must be non-decreasing: depth {} is below depth {i}",
            i + 1
        )));
    }
    Ok(())
}

/// Pieces of one trace segment between consecutive distinct exact crossing
/// parameters. A piece inside a coplanar range lies in a face and is inside;
/// any other piece is open and holds no boundary point, so its float midpoint
/// decides its label.
fn segment_pieces(index: &SolidIndex, line: Line, d0: f64, d1: f64, out: &mut Vec<Interval>) {
    let hits = index.surface.raw_hits(&line);
    let mut breaks: Vec<(Param, Option<u32>)> = vec![(Param::zero(), None), (Param::one(), None)];
    for h in &hits {
        breaks.push((h.lo.clone(), Some(h.face)));
        breaks.push((h.hi.clone(), Some(h.face)));
    }
    breaks.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.unwrap_or(u32::MAX).cmp(&b.1.unwrap_or(u32::MAX)))
    });
    breaks.dedup_by(|b, a| a.0.cmp(&b.0).is_eq());
    let last = breaks.len() - 1;
    let mut depth = d0;
    let depths: Vec<f64> = breaks
        .iter()
        .enumerate()
        .map(|(i, (t, _))| {
            depth = match i {
                0 => d0,
                _ if i == last => d1,
                _ => (d0 + t.value().clamp(0.0, 1.0) * (d1 - d0)).clamp(depth, d1),
            };
            depth
        })
        .collect();
    for i in 0..last {
        let (a, b) = (&breaks[i], &breaks[i + 1]);
        let in_face = hits.iter().any(|h| {
            h.kind == HitKind::Coplanar && h.lo.cmp(&a.0).is_le() && b.0.cmp(&h.hi).is_le()
        });
        let inside = in_face || {
            let mid = line.at(0.5 * (a.0.value() + b.0.value()));
            index.locate(mid) != Location::Outside
        };
        out.push(Interval {
            from: depths[i],
            to: depths[i + 1],
            inside,
            entry_face: a.1,
            exit_face: b.1,
        });
    }
}

pub(crate) fn intervals(index: &SolidIndex, trace: &Trace) -> Result<Vec<Interval>> {
    check(trace)?;
    let mut pieces = vec![];
    for k in 0..trace.xyz.len() - 1 {
        let (p, q) = (trace.xyz[k], trace.xyz[k + 1]);
        let (d0, d1) = (trace.depth[k], trace.depth[k + 1]);
        if p == q {
            pieces.push(Interval {
                from: d0,
                to: d1,
                inside: index.locate(p) != Location::Outside,
                entry_face: None,
                exit_face: None,
            });
        } else {
            segment_pieces(index, Line::Segment(p, q), d0, d1, &mut pieces);
        }
    }
    let mut out: Vec<Interval> = vec![];
    for piece in pieces.into_iter().filter(|x| x.to > x.from) {
        match out.last_mut() {
            Some(prev) if prev.inside == piece.inside => {
                prev.to = piece.to;
                prev.exit_face = piece.exit_face;
            }
            _ => out.push(piece),
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use drusa_core::fixtures;

    fn straight(a: [f64; 3], b: [f64; 3], d: [f64; 2]) -> Trace {
        Trace {
            xyz: vec![a, b],
            depth: d.to_vec(),
        }
    }

    #[test]
    fn vertical_hole_through_a_cube() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let iv = s
            .trace_intervals(&straight([0.25, 0.75, 2.0], [0.25, 0.75, -1.0], [0.0, 3.0]))
            .unwrap();
        assert_eq!(iv.len(), 3);
        assert_eq!(
            iv.iter()
                .map(|x| (x.from, x.to, x.inside))
                .collect::<Vec<_>>(),
            [(0.0, 1.0, false), (1.0, 2.0, true), (2.0, 3.0, false)]
        );
        assert!(iv[1].entry_face.is_some() && iv[1].exit_face.is_some());
        assert_eq!(iv[0].entry_face, None);
        assert_eq!(iv[2].exit_face, None);
    }

    #[test]
    fn runs_in_a_face_count_as_inside() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let iv = s
            .trace_intervals(&straight([-1.0, 0.5, 1.0], [2.0, 0.5, 1.0], [0.0, 3.0]))
            .unwrap();
        assert_eq!(
            iv.iter()
                .map(|x| (x.from, x.to, x.inside))
                .collect::<Vec<_>>(),
            [(0.0, 1.0, false), (1.0, 2.0, true), (2.0, 3.0, false)]
        );
        let iv = s
            .trace_intervals(&straight([-1.0, 2.0, 1.0], [2.0, -1.0, 1.0], [0.0, 3.0]))
            .unwrap();
        assert_eq!(iv.len(), 3);
        assert!(iv[1].inside);
        assert!((iv[1].to - iv[1].from - 1.0).abs() < 1e-15);
        let iv = s
            .trace_intervals(&straight([0.5, 0.0, 2.0], [0.5, 2.0, 0.0], [0.0, 3.0]))
            .unwrap();
        assert_eq!(
            iv.len(),
            1,
            "touching the top edge at one point leaves no inside run"
        );
        assert!(!iv[0].inside);
    }

    #[test]
    fn collar_inside_and_bent_trace() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let t = Trace {
            xyz: vec![
                [0.5, 0.5, 0.5],
                [0.5, 0.5, 0.5],
                [0.5, 0.5, -0.5],
                [2.0, 0.5, -0.5],
            ],
            depth: vec![0.0, 0.0, 1.0, 2.5],
        };
        let iv = s.trace_intervals(&t).unwrap();
        assert_eq!(
            iv.iter()
                .map(|x| (x.from, x.to, x.inside))
                .collect::<Vec<_>>(),
            [(0.0, 0.5, true), (0.5, 2.5, false)]
        );
        let bad = Trace {
            xyz: vec![[0.0; 3], [1.0; 3]],
            depth: vec![1.0, 0.0],
        };
        assert!(s.trace_intervals(&bad).is_err());
        assert!(
            s.trace_intervals(&straight([0.0; 3], [1.0; 3], [0.0, 1.0]))
                .is_ok()
        );
        let one = Trace {
            xyz: vec![[0.0; 3]],
            depth: vec![0.0],
        };
        assert!(s.trace_intervals(&one).is_err());
    }
}
