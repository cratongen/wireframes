//! Results that rounding to `f64` used to leave open or self-intersecting, from
//! the underground acceptance scenarios, reduced to their inputs.

use num_rational::BigRational as Q;
use wireframes_boolean::{
    BooleanOptions, Op, boolean, clip_by_surface, difference, intersection, split_by_plane,
    split_by_plane_report,
};
use wireframes_core::{Error, Mesh, Result, fixtures};
use wireframes_ops::volume;

fn solid(m: &Mesh) -> bool {
    m.validate(0.0).unwrap().summary.is_solid()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

#[test]
fn utm_stope_clipped_by_a_planar_grid_fault_is_two_solids() {
    let stope = utm_stope();
    let r = clip_by_surface(&stope, &planar_grid_fault(), &BooleanOptions::default()).unwrap();
    assert!(solid(&r.below) && solid(&r.above));
    let (lo, hi) = (volume(&r.below).unwrap(), volume(&r.above).unwrap());
    assert!(close(lo + hi, volume(&stope).unwrap()), "{lo} + {hi}");
}

#[test]
fn thin_vein_minus_an_arched_drive_is_a_solid() {
    let (v, d) = (thin_vein(), arched_drive());
    let o = BooleanOptions::default();
    let (cut, common) = (
        difference(&v, &d, &o).unwrap().mesh,
        intersection(&v, &d, &o).unwrap().mesh,
    );
    assert!(solid(&cut) && solid(&common));
    let total = volume(&cut).unwrap() + volume(&common).unwrap();
    assert!(close(total, volume(&v).unwrap()));
}

#[test]
fn lofted_lens_clipped_by_a_tilted_grid_fault_is_two_solids() {
    let lens = lofted_lens();
    let r = clip_by_surface(&lens, &tilted_grid_fault(), &BooleanOptions::default()).unwrap();
    assert!(solid(&r.below) && solid(&r.above));
    let (lo, hi) = (volume(&r.below).unwrap(), volume(&r.above).unwrap());
    assert!(close(lo + hi, volume(&lens).unwrap()), "{lo} + {hi}");
}

/// The two inputs of a reduced case in `tests/data/rounding/<name>.txt`: per
/// mesh a line `tag vertices faces`, then one vertex and one face per line,
/// coordinates written so that they parse to the same `f64`.
fn pair(name: &str) -> (Mesh, Mesh) {
    let path = format!(
        "{}/tests/data/rounding/{name}.txt",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap();
    let mut lines = text.lines();
    let mut read = || {
        let head: Vec<usize> = lines
            .next()
            .unwrap()
            .split(' ')
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        let v: Vec<[f64; 3]> = (0..head[0])
            .map(|_| {
                let p: Vec<f64> = lines
                    .next()
                    .unwrap()
                    .split(' ')
                    .map(|x| x.parse().unwrap())
                    .collect();
                [p[0], p[1], p[2]]
            })
            .collect();
        let t: Vec<[u32; 3]> = (0..head[1])
            .map(|_| {
                let p: Vec<u32> = lines
                    .next()
                    .unwrap()
                    .split(' ')
                    .map(|x| x.parse().unwrap())
                    .collect();
                [p[0], p[1], p[2]]
            })
            .collect();
        Mesh::new(v, t).unwrap()
    };
    let a = read();
    (a, read())
}

const REDUCED: [&str; 9] = [
    "clipped_boxes",
    "pit_pushback",
    "touching_slivers",
    "runaway_union",
    "runaway_intersection",
    "runaway_difference",
    "runaway_difference_reversed",
    "rotated_stopes",
    "crosscut_chain",
];

/// Union, intersection and difference of `a` and `b` are solids (or empty)
/// and satisfy `vol(a | b) + vol(a & b) = vol(a) + vol(b)` and
/// `vol(a - b) + vol(a & b) = vol(a)`.
fn booleans_are_solids_with_exact_volumes(name: &str) {
    let (a, b) = pair(name);
    assert!(solid(&a) && solid(&b), "{name}: inputs");
    let o = BooleanOptions::default();
    let run = |op| {
        let m = boolean(&a, &b, op, &o)
            .unwrap_or_else(|e| panic!("{name} {op:?}: {e}"))
            .mesh;
        assert!(m.is_empty() || solid(&m), "{name} {op:?}");
        if m.is_empty() {
            0.0
        } else {
            volume(&m).unwrap()
        }
    };
    let (u, i, d) = (run(Op::Union), run(Op::Intersection), run(Op::Difference));
    let (va, vb) = (volume(&a).unwrap(), volume(&b).unwrap());
    assert!(close(u + i, va + vb), "{name}: {u} + {i} vs {va} + {vb}");
    assert!(close(d + i, va), "{name}: {d} + {i} vs {va}");
}

#[test]
fn difference_of_prisms_whose_caps_differ_by_one_ulp_is_empty_or_a_solid() {
    let (a, b) = ulp_prisms();
    let m = difference(&b, &a, &BooleanOptions::default()).unwrap().mesh;
    assert!(m.is_empty() || solid(&m));
    let v = if m.is_empty() {
        0.0
    } else {
        volume(&m).unwrap()
    };
    assert!(v.abs() <= 1e-12 * volume(&a).unwrap(), "{v}");
}

#[test]
fn booleans_of_boxes_clipped_by_the_same_topography() {
    booleans_are_solids_with_exact_volumes("clipped_boxes");
}

#[test]
fn booleans_of_a_triangular_pit_and_its_pushback_clipped_by_the_same_ground() {
    booleans_are_solids_with_exact_volumes("pit_pushback");
}

#[test]
fn booleans_of_touching_pushback_slivers() {
    booleans_are_solids_with_exact_volumes("touching_slivers");
}

#[test]
fn booleans_of_the_runaway_pushback_pair_reduced_for_union() {
    booleans_are_solids_with_exact_volumes("runaway_union");
}

#[test]
fn booleans_of_the_runaway_pushback_pair_reduced_for_intersection() {
    booleans_are_solids_with_exact_volumes("runaway_intersection");
}

#[test]
fn booleans_of_the_runaway_pushback_pair_reduced_for_difference() {
    booleans_are_solids_with_exact_volumes("runaway_difference");
}

#[test]
fn booleans_of_the_runaway_pushback_pair_reduced_for_the_reversed_difference() {
    booleans_are_solids_with_exact_volumes("runaway_difference_reversed");
}

#[test]
#[ignore = "rounding: a tessellated stope wall rotated to strike is planar only to rounding"]
fn union_of_a_tessellated_stope_and_a_box_rotated_to_strike() {
    booleans_are_solids_with_exact_volumes("rotated_stopes");
}

#[test]
fn union_of_a_rounded_drive_with_a_crosscut() {
    booleans_are_solids_with_exact_volumes("crosscut_chain");
}

#[test]
fn square_pit_at_utm_clipped_by_tilted_ground_is_two_solids() {
    let (pit, ground) = pair("square_pit_tilted_ground");
    assert!(solid(&pit));
    let r = clip_by_surface(&pit, &ground, &BooleanOptions::default()).unwrap();
    assert!(solid(&r.below) && solid(&r.above));
    let (lo, hi) = (volume(&r.below).unwrap(), volume(&r.above).unwrap());
    assert!(close(lo + hi, volume(&pit).unwrap()), "{lo} + {hi}");
}

#[test]
fn split_of_an_icosphere_through_its_centre_is_two_solids() {
    let m = fixtures::sphere(1);
    let (a, b) = split_by_plane(&m, m.vertices()[12], [1.0, 1.0, 1.0]).unwrap();
    assert!(solid(&a) && solid(&b));
    let (lo, hi) = (volume(&a).unwrap(), volume(&b).unwrap());
    assert!(close(lo + hi, volume(&m).unwrap()), "{lo} + {hi}");
}

/// Six times the signed volume of `m`, exactly.
fn exact_volume6(m: &Mesh) -> Q {
    let q = |x: f64| Q::from_float(x).unwrap();
    m.triangles().iter().fold(q(0.0), |acc, t| {
        let [a, b, c] = t.map(|i| m.vertices()[i as usize].map(q));
        acc + &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
            + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
            + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
    })
}

/// Every inclined split of `sphere(2)` through one of its vertices gives two
/// solids whose exact volumes add up to the sphere's within the sum of their
/// reports' volume bounds, and each part's vertices lie within the move bound.
#[test]
fn inclined_splits_of_an_icosphere_add_up_within_their_volume_bounds() {
    let m = fixtures::sphere(2);
    let whole = exact_volume6(&m);
    for &o in m.vertices().iter().step_by(3) {
        for n in [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 2.0, 3.0],
        ] {
            let (a, b, r) = split_by_plane_report(&m, o, n).unwrap();
            for part in [&a, &b] {
                assert!(part.is_empty() || solid(part), "{o:?} {n:?}");
            }
            let gap = exact_volume6(&a) + exact_volume6(&b) - &whole;
            let gap = if gap < Q::from_float(0.0).unwrap() {
                -gap
            } else {
                gap
            };
            let bound = Q::from_float(6.0 * (r[0].volume_bound + r[1].volume_bound)).unwrap();
            assert!(gap <= bound, "{o:?} {n:?}");
            for k in 0..3 {
                assert!(r[0].max_move[k] <= r[0].move_bound[k]);
                assert!(r[1].max_move[k] <= r[1].move_bound[k]);
            }
        }
    }
}

/// Each reduced case either returns solids or fails as a kernel limit, never
/// as an input defect or a panic; prints how many fail.
#[test]
fn reduced_rounding_cases_return_solids_or_kernel_errors() {
    let o = BooleanOptions::default();
    let mut refused = vec![];
    let mut check = |name: &str, r: Result<Vec<Mesh>>| match r {
        Ok(ms) => assert!(ms.iter().all(|m| m.is_empty() || solid(m)), "{name}"),
        Err(Error::Kernel { .. }) => refused.push(name.to_string()),
        Err(e) => panic!("{name}: {e}"),
    };
    for name in REDUCED {
        let (a, b) = pair(name);
        for op in [Op::Union, Op::Intersection, Op::Difference] {
            check(
                &format!("{name} {op:?}"),
                boolean(&a, &b, op, &o).map(|r| vec![r.mesh]),
            );
        }
    }
    let (a, b) = ulp_prisms();
    check("ulp prisms", difference(&b, &a, &o).map(|r| vec![r.mesh]));
    let (pit, ground) = pair("square_pit_tilted_ground");
    check(
        "square pit",
        clip_by_surface(&pit, &ground, &o).map(|r| vec![r.below, r.above]),
    );
    let m = fixtures::sphere(1);
    check(
        "icosphere",
        split_by_plane(&m, m.vertices()[12], [1.0, 1.0, 1.0]).map(|(a, b)| vec![a, b]),
    );
    println!("{} refused: {refused:?}", refused.len());
}

fn ulp_prisms() -> (Mesh, Mesh) {
    let prism = |half: f64, bottom: f64| {
        let mut v = vec![];
        for top in [false, true] {
            for y in [-half, half] {
                for x in [-half, half] {
                    v.push([x, y, if top { 5.0 + x / 4.0 + y / 2.0 } else { bottom }]);
                }
            }
        }
        Mesh::new(v, utm_stope().triangles().to_vec()).unwrap()
    };
    let a = prism(10.0, -10.0);
    let mut v = prism(6.0, -5.0).vertices().to_vec();
    v[4][2] = v[4][2].next_up();
    let b = Mesh::new(v, a.triangles().to_vec()).unwrap();
    (a, b)
}

fn utm_stope() -> Mesh {
    Mesh::new(
        vec![
            [700060.0, 7500000.0, -500.0],
            [700080.0, 7500000.0, -500.0],
            [700060.0, 7500016.0, -500.0],
            [700080.0, 7500016.0, -500.0],
            [700060.0, 7500000.0, -476.0],
            [700080.0, 7500000.0, -476.0],
            [700060.0, 7500016.0, -476.0],
            [700080.0, 7500016.0, -476.0],
        ],
        vec![
            [0, 2, 1],
            [1, 2, 3],
            [4, 5, 6],
            [5, 7, 6],
            [0, 1, 4],
            [1, 5, 4],
            [2, 6, 3],
            [3, 6, 7],
            [0, 4, 2],
            [2, 4, 6],
            [1, 3, 5],
            [3, 7, 5],
        ],
    )
    .unwrap()
}

fn planar_grid_fault() -> Mesh {
    Mesh::new(
        vec![
            [700035.0, 7500006.0, -510.0],
            [700051.6666666666, 7500011.0, -510.0],
            [700068.3333333334, 7500016.0, -510.0],
            [700085.0, 7500021.0, -510.0],
            [700036.6666666666, 7500009.333333333, -493.3333333333333],
            [700053.3333333333, 7500014.333333333, -493.3333333333333],
            [700070.0, 7500019.333333333, -493.3333333333333],
            [700086.6666666666, 7500024.333333333, -493.3333333333333],
            [700103.3333333333, 7500029.333333333, -493.3333333333333],
            [700038.3333333334, 7500012.666666667, -476.66666666666663],
            [700055.0, 7500017.666666667, -476.66666666666663],
            [700071.6666666667, 7500022.666666667, -476.66666666666663],
            [700088.3333333334, 7500027.666666667, -476.66666666666663],
            [700105.0, 7500032.666666667, -476.66666666666663],
            [700056.6666666666, 7500021.0, -460.0],
            [700073.3333333334, 7500026.0, -460.0],
            [700090.0, 7500031.0, -460.0],
        ],
        vec![
            [0, 1, 5],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 8, 7],
            [4, 5, 10],
            [5, 6, 11],
            [5, 11, 10],
            [6, 7, 12],
            [6, 12, 11],
            [7, 13, 12],
            [9, 10, 14],
            [10, 11, 15],
            [10, 15, 14],
            [11, 12, 16],
            [11, 16, 15],
        ],
    )
    .unwrap()
}

fn thin_vein() -> Mesh {
    Mesh::new(
        vec![
            [0.0, -20.0, 20.0],
            [150.0, -19.99999999999999, 20.0],
            [
                1.0471346844192482e-15,
                -37.10100716628344,
                -26.984631039295415,
            ],
            [150.0, -37.10100716628343, -26.984631039295415],
            [
                -5.753957801139253e-18,
                -19.90603073792141,
                19.965797985667432,
            ],
            [150.0, -19.906030737921398, 19.965797985667432],
            [
                1.041380726618109e-15,
                -37.00703790420485,
                -27.018833053627983,
            ],
            [150.0, -37.00703790420484, -27.018833053627983],
        ],
        vec![
            [0, 2, 1],
            [1, 2, 3],
            [4, 5, 6],
            [5, 7, 6],
            [0, 1, 4],
            [1, 5, 4],
            [2, 6, 3],
            [3, 6, 7],
            [0, 4, 2],
            [2, 4, 6],
            [1, 3, 5],
            [3, 7, 5],
        ],
    )
    .unwrap()
}

fn arched_drive() -> Mesh {
    Mesh::new(
        vec![
            [17.829392144688853, -38.75965265410792, -2.346587169561796],
            [22.170607855311147, -41.24034734589208, -2.346587169561796],
            [22.170607855311147, -41.24034734589208, 1.1534128304382045],
            [22.116186184543963, -41.20924924831083, 1.487194231372676],
            [21.955650102059785, -41.11751434403416, 1.8042384391145418],
            [21.697049557374665, -40.9697426042141, 2.088647533226305],
            [21.353351861620958, -40.77334392092626, 2.3261600541402494],
            [20.941791452420343, -40.5381665442402, 2.5048661322918333],
            [20.48300568721675, -40.27600324983814, 2.6158046987109405],
            [20.0, -40.0, 2.6534128304382048],
            [19.51699431278325, -39.72399675016186, 2.6158046987109405],
            [19.058208547579657, -39.4618334557598, 2.5048661322918333],
            [18.646648138379042, -39.22665607907374, 2.3261600541402494],
            [18.302950442625335, -39.0302573957859, 2.088647533226305],
            [18.04434989794022, -38.88248565596584, 1.8042384391145418],
            [17.883813815456037, -38.79075075168917, 1.487194231372676],
            [17.829392144688853, -38.75965265410792, 1.1534128304382045],
            [57.82939214468885, 31.240347345892083, -2.346587169561796],
            [62.17060785531115, 28.759652654107917, -2.346587169561796],
            [62.17060785531115, 28.759652654107917, 1.1534128304382045],
            [62.11618618454396, 28.790750751689163, 1.487194231372676],
            [61.955650102059785, 28.882485655965837, 1.8042384391145418],
            [61.69704955737467, 29.030257395785902, 2.088647533226305],
            [61.35335186162096, 29.226656079073738, 2.3261600541402494],
            [60.94179145242035, 29.461833455759802, 2.5048661322918333],
            [60.483005687216746, 29.72399675016186, 2.6158046987109405],
            [60.0, 30.0, 2.6534128304382048],
            [59.516994312783254, 30.27600324983814, 2.6158046987109405],
            [59.05820854757965, 30.538166544240198, 2.5048661322918333],
            [58.64664813837904, 30.773343920926262, 2.3261600541402494],
            [58.30295044262533, 30.969742604214098, 2.088647533226305],
            [58.044349897940215, 31.117514344034163, 1.8042384391145418],
            [57.88381381545604, 31.209249248310837, 1.487194231372676],
            [57.82939214468885, 31.240347345892083, 1.1534128304382045],
        ],
        vec![
            [0, 18, 1],
            [0, 17, 18],
            [1, 19, 2],
            [1, 18, 19],
            [2, 20, 3],
            [2, 19, 20],
            [3, 21, 4],
            [3, 20, 21],
            [4, 22, 5],
            [4, 21, 22],
            [5, 23, 6],
            [5, 22, 23],
            [6, 24, 7],
            [6, 23, 24],
            [7, 25, 8],
            [7, 24, 25],
            [8, 26, 9],
            [8, 25, 26],
            [9, 27, 10],
            [9, 26, 27],
            [10, 28, 11],
            [10, 27, 28],
            [11, 29, 12],
            [11, 28, 29],
            [12, 30, 13],
            [12, 29, 30],
            [13, 31, 14],
            [13, 30, 31],
            [14, 32, 15],
            [14, 31, 32],
            [15, 33, 16],
            [15, 32, 33],
            [16, 17, 0],
            [16, 33, 17],
            [0, 1, 2],
            [0, 2, 3],
            [0, 3, 4],
            [0, 4, 5],
            [0, 5, 6],
            [0, 6, 7],
            [0, 7, 8],
            [0, 8, 9],
            [0, 9, 10],
            [0, 10, 11],
            [0, 11, 12],
            [0, 12, 13],
            [0, 13, 14],
            [0, 14, 15],
            [0, 15, 16],
            [17, 19, 18],
            [17, 20, 19],
            [17, 21, 20],
            [17, 22, 21],
            [17, 23, 22],
            [17, 24, 23],
            [17, 25, 24],
            [17, 26, 25],
            [17, 27, 26],
            [17, 28, 27],
            [17, 29, 28],
            [17, 30, 29],
            [17, 31, 30],
            [17, 32, 31],
            [17, 33, 32],
        ],
    )
    .unwrap()
}

fn lofted_lens() -> Mesh {
    Mesh::new(
        vec![
            [699950.0, 7499973.9524677135, 200.0],
            [699950.0, 7499978.937741537, 194.93320808611836],
            [699950.0, 7499988.027898716, 191.35988012613404],
            [699950.0, 7500000.0, 189.30369430523672],
            [699950.0, 7500013.023766143, 190.60090638933463],
            [699950.0, 7500021.062258463, 194.93320808611836],
            [699950.0, 7500023.944202568, 200.0],
            [699950.0, 7500022.2318539, 205.34815284738164],
            [699950.0, 7500013.023766143, 209.39909361066537],
            [699950.0, 7500000.0, 210.1335838277633],
            [699950.0, 7499988.027898716, 208.64011987386596],
            [699950.0, 7499977.7681461, 205.34815284738164],
            [700000.0, 7499939.386575959, 202.951955497041],
            [700000.0, 7499941.424821112, 197.1473097907795],
            [700000.0, 7499946.059835797, 191.81974408276582],
            [700000.0, 7499952.339115938, 186.9387264256723],
            [700000.0, 7499960.148824594, 182.40009241619344],
            [700000.0, 7499970.0, 178.34936490538905],
            [700000.0, 7499982.238745344, 175.28053245620524],
            [700000.0, 7499996.33952583, 173.8133690255172],
            [700000.0, 7500010.870041859, 174.31372050975423],
            [700000.0, 7500024.171187737, 176.65203682239422],
            [700000.0, 7500035.2167947, 180.28984804281257],
            [700000.0, 7500043.972424743, 184.62614775686953],
            [700000.0, 7500050.977789362, 189.3324966352471],
            [700000.0, 7500056.506308676, 194.41989951182688],
            [700000.0, 7500060.0, 200.0],
            [700000.0, 7500060.259075794, 205.9506930489489],
            [700000.0, 7500056.258127477, 211.7724556452702],
            [700000.0, 7500047.952908432, 216.7655282411965],
            [700000.0, 7500036.442236303, 220.39600768056474],
            [700000.0, 7500023.358384187, 222.56284216640793],
            [700000.0, 7500009.96773946, 223.55410816036465],
            [700000.0, 7499996.6830947, 223.7287769390806],
            [700000.0, 7499983.34486673, 223.18000807197967],
            [700000.0, 7499970.0, 221.65063509461098],
            [700000.0, 7499957.502178862, 218.7687744948459],
            [700000.0, 7499947.402346693, 214.41417532921258],
            [700000.0, 7499941.177049709, 208.92075124904926],
            [700050.0, 7499976.060002866, 202.45861638938732],
            [700050.0, 7499976.929547956, 197.63068511713496],
            [700050.0, 7499981.259126261, 193.08209777319706],
            [700050.0, 7499990.694854902, 189.77681201806752],
            [700050.0, 7500003.044072718, 189.55408568545215],
            [700050.0, 7500013.527416846, 191.83423054976853],
            [700050.0, 7500021.651426485, 195.26518948825753],
            [700050.0, 7500026.0475322865, 200.0],
            [700050.0, 7500022.838544608, 204.99441370121312],
            [700050.0, 7500013.711017404, 208.27659917102972],
            [700050.0, 7500002.891862638, 209.9235964844717],
            [700050.0, 7499990.917413347, 209.97867198621466],
            [700050.0, 7499980.427181683, 207.22500163580577],
        ],
        vec![
            [10, 36, 35],
            [10, 37, 36],
            [10, 11, 37],
            [11, 38, 37],
            [11, 12, 38],
            [11, 0, 12],
            [0, 13, 12],
            [0, 14, 13],
            [0, 1, 14],
            [1, 15, 14],
            [1, 16, 15],
            [1, 17, 16],
            [1, 2, 17],
            [2, 18, 17],
            [2, 19, 18],
            [2, 3, 19],
            [3, 20, 19],
            [3, 21, 20],
            [3, 4, 21],
            [4, 22, 21],
            [4, 23, 22],
            [4, 5, 23],
            [5, 24, 23],
            [5, 25, 24],
            [5, 6, 25],
            [6, 26, 25],
            [6, 27, 26],
            [6, 28, 27],
            [6, 7, 28],
            [7, 29, 28],
            [7, 30, 29],
            [7, 8, 30],
            [8, 31, 30],
            [8, 32, 31],
            [8, 9, 32],
            [9, 33, 32],
            [9, 34, 33],
            [9, 10, 34],
            [10, 35, 34],
            [28, 29, 47],
            [29, 48, 47],
            [29, 30, 48],
            [30, 31, 48],
            [31, 32, 48],
            [32, 49, 48],
            [32, 33, 49],
            [33, 34, 49],
            [34, 50, 49],
            [34, 35, 50],
            [35, 36, 50],
            [36, 51, 50],
            [36, 37, 51],
            [37, 39, 51],
            [37, 38, 39],
            [38, 12, 39],
            [12, 40, 39],
            [12, 13, 40],
            [13, 14, 40],
            [14, 41, 40],
            [14, 15, 41],
            [15, 16, 41],
            [16, 42, 41],
            [16, 17, 42],
            [17, 18, 42],
            [18, 43, 42],
            [18, 19, 43],
            [19, 20, 43],
            [20, 44, 43],
            [20, 21, 44],
            [21, 22, 44],
            [22, 23, 44],
            [23, 45, 44],
            [23, 24, 45],
            [24, 46, 45],
            [24, 25, 46],
            [25, 26, 46],
            [26, 27, 46],
            [27, 47, 46],
            [27, 28, 47],
            [3, 2, 9],
            [10, 9, 2],
            [8, 4, 9],
            [3, 9, 4],
            [5, 4, 8],
            [1, 10, 2],
            [7, 5, 8],
            [11, 10, 1],
            [6, 5, 7],
            [0, 11, 1],
            [43, 49, 42],
            [50, 42, 49],
            [44, 49, 43],
            [48, 49, 44],
            [41, 42, 50],
            [51, 41, 50],
            [45, 48, 44],
            [40, 41, 51],
            [47, 48, 45],
            [39, 40, 51],
            [46, 47, 45],
        ],
    )
    .unwrap()
}

fn tilted_grid_fault() -> Mesh {
    Mesh::new(
        vec![
            [699979.5990953575, 7499920.0, 140.0],
            [699983.8452761157, 7499952.0, 140.0],
            [699988.0914568739, 7499984.0, 140.0],
            [699992.337637632, 7500016.0, 140.0],
            [699996.5838183903, 7500048.0, 140.0],
            [700000.8299991485, 7500080.0, 140.0],
            [699988.5341673776, 7499920.0, 170.0],
            [699992.7803481359, 7499952.0, 170.0],
            [699997.0265288941, 7499984.0, 170.0],
            [700001.2727096522, 7500016.0, 170.0],
            [700005.5188904104, 7500048.0, 170.0],
            [700009.7650711687, 7500080.0, 170.0],
            [699997.4692393978, 7499920.0, 200.0],
            [700001.715420156, 7499952.0, 200.0],
            [700005.9616009142, 7499984.0, 200.0],
            [700010.2077816724, 7500016.0, 200.0],
            [700014.4539624306, 7500048.0, 200.0],
            [700018.7001431888, 7500080.0, 200.0],
            [700006.404311418, 7499920.0, 230.0],
            [700010.6504921762, 7499952.0, 230.0],
            [700014.8966729344, 7499984.0, 230.0],
            [700019.1428536925, 7500016.0, 230.0],
            [700023.3890344508, 7500048.0, 230.0],
            [700027.635215209, 7500080.0, 230.0],
            [700015.3393834381, 7499920.0, 260.0],
            [700019.5855641963, 7499952.0, 260.0],
            [700023.8317449546, 7499984.0, 260.0],
            [700028.0779257127, 7500016.0, 260.0],
            [700032.3241064709, 7500048.0, 260.0],
            [700036.5702872291, 7500080.0, 260.0],
        ],
        vec![
            [0, 1, 7],
            [0, 7, 6],
            [1, 2, 8],
            [1, 8, 7],
            [2, 3, 9],
            [2, 9, 8],
            [3, 4, 10],
            [3, 10, 9],
            [4, 5, 11],
            [4, 11, 10],
            [6, 7, 13],
            [6, 13, 12],
            [7, 8, 14],
            [7, 14, 13],
            [8, 9, 15],
            [8, 15, 14],
            [9, 10, 16],
            [9, 16, 15],
            [10, 11, 17],
            [10, 17, 16],
            [12, 13, 19],
            [12, 19, 18],
            [13, 14, 20],
            [13, 20, 19],
            [14, 15, 21],
            [14, 21, 20],
            [15, 16, 22],
            [15, 22, 21],
            [16, 17, 23],
            [16, 23, 22],
            [18, 19, 25],
            [18, 25, 24],
            [19, 20, 26],
            [19, 26, 25],
            [20, 21, 27],
            [20, 27, 26],
            [21, 22, 28],
            [21, 28, 27],
            [22, 23, 29],
            [22, 29, 28],
        ],
    )
    .unwrap()
}
