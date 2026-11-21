//! Attack 2: seeded mining and geology scenarios against the rounding of
//! exact results to `f64`: stope chains between a common footwall and hanging
//! wall minus drive networks, sublevel layouts united one stope at a time,
//! helical declines united with crosscuts, pits with ramps against orebodies
//! and rough topography, pushback sequences, solids sharing whole walls with
//! and without matching vertices, chains that feed each result into the next,
//! and slots a few ulps wide crossed by development.
//!
//! Every operation is checked the same way: no error, a solid within its step
//! and move bounds, the exact volume identity within the reported volume
//! bounds, and the same result on 1 and 8 threads within a time limit.
//! `ATTACK2_SEEDS` (default 1) sets how many seeds each scenario runs,
//! `ATTACK2_SIZE` (default 0) selects the larger layouts, `ATTACK2_LIMIT` the
//! seconds per operation (default 300), `ATTACK2_SPLITS` the splits of each
//! stope block (default 3, or 12 for the larger layouts), and `ATTACK2_DUMP` a folder that
//! receives the two inputs of every failing operation in the format of
//! `tests/data/attacks`.

use num_rational::BigRational as Q;
use std::f64::consts::TAU;
use std::sync::mpsc;
use std::time::{Duration, Instant};
use wireframes_boolean::{
    BooleanOptions, BooleanResult, Op, RoundingReport, boolean, clip_by_surface,
    split_by_plane_report,
};
use wireframes_core::Mesh;

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x2545_f491_4f6c_dd1d) ^ 0x9e37_79b9_7f4a_7c15)
    }

    fn u(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn f(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * (self.u() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn i(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.u() % (hi - lo + 1) as u64) as i64
    }

    fn utm(&mut self) -> [f64; 3] {
        [
            700_000.0 + self.f(0.0, 1500.0),
            7_500_000.0 + self.f(0.0, 1500.0),
            self.f(-600.0, 400.0),
        ]
    }
}

fn env(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

const BOX_FACES: [[u32; 3]; 12] = [
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
];

/// Corners in the order of `BOX_FACES`: bit 0 is x, bit 1 is y, bit 2 is z.
fn corners(lo: [f64; 3], hi: [f64; 3]) -> Vec<[f64; 3]> {
    (0..8)
        .map(|i| [0, 1, 2].map(|k| if i >> k & 1 == 1 { hi[k] } else { lo[k] }))
        .collect()
}

/// The closed mesh on `v` and `t`, its faces reversed when its signed volume
/// is negative.
fn closed(v: Vec<[f64; 3]>, mut t: Vec<[u32; 3]>) -> Mesh {
    let vol: f64 = t
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|i| v[i as usize]);
            a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        })
        .sum();
    if vol < 0.0 {
        for f in &mut t {
            f.swap(1, 2);
        }
    }
    Mesh::new(v, t).unwrap()
}

fn hexa(v: Vec<[f64; 3]>) -> Mesh {
    closed(v, BOX_FACES.to_vec())
}

/// `p` turned to dip `dip` about the strike axis `x` and to strike `strike`
/// about the vertical (both in degrees), then moved to `o`, each step rounded
/// to `f64`.
fn orient(p: [f64; 3], o: [f64; 3], strike: f64, dip: f64) -> [f64; 3] {
    let (st, ct) = (90.0 - dip).to_radians().sin_cos();
    let (y, z) = (p[1] * ct - p[2] * st, p[1] * st + p[2] * ct);
    let (ss, cs) = strike.to_radians().sin_cos();
    [
        o[0] + p[0] * cs - y * ss,
        o[1] + p[0] * ss + y * cs,
        o[2] + z,
    ]
}

fn oriented_box(lo: [f64; 3], hi: [f64; 3], o: [f64; 3], strike: f64, dip: f64) -> Mesh {
    hexa(
        corners(lo, hi)
            .into_iter()
            .map(|p| orient(p, o, strike, dip))
            .collect(),
    )
}

/// A closed solid through rings of equal length, capped by fans from each
/// end ring's centroid.
fn loft(rings: &[Vec<[f64; 3]>]) -> Mesh {
    let n = rings[0].len();
    let mut v: Vec<[f64; 3]> = rings.iter().flatten().copied().collect();
    let mut t = vec![];
    for r in 0..rings.len() - 1 {
        for j in 0..n {
            let a = (r * n + j) as u32;
            let b = (r * n + (j + 1) % n) as u32;
            let c = a + n as u32;
            let d = b + n as u32;
            t.push([a, b, d]);
            t.push([a, d, c]);
        }
    }
    let centroid = |ring: &[[f64; 3]]| {
        [0, 1, 2].map(|k| ring.iter().map(|p| p[k]).sum::<f64>() / ring.len() as f64)
    };
    let (cb, ct) = (v.len() as u32, v.len() as u32 + 1);
    v.push(centroid(&rings[0]));
    v.push(centroid(&rings[rings.len() - 1]));
    let last = ((rings.len() - 1) * n) as u32;
    for j in 0..n {
        let (a, b) = (j as u32, ((j + 1) % n) as u32);
        t.push([cb, b, a]);
        t.push([ct, last + a, last + b]);
    }
    closed(v, t)
}

/// A height field over `[x0, x0 + nx * s] x [y0, y0 + ny * s]`.
fn ground(x0: f64, y0: f64, s: f64, nx: usize, ny: usize, z: impl Fn(f64, f64) -> f64) -> Mesh {
    let mut v = vec![];
    for j in 0..=ny {
        for i in 0..=nx {
            let (x, y) = (x0 + s * i as f64, y0 + s * j as f64);
            v.push([x, y, z(x, y)]);
        }
    }
    let mut t = vec![];
    let id = |i: usize, j: usize| (j * (nx + 1) + i) as u32;
    for j in 0..ny {
        for i in 0..nx {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j + 1));
            if (i + j) % 2 == 0 {
                t.push([a, b, d]);
                t.push([a, d, c]);
            } else {
                t.push([a, b, c]);
                t.push([b, d, c]);
            }
        }
    }
    Mesh::new(v, t).unwrap()
}

/// Rough ground: two long waves, a short one and seeded jitter on every node.
fn rough_ground(rng: &mut Rng, c: [f64; 3], half: f64, s: f64) -> Mesh {
    let n = (2.0 * half / s).ceil() as usize;
    let (p1, p2, p3) = (rng.f(0.0, TAU), rng.f(0.0, TAU), rng.f(0.0, TAU));
    let jitter: Vec<f64> = (0..(n + 1) * (n + 1)).map(|_| rng.f(-0.8, 0.8)).collect();
    let (x0, y0) = (c[0] - half, c[1] - half);
    ground(x0, y0, s, n, n, |x, y| {
        let (i, j) = (
            ((x - x0) / s).round() as usize,
            ((y - y0) / s).round() as usize,
        );
        c[2] + 6.0 * ((x - c[0]) / 90.0 + p1).sin()
            + 4.0 * ((y - c[1]) / 70.0 + p2).sin()
            + 1.5 * ((x + y - c[0] - c[1]) / 17.0 + p3).sin()
            + jitter[j * (n + 1) + i]
    })
}

fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}

/// Six times the signed volume of `m`, exactly.
fn volume6(m: &Mesh) -> Q {
    m.triangles().iter().fold(q(0.0), |acc, t| {
        let [a, b, c] = t.map(|i| m.vertices()[i as usize].map(q));
        acc + &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
            + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
            + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
    })
}

fn solid(m: &Mesh) -> bool {
    m.is_empty() || m.validate(0.0).unwrap().summary.is_solid()
}

fn describe(m: &Mesh) -> String {
    let s = m.validate(0.0).unwrap().summary;
    format!(
        "closed={} boundary={} nm_edges={} nm_vertices={} crossings={} overlaps={} inward={} degenerate={} duplicate_vertices={}",
        s.is_closed,
        s.boundary_edges,
        s.non_manifold_edges,
        s.non_manifold_vertices,
        s.self_intersections,
        s.coplanar_overlaps,
        s.inward_shells,
        s.degenerate_faces,
        s.duplicate_vertices
    )
}

fn bounded(r: &RoundingReport) -> Result<(), String> {
    if r.steps > r.step_bound {
        return Err(format!("steps {} > bound {}", r.steps, r.step_bound));
    }
    for k in 0..3 {
        if r.max_move[k] > r.move_bound[k] {
            return Err(format!("move {:?} > bound {:?}", r.max_move, r.move_bound));
        }
    }
    Ok(())
}

fn within(err: Q, bound: Q) -> bool {
    err <= bound && -err <= bound
}

/// `f` on a pool of `n` threads with its time in seconds, or why there is
/// no result: a panic, or nothing after `limit`.
fn timed<T: Send + 'static>(
    n: usize,
    limit: Duration,
    f: impl FnOnce() -> T + Send + 'static,
) -> Result<(T, f64), String> {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap();
        let t0 = Instant::now();
        let out = pool.install(f);
        let _ = tx.send((out, t0.elapsed().as_secs_f64()));
    });
    rx.recv_timeout(limit).map_err(|e| match e {
        mpsc::RecvTimeoutError::Timeout => format!("no result within {limit:?}"),
        mpsc::RecvTimeoutError::Disconnected => "panicked".into(),
    })
}

/// Counts of checked operations and the failures, per scenario.
struct Ctx {
    scenario: &'static str,
    seed: u64,
    step: usize,
    passed: usize,
    failures: Vec<String>,
    slowest: f64,
    limit: Duration,
}

impl Ctx {
    fn new(scenario: &'static str) -> Self {
        Ctx {
            scenario,
            seed: 0,
            step: 0,
            passed: 0,
            failures: vec![],
            slowest: 0.0,
            limit: Duration::from_secs(env("ATTACK2_LIMIT", 300)),
        }
    }

    fn dump(&self, what: &str, a: &Mesh, b: &Mesh, plane: Option<[f64; 6]>) {
        let Ok(dir) = std::env::var("ATTACK2_DUMP") else {
            return;
        };
        let mut s = String::new();
        for (tag, m) in [("a", a), ("b", b)] {
            s += &format!("{tag} {} {}\n", m.vertices().len(), m.triangles().len());
            for p in m.vertices() {
                s += &format!("{} {} {}\n", p[0], p[1], p[2]);
            }
            for f in m.triangles() {
                s += &format!("{} {} {}\n", f[0], f[1], f[2]);
            }
        }
        if let Some(p) = plane {
            s += &p.map(|x| x.to_string()).join(" ");
            s += "\n";
        }
        let path = format!(
            "{dir}/{}_{}_{}_{what}.txt",
            self.scenario, self.seed, self.step
        );
        std::fs::write(path, s).unwrap();
    }

    fn fail(&mut self, what: &str, msg: String, a: &Mesh, b: &Mesh, plane: Option<[f64; 6]>) {
        let line = format!(
            "{} seed {} step {} {what}: {msg}",
            self.scenario, self.seed, self.step
        );
        eprintln!("FAIL {line}");
        self.dump(what, a, b, plane);
        self.failures.push(line);
    }

    fn run(&mut self, a: &Mesh, b: &Mesh, op: Op, n: usize) -> Result<BooleanResult, String> {
        let (x, y) = (a.clone(), b.clone());
        let out = timed(n, self.limit, move || {
            boolean(&x, &y, op, &BooleanOptions::default()).map_err(|e| e.to_string())
        });
        let (r, s) = out?;
        self.slowest = self.slowest.max(s);
        r
    }

    /// `op` of `a` and `b`, checked; `None` when any check fails.
    fn op(&mut self, what: &str, a: &Mesh, b: &Mesh, op: Op) -> Option<Mesh> {
        self.step += 1;
        let what = format!("{what}_{op:?}");
        if !solid(a) || !solid(b) {
            panic!(
                "{} seed {}: scenario built a non-solid input for {what}",
                self.scenario, self.seed
            );
        }
        let partner = if op == Op::Intersection {
            Op::Union
        } else {
            Op::Intersection
        };
        let mut got = vec![];
        for o in [op, partner] {
            match self.run(a, b, o, 1) {
                Err(e) => {
                    self.fail(&what, format!("{o:?}: {e}"), a, b, None);
                    return None;
                }
                Ok(r) => {
                    if !solid(&r.mesh) {
                        let d = describe(&r.mesh);
                        self.fail(&what, format!("{o:?}: not a solid: {d}"), a, b, None);
                        return None;
                    }
                    if let Err(e) = bounded(&r.rounding) {
                        self.fail(&what, format!("{o:?}: {e}"), a, b, None);
                        return None;
                    }
                    got.push(r);
                }
            }
        }
        let (r, p) = (&got[0], &got[1]);
        let exact = match op {
            Op::Difference => volume6(a),
            _ => volume6(a) + volume6(b),
        };
        let err = volume6(&r.mesh) + volume6(&p.mesh) - exact;
        let bound = q(6.0 * r.rounding.volume_bound) + q(6.0 * p.rounding.volume_bound);
        if !within(err.clone(), bound.clone()) {
            self.fail(
                &what,
                format!("volume identity: 6 err {err} > 6 bound {bound}"),
                a,
                b,
                None,
            );
            return None;
        }
        match self.run(a, b, op, 8) {
            Ok(r8)
                if r8.mesh == r.mesh
                    && r8.face_source == r.face_source
                    && r8.rounding == r.rounding => {}
            other => {
                let m = match other {
                    Err(e) => e,
                    Ok(_) => "differs".into(),
                };
                self.fail(&what, format!("8 threads: {m}"), a, b, None);
                return None;
            }
        }
        self.passed += 1;
        Some(got.swap_remove(0).mesh)
    }

    /// The two pieces of a split or clip, checked against `m`.
    fn pieces(
        &mut self,
        what: &str,
        m: &Mesh,
        other_mesh: &Mesh,
        plane: Option<[f64; 6]>,
        f: impl Fn() -> Result<(Mesh, Mesh, [RoundingReport; 2]), String>
        + Send
        + Sync
        + Clone
        + 'static,
    ) -> Option<(Mesh, Mesh)> {
        self.step += 1;
        assert!(
            solid(m),
            "{} seed {}: non-solid input for {what}",
            self.scenario,
            self.seed
        );
        let g = f.clone();
        let one = match timed(1, self.limit, f) {
            Err(e) => Err(e),
            Ok((r, s)) => {
                self.slowest = self.slowest.max(s);
                r
            }
        };
        let (lo, hi, rep) = match one {
            Ok(x) => x,
            Err(e) => {
                self.fail(what, e, m, other_mesh, plane);
                return None;
            }
        };
        for x in [&lo, &hi] {
            if !solid(x) {
                let d = describe(x);
                self.fail(what, format!("not a solid: {d}"), m, other_mesh, plane);
                return None;
            }
        }
        for r in &rep {
            if let Err(e) = bounded(r) {
                self.fail(what, e, m, other_mesh, plane);
                return None;
            }
        }
        let err = volume6(&lo) + volume6(&hi) - volume6(m);
        let bound = q(6.0 * rep[0].volume_bound) + q(6.0 * rep[1].volume_bound);
        if !within(err.clone(), bound.clone()) {
            self.fail(
                what,
                format!("volume identity: 6 err {err} > 6 bound {bound}"),
                m,
                other_mesh,
                plane,
            );
            return None;
        }
        match timed(8, self.limit, g) {
            Ok((Ok((l8, h8, r8)), _)) if l8 == lo && h8 == hi && r8 == rep => {}
            other => {
                let e = match other {
                    Err(e) | Ok((Err(e), _)) => e,
                    Ok(_) => "differs".into(),
                };
                self.fail(what, format!("8 threads: {e}"), m, other_mesh, plane);
                return None;
            }
        }
        self.passed += 1;
        Some((lo, hi))
    }

    fn clip(&mut self, what: &str, m: &Mesh, surface: &Mesh) -> Option<(Mesh, Mesh)> {
        let (x, y) = (m.clone(), surface.clone());
        self.pieces(&format!("{what}_clip"), m, surface, None, move || {
            clip_by_surface(&x, &y, &BooleanOptions::default())
                .map(|r| (r.below, r.above, r.rounding))
                .map_err(|e| e.to_string())
        })
    }

    fn split(
        &mut self,
        what: &str,
        m: &Mesh,
        origin: [f64; 3],
        normal: [f64; 3],
    ) -> Option<(Mesh, Mesh)> {
        let x = m.clone();
        let plane = [
            origin[0], origin[1], origin[2], normal[0], normal[1], normal[2],
        ];
        self.pieces(&format!("{what}_split"), m, m, Some(plane), move || {
            split_by_plane_report(&x, origin, normal).map_err(|e| e.to_string())
        })
    }

    fn finish(self) {
        let peak = std::fs::read_to_string("/proc/self/status")
            .unwrap_or_default()
            .lines()
            .find(|l| l.starts_with("VmHWM"))
            .unwrap_or("VmHWM: unknown")
            .to_string();
        eprintln!(
            "SUMMARY {}: {} operations passed, {} failed, slowest {:.2} s, {peak}",
            self.scenario,
            self.passed,
            self.failures.len(),
            self.slowest
        );
        assert!(self.failures.is_empty(), "{:#?}", self.failures);
    }
}

fn seeds() -> std::ops::Range<u64> {
    0..env("ATTACK2_SEEDS", 1)
}

fn large() -> bool {
    env("ATTACK2_SIZE", 0) > 0
}

fn dip(rng: &mut Rng) -> f64 {
    match rng.i(0, 5) {
        0 => 90.0,
        1 => 89.0 + rng.f(0.0, 1.0),
        2 => 45.0,
        _ => rng.f(35.0, 90.0),
    }
}

/// Stopes along strike in two levels, all between the common footwall
/// `y = 0` and hanging wall `y = t`, united one at a time; then a footwall
/// drive, a hanging wall drive and crosscuts, some on the stopes' end walls,
/// taken away one at a time from the block.
fn wall_bound_block(cx: &mut Ctx, rng: &mut Rng) {
    let (o, strike, dip) = (rng.utm(), rng.f(0.0, 360.0), dip(rng));
    let n = if large() { 10 } else { 4 };
    let t = rng.f(4.0, 12.0);
    let mut ends = vec![0.0];
    for _ in 0..n {
        let w = ends.last().unwrap() + rng.f(8.0, 16.0);
        ends.push(w);
    }
    let split = rng.f(15.0, 25.0);
    let top = split + rng.f(15.0, 25.0);
    let mut block: Option<Mesh> = None;
    for k in 0..n {
        for (z0, z1) in [(0.0, split), (split, top)] {
            let (x0, x1) = if z0 == 0.0 {
                (ends[k], ends[k + 1])
            } else {
                let s = (ends[k + 1] - ends[k]) * 0.37;
                (ends[k] + s, ends[k + 1] + s)
            };
            let s = oriented_box([x0, 0.0, z0], [x1, t, z1], o, strike, dip);
            block = Some(match block {
                None => s,
                Some(b) => match cx.op("stope", &b, &s, Op::Union) {
                    Some(m) => m,
                    None => return,
                },
            });
        }
    }
    let mut block = block.unwrap();
    let len = *ends.last().unwrap();
    let mut drives = vec![
        ([-8.0, -3.0, -0.5], [len + 8.0, 1.0, 4.0]),
        (
            [-8.0, t - 1.5, split - 2.0],
            [len + 8.0, t + 3.0, split + 2.5],
        ),
    ];
    for k in 1..n {
        let x = if k % 2 == 0 {
            ends[k]
        } else {
            ends[k] + rng.f(-3.0, 3.0)
        };
        drives.push(([x - 2.5, -15.0, 0.0], [x + 2.5, t + 15.0, 4.5]));
    }
    for (lo, hi) in drives {
        let d = oriented_box(lo, hi, o, strike, dip);
        match cx.op("drive", &block, &d, Op::Difference) {
            Some(m) => block = m,
            None => return,
        }
    }
}

#[test]
fn stope_blocks_between_common_walls_minus_drive_networks() {
    let mut cx = Ctx::new("wall_bound_block");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        wall_bound_block(&mut cx, &mut Rng::new(seed));
    }
    cx.finish();
}

/// A sublevel layout: stopes along strike, across in primaries and
/// secondaries and in levels, each level's stopes offset along strike from
/// the level below, oriented and united one at a time.
fn sublevel_layout(cx: &mut Ctx, rng: &mut Rng) -> Option<Mesh> {
    let (o, strike, dip) = (rng.utm(), rng.f(0.0, 360.0), dip(rng));
    let (along, across, levels) = if large() { (10, 3, 4) } else { (3, 2, 2) };
    let ys: Vec<f64> = (0..=across).map(|j| 9.0 * j as f64).collect();
    let zs: Vec<f64> = (0..=levels).map(|k| 25.0 * k as f64).collect();
    let mut block: Option<Mesh> = None;
    for k in 0..levels {
        let shift = rng.f(0.0, 6.0);
        for j in 0..across {
            for i in 0..along {
                let x0 = shift + 15.0 * i as f64;
                let s = oriented_box(
                    [x0, ys[j], zs[k]],
                    [x0 + 15.0, ys[j + 1], zs[k + 1]],
                    o,
                    strike,
                    dip,
                );
                block = Some(match block {
                    None => s,
                    Some(b) => cx.op("stope", &b, &s, Op::Union)?,
                });
            }
        }
    }
    block
}

#[test]
fn sublevel_layouts_united_one_stope_at_a_time() {
    let mut cx = Ctx::new("sublevel_layout");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        sublevel_layout(&mut cx, &mut Rng::new(100 + seed));
    }
    cx.finish();
}

/// Hexahedra along a helix of radius `r` about `c`, from angle `a0` over
/// `turns` with fall `grade` per metre, `w` wide and `h` high, in segments of
/// `step` radians; consecutive segments share their end faces.
fn helix(
    c: [f64; 3],
    r: f64,
    a0: f64,
    turns: f64,
    grade: f64,
    (w, h): (f64, f64),
    step: f64,
) -> Vec<Mesh> {
    let n = (turns * TAU / step).round() as usize;
    let station = |i: usize| {
        let a = a0 + step * i as f64;
        let z = c[2] - grade * r * step * i as f64;
        let (s, k) = a.sin_cos();
        let at = |rr: f64, zz: f64| [c[0] + rr * k, c[1] + rr * s, zz];
        [
            at(r - w / 2.0, z),
            at(r + w / 2.0, z),
            at(r - w / 2.0, z + h),
            at(r + w / 2.0, z + h),
        ]
    };
    (0..n)
        .map(|i| {
            let (p, q) = (station(i), station(i + 1));
            hexa(vec![p[0], q[0], p[1], q[1], p[2], q[2], p[3], q[3]])
        })
        .collect()
}

/// A helical decline united segment by segment, with crosscuts and a level
/// drive at chosen stations, each union fed with the previous result.
fn decline_with_crosscuts(cx: &mut Ctx, rng: &mut Rng) {
    let c = rng.utm();
    let turns = if large() { 2.0 } else { 0.5 };
    let step = rng.f(8.0, 15.0).to_radians();
    let segs = helix(
        c,
        rng.f(25.0, 40.0),
        rng.f(0.0, TAU),
        turns,
        1.0 / 7.0,
        (5.0, 5.5),
        step,
    );
    let mut dec = segs[0].clone();
    for (i, s) in segs.iter().enumerate().skip(1) {
        dec = match cx.op("segment", &dec, s, Op::Union) {
            Some(m) => m,
            None => return,
        };
        if i % 5 == 2 {
            let strike = rng.f(0.0, 360.0);
            let p = s.vertices()[0];
            let cut = oriented_box([-3.0, -2.5, -0.3], [40.0, 2.5, 4.5], p, strike, 90.0);
            dec = match cx.op("crosscut", &dec, &cut, Op::Union) {
                Some(m) => m,
                None => return,
            };
        }
    }
}

#[test]
fn helical_declines_united_with_crosscuts() {
    let mut cx = Ctx::new("decline");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        decline_with_crosscuts(&mut cx, &mut Rng::new(200 + seed));
    }
    cx.finish();
}

/// A benched pit about `c`: `nb` benches of height `h` from `c[2]`, the toe
/// of bench `j` at radius `r(j, angle)`, faces at 70 degrees, berms of 6 m,
/// walls through `n` angles from `phase`, and the top ring at `top`.
fn pit(
    c: [f64; 3],
    n: usize,
    phase: f64,
    nb: usize,
    h: f64,
    top: f64,
    r: impl Fn(usize, f64) -> f64,
) -> Mesh {
    let run = h / 70f64.to_radians().tan();
    let ring = |rad: &dyn Fn(f64) -> f64, z: f64| -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let a = phase + TAU * i as f64 / n as f64;
                let rr = rad(a);
                [c[0] + rr * a.cos(), c[1] + rr * a.sin(), z]
            })
            .collect()
    };
    let mut rings = vec![];
    for j in 0..nb {
        let z = c[2] + h * j as f64;
        rings.push(ring(&|a| r(j, a), z));
        rings.push(ring(&|a| r(j, a) + run, z + h));
    }
    rings.push(ring(
        &|a| r(nb - 1, a) + run + 2.0 * (top - c[2] - h * nb as f64),
        top,
    ));
    loft(&rings)
}

/// A pit with a ramp, the ramp's segments united one at a time, clipped by
/// rough ground; a dipping orebody clipped by the same ground, minus the pit,
/// and the ore inside the pit.
fn pit_with_ramp(cx: &mut Ctx, rng: &mut Rng) {
    let c = rng.utm();
    let (nb, h) = (rng.i(4, 6) as usize, 10.0);
    let n = rng.i(8, 24) as usize;
    let r0 = rng.f(25.0, 40.0);
    let wob = rng.f(0.0, 0.15);
    let top = c[2] + h * nb as f64 + 15.0;
    let radius = move |j: usize, a: f64| (r0 + 9.6 * j as f64) * (1.0 + wob * (3.0 * a).cos());
    let mut p = pit(c, n, rng.f(0.0, TAU), nb, h, top, radius);
    let a0 = rng.f(0.0, TAU);
    let rr = r0 + 9.6 * (nb as f64) / 2.0;
    let segs = helix(
        [c[0], c[1], c[2] + h * nb as f64],
        rr,
        a0,
        0.6,
        0.1,
        (12.0, 12.0),
        rng.f(10.0, 20.0).to_radians(),
    );
    for s in &segs {
        p = match cx.op("ramp", &p, s, Op::Union) {
            Some(m) => m,
            None => return,
        };
    }
    let g = rough_ground(rng, [c[0], c[1], c[2] + h * nb as f64 - 3.0], 260.0, 15.0);
    let Some((pit_below, _)) = cx.clip("pit", &p, &g) else {
        return;
    };
    let ore = oriented_box(
        [-120.0, -8.0, -60.0],
        [120.0, 8.0, 120.0],
        [c[0] + rng.f(-20.0, 20.0), c[1] + rng.f(-20.0, 20.0), c[2]],
        rng.f(0.0, 360.0),
        rng.f(40.0, 80.0),
    );
    let Some((ore_below, _)) = cx.clip("ore", &ore, &g) else {
        return;
    };
    cx.op("ore_minus_pit", &ore_below, &pit_below, Op::Difference);
    cx.op("ore_in_pit", &ore_below, &pit_below, Op::Intersection);
}

#[test]
fn pits_with_ramps_against_orebodies_and_rough_ground() {
    let mut cx = Ctx::new("pit_with_ramp");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        pit_with_ramp(&mut cx, &mut Rng::new(300 + seed));
    }
    cx.finish();
}

/// Pushback stages of one pit growing in one sector, each clipped by the
/// same rough ground; each stage minus the one before, and the running union
/// of the stages, fed forward.
fn pushbacks(cx: &mut Ctx, rng: &mut Rng) {
    let c = rng.utm();
    let stages = if large() { 12 } else { 4 };
    let (nb, h) = (4, 12.0);
    let n = rng.i(12, 32) as usize;
    let phase = rng.f(0.0, TAU);
    let toward = rng.f(0.0, TAU);
    let r0 = rng.f(30.0, 50.0);
    let grow = rng.f(4.0, 12.0);
    let top = c[2] + h * nb as f64 + 12.0;
    let g = rough_ground(
        rng,
        [c[0], c[1], c[2] + h * nb as f64 - 4.0],
        300.0 + grow * stages as f64,
        20.0,
    );
    let mut prev: Option<Mesh> = None;
    let mut all: Option<Mesh> = None;
    for k in 0..stages {
        let r = move |j: usize, a: f64| {
            let d = (a - toward).cos();
            r0 + 11.0 * j as f64
                + if d > 0.5 {
                    grow * k as f64 * (d - 0.5) * 2.0
                } else {
                    0.0
                }
        };
        let p = pit(c, n, phase, nb, h, top, r);
        let Some((stage, _)) = cx.clip("stage", &p, &g) else {
            return;
        };
        if let Some(q) = &prev
            && cx.op("increment", &stage, q, Op::Difference).is_none()
        {
            return;
        }
        all = Some(match all {
            None => stage.clone(),
            Some(u) => match cx.op("cumulative", &u, &stage, Op::Union) {
                Some(m) => m,
                None => return,
            },
        });
        prev = Some(stage);
    }
}

#[test]
fn pushback_sequences_clipped_by_the_same_ground() {
    let mut cx = Ctx::new("pushbacks");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        pushbacks(&mut cx, &mut Rng::new(400 + seed));
    }
    cx.finish();
}

/// A box whose faces are split into `m x m` grids, corners in the order of
/// `BOX_FACES`, each point mapped by `f`.
fn gridded_box(lo: [f64; 3], hi: [f64; 3], m: usize, f: impl Fn([f64; 3]) -> [f64; 3]) -> Mesh {
    let mut v: Vec<[f64; 3]> = vec![];
    let mut index = std::collections::HashMap::new();
    let mut t = vec![];
    let lerp = |k: usize, i: usize| lo[k] + (hi[k] - lo[k]) * i as f64 / m as f64;
    let mut id = |p: [usize; 3]| -> u32 {
        *index.entry(p).or_insert_with(|| {
            v.push(f([lerp(0, p[0]), lerp(1, p[1]), lerp(2, p[2])]));
            (v.len() - 1) as u32
        })
    };
    for axis in 0..3 {
        let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
        for side in [0, m] {
            for i in 0..m {
                for j in 0..m {
                    let at = |a: usize, b: usize| {
                        let mut p = [0; 3];
                        p[axis] = side;
                        p[u] = a;
                        p[w] = b;
                        p
                    };
                    let (a, b, c, d) = (
                        id(at(i, j)),
                        id(at(i + 1, j)),
                        id(at(i, j + 1)),
                        id(at(i + 1, j + 1)),
                    );
                    if side == 0 {
                        t.push([a, c, b]);
                        t.push([b, c, d]);
                    } else {
                        t.push([a, b, c]);
                        t.push([b, d, c]);
                    }
                }
            }
        }
    }
    closed(v, t)
}

/// A row of stopes oriented to strike and dip, each sharing its whole end
/// wall with the next: alternately two triangles against a gridded wall
/// (vertices that do not match) and gridded against gridded (vertices that
/// match), united one at a time; then each stope taken away from the row
/// again, in another order.
fn shared_walls(cx: &mut Ctx, rng: &mut Rng) {
    let (o, strike, dip) = (rng.utm(), rng.f(0.0, 360.0), dip(rng));
    let n = if large() { 12 } else { 5 };
    let (t, h) = (rng.f(4.0, 10.0), rng.f(15.0, 30.0));
    let stopes: Vec<Mesh> = (0..n)
        .map(|k| {
            let (lo, hi) = ([12.0 * k as f64, 0.0, 0.0], [12.0 * (k + 1) as f64, t, h]);
            let m = [1, 3, 3, 4][(k + rng.i(0, 1) as usize) % 4];
            gridded_box(lo, hi, m, |p| orient(p, o, strike, dip))
        })
        .collect();
    let mut row = stopes[0].clone();
    for s in &stopes[1..] {
        row = match cx.op("stope", &row, s, Op::Union) {
            Some(m) => m,
            None => return,
        };
    }
    for k in (1..n).step_by(2) {
        row = match cx.op("take", &row, &stopes[k], Op::Difference) {
            Some(m) => m,
            None => return,
        };
    }
}

#[test]
fn stopes_sharing_whole_walls_with_and_without_matching_vertices() {
    let mut cx = Ctx::new("shared_walls");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        shared_walls(&mut cx, &mut Rng::new(500 + seed));
    }
    cx.finish();
}

/// A long chain on one orebody: drives, raises and stopes taken away and
/// backfill put back, each result fed into the next, then `x = (x | b) - b`
/// three times for the last body.
fn long_chain(cx: &mut Ctx, rng: &mut Rng) {
    let (o, strike, dip) = (rng.utm(), rng.f(0.0, 360.0), rng.f(50.0, 90.0));
    let mut x = oriented_box([0.0, 0.0, 0.0], [120.0, 20.0, 80.0], o, strike, dip);
    let steps = if large() { 40 } else { 10 };
    let mut last = None;
    for k in 0..steps {
        let (lo, hi) = match k % 4 {
            0 => {
                let z = rng.f(0.0, 70.0);
                ([-5.0, rng.f(-2.0, 18.0), z], [125.0, 0.0, z + 4.5])
            }
            1 => {
                let x0 = rng.f(0.0, 110.0);
                ([x0, -5.0, rng.f(0.0, 40.0)], [x0 + 3.0, 25.0, 0.0])
            }
            2 => {
                let x0 = rng.f(0.0, 100.0);
                (
                    [x0, rng.f(-1.0, 1.0), rng.f(-5.0, 30.0)],
                    [x0 + 15.0, 20.0, 0.0],
                )
            }
            _ => {
                let x0 = rng.f(0.0, 100.0);
                ([x0, 0.0, rng.f(0.0, 40.0)], [x0 + 12.0, 20.0, 0.0])
            }
        };
        let hi = [hi[0], hi[1].max(lo[1] + 3.0), hi[2].max(lo[2] + 4.0)];
        let b = oriented_box(lo, hi, o, strike + rng.f(-0.5, 0.5), dip + rng.f(-0.5, 0.5));
        let op = if k % 4 == 3 {
            Op::Union
        } else {
            Op::Difference
        };
        x = match cx.op("chain", &x, &b, op) {
            Some(m) => m,
            None => return,
        };
        last = Some(b);
    }
    let b = last.unwrap();
    for _ in 0..3 {
        let Some(u) = cx.op("again", &x, &b, Op::Union) else {
            return;
        };
        x = match cx.op("again", &u, &b, Op::Difference) {
            Some(m) => m,
            None => return,
        };
    }
}

#[test]
fn long_chains_on_one_orebody_feed_each_result_forward() {
    let mut cx = Ctx::new("long_chain");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        long_chain(&mut cx, &mut Rng::new(600 + seed));
    }
    cx.finish();
}

fn nudge(mut x: f64, n: i64) -> f64 {
    for _ in 0..n.abs() {
        x = if n > 0 { x.next_up() } else { x.next_down() };
    }
    x
}

/// Stopes along a northing at UTM coordinates whose end walls overlap, touch
/// or leave a gap of up to three ulps, optionally turned to strike by a few
/// degrees or to about 45 degrees, united one at a time; then development
/// crossing every slot, a crosscut whose wall lies on a slot wall, a split
/// through a slot and a split across all of them.
fn slots(cx: &mut Ctx, rng: &mut Rng) {
    let o = rng.utm();
    let turn = match rng.i(0, 2) {
        0 => 0.0,
        1 => rng.f(-5.0, 5.0),
        _ => 45.0 + rng.f(-0.01, 0.01),
    };
    let dip = if rng.i(0, 1) == 0 {
        90.0
    } else {
        rng.f(60.0, 89.9)
    };
    let n = if large() { 10 } else { 4 };
    let (t, h) = (rng.f(5.0, 12.0), rng.f(20.0, 30.0));
    let exact = turn == 0.0 && dip == 90.0;
    let place = |p: [f64; 3]| {
        if exact {
            p
        } else {
            orient([p[0] - o[0], p[1] - o[1], p[2] - o[2]], o, turn, dip)
        }
    };
    let mut y = o[1];
    let mut stopes = vec![];
    let mut walls = vec![];
    let mut global = vec![];
    for _ in 0..n {
        let y1 = y + rng.f(8.0, 15.0);
        let lo = [o[0], y, o[2]];
        let hi = [o[0] + t, y1, o[2] + h];
        stopes.push(hexa(corners(lo, hi).into_iter().map(place).collect()));
        walls.push(y1 - o[1]);
        global.push(y1);
        y = nudge(y1, rng.i(-3, 3));
    }
    let y = y - o[1];
    let mut x = stopes[0].clone();
    for s in &stopes[1..] {
        x = match cx.op("stope", &x, s, Op::Union) {
            Some(m) => m,
            None => return,
        };
    }
    let len = y + 10.0;
    for (lo, hi) in [
        ([-3.0, -5.0, 0.5], [2.0, len, 5.0]),
        ([t - 2.0, -5.0, h - 6.0], [t + 3.0, len, h - 1.5]),
        ([t / 2.0 - 1.5, -5.0, -2.0], [t / 2.0 + 1.5, len, h + 2.0]),
    ] {
        let d = oriented_box(lo, hi, o, turn, dip);
        let Some(m) = cx.op("drive", &x, &d, Op::Difference) else {
            return;
        };
        x = m;
    }
    let k = rng.i(0, n as i64 - 2) as usize;
    let (zc, gw) = (o[2] + h / 2.0, global[k]);
    let cut = corners(
        [o[0] - 3.0, gw, zc - 2.0],
        [o[0] + t + 3.0, gw + 6.0, zc + 2.0],
    );
    let cut = hexa(cut.into_iter().map(place).collect());
    let Some(m) = cx.op("on_wall", &x, &cut, Op::Difference) else {
        return;
    };
    x = m;
    let w = walls[k];
    let at = orient([t / 2.0, w, h / 2.0], o, turn, dip);
    let across = orient([0.0, 1.0, 0.0], [0.0; 3], turn, dip);
    cx.split("on_slot", &x, at, across);
    let along = orient([1.0, 0.0, 0.3], [0.0; 3], turn, dip);
    cx.split("along", &x, at, along);
}

#[test]
fn slots_a_few_ulps_wide_crossed_by_development() {
    let mut cx = Ctx::new("slots");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        slots(&mut cx, &mut Rng::new(700 + seed));
    }
    cx.finish();
}

/// A stope block split by bench planes through its own vertices and by
/// inclined planes through the corners of its faces, intersected with a box
/// whose top is each bench plane, and clipped by a rough fault surface.
fn splits_and_faults(cx: &mut Ctx, rng: &mut Rng) {
    let Some(block) = sublevel_layout(cx, rng) else {
        return;
    };
    let count = env("ATTACK2_SPLITS", if large() { 12 } else { 3 });
    for k in 0..count {
        let v = block.vertices()[rng.i(0, block.vertices().len() as i64 - 1) as usize];
        let normal = if k % 3 == 0 {
            [0.0, 0.0, 1.0]
        } else {
            let a = rng.f(0.0, TAU);
            let d = rng.f(20.0, 80.0).to_radians();
            [d.sin() * a.cos(), d.sin() * a.sin(), d.cos()]
        };
        cx.split("block", &block, v, normal);
        if k % 3 == 0 {
            let lo = [v[0] - 500.0, v[1] - 500.0, v[2] - 500.0];
            let below = hexa(corners(lo, [v[0] + 500.0, v[1] + 500.0, v[2]]));
            cx.op("level_box", &block, &below, Op::Intersection);
        }
    }
    let c = block.vertices()[0];
    let lift = rng.f(10.0, 40.0);
    let g = rough_ground(rng, [c[0], c[1], c[2] + lift], 150.0, 9.0);
    cx.clip("fault", &block, &g);
}

#[test]
fn stope_blocks_split_by_planes_through_their_vertices_and_clipped_by_faults() {
    let mut cx = Ctx::new("splits");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        splits_and_faults(&mut cx, &mut Rng::new(800 + seed));
    }
    cx.finish();
}

/// Stopes stacked in levels at UTM coordinates without turning, each level's
/// floor a few ulps above, on or below the roof of the level under it and its
/// walls a few ulps off those below, united one at a time; then a raise
/// through every level, a drive whose floor lies on a level's floor, a split
/// on a level's floor and a clip by flat ground at that elevation.
fn stacked_slots(cx: &mut Ctx, rng: &mut Rng) {
    let o = rng.utm();
    let levels = if large() { 8 } else { 3 };
    let (l, t) = (rng.f(10.0, 25.0), rng.f(5.0, 12.0));
    let mut z = o[2];
    let mut floors = vec![];
    let mut x: Option<Mesh> = None;
    for _ in 0..levels {
        let z1 = z + rng.f(15.0, 30.0);
        let (x0, y0) = (nudge(o[0], rng.i(-3, 3)), nudge(o[1], rng.i(-3, 3)));
        let s = hexa(corners([x0, y0, z], [x0 + l, y0 + t, z1]));
        floors.push(z);
        x = Some(match x {
            None => s,
            Some(b) => match cx.op("level", &b, &s, Op::Union) {
                Some(m) => m,
                None => return,
            },
        });
        z = nudge(z1, rng.i(-3, 3));
    }
    let mut x = x.unwrap();
    let (cx0, cy0) = (o[0] + l / 2.0, o[1] + t / 2.0);
    let f = floors[rng.i(1, levels as i64 - 1) as usize];
    for (lo, hi) in [
        (
            [cx0 - 1.5, cy0 - 1.5, o[2] - 5.0],
            [cx0 + 1.5, cy0 + 1.5, z + 5.0],
        ),
        (
            [o[0] - 5.0, cy0 - 2.0, f],
            [o[0] + l + 5.0, cy0 + 2.0, f + 4.0],
        ),
    ] {
        x = match cx.op("development", &x, &hexa(corners(lo, hi)), Op::Difference) {
            Some(m) => m,
            None => return,
        };
    }
    cx.split("floor", &x, [cx0, cy0, f], [0.0, 0.0, 1.0]);
    let g = ground(o[0] - 20.0, o[1] - 20.0, 8.0, 8, 6, |_, _| f);
    cx.clip("flat_ground", &x, &g);
}

#[test]
fn stacked_levels_a_few_ulps_apart_crossed_by_raises() {
    let mut cx = Ctx::new("stacked_slots");
    for seed in seeds() {
        cx.seed = seed;
        cx.step = 0;
        stacked_slots(&mut cx, &mut Rng::new(900 + seed));
    }
    cx.finish();
}
