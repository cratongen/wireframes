# Rounding exact results to f64

Status: method C is implemented in `crates/boolean/src/round.rs`; section 7
records where the implementation departs from this design and why, and section
10 traces the two refusals left after section 9 and decides the next change. The API
contract stays in `spec.md`; this file records why booleans, `clip_by_surface`
and `split_by_plane` refused results, which materialization method replaced the
refusal, and how we prove it.

## 1. The refusals, reproduced and reduced

Every case from README (M4 known gaps), CHANGELOG (Known issues) and spec
tolerance rule 5 was rerun on 2026-11-19 at commit `ee20f54`. Each one that still
raises `Error::Kernel` was cut down to the smallest input found that still fails
and committed as an ignored test in `crates/boolean/tests/rounding.rs`, with its
inputs in `crates/boolean/tests/data/rounding/*.txt` (decimal text that parses back
to the same `f64`). The running test
`reduced_rounding_cases_return_solids_or_kernel_errors` checks that every reduced
case fails as a kernel error, never as an input error or a panic.

| Listed case | Reproduced | Reduced test (faces of a + b) |
|---|---|---|
| Prisms whose caps differ by one ulp | `b - a` refuses; union, intersection and `a - b` pass | `difference_of_prisms_whose_caps_differ_by_one_ulp_is_empty_or_a_solid` (12 + 12) |
| Pits and pushbacks clipped by the same topography | boxes clipped by one DTM: difference refuses; square pit and pushback, union refuses | `booleans_of_boxes_clipped_by_the_same_topography` (20 + 20), `booleans_of_a_triangular_pit_and_its_pushback_clipped_by_the_same_ground` (32 + 42, 3 ring directions, 1 bench) |
| Touching pushback slivers (`pushback_pair_touching.npz`) | union, intersection and `pit - pushback` refuse | `booleans_of_touching_pushback_slivers` (18 + 28) |
| Runaway pushback pair (`pushback_pair_runaway.npz`) | all four refuse; the union is not even closed (4 boundary edges) | four tests `booleans_of_the_runaway_pushback_pair_reduced_for_*` (124 + 136 to 164) |
| Face-sharing stopes rotated to strike | union refuses | `union_of_a_tessellated_stope_and_a_box_rotated_to_strike` (68 + 12: an 8 x 1 x 1 grid box and a plain box) |
| Union 11 of 12 in example 07 | drive united with crosscuts 1 and 2, then with crosscut 11, refuses | `union_of_a_rounded_drive_with_a_crosscut` (156 + 52 after cropping) |
| Octagonal pit clipped by flat ground at UTM | flat ground at 119 elevations from 200.5 to 259.5 m: no refusal. Tilted ground (200 seeded planes, seed 20261119, slopes up to 0.2): octagon 9, hexagon 0, square 1 | `square_pit_at_utm_clipped_by_tilted_ground_is_two_solids` (30 + 2) |
| 188 of 648 inclined splits of `sphere(2)` | 188 of 648; `sphere(1)` refuses 10 of 320 | `split_of_an_icosphere_through_its_centre_is_two_solids` (80 faces, plane through vertex 12 with normal (1, 1, 1)) |

The flat-ground octagon clip no longer refuses; README and CHANGELOG still list it.
The tilted sweep is the closest case that still fails.

Reduction worked in two ways. Cropping both inputs with six plane splits around a
problem face shrank the crosscut chain. It could not shrink the compact pairs,
because the caps it adds outnumber the faces it removes. Generator parameters
shrank the rest: fewer ring directions and benches, fewer stope panels, and the
smallest subset of earlier crosscuts.

### Mechanisms

For each self-intersecting or overlapping face pair of the rounded result, a
diagnostic matched both faces to the input face whose plane holds them, then
measured the angle between those planes and how far the nearer vertex lies from
the other plane, in ulps of the largest coordinate. Four mechanisms cover every
pair:

- **M1, sub-ulp shell.** The whole exact shell is thinner than the f64 spacing.
  The difference of the ulp prisms rounds to two faces on three vertices, an
  inward shell of zero volume.
- **M2, two inputs that agree only to rounding.** Faces of `a` and `b` lie at
  angles below 1e-9 rad or with a vertex within one ulp of the other plane. The
  exact result holds slivers and pockets between them thinner than an ulp, and
  rounding the implicit points on their intersection curves folds the slivers
  across each other. Clipped boxes, pit and pushback, both pushback pairs.
- **M3, one input whose faces agree only to rounding.** A tessellated wall
  rotated to strike, or a drive that earlier unions rounded, has neighbouring
  faces that are planar only to rounding. Implicit points where the other input
  crosses them round off the wall, and the sub-triangles of adjacent input faces
  overlap. Rotated stopes (all 22 problem pairs inside one input), crosscut chain.
- **M4, a vertex within ulps of a transverse cut.** An input vertex lies within an
  ulp of a cutting face or plane it does not touch exactly. The crossing point on
  each of its edges lies within an ulp of the vertex, and rounding puts it on the
  wrong side (angles near 1.5 rad, distances 0 to 0.1 ulp). Tilted ground, the
  icosphere split.

## 2. Where validity is lost today

The arrangement (`arrangement.rs`), the implicit points (`hp.rs`) and the
classification (`classify.rs`) are exact. Validity is lost at these places:

1. `ops.rs:115`, in `assemble`: `arr.pts[p].to_f64()` rounds each coordinate of
   each implicit point to the nearest `f64` (`hp.rs:110`, `round_quotient`). Each
   point moves by at most half an ulp per coordinate. Nothing checks that faces
   which were disjoint stay disjoint after the move; a vertex closer than one
   ulp to a face it does not touch can cross it.
2. `ops.rs:132`: `collapse_rounded` runs on one shell at a time. Problems between
   two shells (touching shells of a union, a sliver shell inside another) are
   never seen.
3. `ops.rs:165-417`, `collapse_rounded` itself. It merges vertices with equal
   coordinates, cancels opposite faces on the same three vertices, splits away
   faces that are collinear, folded or turned (only into a neighbour from the same
   input face), and collapses short edges. The last stage (lines 250-385) is a
   heuristic: it calls `validate`, retries splits around problem faces, and at
   line 370 undoes a round that does not lower the problem count, then returns
   the shell as it is. It has no operation for a vertex against a face of
   another input face (M2, M3, M4) and none that removes a shell (M1). Two
   same-orientation copies of a face survive it (the M1 case).
4. `ops.rs:714`, `solid_result`: the only guard. It runs `validate` and turns any
   problem into `Error::Kernel`.
5. `clip.rs`: the same `assemble` and `solid_result` path, so it loses validity at
   the same places.
6. `slice.rs:320`: `split_by_plane` rounds every section point to `f64` before it
   triangulates the cap (`slice.rs:428`), so the cap and the cut side faces are
   built on rounded points, not exact ones. `coincident` (`slice.rs:150`) merges
   section points whose rounded projections are equal. `clean`
   (`slice.rs:172-215`) runs `collapse_rounded` without normals, so the
   validate-driven stage never runs for splits. `slice.rs:495` refuses.

## 3. Methods compared

Notation: `p` is an exact output vertex, `r(p)` its coordinate-wise nearest
`f64`, `u(p)` the ulp of each coordinate of `p`. Rounding moves `p` by
`|r(p) - p| <= u(p) / 2` per coordinate, at most `rho(p) = |u(p)| / 2` in length.
For an output face `f` of area `A_f` and perimeter `P_f` whose corners move by at
most `rho_f`, every point of `f` moves by at most `rho_f` (a convex combination of
corner motions), so `f` sweeps a region inside `f + B(rho_f)`. By Steiner's
formula for a flat triangle,

```text
vol(f + B(rho)) = 2 A_f rho + (pi / 2) P_f rho^2 + (4 / 3) pi rho^3 =: T(f, rho).
```

Moving one closed surface into another changes the winding number of a point
only when a face sweeps across it, so the volume change is bounded:

```text
|vol(out) - vol(exact)| <= sum over moved faces of T(f, rho_f).          (V)
```

Each method below either fits (V) with a documented `rho_f` or does not.

**Separation lemma (S).** If every vertex lies farther than `rho_v + rho_f` from
every face it does not belong to, every edge farther than `rho_e + rho_e'` from
every edge it does not share a vertex with, and every corner of a face farther
than the sum of the motions from the opposite edge, then moving each vertex by at
most its `rho` keeps the surface embedded and keeps every face's orientation. A
crossing after the move needs a vertex-face or edge-edge pair to reach distance
zero, and their distance shrinks by at most the sum of their motions. A face can
flip only by passing through a degenerate state, which needs a corner to reach its
opposite edge. Input vertices have `rho = 0`, so only features near implicit
points can break. (S) is sufficient, not necessary; exact `validate` on the
rounded mesh decides the necessary and sufficient case.

### A. Iterated snap rounding on a grid derived from the input extent

Pick a grid of spacing `g` per axis, a power of two no finer than the ulp of the
largest coordinate in the bounding box, so every grid point is an `f64`. Every
vertex snaps to its grid point; every face that passes through a hot cell (a
cell holding a vertex) is rerouted through that cell's centre; iterated snap
rounding repeats until no face passes within half a cell of a vertex it does not
hold.

- **Correctness.** In 2D, snap rounding is a continuous deformation that only
  collapses features (Guibas and Marimont), and iterated snap rounding keeps
  vertices at least half a cell from foreign edges (Halperin and Packer). In 3D,
  rerouting a face through hot cells needs a retriangulation that does not cross
  itself; the only general construction with a proof we know is Devillers,
  Lazard and Lenhart, "3D snap rounding" (SoCG 2018). We have not checked their
  proof or their bounds line by line; as stated in their abstract, faces move by
  at most 3/2 cell in L-infinity and the output has `O(n^15)` faces in the worst
  case.
- **Termination.** Rerouting only adds vertices at existing hot-cell centres, so
  the number of (cell, face) insertions is finite: at most hot cells times faces.
- **Vertex movement.** Up to `3 g / 2` per coordinate for every vertex, inputs
  included. An input coordinate that is not a multiple of `g` moves: with a UTM
  northing near 7.5e6, `g >= 2^-30` m (about 9.3e-10 m), so a vertex at
  `z = 200.1` moves even though it is an exact input.
- **Volume.** (V) with `rho_f = (3 sqrt(3) / 2) g` on every face.
- **Coincident input vertices and shared walls.** Equal vertices snap together.
  Walls that differ by less than a cell become coincident faces, which then must
  cancel in opposite pairs. Every other input vertex not on the grid moves and must
  be reported.
- **Cost.** Every face is tested against every hot cell near it, and in 3D the
  rerouted faces must be retriangulated; the worst case above is out of reach for
  meshes of 1e5 faces on a 30 GB machine. Not measured.
- **Reduced cases resolved.** All of them, if the 3D construction is
  implemented in full. A simplified 3D rerouting has no proof.

### B. Round, then exactly re-arrange the faces that rounding made cross

Round every vertex. Find the faces whose rounded triangles intersect (exact
predicates on `f64`). Build the exact arrangement of those faces with the
existing kernel, keep the cells of winding number one or more, round the new
implicit points, and repeat until nothing crosses.

- **Correctness of one round.** The exact self-union of a closed oriented mesh is
  a valid solid, so each round's exact result is valid; only its new points are
  not `f64`.
- **Termination.** We cannot bound it. Old vertices are `f64` and do not move, so
  round `k + 1` can only break faces around the points created in round `k`, but
  no quantity decreases. Two faces that cross at an angle `theta` put their
  crossing line anywhere within about `ulp / theta` of where rounding puts it, so
  a near-coplanar pair (M2, M3) can produce a fresh crossing each round. The
  number of `f64` points in a bounding box is finite, so the loop ends in theory,
  but the count is near `2^150`. This repository already had a rounding cleanup
  that looped until memory ran out. With an iteration cap, B becomes a heuristic
  that still raises `KernelError`.
- **Vertex movement.** Each round moves new points by half an ulp; points can
  drift by the sum over rounds.
- **Volume.** (V) per round; the total grows with the number of rounds and has no
  bound without one on the rounds.
- **Coincident input vertices and shared walls.** Untouched; inputs never move.
- **Cost.** One exact arrangement per round on the crossing faces; the number of
  rounds is unbounded.
- **Reduced cases resolved.** Unknown for each case until run; no proof for any.

### C. Collapse sub-ulp features before materialization, by exact tests

Keep the exact result. Decide, with exact tests on the exact points and their
roundings, which features rounding cannot place, collapse them in a fixed order,
then round and certify with exact `validate`. Every coordinate the output uses is
`r(p)` of some exact result vertex `p`; no operation invents a coordinate. Call
`H` that set of hot points (input vertices are their own hot points).

Operations, each tested exactly:

- **C1 merge.** Vertices `p != q` with `r(p) = r(q)` become one (today's merge,
  extended across shells, never merging copies of one exact point). Motion: the
  rounding only.
- **C2 edge collapse.** On an edge `pq` whose rounded neighbourhood is degenerate
  or turned against its exact plane normal, with `|r(p) - r(q)| <= u` per
  coordinate, link condition holding and no moved face turning: `p` moves to
  `r(q)`. An input vertex is kept when one end is input. Motion: at most
  `3 u / 2` per coordinate.
- **C3 sub-ulp shell removal.** A rounded shell whose exact signed volume
  (exact sum on the `f64` corners) is zero or of the wrong sign is removed. Its
  exact volume is then at most `sum T(f, rho_f)` over its faces by (V), because
  its rounded volume does not exceed zero. Removing an inward shell fills a void
  by the same bound.
- **C4 cancellation.** Opposite faces on the same three vertices cancel
  (today's rule, across shells). Two faces with the same orientation on the same
  three vertices cannot come from a valid exact solid after C1 to C3: winding
  numbers on both sides would differ by two. If one appears, certification fails.
- **C5 hot-point insertion.** A hot point `h` whose rounding box
  `{x : r(x) = h}` meets an output face `g` that does not have `h` as a corner
  (exact box-triangle test with `orient3d` on `f64` corners) is inserted into `g`:
  split `g` at an exact point of `g` inside the box, then give that vertex the
  coordinates `h`. This is snap rounding's rerouting, with hot cells equal to the
  rounding boxes of the exact result's own vertices instead of a uniform grid.
  Motion: at most one box width, `u(h)` per coordinate.

- **Correctness.** The output is valid when exact `validate` reports a solid;
  every operation keeps each edge used once in each direction (C2 checks the link
  condition, C3 removes whole shells, C4 removes pairs, C5 splits a face and its
  edge neighbours conformally). Where (S) holds after C1 to C5 with the motions
  above, validity follows from the separation lemma without running `validate`.
  C has no general 3D theorem: two edges from different sheets can pass within
  an ulp of each other with no hot point near either, and C5 adds no point there.
- **Termination.** Let `F0` be the faces of the exact result and `H` the hot
  points; neither grows. Keep a set of (hot point, original face) pairs already
  inserted; C5 adds one pair and never removes one, so it runs at most
  `|H| * |F0|` times. C1 and C2 each remove a vertex; C3 and C4 each remove faces;
  C5 adds at most two vertices' worth of faces. Order the potential as
  `(|H| * |F0| - inserted pairs, vertices + faces)`: C5 lowers the first part,
  C1 to C4 keep it and lower the second. The total number of steps is at most
  `(|H| |F0| + 1) * (V0 + F0 + 2 |H| |F0|)`, a polynomial bound proven before the
  loop starts, and the loop asserts it.
- **Vertex movement.** Rounding: half an ulp per coordinate, not reported. Any
  vertex that ends farther than half an ulp from its exact position (C2, C5, or an
  input vertex merged into another) is counted. A vertex may take part in several
  operations, so the bound is the sum of its moves; the implementation refuses any
  operation that would push one vertex past `2 u` per coordinate from its exact
  position, and the report states that bound.
- **Volume.** (V) with `rho_f` the largest total move of a corner of `f`. Shell
  removal is inside (V) as shown under C3. Cancellation changes nothing.
- **Coincident input vertices and shared walls.** Equal input vertices are one
  exact point in the arrangement already; copies of one exact point that keep
  touching shells apart are never merged. Input vertices move only through C1 or
  C2 when two distinct input vertices round together or lie an ulp apart in a
  conflict; each such move is counted separately. Walls that agree only to
  rounding (M2, M3) become shared vertices through C1 and C5, and the thin pocket
  between them cancels through C4 or disappears through C3. A face contact
  between neighbouring solids that the inputs already had stays invalid, as the
  glossary requires.
- **Cost.** Results that round cleanly pay one `validate`, which runs today. The
  operations touch only faces near problems, found from the BVH `validate` already
  builds. Not measured.
- **Reduced cases resolved.** M1 by C3 (the ulp prisms). M4 by C1 or C2: the
  crossing point lies on an input edge within an ulp of its input end, so it rounds
  onto that end or collapses into it (tilted ground, icosphere). M2 and M3 by C1, C5
  and C4: the two sheets share hot points, and the faces between them cancel
  (clipped boxes, pit and pushback, both pushback pairs, rotated stopes, crosscut
  chain). Each of these is an argument from the measured mechanism; the tests
  decide.

## 4. Choice: C

C keeps every input vertex where it is unless two inputs collide within an ulp,
adds no coordinate the exact result did not round to, ends within a bound proven
before it starts, and costs nothing on results that round cleanly. A moves input
vertices that do not sit on its grid, and the only 3D construction with a proof
carries a worst case far beyond these meshes. B has no termination bound, and this
repository has already met that failure.

C does not reach the owner's target in theory: it has no 3D validity theorem, so
a case it cannot certify is not proven unrepresentable. The implementation keeps
`Error::Kernel` for exactly two reasons and names them in the message:

- `non-finite`: an exact coordinate rounds to infinity. No `f64` mesh represents
  it; this is the only refusal we can prove necessary.
- `uncertified`: the loop ended within its bound and `validate` still finds a
  problem. The error carries the report, the mesh and the rounding report.

Any `uncertified` case found in a test is a defect to reduce and analyze, not a
limit to document. If one survives analysis, the next step is A's 3D construction
applied to that case's conflict cluster only.

## 5. Changes

Functions, all in `crates/boolean/src`:

- `hp.rs`: `Pt::round_box(&self) -> ([f64; 3], [f64; 3])`, the `f64` box of points
  that round to `self.to_f64()`, for C5's exact box test.
- New `round.rs`:
  - `pub struct RoundingReport { rounded_vertices: usize, moved_vertices: usize,
    moved_input_vertices: usize, max_move: [f64; 3], move_bound: [f64; 3],
    merged_vertices: usize, collapsed_edges: usize, inserted_points: usize,
    removed_shells: usize, cancelled_faces: usize, volume_bound: f64, steps: usize,
    step_bound: usize }`. `max_move` and `move_bound` are per axis, in the mesh's
    length unit, rounded up; `volume_bound` is the sum (V), rounded up.
  - `pub(crate) fn materialize(arr: &Arrangement, shells: &[Shell]) ->
    Result<(Assembled, RoundingReport)>`: C1 to C5 on the whole output (all
    shells together), in the order C1, C4, C3, C2, C5, repeated to the bound, then
    exact `validate`.
  - `fn hot_boxes`, `fn insert_hot_point`, `fn remove_flat_shells`,
    `fn step_bound`.
- `ops.rs`: `assemble` calls `materialize`; `collapse_rounded`, `edge_collapse` and
  `rounding_problems` move into `round.rs` as C1, C2 and C4, and the
  validate-count loop with its undo goes away. `solid_result` raises
  `non-finite` or `uncertified`. `BooleanResult` gains `rounding: RoundingReport`.
- `clip.rs`: `ClipResult` gains `rounding: [RoundingReport; 2]`.
- `slice.rs`: `split_by_plane` keeps section points exact (`Pt`) through the cap
  triangulation and hands both parts to `materialize`; `coincident` and `clean` go
  away. A new `split_by_plane_report` returns the two meshes and their reports;
  `split_by_plane` keeps its signature.
- `crates/py`, `_wireframes.pyi`, `python/wireframes`: `wf.RoundingReport` with the
  same fields; booleans, `clip_by_surface` and `split_by_plane` take
  `rounding_report=False` and, when true, return the report after the meshes (and
  after the repair reports when `repair=True`).
- `spec.md`: tolerance rule 5 and section 6.2 item 7 are rewritten for C in the
  commit that changes the behaviour.

## 6. Tests that will prove it

- The 12 ignored tests in `crates/boolean/tests/rounding.rs`, the ignored test in
  `acceptance_open_pit.rs` and the 6 skipped Python acceptance tests run unchanged
  and pass.
- `cutting.rs`: `inclined_splits_through_every_sphere_vertex_cap_the_section`
  bounds refusals at 188 today; it must reach 0, and the bound tightens to 0.
- For every reduced case and every operation: each output vertex lies within
  `move_bound` of its exact source position, compared exactly with
  `num-rational` (already a dev-dependency); `moved_vertices` equals the count of
  vertices farther than half an ulp; and `|vol(out) - vol(exact)| <=
  volume_bound`, with `vol(exact)` summed exactly over the exact result's faces.
- One unit test per operation C1 to C5 on a hand-built shell, including a shell
  that needs C5 then C4.
- A termination test: a constructed worst case (many near-coplanar sheets at UTM
  coordinates) where `steps <= step_bound`, run under the 4 GB scope first.
- A proptest of seeded near-coplanar pairs (two boxes, one rotated or offset by a
  few ulps at UTM coordinates), moderate case count, asserting solids, the
  volume identity within `volume_bound`, and no `uncertified` error.
- Determinism: the reduced cases and the stage scenario give the same bytes with
  `RAYON_NUM_THREADS` 1, 3 and 4.
- Inputs unchanged: every boolean, clip and split test compares the input meshes
  before and after the call.

## 7. Implementation notes and departures

Measured on the reduced cases and on the 648 inclined splits of `sphere(2)`
(`cutting.rs`); counts, not timings.

- **When the operations run.** The per-shell cleanup of `spec.md` section 6.2 item
  7 runs first, as before. The operations over all shells run only when its
  output is not a solid, so a result that rounded cleanly is unchanged. When they
  end without a solid, the error carries the mesh of the per-shell cleanup run
  without the move cap below: the cleanup as it was before this method, which
  keeps the mesh closed, so `error.mesh.volume` still works. Its edge collapses
  can move a vertex far: about 9e11 ulps in the pushback of `docs/guide/booleans.md`.
  That mesh is never returned as a result, and its report raises `move_bound` to
  its `max_move`.
- **Volume bound.** `volume_bound` is `sum T(f_R, 2 rho_f) + |vol(F) - vol(R)|`,
  where `R` is the draft with every exact vertex at its rounding, `rho_f` the
  largest rounding radius of a corner of `f`, and `F` the output. The first term
  is (V) for the motion from the exact result to `R`: every face of the exact
  result lies within `rho_f` of its rounded face, so its sweep lies in
  `f_R + B(2 rho_f)`. The second term is computed exactly on `f64` points over the
  faces in which `R` and `F` differ, so it covers every later merge, flip,
  collapse, insertion and removal without a bound per operation.
- **C3** removes a component whose rounded volume is zero, or negative while
  `validate` reports it as an inward shell, and only when that volume is at most
  the component's own (V) sum.
- **C2 and C5 are driven by the problem pairs `validate` reports.** For each
  self-intersecting or overlapping pair, a corner `h` of one face whose hot box
  meets the other face triggers C2 when the box also holds a corner of that face
  (the two vertices become one), and C5 otherwise. C2 is a vertex identification,
  not only an edge collapse; the link condition is not checked, and the final
  `validate` decides.
- **Hot box.** The design used the rounding box of `h`. The implementation first
  tries one ulp of `h` per axis, rounded outward to `f64` ends so that the exact
  separating-axis test runs on `f64` values. A round where that box finds nothing
  retries with twice the sum of the rounding radii of `h` and of the face's
  farthest-rounding corner on every axis, the reach the separation lemma (S)
  allows (see "Isotropic reach" below). With the narrow box alone, 6 of the 648
  splits stayed uncertified: crossings that round two ulps from the vertex next to
  them. With the wide box alone, 18 did, because it acts on the first candidate it
  finds. Narrow then wide leaves none.
- **Move cap.** The design capped a vertex at two ulps per axis from its exact
  position. Collapsing those crossings onto their vertex needs two and a half, so
  the cap is four ulps of the vertex's largest coordinate on every axis
  (`move_bound`), checked for every collapse in both the per-shell cleanup and the
  operations over all shells.
- **Step bound.** `step_bound` adds the per-shell bound (each flip removes an edge
  the cleanup did not make; at most `3 f0 / 2 + v0 f0` such edges per round and
  `f0^2 + 1` rounds) to the bound over all shells (`v f` insertions, each adding at
  most two faces, plus one removal per vertex or face that ever exists, plus one
  round per action). Both are computed before the first step and asserted.
- **Splits** follow the design: section points stay exact through the cap CDT, and
  both parts go through `materialize`. `coincident` remains, but on exact points it
  only joins section points at the same exact position; `clean` is gone.
- **Isotropic reach.** The wide hot box and the move cap are the same on every
  axis: the rounding radius of a point is the length of its half-ulp vector, set by
  its largest coordinate. At UTM northings an x or y ulp is about 9.3e-10 m and a z
  ulp at 300 m about 5.7e-14 m; the runaway pair holds vertices with equal x and y
  whose z values spread over about 3.5e-11 m, which a per-axis box never joins
  while a face rounding by half an x ulp can still tilt across them.
- **Flat faces and copies.** The loop also takes `validate`'s degenerate faces: the
  middle corner of a flat face is inserted on that face's longest edge (C5 on the
  face across it), which is the cleanup's flip. No insertion or collapse may put two
  copies of one exact point on one face. Two opposite faces on the same three
  coordinates whose corners differ only by such copies have the copies joined, so
  they cancel (C4). No insertion may repeat a live face with the same
  orientation. The loop validates a mesh of the vertices its live faces use.
- **Flat runs.** The property test of chains of twelve unions at UTM coordinates
  found lines of rounded points: a crosscut wall crossing a drive wall that
  earlier unions left planar only to rounding meets it along a vertical line in
  exact arithmetic only to within sub-ulp kinks, so its crossings round to points
  whose x differ by up to two ulps either side. Faces of the crosscut wall across
  three such points are flat or nearly so, and inserting the middle corner of one
  flat face into another flat face makes no progress. A run of such faces has no
  area to rounding, so it is removed whole and the faces around it are split at
  its vertices along the line (ordered by the axis the run spreads most on), so
  the two sides meet edge for edge. A run is taken only when a face of it is
  collinear or a problem face, and only when every vertex lies within the wide
  reach of the segment between the run's extreme vertices; a run along a curve,
  such as thin faces around an icosphere's cap, would otherwise be split in the
  wrong order. It is put back when the edges along the line do not pair, or
  when the faces it makes hold as many flat faces as the run (section 9). The
  loop also counts removals and takes none past `f + 4 v f` (`zip_bound`),
  which is part of `step_bound`.

## 8. Conforming inputs that agree to rounding

C2 and C5 work from hot points, so they cannot help where surfaces of the two
inputs agree only to an ulp with different triangulations (M2, M3): their edges
pass within an ulp of each other with no vertex near, the gap section 3.C names.
The reduced cases for clipped boxes, the pit and pushback, the rotated stopes and
the runaway pair all stayed uncertified after section 7, with nearly every problem
pair a face of `a` and a face of `b`, parallel and opposite.

The fix moves the snapping before the arrangement. When a boolean or a clip is
uncertified, `conform.rs` makes copies of both inputs: every vertex of one input
whose narrow hot box meets a face of the other, not at a corner, is inserted into
that face; then every face of one input whose three corners are now vertices of the
other has its three edges recovered in the other by flips (Sloan 1993, in the
projection along the face's dominant normal axis, only through strictly convex
quads within four ulps of the face's plane). Both triangulations then hold the same
faces wherever the surfaces agree to rounding, and the exact arrangement sees those
faces as coincident and classifies them exactly, so no sliver is left to round.
Every insertion and flip is kept only when its new faces are not flat and cross or
overlap no face near them (exact pair tests, found through a BVH of the input faces).

- **No input vertex moves.** Insertions add vertices of the other input at their own
  coordinates.
- **Volume.** An interior insertion replaces a face by three across the tetrahedron
  of its corners and the point; an edge insertion does so on both faces of the edge;
  a flip replaces two faces across the tetrahedron of the quad. The symmetric
  difference between an input and its copy lies in the union of those tetrahedra,
  so `vol(a xor a') + vol(b xor b')` is at most the sum of their exact volumes,
  which bounds `|vol(op(a', b')) - vol(op(a, b))|` and is added to `volume_bound`.
  The rational check `public_booleans_of_reduced_cases_stay_within_their_bounds`
  compares the public result with the exact result of the original inputs.
- **Termination.** Each vertex enters the other input at most once, each adding two
  live faces, so live faces stay below `F = f + 2 v`; flips keep that count. There
  are at most `6 F` edge recoveries, each crossing at most `3 F` edges and allowed
  `c^2 + 4` attempts for its `c` crossing edges. The bound is computed before the
  first step and asserted, and its steps join the report's.
- **Result.** The result of the conformed inputs is returned only when certified;
  otherwise the first, uncertified result is the error.

## 9. Adversarial cases

An adversarial tester built mining scenarios against sections 7 and 8; each
failure is reduced to the two inputs of the operation that failed, in
`crates/boolean/tests/data/attacks/`, and runs in
`crates/boolean/tests/attacks.rs` with the volume identity checked in rationals
against the sum of the volume bounds. Two root causes and one missing operation
explain all eight.

### Flat runs that come back

The union of two pit shells clipped by one topography never finished. The loop
removed a one-face flat run, split the face across it into two faces, one of
which was flat by the hot-box test, and on the next round removed that face and
recreated the first: two states, alternating, with the live face count fixed.
The `zip_bound` cap ended the loop after about `4 v f` rounds, a validate each,
so the bound proved termination and no useful time.

A removal is now kept only when the faces it makes hold fewer flat faces than
the run it removed. Faces next to a run are not flat (flat faces joined by an
edge belong to one run), so the number of flat live faces drops by at least one
with each kept removal. Between two other operations the loop can therefore
take at most as many removals as there are live faces, and a removal can never
undo the one before it. With this change the pits union returns a solid; the
three chains (lids, caps, the turned body) finish too, then needed the descent
operations themselves, not of the split.

### Sub-ulp parts with crossing faces

`split_by_plane` along a rotated stope's own face plane (origin a face vertex,
normal the `f64` cross product of the face's edges) puts the plane within about
1e-15 m of the face, a fraction of the 9.3e-10 m ulp at UTM northings, at an
angle that lets it cross the nearly parallel edges anywhere along their length.
One part is a wedge thinner than an ulp whose crossing points round off its
faces by half an ulp, so its cap and its wall cross. C3 removed a component only
when its rounded volume was zero or wrong-signed; this one has a small positive
volume and two self-intersections. C3 now also takes a component that holds a
face of a self-intersecting, overlapping or flat problem, under the same test:
its rounded volume is at most its own (V) sum, so its exact volume is at most
twice that sum, and `volume_bound` counts the removal through the exact volume
change. The upper part of the split comes back empty, which the test accepts.

A second round of attacks split rotated grid stopes the same way and found 9 of
218 face splits that failed with one non-manifold edge instead
(`split_along_own_face_a` to `_i` in `attacks.rs`). The plane passes through a
vertex of the wall with the wall's `f64` normal, so every other wall vertex lies
within a fraction of an ulp of it, on either side. One part is the sliver between
the plane and the wall: the wall's faces and cap faces a fraction of an ulp away.
Where a wall edge meets the cap, rounding leaves the edge with two wall faces and
two cap faces in an order that does not alternate, which `validate` reports as a
non-manifold edge, not as a crossing. C3 now also takes a component that holds a
face of a non-manifold edge or vertex, under the same volume test, and the sliver
part comes back empty.
The wider test runs only where descents run (below). Run on every boolean and
clip, it certified the first run of a pit clipped twice by one topography with
53 sliver shells above it, where the conformed inputs give the empty part the
Python test `test_rounding_report_comes_last_and_its_bounds_hold` expects.

### Descent

The other five (repeated union then difference, the stope block minus a raise,
development turned to strike, the lid and cap chains) stalled with every
operation of section 7 refusing or making no progress: insertions created
non-manifold edges and further crossings, and nothing undid them. Their inputs
are themselves rounded results, holding fins and walls a few ulps thick that
earlier operations left (M3 compounded over a chain).

The descent replaces the heuristic order with a measured one. When neither C2
nor C5 finds anything, the loop collects candidate changes from the problems
`validate` reports: for each self-intersecting, overlapping or flat face, the
flip of each of its edges with the face across; for each pair, the join of
opposite copies (C4) and, for each corner `h` of one face whose wide hot box
meets the other, the collapse of `h` into each corner within the box (C2) and
its insertion into the face, inside or on each edge within the box (C5). Each
candidate is applied to a copy, opposite pairs cancel, and the loop counts the
problems `validate` finds among the faces whose bounding boxes meet the box of
every face the change made or removed. The change that lowers that count the
most is kept, the first in candidate order on ties; when none lowers it, the
loop ends.

- **The local count is the global count.** A face the change made or removed
  lies inside the box, so every pair it is in, and every edge it has, is among
  the faces counted; faces outside the box and their pairs are the same before
  and after. Each kept descent therefore lowers the problem count of the whole
  mesh by at least one, and no sequence of descents repeats a state.
- **Only descents after the first.** Once the loop descends it takes nothing
  else, so the insertions that the descent undoes cannot be made again; with
  them interleaved the problem count of a reduced case alternated between two
  values.
- **Termination.** Each descent lowers a non-negative count, and the loop also
  counts descents and takes none past `4 (v + f)` (`descent_bound`); each
  descent adds its change, its round and the removal of the at most two faces
  it adds, `4 * descent_bound`, to `step_bound`. A descent tries at most
  49 candidates per problem pair, 3 per flat face, and for a non-manifold edge
  of `k` faces `3 k` flips and 43 candidates for each of its `k (k - 1) / 2`
  pairs, each with one copy of the mesh and one `validate` of the faces in its
  box. Time: not measured.
- **Bounds.** Flips and insertions move no vertex. Collapses go through C2 and
  keep the move cap of section 7, so `move_bound` stands. Every change replaces
  faces, and `volume_bound` already adds the exact volume change from the draft
  to the output over the faces that differ.
- **When it runs.** Descents cost a copy and a local validate per candidate, so
  a boolean or clip first runs as before: the inputs, then the conformed inputs
  (section 8). Only when neither is certified do both run again with descents.
  A split has no conformed retry and descends on its first run. The descent
  leaves results that sections 7 and 8 certify as they were.

### Faces around a non-manifold edge

The second round of attacks chained union then difference of a body turned by
under a degree at UTM coordinates (`repeat_fifth_union`,
`repeat_turned_third_difference`). Their first inputs touch themselves along 2
to 6 edges, so the exact result holds several copies of one exact point. After
rounding, two faces of one input landed on the same coordinates with opposite
orientation, differing only by such copies, and `validate` reported the edge
they share as non-manifold: four faces whose order does not alternate. C4 joined
copies only for crossing or overlapping pairs, and the descent took candidates
only from those pairs and from flat faces, so nothing acted on the edge. Now C4
also joins copies among the faces around each non-manifold edge, and the
descent treats every two faces around such an edge as a pair, with the same
joins, collapses, insertions and flips. Both chains then return solids within
their bounds.

### Thin slots inside one input

Stopes rotated to strike 118.4 and dip 62.5, united one at a time, then minus a
drive along strike (`stopes_dipping_62_minus_drive`) failed with 231
self-intersections; the descent left about ten. Of the crossing pairs, 106 joined
two faces of the stope block, 86 a face of the block and one of the drive, and
44 two faces of the drive. No face of the drive is within 1e-6 rad of parallel to
a face of the block, so the conforming of section 8 has nothing to match. The
block holds 321 pairs of opposite faces closer than an ulp: where two stopes
agreed only to rounding when they were united, the union left slots between their
walls a fraction of an ulp wide. Rounding the drive's crossing points on both
walls of a slot swaps the walls, and no single local change undoes it: moving a
merged crossing vertex by one ulp put two walls back in order on one edge of the
stope chain case and made eleven new problems nearby.

`conform.rs` now closes such slots before the arrangement, as a last retry after
sections 7, 8 and the descent. Faces of one input that face opposite ways, lie
within four ulps of each other's planes and overlap in the projection along the
first face's dominant axis are paired, and the pairs are two-coloured (a slot
has two walls); an odd cycle drops the retry. One sheet holds the whole input with
a colour per face, and edges and clashes are kept per colour, so the two walls of
a slot can share vertices and edges. Each vertex of one colour is inserted, at
most once, into a face of the other colour that its wide hot box meets away from
the face's corners; then each face whose corners are all vertices of the other
colour has its edges recovered on the other colour by section 8's flips. Faces
that then hold the same three vertices in opposite orientation cancel, which
closes the slot where its walls agree, and an inward shell that seals off is
removed when it holds at most the volume of a sheet `4 rho` thick over its faces.
The pass repeats on its own output while it cancels faces, at most once per face
of the input. The closed copy is kept only when exact `validate` passes it as a
solid, and the boolean then runs on it with section 8's conforming and with
descents.

- **No vertex moves.** Inserted points are vertices of the same input.
- **Volume.** Each insertion and flip replaces faces across tetrahedra; their
  exact volumes, summed, bound `vol(a xor a')`, as in section 8, and join
  `volume_bound`. Cancelling two faces on the same vertices changes no volume.
  Removing a sealed-off void adds its exact volume, rounded up. The unit test
  `closing_the_slots_of_rotated_stopes_ends_within_its_bounds` checks the bound
  on the closed input against rationals.
- **Termination.** One pass makes at most `2 v` insertions (each vertex once per
  colour) and the recoveries of section 8 on at most `f + 4 v` faces; its bound,
  `step_bound(2 v, f)`, is computed before the first step and asserted. Passes
  stop when one cancels nothing, and a counter allows at most one per input face.
- **Reach.** Pairing and insertion use the wide hot box, about 3.5 ulps of the
  largest coordinate at UTM northings. Slots wider than that stay open.

With this, `stopes_dipping_62_minus_drive` returns a solid within its bounds.

Determinism: candidate order follows `validate`'s problem order and fixed loops
over corners and edges, so `attack_results_do_not_depend_on_thread_count`
compares four attack cases on one and four threads, mesh, face sources and
report.

The descent has no more of a 3D theorem than the rest of method C: it ends when
no single local change lowers the count, which can leave problems a sequence of
two changes would remove. Such a case still raises `uncertified`; none in the
tests does.

## 10. Rounding 2: the two refusals left after section 9

Two attack cases still raise `uncertified`:
`stope_chain_minus_a_drive_along_strike_is_a_solid` and
`stopes_dipping_62_minus_a_footwall_drive_is_a_solid` in `attacks.rs`. This
section traces both, reduces them, and decides the next change. It changes no
kernel code.

### Reductions

Cropping by axis-aligned planes through both inputs, as in section 1, changed the
mechanism: the caps of `a` and `b` on one plane agree only to rounding (M2), and
conforming (section 8) then acted where the original never did. The reduction
that keeps the mechanism crops `a` with six planes (along the drive's horizontal
axis, across it, and in `z`) about the slot that fails, crops the drive only
along its axis, 0.5 m beyond the box, and accepts a crop only when the
difference still raises `uncertified` and its trace has the same signature as
the original: conforming does nothing, closing ends the same way, and every
problem pair of the last run joins two faces of `a`. Bisection ran on a cube
half-width, then on each of the six sides.

| Case | Original (faces of a + b) | Reduced (faces of a + b) | Test, ignored as `rounding 2: ...` |
|---|---|---|---|
| Stope chain minus drive | 934 + 64 | 232 + 166 (box 3.7 x 9.9 x 6.6 m, drive cropped) | `stope_chain_minus_a_drive_reduced_is_a_solid` |
| Stopes dipping 62.5 minus footwall drive | 590 + 64 | 282 + 64 (cube 54.5 m, whole drive) | `stopes_dipping_62_minus_a_footwall_drive_reduced_is_a_solid` |

The cropped drive of the chain has more faces than the whole one, because each
cut splits its long side faces; with the whole drive the cropped chain
certifies, so the crop of the drive stays. The footwall case stops at a cube of
half-width 27.25 m: at 27.19 m it certifies, and at 24 m closing fails
differently. In a release build each of the four inputs gives an `uncertified`
error, no other error and no panic; the check stays out of the default test run
because it is slow in the debug build the gate uses (not measured as a benchmark).

### Where each one refuses

Both run every stage of `ops::boolean`: the plain run, conforming (section 8,
which finds nothing to match in either case), the run with descents, and the
closed inputs (section 9) with and without descents. Every problem left at the
end of each run joins two faces of `a` that face opposite ways a fraction of an
ulp apart: two walls of a slot in the stope block, cut by the drive. The two
walls hold distinct exact crossing points with each drive face, closer than an
ulp; each rounds to its own nearest `f64`, and when the two roundings land in
the wrong order the walls cross.

**Footwall drive.** The plain run ends with 153 self-intersections inside `a`,
13 non-manifold edges and 8 non-manifold vertices after C1 to C5; the run with
descents gets down to 6 self-intersections, then no single change lowers the
count. Closing pairs 55 faces in its first pass, inserts 20 points, flips 23
edges and cancels 74 faces; its second pass pairs 10 faces and changes nothing.
The closed run still ends with 3 self-intersections between faces 550 and 137
of the closed `a` (anti-parallel, 0 to 0.48 ulp apart), 1 non-manifold edge and
4 non-manifold vertices. The slot left open is three faces: 550 = (76, 48, 45)
on one side, 137 = (45, 46, 76) and 551 = (76, 46, 48) on the other, with
vertex 46 within reach of edge 48-45. Inserting 46 on that edge would split 550
into two faces that cancel 137 and 551 exactly. `Sheet::commit` refuses it:
one face it makes from the face across the edge, (45, 46, 27), overlaps 137,
and the commit tests that overlap before 137 cancels against the made face
(46, 45, 76). A prototype that cancels such pairs inside the commit closes
more (80 faces cancelled), but 8 faces stay paired. In the reduced case the
reason shows: the face across the slot's rim, (42, 45, 24), is itself thinner
than rounding (its corner 42 lies within reach of its edge 24-45), so every
single insertion leaves a flat or crossing face that the commit refuses, and no
order of single insertions reaches the state where the slot and the sliver
vanish together.

**Stope chain.** The plain run ends with 20 self-intersections inside `a`,
8 non-manifold edges and 8 non-manifold vertices; descents leave 5
self-intersections between faces 463 and 511 of `a` (anti-parallel, 0 to
0.81 ulp apart), 2 non-manifold edges and 8 non-manifold vertices. Closing
pairs 34 faces, inserts 8 points, flips 6 edges and cancels 24 faces, then
`validate` rejects the closed copy (2 non-manifold edges, 2
self-intersections), so `close_and_conform` returns nothing and the boolean
never runs on a closed input. The crossing joins a face the closing made on one
wall (colour true) with two faces next to the slot that no pair coloured,
(294, 295, 293) and (294, 253, 295); vertices 293 and 295 differ only in `z`,
by about 3e-10 m, under the 9.3e-10 m ulp of the northing. `Sheet::commit`
skips the overlap test between faces of different sides, and faces outside any
slot sit on side false, so nothing tested the made face against them. With the
test kept for every pair except two walls of opposite sides (prototype, not
committed), closing runs three passes, the closed copy is a solid, and the
difference certifies. The reduced chain still fails with that prototype: its
closed copy has 4 non-manifold edges, 5 non-manifold vertices, 4
self-intersections and 15 overlaps.

### Why C1 to C5 cannot repair it on the output

Let `p1` and `p2` be the crossings of one drive face with the two walls, closer
than an ulp, rounding to different `f64` points in swapped order. C1 needs
`r(p1) = r(p2)`. C2 joins `p1` to `p2` when one hot box holds the other, but the
two wall faces then share one vertex and nothing along the cut, which is the
rise in non-manifold edges and vertices the descents show. C5 puts `p1` into
the other wall, on the wrong side of it. C4 needs both walls to hold the same
three vertices, which takes every crossing along the cut joined at once, a
sequence of C2 steps each of which raises the local problem count, so the
descent stops before it. C3 does not apply: the slot belongs to the large shell.
Removing the slot needs one change over a cluster of faces, before the
arrangement or after it.

### Decision: snap the failing cluster as one change

Two causes, one route.

1. **Fix inside closing (chain).** `Sheet::commit` tests every made face against
   every live face near it, except a made face of one wall against a face of the
   opposite wall of the same slot; before the test it cancels each made face
   against a live face on the same three vertices in opposite orientation, and
   the pass counts those cancellations with the rest.
2. **Escalation: cluster snap (footwall).** When a closing pass leaves paired
   faces it cannot close, and when the output loop of section 7 ends uncertified
   after descents, the loop takes a cluster `K`: the faces of the problems (the
   paired faces left, or the problem faces `validate` reports), grown by every
   face whose box meets the wide hot box of a vertex of `K`, until no face joins.
   Let `H` be the vertices of `K`. For every face `f` of `K` and every `h` in `H`
   that is not a corner of `f` and whose wide hot box meets `f` (exact test of
   section 7), `h` is placed on `f`: on the edge the box meets, for both faces on
   that edge, else inside. Each face of `K` is then triangulated at once over its
   corners and placed points (points on an edge ordered along it, the
   triangulation built by exact `orient2d` in the projection along the face's
   dominant axis), faces whose corners are collinear go as flat runs (section
   7), opposite pairs cancel (C4), and components of zero exact volume go (C3).
   One exact `validate` of the faces whose boxes meet the cluster's box accepts
   the whole change or rejects it. This is snap rounding with the rounding cells
   of the result's own vertices as hot cells, applied to one cluster: the
   rerouting of section 3.A without a grid, and without its proof.

Both changes run only where today's pipeline ends `uncertified`: the corrected
closing as a new retry after the current one, the cluster snap after that and
as the last step of the output loop. A result certified today runs the same
code and gives the same bytes; its timing is not measured here.

**Guarantee, plainly.** This route gives no validity guarantee for all inputs,
and none proven for a named class. Exact `validate` certifies each result, as
for the rest of method C, and `uncertified` stays possible. The route aims at
one class: inputs holding slots and slivers thinner than rounding (left by
earlier unions of bodies that agreed only to rounding) crossed by the other
input. The traces above argue that it closes the two refusals; that is
evidence, and the tests will decide. A guarantee for all inputs would need the
3D snap rounding of section 3.A (Devillers, Lazard and Lenhart) run on the
cluster, with its worst-case output size and its proof checked; this route does
not implement it.

**Correctness of what is proven.** Placing `h` on an edge splits both faces on
that edge, and each face of `K` is triangulated as the same polygon its
neighbours see (the same points on each shared edge, in the same order), so
every edge of the cluster stays used once in each direction; cancellation
removes pairs and C3 removes whole components, so the mesh stays closed. The
commit or `validate` decides the rest.

**Termination.** Cluster growth adds faces and stops at the live face count `F`.
One snap places at most `|H| |K| <= V F` points and triangulates each face once;
a face holding `m` placed points becomes at most `2 m + 1` faces, so one snap
makes at most `F + 2 V F` faces. Each cluster is seeded by problem faces, and a
face seeds at most one cluster in a run, accepted or rejected, so a run takes at
most `F` snaps. Closing passes stay bounded as in section 9 (one per input face)
and each adds at most `F` snaps. The bound `F (V F + F + 2 V F)` per run joins
`step_bound` before the first step and is asserted, as every bound in this
file.

**Movement.** Placed points are vertices that already exist: in closing,
vertices of the same input; on the output, rounded vertices of the result. No
vertex moves, so `moved_vertices`, `max_move` and `move_bound` keep their
meaning. A face moves: every point of a face of `K` ends within the wide reach
(twice the sum of the rounding radii of its corner and the placed point) of
where it was, because each new face has corners on the old face or within that
reach of it.

**Volume.** Placements are applied in a fixed order (faces, then points by
index), each replacing faces across the tetrahedron of the face and the point,
as `Sheet::place` does; the exact volumes of those tetrahedra, summed, bound
`vol(a xor a')` in closing and join `volume_bound` as in sections 8 and 9. On
the output, `volume_bound` already adds the exact volume change between the
draft and the output over the faces that differ.

**Functions to change.**

- `conform.rs`: `Sheet::commit` (overlap test across sides except between two
  walls, cancellation of opposite pairs before the test, a count of them);
  `Sheet` gains the wall flag per input face; `close_slots` counts
  cancellations made in commits and calls a new `snap_cluster` when a pass
  leaves paired faces; new `fn snap_cluster(&mut Sheet, seeds: &[usize]) ->
  Option<f64>` (the volume of its tetrahedra) and `fn retriangulate(&mut Sheet,
  f: usize, points: &[(u32, Option<usize>)])`.
- `round.rs`: `Global::snap_cluster` with the same steps over the rounded
  result, called in `run` after the descents end without a solid;
  `global_step_bound` adds the bound above.
- `ops.rs` (`boolean`) and `clip.rs`: a third retry with the corrected closing
  after the current stages.
- `RoundingReport` gains `snapped_clusters: usize`, with its `pyo3` binding,
  `.pyi` stub and `spec.md` tolerance rule 5, in the commit that changes the
  behaviour.

**Tests that will prove it.** The four tests of this section run unchanged and
pass, with the volume identity checked in rationals; a unit test snaps the
footwall residue (the three slot faces, the sliver on the rim and their
neighbours) by hand; a termination test asserts `steps <= step_bound` on a
cluster built to grow to the whole mesh, run under the 4 GB scope first; the
determinism test of section 9 takes the two cases; every attack, rounding and
acceptance test passes as before.

### Step 1 landed: the strict commit

`Sheet` carries a `strict` flag and a wall flag per input face; only a third
closing retry, after the closing of section 9 and its runs, sets them, so a
result certified before runs the same code (`results_certified_before_rounding_2_are_unchanged`
pins five of them by a hash of their bytes). With the flag, a commit tests a made
face against every live face near it except a face of the opposite colour when
both are walls, and leaves out of the test each made face that has a live face
on the same three vertices in opposite orientation, together with that face;
the pair cancels at the end of the pass as before. This departs from the
decision above in one point: the commit does not cancel such pairs itself, it
only exempts them, so the undo log and the edge maps stay as they were and the
cancellation count is the one the pass already makes. Termination and bounds
are those of section 9: the commit decides which steps are kept, not how many
are tried.

`stope_chain_minus_a_drive_along_strike_is_a_solid` now passes unchanged,
difference and intersection, with the volume identity in rationals. The
reduced chain and both footwall cases still raise `uncertified`.

### Step 2 landed: absorbing slivers and coplanar walls; no cluster snap

The escalation the decision names, the cluster snap, turned out not to be
needed for these inputs, and it is not implemented. Two smaller changes to the
strict closing of step 1 certify all three cases left, and this section
replaces the cluster snap part of the decision above.

1. **Absorbing the face across (footwall).** With `strict`, `Sheet::place` of a
   point `h` on edge `a b` of face `f = (a, b, x)` checks whether `h` is the far
   corner `d` of the face across, `g = (b, a, d)`. Before, that made two flat
   faces, which the commit refused, so the rim sliver of the reduced footwall
   slot could never go. Now `f` and `g` are both removed and `(a, h, x)`,
   `(h, b, x)` take their place. The edges `a h` and `h b` of `g` stay used in
   the same direction by the new faces, and `h x` is used once each way, so
   the shell stays closed (`a_point_on_the_far_corner_absorbs_the_face_across`).
   The faces before and after bound the tetrahedron `(a, b, x, h)`, whose exact
   volume, rounded up, is returned and joins `volume_bound`. No vertex moves.
2. **Coplanar walls (reduced chain).** The reduced chain pairs two walls that
   lie exactly in one plane (their vertices sit on two vertical lines) and do
   not overlap in the input. The strict commit skips the test between opposite
   walls, so insertions and flips made them overlap in their plane, and they
   could not cancel because their triangulations differ around two vertices
   0.16 m apart. A second strict variant (`Sheet::coplanar`) also refuses a
   made face that overlaps a face of the opposite wall exactly in a common
   plane without holding the same three corners. It runs only when the plain
   strict pass gives nothing: run first, it leaves the full footwall case with
   one non-manifold edge, because there the walls must overlap in their plane
   for a moment before they cancel.

`close_all` with `strict` therefore tries, per pass, the strict pass and then
the coplanar one. Each pass keeps the step bound of section 9, asserted, and the
number of passes stays at most the face count of the input, so the closing does
at most twice the work of section 9 per input face. No new loop was added.

**Results.** `stope_chain_minus_a_drive_reduced_is_a_solid`,
`stopes_dipping_62_minus_a_footwall_drive_is_a_solid` and
`stopes_dipping_62_minus_a_footwall_drive_reduced_is_a_solid` pass unchanged,
with the volume identity in rationals. The unit test
`strict_closing_of_the_reduced_cases_ends_within_its_bounds` checks the step
bound and, in rationals, the volume bound of the closed inputs. The thread-count
test takes both reduced cases. `RoundingReport` gains no field: absorbed faces
count among `cancelled_faces` through the pass, and its insertions among
`inserted_points`.

**Guarantee, unchanged.** This gives no validity guarantee for all inputs and
none proven for a named class. Exact `validate` certifies each result, and
`uncertified` stays possible. The cluster snap stays the documented escalation
for a slot that neither strict variant closes.

### The class as a property test

`stope_chains_minus_a_drive_along_strike` in `rounding_props.rs` builds the
class these cases stand for: eight stopes, four along strike in two levels,
each upper stope standing on the top of the one below, all with their footwall
on one common wall, each with its own height (18 to 30 m) and thickness (4 to
8 m), turned to a dip from 45 to 90 degrees (90 drawn on its own) and to a
strike, and moved to UTM coordinates, every vertex rounded to `f64` at each
step. The stopes are united one at a time, every union checked to be a solid,
and a drive along strike is taken away, on the floor of the lower level or
through the common wall. The difference and the intersection must be solids
within their step and move bounds, the difference the same on 1 and 4 threads,
and `vol(a - b) + vol(a & b) = vol(a)` must hold in rationals within the summed
volume bounds. It runs 6 cases in the gate.

### Sliver shells nested in the block

With 300 cases on two other seeds in a release build (not part of the gate),
that test found a union that the stages above left with one problem: an inward
shell. Eight vertical stopes in two levels against
a common wall, turned to strike at UTM coordinates, united one at a time; the
union with the seventh stope (`seventh_union_of_a_stacked_stope_chain_is_a_solid`)
rounded to a block and an eight-face shell about 5 m long whose volume is about
2e-10 m3. The shell lies inside the block and touches it along one edge and at
four vertices, so `validate` reports it as wound against its nesting parity.
C3 missed it twice: it groups faces by edge connectivity, which joins the shell
to the block through the edge they share, and it took inward components only
when their volume is negative, while this one is positive.

With `crossing` (the runs with descents), when no component goes, C3 now also
takes each shell, as `validate` pairs edges into wedges, that holds the face of
an `InwardShell` problem, of either sign, under the same test: its rounded
volume at most its own (V) sum, so its exact volume is at most twice that sum.
`volume_bound` counts the removal through the exact volume change from the
draft to the output, as for every C3 removal. The rule runs only in runs with
descents, only when no component went and `validate` reports an inward shell;
the five results that `results_certified_before_rounding_2_are_unchanged` pins
keep their bytes. It adds no loop: it is one more test inside C3, whose
removals `step_bound` already counts.

### A refusal the property test found

Run with 300 cases on two other seeds in a release build (not part of the
gate), it found two refusals. The first, a union of the chain, is fixed above.
The second, on seed 11 after 90 passing cases, still raises `uncertified` and runs as the ignored test
`stacked_stope_chain_minus_a_floor_drive_is_a_solid` (inputs of 236 and 12
faces, dip 89.6, the drive on the floor). It is not a slot: closing finds no
paired walls in either input and conforming finds nothing to match. The
descent ends after five descents with one edge of three faces, five
non-manifold vertices and three crossings, all between faces of the stope block
where the drive crosses a near-vertical wall along a line of rounded points a
few millimetres apart in `x` and `y`. No single change lowers the count from
there. It is not reduced and not traced further in this pass.

## 11. Attack 2: splits through a block's own vertices

The mining scenarios of `attacks_2.rs` split twelve- and 120-stope blocks,
turned to strike and dip at UTM coordinates, by planes through their own
vertices. Seven splits failed: two panicked and five raised `uncertified`.
Two defects caused four of them; neither is in the rounding operations.

**The panics.** `radial_pairs_at` required every face around an edge to hold
one corner at each end and a third off their line, and panicked otherwise
(`third corner off the edge`, or an index of 3 when both ends of the edge share
coordinates). `validate` keeps that precondition by building its edge map
without degenerate faces, but C3 builds `Shells` over the live mesh of the
rounding loop, where flat faces and faces with two corners on one point are
ordinary states. Such a face bounds no wedge, so the edge is now
`Radial::Unbalanced` and `Shells` reports it as non-manifold, which is what it
is (`faces_that_bound_no_wedge_leave_an_edge_unpaired`).

**The open draft.** With the panic gone, both cases and two of the level splits
raised `uncertified` with boundary edges. Instrumenting `split_by_plane_report`
showed that each part's faces, counted by section point, were closed, while the
draft handed to `materialize` had 4 to 12 edges used unevenly. Section points at
one exact position (shells of the block touching on the plane) get one copy
per output shell, and `shells_of` labelled shells on the rounded section points
with degenerate faces skipped. A face that is a triangle exactly but whose
corners round onto one line (in the level split a, a wall face whose two
section points differ by 9e-11 m in `z`) got no label, hence a copy of its own,
and came apart from its neighbours: an open draft no rounding can certify.
`shells_of` now keeps every face. An edge with two faces joins them whatever
their shape, and an edge with more faces where one bounds no wedge joins them
all, which keeps fewer copies and never opens a seam. Level splits a and d, the
inclined split and the 120-stope split now pass unchanged.

The level splits b and c and the split through a slot reached the rounding
closed and still raised `uncertified`: they are failures of the rounding
operations themselves, not of the split.

**The split's draft against the boolean's.** Intersecting the block of level
split b or c with a box whose top lies on the level certifies in stage 0 of
`ops::boolean` (the plain run or the run with descents), with no conforming and
no slot closing. So the gap lies in the draft the split builds. Two things
differ: the boolean gives each exact point one copy per umbrella of faces about
it and hands `materialize` one group per output shell, while the split gave a
copy only where section points coincide and put all faces in one group. A part
the first draft does not certify is now drafted the boolean's way and rounded
again; the second result is kept only when it certifies, so every part
certified before keeps its bytes. That certifies the split through a slot (each
of its walls is cleaned up with its own shell). It does not certify b or c.

**Level splits b and c.** In c the rounding loop ends with no crossings, one
edge of four faces and six non-manifold vertices around it, and no descent
lowers the count. The edge joins two rounded section points on the level, two
cap faces and two faces above: a wall rising 25 m, and the upper part of an
input face whose third corner lies 2.66e-11 m above the level and, in plan,
almost on the edge's line. The input block holds a step 1.9e-11 m high there
(two input vertices with equal `x` and `y`), left by the unions that built it.
Rounding moves the two section points by about 1e-10 m in plan, more than the
height of that face, so its half-plane about the edge swings past the wall and
the faces no longer alternate in direction around the edge. Every operation of
method C places the corner of that face at an existing point, and none of them
moves it below the wall without inventing a coordinate. This is the gap section
3.C names (a feature thinner than rounding with no hot point that repairs it),
met inside one input. Not fixed by any operation of method C.

**Coordinate planes.** The boolean of the block with a box whose top lies on the
level certifies b and c, as it did for the 120 tries of the attack sweep. Its
exact result is the same set as the split's part, so a part that neither split
draft certifies is now, for a coordinate plane, computed that way and kept when
the boolean certifies (`boxed` in `slice.rs`): a box face lies exactly on
`x_k = origin[k]` because that is an `f64`. Level splits b and c pass unchanged.
Why the boolean's rounding certifies them was not traced: its arrangement
triangulates the cap through the box's top face, which gives the descent
other faces to act on. This is evidence, not a proof. An inclined plane has no
box face exactly on it, so an inclined split of a block holding a step thinner
than rounding at the plane can still raise `uncertified`; no attack case does
today.

## 12. Time

The attack cases that certify but take minutes spend their time in two places
(release build, one profiled run each, not a benchmark). In the union of the
lid chain, `insert_hot_points` took 8.2 of 13.4 s of rounding: every insertion,
collapse, flip and flat-face split looked for faces on an edge or at a vertex by
scanning every face ever made. `Global` now keeps, per vertex, every face that
held it at some write; a lookup filters that list by liveness and current
corners and sorts it, so it returns the faces a scan returns, in the same
order, and results keep their bytes
(`the_vertex_index_finds_what_a_scan_finds`). Faces are written only by `push`,
`merge_equal`, `join_opposite` and `collapse`, and each lists the face it
writes; a face that dies and comes back keeps its corners, so it is still
listed.

In the chain difference and the pushback union of `attacks.rs` the descent
takes most of the time (34 of 38 s, and 87 of 104 s): for each candidate
change it validates the faces whose boxes meet the changed region twice, before
and after, about 5 ms per candidate over 8693 candidates in the chain
difference. Counting only the problems that involve a changed face would make
that local, but the count must equal `validate`'s exactly for the descent to
choose the same change, which needs `validate`'s pair tests and edge rules
exposed per face. Not done in this pass; the tests that ask these cases to
finish within a minute in the debug build stay ignored.
