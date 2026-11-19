use wireframes_core::{Mesh, fixtures};
use wireframes_io::{WriteOptions, read_mesh, write_mesh};

fn volume(m: &Mesh) -> f64 {
    (0..m.n_faces() as u32)
        .map(|f| {
            let [a, b, c] = m.corners(f);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

#[test]
fn touching_shell_solids_round_trip_as_solids() {
    let dir = std::env::temp_dir().join(format!("wireframes-io-{}-touching", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for (name, m) in [
        ("apart", fixtures::edge_touching_cubes(false)),
        ("welded", fixtures::edge_touching_cubes(true)),
        ("vertex", fixtures::vertex_touching_cubes()),
        ("tetrahedra", fixtures::touching_tetrahedra()),
    ] {
        let want = volume(&m);
        for ext in ["obj", "ply", "off", "stl", "vtk", "vtu", "dxf", "_tr.csv"] {
            let sep = if ext.starts_with('_') { "" } else { "." };
            let p = dir.join(format!("{name}{sep}{ext}"));
            write_mesh(&p, &m, &WriteOptions::default()).unwrap();
            let back = read_mesh(&p).unwrap_or_else(|e| panic!("{name} {ext}: {e}"));
            let s = back.validate(0.0).unwrap().summary;
            assert!(s.is_solid(), "{name} {ext}: {s:?}");
            assert_eq!(s.shells, 2, "{name} {ext}");
            assert!(s.touching_shell_vertices > 0, "{name} {ext}: {s:?}");
            assert!((volume(&back) - want).abs() < 1e-12, "{name} {ext}");
        }
    }
    std::fs::remove_dir_all(dir).unwrap();
}
