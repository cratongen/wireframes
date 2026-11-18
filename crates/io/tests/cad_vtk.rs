mod cad_common;

use cad_common::*;
use drusa_core::{Column, ColumnData, Error, Mesh, PointSet, Polylines, Table};
use drusa_io::{Data, UnstructuredGrid, WriteOptions, vtk};
use proptest::prelude::*;

pub fn volume_grid() -> UnstructuredGrid {
    let mut points: Vec<[f64; 3]> = Vec::new();
    for z in [0.0, 1.0] {
        for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
            points.push([x, y, z]);
        }
    }
    points.push([0.25, 0.25, -1.0]);
    UnstructuredGrid {
        points,
        cell_types: vec![12, 10, 5],
        offsets: vec![0, 8, 12, 15],
        connectivity: vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 3, 8, 4, 5, 6],
        faces: None,
        point_data: Table::new(9)
            .with_column("p", Column::f64((0..9).map(|i| i as f64 * 0.5).collect()))
            .unwrap(),
        cell_data: Table::new(3)
            .with_column("perm", Column::f64(vec![1e-3, 2.5, -7.0]))
            .unwrap()
            .with_column(
                "zone",
                Column::i64(vec![3, 4, 5])
                    .with_valid(vec![true, false, true])
                    .unwrap(),
            )
            .unwrap(),
    }
}

fn round(data: &Data, opts: &WriteOptions) -> Data {
    let bytes = vtk::write(data, opts).unwrap();
    vtk::read(&bytes).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn mesh_round_trips(m in attributed_mesh(), ascii_mode in any::<bool>()) {
        let opts = WriteOptions { ascii: ascii_mode, ..WriteOptions::default() };
        let back = round(&mesh_data(m.clone()), &opts);
        prop_assert_eq!(back.mesh, Some(m));
        prop_assert!(back.polylines.is_none() && back.points.is_none() && back.grid.is_none());
    }

    /// Mixed files drop columns that are null on every row of one kind, and a
    /// mesh without faces is not told apart from loose points there.
    #[test]
    fn mixed_round_trips(m in attributed_mesh(), pl in polylines(true), ps in points(), ascii_mode in any::<bool>()) {
        prop_assume!(m.n_faces() > 0);
        prop_assume!(m.vertex_attributes().columns().iter().all(|c| c.null_count() < c.len()));
        let ps = ps.clone().with_column("w", Column::f64(vec![2.0; ps.len()])).unwrap();
        let data = Data { mesh: Some(m), polylines: Some(pl), points: Some(ps), ..Data::default() };
        let opts = WriteOptions { ascii: ascii_mode, ..WriteOptions::default() };
        let back = round(&data, &opts);
        prop_assert_eq!(back.polylines, data.polylines);
        prop_assert_eq!(back.points, data.points);
        prop_assert_eq!(back.mesh, data.mesh);
    }
}

#[test]
fn multipart_features_and_closed_parts_round_trip() {
    let pl = Polylines::from_parts(
        vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            vec![[0.2, 0.2, 0.0], [0.4, 0.2, 0.0], [0.4, 0.4, 0.0]],
            vec![[5.0, 5.0, 1.0], [6.0, 5.0, 1.0]],
        ],
        vec![true, true, false],
        Some(vec![0, 0, 1]),
    )
    .unwrap()
    .with_column("name_id", Column::f64(vec![1.5, 2.5]))
    .unwrap();
    let data = Data {
        polylines: Some(pl.clone()),
        ..Data::default()
    };
    for opts in [ascii(), WriteOptions::default()] {
        assert_eq!(round(&data, &opts).polylines, Some(pl.clone()));
    }
}

#[test]
fn volume_grid_round_trips() {
    let g = volume_grid();
    let data = Data {
        grid: Some(g.clone()),
        ..Data::default()
    };
    for opts in [ascii(), WriteOptions::default()] {
        let back = round(&data, &opts);
        assert_eq!(back.grid, Some(g.clone()));
        assert!(back.mesh.is_none());
    }
}

#[test]
fn utf8_columns_need_the_skip_flag() {
    let m = Mesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
    .with_face_column("rock", Column::utf8(vec!["ox".into()]))
    .unwrap();
    let data = mesh_data(m.clone());
    assert!(matches!(
        vtk::write(&data, &ascii()),
        Err(Error::Unsupported(_))
    ));
    let opts = WriteOptions {
        skip_unsupported_columns: true,
        ..ascii()
    };
    let back = vtk::read(&vtk::write(&data, &opts).unwrap()).unwrap();
    assert_eq!(back.mesh.unwrap().face_attributes().num_columns(), 0);
}

#[test]
fn grid_with_surface_data_is_refused() {
    let data = Data {
        grid: Some(volume_grid()),
        points: Some(PointSet::new(vec![[0.0; 3]]).unwrap()),
        ..Data::default()
    };
    assert!(matches!(
        vtk::write(&data, &ascii()),
        Err(Error::Unsupported(_))
    ));
}

const LEGACY_42: &str = "# vtk DataFile Version 4.2
a cube face and a strip
ASCII
DATASET POLYDATA
POINTS 6 float
0 0 0 1 0 0 1 1 0
0 1 0 2 0 0 2 1 0
METADATA
INFORMATION 0

POLYGONS 1 5
4 0 1 2 3
TRIANGLE_STRIPS 1 5
4 1 4 2 5
POINT_DATA 6
SCALARS my%20temp double
LOOKUP_TABLE default
1 2 3 4 5 6
VECTORS vel float
1 0 0 0 1 0 0 0 1 1 1 1 0 0 0 2 2 2
CELL_DATA 2
FIELD FieldData 1
rock 1 2 int
7 8
";

#[test]
fn reads_legacy_polygons_strips_and_attributes() {
    let d = vtk::read(LEGACY_42.as_bytes()).unwrap();
    let m = d.mesh.unwrap();
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3], [1, 4, 2], [2, 4, 5]]);
    let fa = m.face_attributes();
    assert_eq!(fa.names(), &["rock", "polygon"]);
    assert_eq!(
        fa.column("rock").unwrap().data,
        ColumnData::Int64(vec![7, 7, 8, 8])
    );
    assert_eq!(
        fa.column("polygon").unwrap().data,
        ColumnData::Int64(vec![0, 0, 1, 2])
    );
    let va = m.vertex_attributes();
    assert_eq!(va.names(), &["my temp", "vel_0", "vel_1", "vel_2"]);
    assert_eq!(
        va.column("vel_2").unwrap().data,
        ColumnData::Float64(vec![0.0, 0.0, 1.0, 1.0, 0.0, 2.0])
    );
    let again =
        vtk::read(&vtk::write(&mesh_data(m.clone()), &Default::default()).unwrap()).unwrap();
    assert_eq!(again.mesh, Some(m));
}

fn be_i32(v: &[i32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_be_bytes()).collect()
}

#[test]
fn reads_old_binary_unstructured_grid() {
    let mut b =
        b"# vtk DataFile Version 3.0\nt\nBINARY\nDATASET UNSTRUCTURED_GRID\nPOINTS 4 float\n"
            .to_vec();
    for p in [[0f32, 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]] {
        for x in p {
            b.extend(x.to_be_bytes());
        }
    }
    b.extend(b"\nCELLS 1 5\n");
    b.extend(be_i32(&[4, 0, 1, 2, 3]));
    b.extend(b"\nCELL_TYPES 1\n");
    b.extend(be_i32(&[10]));
    b.extend(b"\nCELL_DATA 1\nSCALARS id int 1\nLOOKUP_TABLE default\n");
    b.extend(be_i32(&[42]));
    b.push(b'\n');
    let g = vtk::read(&b).unwrap().grid.unwrap();
    assert_eq!(g.cell_types, vec![10]);
    assert_eq!(g.offsets, vec![0, 4]);
    assert_eq!(g.connectivity, vec![0, 1, 2, 3]);
    assert_eq!(g.points[3], [0.0, 0.0, 1.0]);
    assert_eq!(
        g.cell_data.column("id").unwrap().data,
        ColumnData::Int64(vec![42])
    );
}

#[test]
fn surface_only_unstructured_grid_reads_as_mesh() {
    let src = "# vtk DataFile Version 5.1\nx\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 4 double\n0 0 0 1 0 0 1 1 0 0 1 0\nCELLS 2 4\nOFFSETS vtktypeint64\n0 4\nCONNECTIVITY vtktypeint64\n0 1 2 3\nCELL_TYPES 1\n9\n";
    let m = vtk::read(src.as_bytes()).unwrap().mesh.unwrap();
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3]]);
}

#[test]
fn malformed_files_error() {
    let cases: &[&str] = &[
        "",
        "hello\n",
        "# vtk DataFile Version 2.0\nt\nXML\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET POLYDATA\nPOINTS 2 float\n0 0 0\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET POLYDATA\nPOINTS 1 float\n0 0 0\nPOLYGONS 1 4\n3 0 1 2\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET POLYDATA\nPOINTS 1 float\n0 0 0\nPOLYGONS 1 4\n9 0 0 0\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET STRUCTURED_POINTS\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET POLYDATA\nPOINTS 1 float\n0 0 0\nPOINT_DATA 2\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET POLYDATA\nPOINTS 1 bogus\n0 0 0\n",
        "# vtk DataFile Version 2.0\nt\nASCII\nDATASET POLYDATA\nSCALARS a float\n",
    ];
    for c in cases {
        assert!(vtk::read(c.as_bytes()).is_err(), "{c:?}");
    }
}

#[test]
fn damaged_files_never_panic() {
    let data = Data {
        mesh: Some(
            Mesh::new(
                vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
                vec![[0, 1, 2], [1, 3, 2]],
            )
            .unwrap(),
        ),
        polylines: Some(
            Polylines::from_parts(vec![vec![[0.0; 3], [2.0, 2.0, 2.0]]], vec![false], None)
                .unwrap(),
        ),
        ..Data::default()
    };
    for opts in [ascii(), WriteOptions::default()] {
        let bytes = vtk::write(&data, &opts).unwrap();
        never_panics(&bytes, |b| {
            let _ = vtk::read(b);
        });
    }
    let grid = Data {
        grid: Some(volume_grid()),
        ..Data::default()
    };
    never_panics(&vtk::write(&grid, &Default::default()).unwrap(), |b| {
        let _ = vtk::read(b);
    });
    never_panics(LEGACY_42.as_bytes(), |b| {
        let _ = vtk::read(b);
    });
}
