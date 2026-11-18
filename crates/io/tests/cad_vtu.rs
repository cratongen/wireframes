mod cad_common;

use std::io::Write;

use cad_common::*;
use drusa_core::{Column, ColumnData, Error, Mesh, PointSet, Polylines, Table};
use drusa_io::{Data, UnstructuredGrid, WriteOptions, vtu};
use flate2::{Compression, write::ZlibEncoder};
use proptest::prelude::*;

fn volume_grid() -> UnstructuredGrid {
    let mut points: Vec<[f64; 3]> = Vec::new();
    for z in [0.0, 1.0] {
        for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
            points.push([x, y, z]);
        }
    }
    points.push([0.25, 0.25, -1.0]);
    UnstructuredGrid {
        points,
        cell_types: vec![12, 10],
        offsets: vec![0, 8, 12],
        connectivity: vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 3, 8],
        faces: None,
        point_data: Table::new(9)
            .with_column("p", Column::f64((0..9).map(|i| i as f64).collect()))
            .unwrap(),
        cell_data: Table::new(2)
            .with_column("zone", Column::i64(vec![7, -2]))
            .unwrap()
            .with_column("active", Column::bool(vec![true, false]))
            .unwrap(),
    }
}

/// A tetrahedron stored as a polyhedron (type 42) next to a plain tetra.
fn polyhedron_grid() -> UnstructuredGrid {
    let mut g = volume_grid();
    g.cell_types = vec![42, 10];
    g.offsets = vec![0, 4, 8];
    g.connectivity = vec![0, 1, 3, 8, 0, 1, 3, 8];
    g.cell_data = Table::new(2);
    let stream = vec![4, 3, 0, 1, 3, 3, 0, 1, 8, 3, 1, 3, 8, 3, 3, 0, 8];
    g.faces = Some((vec![stream.len() as i64, -1], stream));
    g
}

fn round(data: &Data, opts: &WriteOptions, vtp: bool) -> Data {
    let bytes = if vtp {
        vtu::write_vtp(data, opts).unwrap()
    } else {
        vtu::write(data, opts).unwrap()
    };
    vtu::read(&bytes).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn mesh_round_trips(m in attributed_mesh(), z in any::<bool>(), vtp in any::<bool>()) {
        let opts = WriteOptions { zlib: z, ..WriteOptions::default() };
        let back = round(&mesh_data(m.clone()), &opts, vtp);
        prop_assert_eq!(back.mesh, Some(m));
    }

    #[test]
    fn vtp_mixed_round_trips(m in attributed_mesh(), pl in polylines(true), ps in points(), z in any::<bool>()) {
        prop_assume!(m.n_faces() > 0);
        prop_assume!(m.vertex_attributes().columns().iter().all(|c| c.null_count() < c.len()));
        let data = Data { mesh: Some(m), polylines: Some(pl), points: Some(ps), ..Data::default() };
        let opts = WriteOptions { zlib: z, ..WriteOptions::default() };
        let back = round(&data, &opts, true);
        prop_assert_eq!(back.mesh, data.mesh);
        prop_assert_eq!(back.polylines, data.polylines);
        prop_assert_eq!(back.points, data.points);
    }
}

#[test]
fn volume_and_polyhedron_grids_round_trip() {
    for g in [volume_grid(), polyhedron_grid()] {
        for opts in [WriteOptions::default(), zlib()] {
            let data = Data {
                grid: Some(g.clone()),
                ..Data::default()
            };
            let back = round(&data, &opts, false);
            let mut want = g.clone();
            if let Some(c) = want.cell_data.column("active").cloned() {
                let ColumnData::Bool(b) = c.data else {
                    unreachable!()
                };
                want.cell_data = want
                    .cell_data
                    .with_column("active", Column::i64(b.iter().map(|&x| x as i64).collect()))
                    .unwrap();
            }
            assert_eq!(back.grid, Some(want));
        }
    }
}

#[test]
fn wrong_containers_are_refused() {
    let pts = Data {
        points: Some(PointSet::new(vec![[0.0; 3]]).unwrap()),
        ..Data::default()
    };
    assert!(matches!(
        vtu::write(&pts, &Default::default()),
        Err(Error::Unsupported(_))
    ));
    let grid = Data {
        grid: Some(volume_grid()),
        ..Data::default()
    };
    assert!(matches!(
        vtu::write_vtp(&grid, &Default::default()),
        Err(Error::Unsupported(_))
    ));
}

fn b64(b: &[u8]) -> String {
    const T: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut s = String::new();
    for c in b.chunks(3) {
        let n = (c[0] as u32) << 16
            | (*c.get(1).unwrap_or(&0) as u32) << 8
            | *c.get(2).unwrap_or(&0) as u32;
        for k in 0..4 {
            if k <= c.len() {
                s.push(T[(n >> (18 - 6 * k) & 63) as usize] as char);
            } else {
                s.push('=');
            }
        }
    }
    s
}

fn le_f32(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

fn le_i32(v: &[i32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

const TETRA_POINTS: [f32; 12] = [0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.];

fn tetra_vtu(points: &str, conn: &str, offsets: &str, types: &str, extra: &str) -> String {
    format!(
        r#"<?xml version="1.0"?>
<!-- a comment -->
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"{extra}>
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points>{points}</Points>
      <Cells>{conn}{offsets}{types}</Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
"#
    )
}

fn check_tetra(d: Data) {
    let g = d.grid.expect("a grid");
    assert_eq!(g.cell_types, vec![10]);
    assert_eq!(g.offsets, vec![0, 4]);
    assert_eq!(g.connectivity, vec![0, 1, 2, 3]);
    assert_eq!(g.points[1], [1.0, 0.0, 0.0]);
}

#[test]
fn reads_ascii_inline() {
    let x = tetra_vtu(
        r#"<DataArray type="Float32" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 0 1 0 0 0 1</DataArray>"#,
        r#"<DataArray type="Int32" Name="connectivity" format="ascii">0 1 2 3</DataArray>"#,
        r#"<DataArray type="Int32" Name="offsets" format="ascii">4</DataArray>"#,
        r#"<DataArray type="UInt8" Name="types" format="ascii">10</DataArray>"#,
        "",
    );
    check_tetra(vtu::read(x.as_bytes()).unwrap());
}

fn inline(ty: &str, name: &str, ncomp: usize, data: &[u8], separate: bool) -> String {
    let header = (data.len() as u32).to_le_bytes();
    let text = if separate {
        b64(&header) + &b64(data)
    } else {
        b64(&[&header[..], data].concat())
    };
    format!(
        r#"<DataArray type="{ty}" Name="{name}" NumberOfComponents="{ncomp}" format="binary">{text}</DataArray>"#
    )
}

#[test]
fn reads_inline_base64_joint_and_separate() {
    for sep in [false, true] {
        let x = tetra_vtu(
            &inline("Float32", "Points", 3, &le_f32(&TETRA_POINTS), sep),
            &inline("Int32", "connectivity", 1, &le_i32(&[0, 1, 2, 3]), sep),
            &inline("Int32", "offsets", 1, &le_i32(&[4]), sep),
            &inline("UInt8", "types", 1, &[10], sep),
            "",
        );
        check_tetra(vtu::read(x.as_bytes()).unwrap());
    }
}

fn zblock(data: &[u8], hdr8: bool) -> Vec<u8> {
    let mut e = ZlibEncoder::new(Vec::new(), Compression::default());
    e.write_all(data).unwrap();
    let c = e.finish().unwrap();
    let h = [1u64, data.len() as u64, data.len() as u64, c.len() as u64];
    let mut out: Vec<u8> = if hdr8 {
        h.iter().flat_map(|x| x.to_le_bytes()).collect()
    } else {
        h.iter().flat_map(|&x| (x as u32).to_le_bytes()).collect()
    };
    out.extend(c);
    out
}

#[test]
fn reads_inline_zlib_with_uint32_headers() {
    let arr = |ty: &str, name: &str, n: usize, d: &[u8]| {
        let z = zblock(d, false);
        let (h, body) = z.split_at(16);
        format!(
            r#"<DataArray type="{ty}" Name="{name}" NumberOfComponents="{n}" format="binary">{}{}</DataArray>"#,
            b64(h),
            b64(body)
        )
    };
    let x = tetra_vtu(
        &arr("Float32", "Points", 3, &le_f32(&TETRA_POINTS)),
        &arr("Int32", "connectivity", 1, &le_i32(&[0, 1, 2, 3])),
        &arr("Int32", "offsets", 1, &le_i32(&[4])),
        &arr("UInt8", "types", 1, &[10]),
        r#" compressor="vtkZLibDataCompressor""#,
    );
    check_tetra(vtu::read(x.as_bytes()).unwrap());
}

#[test]
fn reads_appended_base64_and_big_endian_raw() {
    let parts: Vec<Vec<u8>> = vec![
        le_f32(&TETRA_POINTS),
        le_i32(&[0, 1, 2, 3]),
        le_i32(&[4]),
        vec![10],
    ];
    let mut text = String::new();
    let mut offs = Vec::new();
    for p in &parts {
        offs.push(text.len());
        text += &b64(&(p.len() as u32).to_le_bytes());
        text += &b64(p);
    }
    let arr = |ty: &str, name: &str, n: usize, k: usize| {
        format!(
            r#"<DataArray type="{ty}" Name="{name}" NumberOfComponents="{n}" format="appended" offset="{}"/>"#,
            offs[k]
        )
    };
    let x = tetra_vtu(
        &arr("Float32", "Points", 3, 0),
        &arr("Int32", "connectivity", 1, 1),
        &arr("Int32", "offsets", 1, 2),
        &arr("UInt8", "types", 1, 3),
        "",
    )
    .replace(
        "</VTKFile>",
        &format!("<AppendedData encoding=\"base64\">\n _{text}\n</AppendedData>\n</VTKFile>"),
    );
    check_tetra(vtu::read(x.as_bytes()).unwrap());

    let be: Vec<Vec<u8>> = vec![
        TETRA_POINTS.iter().flat_map(|x| x.to_be_bytes()).collect(),
        [0i64, 1, 2, 3]
            .iter()
            .flat_map(|x| x.to_be_bytes())
            .collect(),
        4i64.to_be_bytes().to_vec(),
        vec![10],
    ];
    let mut raw = Vec::new();
    let mut offs = Vec::new();
    for p in &be {
        offs.push(raw.len());
        raw.extend((p.len() as u64).to_be_bytes());
        raw.extend(p);
    }
    let arr = |ty: &str, name: &str, n: usize, k: usize| {
        format!(
            r#"<DataArray type="{ty}" Name="{name}" NumberOfComponents="{n}" format="appended" offset="{}"/>"#,
            offs[k]
        )
    };
    let head = tetra_vtu(
        &arr("Float32", "Points", 3, 0),
        &arr("Int64", "connectivity", 1, 1),
        &arr("Int64", "offsets", 1, 2),
        &arr("UInt8", "types", 1, 3),
        r#" header_type="UInt64""#,
    )
    .replace("LittleEndian", "BigEndian")
    .replace("</VTKFile>\n", "<AppendedData encoding=\"raw\">_");
    let mut bytes = head.into_bytes();
    bytes.extend(raw);
    bytes.extend(b"</AppendedData></VTKFile>");
    check_tetra(vtu::read(&bytes).unwrap());
}

fn ascii_arr(ty: &str, name: &str, n: usize, v: &str) -> String {
    format!(
        r#"<DataArray type="{ty}" Name="{name}" NumberOfComponents="{n}" format="ascii">{v}</DataArray>"#
    )
}

#[test]
fn reads_several_pieces_with_id_offsets() {
    let piece = |pts: &str, cd: &str| {
        format!(
            r#"<Piece NumberOfPoints="4" NumberOfCells="1"><CellData>{}</CellData><Points>{}</Points><Cells>{}{}{}</Cells></Piece>"#,
            ascii_arr("Float64", "k", 1, cd),
            ascii_arr("Float64", "", 3, pts),
            ascii_arr("Int64", "connectivity", 1, "0 1 2 3"),
            ascii_arr("Int64", "offsets", 1, "4"),
            ascii_arr("UInt8", "types", 1, "10"),
        )
    };
    let x = format!(
        r#"<VTKFile type="UnstructuredGrid"><UnstructuredGrid>{}{}</UnstructuredGrid></VTKFile>"#,
        piece("0 0 0 1 0 0 0 1 0 0 0 1", "1.5"),
        piece("5 0 0 6 0 0 5 1 0 5 0 1", "2.5"),
    );
    let g = vtu::read(x.as_bytes()).unwrap().grid.unwrap();
    assert_eq!(g.points.len(), 8);
    assert_eq!(g.offsets, vec![0, 4, 8]);
    assert_eq!(g.connectivity, vec![0, 1, 2, 3, 4, 5, 6, 7]);
    assert_eq!(
        g.cell_data.column("k").unwrap().data,
        ColumnData::Float64(vec![1.5, 2.5])
    );
}

#[test]
fn reads_vtp_polys_and_strips_with_cell_data_in_vtk_order() {
    let x = format!(
        r#"<VTKFile type="PolyData" version="1.0"><PolyData>
<Piece NumberOfPoints="5" NumberOfVerts="0" NumberOfLines="1" NumberOfStrips="1" NumberOfPolys="1">
<CellData>{}</CellData>
<Points>{}</Points>
<Lines>{}{}</Lines>
<Strips>{}{}</Strips>
<Polys>{}{}</Polys>
</Piece></PolyData></VTKFile>"#,
        ascii_arr("Int32", "tag", 1, "10 20 30"),
        ascii_arr("Float32", "Points", 3, "0 0 0 1 0 0 1 1 0 0 1 0 2 0 0"),
        ascii_arr("Int32", "connectivity", 1, "0 4"),
        ascii_arr("Int32", "offsets", 1, "2"),
        ascii_arr("Int32", "connectivity", 1, "1 4 2"),
        ascii_arr("Int32", "offsets", 1, "3"),
        ascii_arr("Int32", "connectivity", 1, "0 1 2 3"),
        ascii_arr("Int32", "offsets", 1, "4"),
    );
    let d = vtu::read(x.as_bytes()).unwrap();
    let m = d.mesh.unwrap();
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3], [1, 4, 2]]);
    assert_eq!(
        m.face_attributes().column("tag").unwrap().data,
        ColumnData::Int64(vec![20, 20, 30])
    );
    let pl = d.polylines.unwrap();
    assert_eq!(
        pl.attributes().column("tag").unwrap().data,
        ColumnData::Int64(vec![10])
    );
}

#[test]
fn reads_vtk94_polyhedron_arrays_into_legacy_stream() {
    let x = format!(
        r#"<VTKFile type="UnstructuredGrid" version="2.3"><UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="1">
<Points>{}</Points><Cells>{}{}{}{}{}{}{}</Cells></Piece></UnstructuredGrid></VTKFile>"#,
        ascii_arr("Float32", "Points", 3, "0 0 0 1 0 0 0 1 0 0 0 1"),
        ascii_arr("Int64", "connectivity", 1, "0 1 2 3"),
        ascii_arr("Int64", "offsets", 1, "4"),
        ascii_arr("UInt8", "types", 1, "42"),
        ascii_arr("Int64", "face_connectivity", 1, "0 2 1 0 1 3 1 2 3 2 0 3"),
        ascii_arr("Int64", "face_offsets", 1, "3 6 9 12"),
        ascii_arr("Int64", "polyhedron_to_faces", 1, "0 1 2 3"),
        ascii_arr("Int64", "polyhedron_offsets", 1, "4"),
    );
    let g = vtu::read(x.as_bytes()).unwrap().grid.unwrap();
    let (fo, fs) = g.faces.unwrap();
    assert_eq!(fs, vec![4, 3, 0, 2, 1, 3, 0, 1, 3, 3, 1, 2, 3, 3, 2, 0, 3]);
    assert_eq!(fo, vec![17]);
}

#[test]
fn other_compressors_are_refused() {
    let x = r#"<VTKFile type="UnstructuredGrid" compressor="vtkLZ4DataCompressor"><UnstructuredGrid/></VTKFile>"#;
    assert!(matches!(
        vtu::read(x.as_bytes()),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn damaged_files_never_panic() {
    let m = Mesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
        vec![[0, 1, 2], [1, 3, 2]],
    )
    .unwrap();
    let pl =
        Polylines::from_parts(vec![vec![[0.0; 3], [2.0, 2.0, 2.0]]], vec![false], None).unwrap();
    let data = Data {
        mesh: Some(m),
        polylines: Some(pl),
        ..Data::default()
    };
    for opts in [WriteOptions::default(), zlib()] {
        never_panics(&vtu::write_vtp(&data, &opts).unwrap(), |b| {
            let _ = vtu::read(b);
        });
        let g = Data {
            grid: Some(polyhedron_grid()),
            ..Data::default()
        };
        never_panics(&vtu::write(&g, &opts).unwrap(), |b| {
            let _ = vtu::read(b);
        });
    }
}
