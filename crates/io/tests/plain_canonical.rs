mod plain_support;

use plain_support::*;
use proptest::prelude::*;
use wireframes_core::{Column, Mesh, Table};
use wireframes_io::{Data, Format, WriteOptions, canonicalize, obj, off, ply, stl};

const PLAIN: [Format; 4] = [Format::Obj, Format::Ply, Format::Off, Format::Stl];

fn write(f: Format, d: &Data, o: &WriteOptions) -> wireframes_core::Result<Vec<u8>> {
    match f {
        Format::Obj => obj::write(d, o),
        Format::Ply => ply::write(d, o),
        Format::Off => off::write(d, o),
        Format::Stl => stl::write(d, o),
        _ => unreachable!(),
    }
}

fn read(f: Format, b: &[u8]) -> Mesh {
    let d = match f {
        Format::Obj => obj::read(b),
        Format::Ply => ply::read(b),
        Format::Off => off::read(b),
        Format::Stl => stl::read(b),
        _ => unreachable!(),
    };
    d.unwrap().mesh.unwrap_or_else(Mesh::empty)
}

/// STL ASCII stores `solid` and binary `stl_attribute`, never both.
fn for_mode(m: Mesh, f: Format, ascii: bool) -> Mesh {
    if f != Format::Stl {
        return m;
    }
    let drop = if ascii { "stl_attribute" } else { "solid" };
    let vt = m.vertex_attributes().clone();
    let ft = m.face_attributes().clone().without_column(drop);
    m.with_attributes(vt, ft).unwrap()
}

fn column(n: usize) -> impl Strategy<Value = Column> {
    let nulls = prop::collection::vec(prop::bool::weighted(0.85), n);
    let ints = prop_oneof![
        -5i64..5,
        any::<i32>().prop_map(i64::from),
        Just(1i64 << 53),
        Just((1i64 << 53) + 1),
        Just(i64::MIN),
        0i64..70_000,
    ];
    let names = prop::sample::select(vec![
        "",
        "a",
        "wireframes",
        " pit ",
        "zone#2",
        "a  b",
        "%null",
    ]);
    prop_oneof![
        (
            prop::collection::vec(prop_oneof![coord(), Just(f64::NAN)], n),
            nulls.clone()
        )
            .prop_map(|(v, ok)| Column::f64(v).with_valid(ok).unwrap()),
        (prop::collection::vec(ints, n), nulls.clone(), any::<bool>()).prop_map(|(v, ok, keep)| {
            let c = Column::i64(v);
            if keep { c } else { c.with_valid(ok).unwrap() }
        }),
        (prop::collection::vec(any::<bool>(), n), nulls.clone()).prop_map(|(v, ok)| Column::bool(
            v
        )
        .with_valid(ok)
        .unwrap()),
        (prop::collection::vec(names, n), nulls).prop_map(|(v, ok)| {
            Column::utf8(v.into_iter().map(String::from).collect())
                .with_valid(ok)
                .unwrap()
        }),
    ]
}

fn table(n: usize, names: &'static [&'static str]) -> impl Strategy<Value = Table> {
    prop::collection::vec(prop::option::of(column(n)), names.len()).prop_map(move |cols| {
        let mut t = Table::new(n);
        for (name, c) in names.iter().zip(cols) {
            if let Some(c) = c {
                t = t.with_column(name, c).unwrap();
            }
        }
        t
    })
}

const VERTEX: &[&str] = &["grade", "x", "bad name", "flag"];
const FACE: &[&str] = &[
    "object",
    "rock",
    "group",
    "solid",
    "stl_attribute",
    "vertex_index",
];

fn any_mesh() -> impl Strategy<Value = Mesh> {
    (polygon_mesh(10, 6), any::<bool>(), any::<bool>()).prop_flat_map(|((m, _), first, crs)| {
        let (nv, nf) = (m.n_vertices(), m.n_faces());
        (Just(m), table(nv, VERTEX), table(nf, FACE)).prop_map(move |(m, vt, ft)| {
            let poly = m.face_attributes().column("polygon").cloned();
            let mut ft = ft;
            if let Some(p) = poly {
                ft = if first {
                    let mut t = Table::new(nf).with_column("polygon", p).unwrap();
                    for (name, c) in ft.names().iter().zip(ft.columns()) {
                        t = t.with_column(name, c.clone()).unwrap();
                    }
                    t
                } else {
                    ft.with_column("polygon", p).unwrap()
                };
            }
            let m = m.with_attributes(vt, ft).unwrap();
            if crs {
                m.with_crs(Some("EPSG:31983".into()))
            } else {
                m
            }
        })
    })
}

fn skip(ascii: bool) -> WriteOptions {
    skipping(ascii)
}

#[test]
fn canonicalize_does_not_panic() {
    let m = Mesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    for f in PLAIN {
        assert_eq!(canonicalize(&m, f), m, "{f:?}");
    }
}

#[test]
fn stl_canonical_welds_rounds_and_drops_unused() {
    let m = Mesh::new(
        vec![
            [9.0, 9.0, 9.0],
            [0.1, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-0.0, 1.0, 0.0],
            [0.1, 0.0, 0.0],
        ],
        vec![[1, 2, 3], [5, 4, 2]],
    )
    .unwrap();
    let c = canonicalize(&m, Format::Stl);
    assert_eq!(
        c.vertices(),
        &[
            [f64::from(0.1f32), 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
    );
    assert_eq!(c.triangles(), &[[0, 1, 2], [0, 2, 1]]);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn write_then_read_is_identity_for_canonical_meshes(m in any_mesh()) {
        for f in PLAIN {
            let c = canonicalize(&m, f);
            prop_assert_eq!(&canonicalize(&c, f), &c, "{:?} is not idempotent", f);
            for ascii in [true, false] {
                let c = for_mode(c.clone(), f, ascii);
                let d = mesh_data(c.clone());
                let opts = WriteOptions { ascii, ..WriteOptions::default() };
                let back = read(f, &write(f, &d, &opts).unwrap());
                prop_assert_eq!(back, c, "{:?} ascii={}", f, ascii);
            }
        }
    }

    #[test]
    fn canonicalize_predicts_a_skipping_round_trip(m in any_mesh()) {
        for f in PLAIN {
            for ascii in [true, false] {
                if f == Format::Stl && ascii {
                    continue;
                }
                let src = for_mode(m.clone(), f, ascii);
                let Ok(bytes) = write(f, &mesh_data(src.clone()), &skip(ascii)) else {
                    continue;
                };
                let want = for_mode(canonicalize(&src, f), f, ascii);
                prop_assert_eq!(read(f, &bytes), want, "{:?} ascii={}", f, ascii);
            }
        }
    }
}
