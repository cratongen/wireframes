use std::path::PathBuf;

use crate::MeshReport;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0}")]
    InvalidInput(String),
    #[error(
        "invalid mesh: {context}: {}; the error's report lists each problem (as mesh.validate() does), \
         and mesh.repair() or repair=True fixes what repair can",
        report.summary.defects()
    )]
    InvalidMesh {
        context: String,
        report: Box<MeshReport>,
    },
    #[error(
        "kernel limitation, not an input defect: {message}{}",
        report.as_ref().map_or(String::new(), |r| format!(
            ": {}; the error's report lists each problem and its mesh holds the result",
            r.summary.defects()
        ))
    )]
    Kernel {
        message: String,
        report: Option<Box<MeshReport>>,
        mesh: Option<Box<crate::Mesh>>,
    },
    #[error("{}", format_location(path, line, message))]
    Format {
        path: Option<PathBuf>,
        line: Option<u64>,
        message: String,
    },
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Unsupported(String),
    #[error("crs mismatch: {left} vs {right}")]
    CrsMismatch { left: String, right: String },
    #[error("length unit mismatch: {left} vs {right}")]
    UnitMismatch { left: String, right: String },
    #[error("missing column {name:?}; columns are {available:?}")]
    MissingColumn {
        name: String,
        available: Vec<String>,
    },
}

pub type Result<T> = std::result::Result<T, Error>;

fn format_location(path: &Option<PathBuf>, line: &Option<u64>, message: &str) -> String {
    match (path, line) {
        (Some(p), Some(l)) => format!("{}:{l}: {message}", p.display()),
        (Some(p), None) => format!("{}: {message}", p.display()),
        (None, Some(l)) => format!("line {l}: {message}"),
        (None, None) => message.to_string(),
    }
}

impl Error {
    pub fn invalid(message: impl Into<String>) -> Error {
        Error::InvalidInput(message.into())
    }

    pub fn kernel(message: impl Into<String>) -> Error {
        Error::Kernel {
            message: message.into(),
            report: None,
            mesh: None,
        }
    }

    /// A result that `validate` does not pass as a solid, with its report.
    pub fn kernel_result(message: impl Into<String>, mesh: crate::Mesh) -> Result<Error> {
        let report = mesh.validate(0.0)?;
        Ok(Error::Kernel {
            message: message.into(),
            report: Some(Box::new(report)),
            mesh: Some(Box::new(mesh)),
        })
    }

    pub fn format(message: impl Into<String>) -> Error {
        Error::Format {
            path: None,
            line: None,
            message: message.into(),
        }
    }

    pub fn format_at(line: u64, message: impl Into<String>) -> Error {
        Error::Format {
            path: None,
            line: Some(line),
            message: message.into(),
        }
    }

    /// Attaches `path` to a `Format` error that has none; other variants pass through.
    pub fn with_path(self, path: impl Into<PathBuf>) -> Error {
        match self {
            Error::Format {
                path: None,
                line,
                message,
            } => Error::Format {
                path: Some(path.into()),
                line,
                message,
            },
            e => e,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_messages_name_their_location() {
        assert_eq!(Error::format("bad").to_string(), "bad");
        assert_eq!(Error::format_at(3, "bad").to_string(), "line 3: bad");
        let e = Error::format_at(3, "bad").with_path("a.obj");
        assert_eq!(e.to_string(), "a.obj:3: bad");
        let e = Error::format("bad").with_path("a.obj");
        assert_eq!(e.to_string(), "a.obj: bad");
        assert!(matches!(
            Error::invalid("x").with_path("p"),
            Error::InvalidInput(_)
        ));
    }

    #[test]
    fn invalid_mesh_messages_name_the_defects_and_the_way_out() {
        let open = crate::Mesh::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let e = Error::InvalidMesh {
            context: "a".into(),
            report: Box::new(open.validate(0.0).unwrap()),
        };
        let m = e.to_string();
        assert!(m.starts_with("invalid mesh: a: boundary_edges=3; "), "{m}");
        assert!(
            m.contains("mesh.validate()") && m.contains("repair=True"),
            "{m}"
        );
    }

    #[test]
    fn kernel_messages_say_the_input_is_not_at_fault() {
        assert_eq!(
            Error::kernel("no ear").to_string(),
            "kernel limitation, not an input defect: no ear"
        );
        let open = crate::Mesh::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let e = Error::kernel_result("rounding", open).unwrap();
        assert_eq!(
            e.to_string(),
            "kernel limitation, not an input defect: rounding: boundary_edges=3; \
             the error's report lists each problem and its mesh holds the result"
        );
    }
}
