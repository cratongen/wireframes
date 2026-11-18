use std::path::PathBuf;

use crate::MeshReport;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0}")]
    InvalidInput(String),
    #[error("invalid mesh: {context}")]
    InvalidMesh {
        context: String,
        report: Box<MeshReport>,
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
}
