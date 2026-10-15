// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Minimal .xyz reader: atom count, a comment line, then `label x y z` lines in Angstrom.
//!
//! Extra columns after `x y z` (as in extended xyz) are ignored.

use std::error::Error;
use std::fs;
use std::path::Path;

/// One atom of an .xyz file.
#[derive(Debug, Clone)]
pub struct XyzAtom {
    /// Label as given in the file, e.g. "Y" or "Y1".
    pub label: String,
    /// Element read from the label: leading letters, e.g. "Y1" -> "Y", "OX" -> "O"... see
    /// [`element_from_label`].
    pub element: String,
    /// Cartesian position in Angstrom.
    pub position: [f64; 3],
}

pub fn read_file(path: impl AsRef<Path>) -> Result<Vec<XyzAtom>, Box<dyn Error>> {
    let path = path.as_ref();
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&text).map_err(|e| format!("{}: {e}", path.display()).into())
}

pub fn parse(text: &str) -> Result<Vec<XyzAtom>, Box<dyn Error>> {
    let mut lines = text.lines();
    let count: usize = lines.next()
        .ok_or("empty file")?
        .trim()
        .parse()
        .map_err(|_| "first line is not an atom count")?;
    lines.next().ok_or("missing comment line")?;
    let mut atoms = Vec::with_capacity(count);
    for line in lines.filter(|l| !l.trim().is_empty()).take(count) {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(format!("expected 'label x y z', got '{line}'").into());
        }
        let coordinate = |s: &str| s.parse::<f64>().map_err(|_| format!("bad coordinate '{s}' in '{line}'"));
        let label = fields[0].to_string();
        atoms.push(XyzAtom {
            element: element_from_label(&label).ok_or_else(|| format!("cannot tell the element of '{label}'"))?,
            label,
            position: [coordinate(fields[1])?, coordinate(fields[2])?, coordinate(fields[3])?],
        });
    }
    if atoms.len() != count {
        return Err(format!("expected {count} atoms, found {}", atoms.len()).into());
    }
    Ok(atoms)
}

/// "Y1" -> "Y", "pr" -> "Pr", "O" -> "O": first letter plus following lower-case letters.
fn element_from_label(label: &str) -> Option<String> {
    let mut chars = label.chars();
    let first = chars.next().filter(|c| c.is_ascii_alphabetic())?;
    let mut element = first.to_ascii_uppercase().to_string();
    element.extend(chars.take_while(|c| c.is_ascii_lowercase()));
    Some(element)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extended_xyz() {
        let atoms = parse("2\nLattice=\"1 0 0 0 1 0 0 0 1\"\n Y    0.0  2.65  2.99\nO1 1 2 3 extra\n").unwrap();
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].element, "Y");
        assert_eq!(atoms[1].label, "O1");
        assert_eq!(atoms[1].element, "O");
        assert_eq!(atoms[1].position, [1.0, 2.0, 3.0]);
        assert!(parse("3\n\nY 0 0 0\n").is_err());
    }
}
