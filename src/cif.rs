// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Minimal CIF 1.1 reader.
//!
//! Turns a file into data blocks made of plain items (`_tag value`) and
//! loops (`loop_` + tags + rows). It knows nothing about crystallography;
//! that lives in `crystal.rs`. Lookups are by tag name, so the order of
//! lines and columns in the file does not matter.

use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::path::Path;

#[derive(Debug)]
pub struct CifError(pub String);

impl fmt::Display for CifError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for CifError {}

pub type Result<T> = std::result::Result<T, CifError>;

pub fn err<T>(msg: impl Into<String>) -> Result<T> {
    Err(CifError(msg.into()))
}

/// A loop: column names plus rows of raw string values.
#[derive(Debug, Clone)]
pub struct Loop {
    pub tags: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl Loop {
    /// Index of the column with this tag, if present.
    pub fn column(&self, tag: &str) -> Option<usize> {
        let tag = normalize_tag(tag);
        self.tags.iter().position(|t| *t == tag)
    }
}

#[derive(Debug, Clone)]
pub struct DataBlock {
    pub name: String,
    items: HashMap<String, String>,
    loops: Vec<Loop>,
}

impl DataBlock {
    /// Value of a plain (non-loop) item.
    pub fn get(&self, tag: &str) -> Option<&str> {
        self.items.get(&normalize_tag(tag)).map(|s| s.as_str())
    }

    /// The loop that contains this tag as one of its columns.
    pub fn loop_with(&self, tag: &str) -> Option<&Loop> {
        self.loops.iter().find(|l| l.column(tag).is_some())
    }
}

/// CIF tags are case-insensitive, and mmCIF/CIF2 write `_a.b` for `_a_b`.
pub fn normalize_tag(tag: &str) -> String {
    tag.to_ascii_lowercase().replace('.', "_")
}

/// Parses a CIF number, dropping a standard uncertainty: "10.6078(3)" -> 10.6078.
/// "?" (unknown) and "." (not applicable) give None.
pub fn parse_number(s: &str) -> Option<f64> {
    let s = s.trim();
    let s = match s.find('(') {
        Some(i) => &s[..i],
        None => s,
    };
    s.parse().ok()
}

pub fn read_file(path: impl AsRef<Path>) -> Result<Vec<DataBlock>> {
    let path = path.as_ref();
    let text = fs::read_to_string(path)
        .map_err(|e| CifError(format!("cannot read {}: {e}", path.display())))?;
    parse(&text)
}

pub fn parse(text: &str) -> Result<Vec<DataBlock>> {
    let mut blocks: Vec<DataBlock> = Vec::new();
    let mut tokens = tokenize(text).into_iter().peekable();

    while let Some(token) = tokens.next() {
        let word = match token {
            Token::Bare(w) => w,
            Token::Quoted(v) => return err(format!("value '{v}' without a tag")),
        };
        let lower = word.to_ascii_lowercase();

        if lower.starts_with("data_") {
            blocks.push(DataBlock {
                name: word[5..].to_string(),
                items: HashMap::new(),
                loops: Vec::new(),
            });
        } else if lower == "loop_" {
            let block = current_block(&mut blocks)?;
            let mut tags = Vec::new();
            while let Some(Token::Bare(t)) = tokens.peek() {
                if !t.starts_with('_') {
                    break;
                }
                tags.push(normalize_tag(t));
                tokens.next();
            }
            if tags.is_empty() {
                return err("loop_ without any tags");
            }
            let mut values = Vec::new();
            while let Some(t) = tokens.peek() {
                if t.is_keyword() {
                    break;
                }
                values.push(tokens.next().unwrap().into_string());
            }
            if values.len() % tags.len() != 0 {
                return err(format!(
                    "loop starting with {} has {} values, not a multiple of its {} columns",
                    tags[0],
                    values.len(),
                    tags.len()
                ));
            }
            let rows = values.chunks(tags.len()).map(|r| r.to_vec()).collect();
            block.loops.push(Loop { tags, rows });
        } else if word.starts_with('_') {
            let block = current_block(&mut blocks)?;
            match tokens.next() {
                Some(t) if !t.is_keyword() => {
                    block.items.insert(normalize_tag(&word), t.into_string());
                }
                _ => return err(format!("tag {word} has no value")),
            }
        } else if lower.starts_with("save_") || lower == "global_" {
            return err(format!("'{word}' (save frames / global blocks) is not supported"));
        } else {
            return err(format!("value '{word}' without a tag"));
        }
    }
    Ok(blocks)
}

fn current_block(blocks: &mut [DataBlock]) -> Result<&mut DataBlock> {
    match blocks.last_mut() {
        Some(b) => Ok(b),
        None => err("data found before the first data_ line"),
    }
}

enum Token {
    Bare(String),
    /// Quoted strings and semicolon text fields: always values, never tags.
    Quoted(String),
}

impl Token {
    fn is_keyword(&self) -> bool {
        match self {
            Token::Bare(w) => {
                let w = w.to_ascii_lowercase();
                w.starts_with('_')
                    || w == "loop_"
                    || w == "global_"
                    || w.starts_with("data_")
                    || w.starts_with("save_")
            }
            Token::Quoted(_) => false,
        }
    }

    fn into_string(self) -> String {
        match self {
            Token::Bare(s) | Token::Quoted(s) => s,
        }
    }
}

fn tokenize(text: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        // Semicolon text field: from a line starting with ';' to the next one.
        if let Some(first) = line.strip_prefix(';') {
            let mut value = first.to_string();
            let mut rest_of_closing_line = "";
            for l in lines.by_ref() {
                if let Some(rest) = l.strip_prefix(';') {
                    rest_of_closing_line = rest;
                    break;
                }
                value.push('\n');
                value.push_str(l);
            }
            tokens.push(Token::Quoted(value));
            tokenize_line(rest_of_closing_line, &mut tokens);
        } else {
            tokenize_line(line, &mut tokens);
        }
    }
    tokens
}

fn tokenize_line(line: &str, tokens: &mut Vec<Token>) {
    let chars: Vec<char> = line.chars().collect();
    let n = chars.len();
    let mut i = 0;
    while i < n {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
        } else if c == '#' {
            break; // comment to end of line
        } else if c == '\'' || c == '"' {
            // A quote only closes if followed by whitespace or end of line,
            // so 'O'Brien' style values survive.
            let start = i + 1;
            let mut j = start;
            while j < n && !(chars[j] == c && (j + 1 == n || chars[j + 1].is_whitespace())) {
                j += 1;
            }
            tokens.push(Token::Quoted(chars[start..j.min(n)].iter().collect()));
            i = j + 1;
        } else {
            let start = i;
            while i < n && !chars[i].is_whitespace() {
                i += 1;
            }
            tokens.push(Token::Bare(chars[start..i].iter().collect()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn items_loops_quotes_comments() {
        let text = "\
# comment
data_test
_cell_length_a 10.6078(3)  # trailing comment
_symmetry_space_group_name_H-M 'I a -3'
_note
;
multi-line
text
;
loop_
_symmetry_equiv_pos_as_xyz
'x,y,z'
\"-x, -y, -z\"
";
        let blocks = parse(text).unwrap();
        let b = &blocks[0];
        assert_eq!(b.name, "test");
        assert_eq!(parse_number(b.get("_CELL_LENGTH_A").unwrap()), Some(10.6078));
        assert_eq!(b.get("_symmetry_space_group_name_h-m"), Some("I a -3"));
        assert_eq!(b.get("_note"), Some("\nmulti-line\ntext"));
        let l = b.loop_with("_symmetry_equiv_pos_as_xyz").unwrap();
        assert_eq!(l.rows, vec![vec!["x,y,z".to_string()], vec!["-x, -y, -z".to_string()]]);
    }

    #[test]
    fn ragged_loop_is_an_error() {
        assert!(parse("data_x\nloop_\n_a\n_b\n1 2 3\n").is_err());
    }
}
