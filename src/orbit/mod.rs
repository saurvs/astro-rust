//! Elliptic, parabolic and near-parabolic orbits

pub mod elliptic;
pub mod parabolic;
pub mod near_parabolic;

/// Represents an orbital node
pub enum Node {
    /// Ascending node
    Ascend,
    /// Descending node
    Descend
}
