//! Private helpers shared by the lunar event modules (`phase`, `apsis`,
//! `node`, `declination`).
//!
//! Each of those chapters of Meeus describes a sequence of events whose
//! *mean* instants are equally spaced. Numbering every event of a sequence
//! with a single integer "slot" (for instance four slots per lunation for
//! the phases) lets one iterator type walk through any of them in
//! chronological order.

/// An event that happens at a definite instant.
pub(super) trait Timed {
    /// Instant of the event, as a Julian Ephemeris Day.
    fn jde(&self) -> f64;
}

/// How many slots the iterator may skip while looking for the first event
/// of a range. The starting guess is only a couple of slots early, so this
/// is just a guard against pathological input (huge or rounded-off JDs).
const MAX_SKIP: u32 = 64;

/// Iterator over the events of a sequence whose instants lie in the
/// half-open interval `[start, end)`, in chronological order.
///
/// The true instants of successive events in all four sequences differ from
/// their mean instants by much less than the spacing of the events, so the
/// true instants are monotonic in the slot number and iteration can stop
/// at the first event at or after `end`.
#[derive(Debug, Clone)]
pub(super) struct SlotIter<E> {
    slot: i64,
    end: f64,
    event: fn(i64) -> E,
    done: bool,
}

impl<E: Timed> SlotIter<E> {
    /// Creates the iterator. `first_slot` is an estimate of the slot of the
    /// first event in the range; it must not be later than that slot, so
    /// callers subtract a margin from their estimate.
    pub(super) fn new(start: f64, end: f64, first_slot: f64, event: fn(i64) -> E) -> Self {
        let mut iter = SlotIter {
            slot: 0,
            end,
            event,
            done: true,
        };
        // Also rejects NaN bounds.
        if !(start.is_finite() && first_slot.is_finite() && start < end) {
            return iter;
        }
        // Saturating float-to-int conversion: never panics.
        iter.slot = first_slot.floor() as i64;
        for _ in 0..MAX_SKIP {
            if (iter.event)(iter.slot).jde() >= start {
                iter.done = false;
                return iter;
            }
            match iter.slot.checked_add(1) {
                Some(next) => iter.slot = next,
                None => return iter,
            }
        }
        iter
    }
}

impl<E: Timed> Iterator for SlotIter<E> {
    type Item = E;

    fn next(&mut self) -> Option<E> {
        if self.done {
            return None;
        }
        let event = (self.event)(self.slot);
        if event.jde().is_nan() || event.jde() >= self.end {
            self.done = true;
            return None;
        }
        match self.slot.checked_add(1) {
            Some(next) => self.slot = next,
            None => self.done = true,
        }
        Some(event)
    }
}

impl<E: Timed> std::iter::FusedIterator for SlotIter<E> {}

/// Returns the event of a sequence that is nearest in time to `jd`.
///
/// `mean_number` is the (non-integer) event number whose *mean* instant is
/// `jd`; `event(n)` computes the true event number `n`. The nearest mean
/// event and its neighbour on the other side of `jd` are compared, which is
/// enough because true and mean instants differ by less than half the
/// spacing of the events. `jd` must be finite (callers check).
pub(super) fn nearest<E: Timed>(jd: f64, mean_number: f64, event: impl Fn(i64) -> E) -> E {
    let n0 = mean_number.round() as i64;
    let e0 = event(n0);
    let d0 = jd - e0.jde();
    let n1 = if d0 >= 0.0 {
        n0.checked_add(1)
    } else {
        n0.checked_sub(1)
    };
    if let Some(n1) = n1 {
        let e1 = event(n1);
        if (e1.jde() - jd).abs() < d0.abs() {
            return e1;
        }
    }
    e0
}

/// Reduces an angle in degrees to `[0, 360)` and converts it to radians.
#[inline]
pub(super) fn reduced_rad(deg: f64) -> f64 {
    deg.rem_euclid(360.0).to_radians()
}

/// The factor E (Meeus, eq. 47.6) that accounts for the decreasing
/// eccentricity of the Earth's orbit in terms containing the Sun's mean
/// anomaly. `T` is in Julian centuries from J2000.0.
#[inline]
pub(super) fn eccentricity_factor(T: f64) -> f64 {
    1.0 - T * (0.002_516 + T * 0.000_007_4)
}
