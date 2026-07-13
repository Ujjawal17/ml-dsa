//! Cachegrind target for `SampleInBall`'s data-dependent memory behaviour.
//!
//! Runs `sample_in_ball` on a seed byte chosen from argv, so two runs with different
//! seeds can be compared under `valgrind --tool=cachegrind`. `cg_annotate --auto=yes`
//! then shows, per source line in `sample.rs`:
//!   * the rejection-sampling line (`squeeze` / `j = byte`) — data accesses that VARY
//!     with the seed (the variable loop-count leak), and
//!   * the write `c.coeffs[j] = c.coeffs[i]` — executed a FIXED τ = 49 times, but at
//!     a secret-derived index `j` (the data-address / cache-line leak).

use std::hint::black_box;

use ml_dsa::params::MlDsa65;
use ml_dsa::sample::sample_in_ball;

fn main() {
    let b: u8 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    let seed = vec![b; 48]; // c~ is 48 bytes for ML-DSA-65
    black_box(sample_in_ball::<MlDsa65>(black_box(&seed)));
}
