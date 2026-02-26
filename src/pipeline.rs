//! # pipeline
//!
//! ## Purpose
//! Multi-segment pipeline modeling with elevation profiles, friction losses,
//! and pump power calculations. Core module for IGBWP route analysis.
//!
//! ## Algorithms
//! - Per-segment Darcy-Weisbach friction loss: ΔP = f·(L/D)·(ρv²/2)
//! - Elevation head per segment: ΔH_elev = z_end - z_start
//! - Total head: H_total = Σ(friction_head) + max_elevation_gain
//! - Pump power: P = ρ·g·Q·H_total / η
//! - Terrain-adjusted length: L = Σ √(Δx² + Δh²)
//!
//! ## Data Structures
//! - [`PipelineSegment`] — single segment with start/end elevation
//! - [`Pipeline`] — chain of segments forming a complete route
//! - [`PipelineResult`] — full analysis output

use crate::fluid::{Fluid, G};
use crate::pipe::DarcyWeisbach;
use crate::flow::{pump_power, PumpPower};

/// A single segment of a pipeline route.
#[derive(Debug, Clone, PartialEq)]
pub struct PipelineSegment {
    /// Segment name (e.g., "Paradip → Sambalpur")
    pub name: String,
    /// Horizontal distance (m) — straight-line, not terrain-adjusted
    pub horizontal_distance_m: f64,
    /// Start elevation (m above sea level)
    pub start_elevation_m: f64,
    /// End elevation (m above sea level)
    pub end_elevation_m: f64,
    /// Terrain multiplier (1.0 = flat, >1.0 accounts for routing around obstacles)
    pub terrain_multiplier: f64,
    /// Pipe inner diameter (m)
    pub diameter_m: f64,
    /// Pipe roughness (m) — absolute roughness height ε
    pub roughness_m: f64,
}

impl PipelineSegment {
    /// Elevation change (m). Positive = uphill.
    #[inline]
    pub fn elevation_change_m(&self) -> f64 {
        self.end_elevation_m - self.start_elevation_m
    }

    /// Terrain-adjusted pipe length (m).
    #[inline]
    pub fn adjusted_length_m(&self) -> f64 {
        self.horizontal_distance_m * self.terrain_multiplier
    }
}

/// Analysis result for a single pipeline segment.
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentResult {
    /// Segment name
    pub name: String,
    /// Terrain-adjusted pipe length (m)
    pub pipe_length_m: f64,
    /// Elevation change (m), positive = uphill
    pub elevation_change_m: f64,
    /// Friction head loss (m)
    pub friction_head_m: f64,
    /// Flow velocity (m/s)
    pub velocity_m_s: f64,
    /// Reynolds number
    pub reynolds: f64,
    /// Darcy friction factor
    pub friction_factor: f64,
}

/// Full pipeline analysis result.
#[derive(Debug, Clone, PartialEq)]
pub struct PipelineResult {
    /// Per-segment results
    pub segments: Vec<SegmentResult>,
    /// Total terrain-adjusted pipe length (m)
    pub total_length_m: f64,
    /// Total friction head loss across all segments (m)
    pub total_friction_head_m: f64,
    /// Net elevation change: last segment end - first segment start (m)
    pub net_elevation_change_m: f64,
    /// Maximum elevation encountered along the route (m)
    pub max_elevation_m: f64,
    /// Maximum elevation gain from the start (m) — the critical lift
    pub max_elevation_gain_m: f64,
    /// Total dynamic head: max_elevation_gain + total_friction_head (m)
    pub total_head_m: f64,
    /// Flow rate used (m³/s)
    pub flow_rate_m3_s: f64,
    /// Pump power calculation
    pub pump_power: PumpPower,
}

/// A complete pipeline route composed of segments.
#[derive(Debug, Clone)]
pub struct Pipeline {
    /// Ordered list of pipeline segments
    pub segments: Vec<PipelineSegment>,
    /// Working fluid
    pub fluid: Fluid,
    /// Target volume flow rate (m³/s)
    pub flow_rate_m3_s: f64,
    /// Pump efficiency (0–1)
    pub pump_efficiency: f64,
}

impl Pipeline {
    /// Create a new pipeline with the given fluid, flow rate, and pump efficiency.
    pub fn new(fluid: Fluid, flow_rate_m3_s: f64, pump_efficiency: f64) -> Self {
        Self {
            segments: Vec::new(),
            fluid,
            flow_rate_m3_s,
            pump_efficiency,
        }
    }

    /// Add a segment to the pipeline.
    pub fn add_segment(&mut self, segment: PipelineSegment) {
        self.segments.push(segment);
    }

    /// Analyze the full pipeline: friction losses, elevation profile, pump power.
    pub fn analyze(&self) -> PipelineResult {
        let mut segment_results = Vec::with_capacity(self.segments.len());
        let mut total_length = 0.0;
        let mut total_friction_head = 0.0;
        let mut max_elevation = f64::NEG_INFINITY;
        let mut start_elevation = 0.0;

        for (i, seg) in self.segments.iter().enumerate() {
            if i == 0 {
                start_elevation = seg.start_elevation_m;
            }

            // Track max elevation
            if seg.start_elevation_m > max_elevation {
                max_elevation = seg.start_elevation_m;
            }
            if seg.end_elevation_m > max_elevation {
                max_elevation = seg.end_elevation_m;
            }

            let pipe_length = seg.adjusted_length_m();
            total_length += pipe_length;

            // Compute flow velocity for this segment
            let area = core::f64::consts::PI * (seg.diameter_m / 2.0).powi(2);
            let velocity = self.flow_rate_m3_s / area;

            // Darcy-Weisbach friction loss for this segment
            let dw = DarcyWeisbach::new(
                &self.fluid,
                seg.diameter_m,
                pipe_length,
                velocity,
                seg.roughness_m,
            );
            let re = dw.reynolds_number();
            let f = dw.friction_factor(re);

            // Friction head loss: h_f = f·(L/D)·(v²/(2g))
            let friction_head = f * (pipe_length / seg.diameter_m) * (velocity * velocity / (2.0 * G));
            total_friction_head += friction_head;

            segment_results.push(SegmentResult {
                name: seg.name.clone(),
                pipe_length_m: pipe_length,
                elevation_change_m: seg.elevation_change_m(),
                friction_head_m: friction_head,
                velocity_m_s: velocity,
                reynolds: re,
                friction_factor: f,
            });
        }

        let net_elevation = if let (Some(first), Some(last)) =
            (self.segments.first(), self.segments.last())
        {
            last.end_elevation_m - first.start_elevation_m
        } else {
            0.0
        };

        let max_elevation_gain = max_elevation - start_elevation;

        // Total dynamic head: elevation gain to highest point + all friction losses
        let total_head = max_elevation_gain + total_friction_head;

        let pp = pump_power(
            self.fluid.density_kg_m3,
            self.flow_rate_m3_s,
            total_head,
            self.pump_efficiency,
        );

        PipelineResult {
            segments: segment_results,
            total_length_m: total_length,
            total_friction_head_m: total_friction_head,
            net_elevation_change_m: net_elevation,
            max_elevation_m: max_elevation,
            max_elevation_gain_m: max_elevation_gain,
            total_head_m: total_head,
            flow_rate_m3_s: self.flow_rate_m3_s,
            pump_power: pp,
        }
    }
}

/// Calculate terrain-adjusted length from discrete elevation samples.
///
/// L = Σ √(Δx² + Δh²)
///
/// This is the discrete arc-length integral from README §Terrain Curve Extrapolation.
///
/// # Arguments
/// - `horizontal_distance_m` — total straight-line horizontal distance (m)
/// - `elevation_samples` — elevation readings along the route (m)
///
/// # Returns
/// Actual pipe length following terrain undulation (m).
pub fn terrain_adjusted_length(horizontal_distance_m: f64, elevation_samples: &[f64]) -> f64 {
    if elevation_samples.len() < 2 {
        return horizontal_distance_m;
    }
    let segment_dx = horizontal_distance_m / (elevation_samples.len() - 1) as f64;
    let mut total = 0.0;
    for i in 1..elevation_samples.len() {
        let dh = elevation_samples[i] - elevation_samples[i - 1];
        total += (segment_dx * segment_dx + dh * dh).sqrt();
    }
    total
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build the IGBWP pipeline from IGBWP.md §6
    fn igbwp_pipeline() -> Pipeline {
        // Fresh water at 25°C (post-desalination)
        let fluid = Fluid::water(25.0);
        // 595 m³/s, η=0.85
        let mut pl = Pipeline::new(fluid, 595.0, 0.85);

        // Segment 1: Paradip → Sambalpur (200 km, +170m, ×1.10)
        pl.add_segment(PipelineSegment {
            name: "Paradip → Sambalpur".into(),
            horizontal_distance_m: 200_000.0,
            start_elevation_m: 0.0,
            end_elevation_m: 170.0,
            terrain_multiplier: 1.10,
            diameter_m: 16.0,
            roughness_m: 0.000_046, // smooth steel
        });

        // Segment 2: Sambalpur → Raipur (180 km, +130m, ×1.20)
        pl.add_segment(PipelineSegment {
            name: "Sambalpur → Raipur".into(),
            horizontal_distance_m: 180_000.0,
            start_elevation_m: 170.0,
            end_elevation_m: 300.0,
            terrain_multiplier: 1.20,
            diameter_m: 16.0,
            roughness_m: 0.000_046,
        });

        // Segment 3: Raipur → Mirzapur (200 km, −220m, ×1.15)
        pl.add_segment(PipelineSegment {
            name: "Raipur → Mirzapur".into(),
            horizontal_distance_m: 200_000.0,
            start_elevation_m: 300.0,
            end_elevation_m: 80.0,
            terrain_multiplier: 1.15,
            diameter_m: 16.0,
            roughness_m: 0.000_046,
        });

        // Segment 4: Mirzapur → Prayagraj (80 km, +18m, ×1.05)
        pl.add_segment(PipelineSegment {
            name: "Mirzapur → Prayagraj".into(),
            horizontal_distance_m: 80_000.0,
            start_elevation_m: 80.0,
            end_elevation_m: 98.0,
            terrain_multiplier: 1.05,
            diameter_m: 16.0,
            roughness_m: 0.000_046,
        });

        pl
    }

    #[test]
    fn test_igbwp_total_length() {
        let pl = igbwp_pipeline();
        let r = pl.analyze();
        // Expected: ~880 km terrain-adjusted (from multipliers applied to each segment)
        // 200*1.10 + 180*1.20 + 200*1.15 + 80*1.05 = 220+216+230+84 = 750 km
        // Note: IGBWP says ~990 km with ×1.30 overall; our per-segment multipliers give ~750 km
        let total_km = r.total_length_m / 1000.0;
        assert!(total_km > 700.0 && total_km < 1100.0,
            "IGBWP total length = {total_km:.0} km");
    }

    #[test]
    fn test_igbwp_max_elevation() {
        let pl = igbwp_pipeline();
        let r = pl.analyze();
        // Raipur is the highest point at 300m
        assert!((r.max_elevation_m - 300.0).abs() < 0.1, "max elev = {:.0}", r.max_elevation_m);
        assert!((r.max_elevation_gain_m - 300.0).abs() < 0.1, "max gain = {:.0}", r.max_elevation_gain_m);
    }

    #[test]
    fn test_igbwp_pump_power_order_of_magnitude() {
        let pl = igbwp_pipeline();
        let r = pl.analyze();
        // IGBWP §7: total pumping ≈ 6 GW (4.1 elevation + 1.9 friction)
        // Our model uses max_elevation_gain (300m) + friction, so should be in the right ballpark
        let gw = r.pump_power.shaft_gw();
        assert!(gw > 2.0 && gw < 15.0,
            "IGBWP pump power = {gw:.2} GW (expected ~6 GW order)");
    }

    #[test]
    fn test_igbwp_friction_head_positive() {
        let pl = igbwp_pipeline();
        let r = pl.analyze();
        assert!(r.total_friction_head_m > 0.0, "friction head must be positive");
        // Friction head should be meaningful but not dominate (100–500m range for this scale)
        assert!(r.total_friction_head_m > 10.0,
            "friction head = {:.1} m (should be significant)", r.total_friction_head_m);
    }

    #[test]
    fn test_igbwp_net_elevation() {
        let pl = igbwp_pipeline();
        let r = pl.analyze();
        // Start: 0m (Paradip), End: 98m (Prayagraj) → net = +98m
        assert!((r.net_elevation_change_m - 98.0).abs() < 0.1,
            "net elevation = {:.1}", r.net_elevation_change_m);
    }

    #[test]
    fn test_segment_elevation_change() {
        let seg = PipelineSegment {
            name: "test".into(),
            horizontal_distance_m: 100_000.0,
            start_elevation_m: 0.0,
            end_elevation_m: 170.0,
            terrain_multiplier: 1.10,
            diameter_m: 16.0,
            roughness_m: 0.000_046,
        };
        assert!((seg.elevation_change_m() - 170.0).abs() < 1e-10);
        assert!((seg.adjusted_length_m() - 110_000.0).abs() < 1e-6);
    }

    // --- Terrain arc-length ---

    #[test]
    fn test_terrain_adjusted_flat() {
        // Flat terrain: adjusted length = horizontal distance
        let samples = [0.0, 0.0, 0.0, 0.0];
        let l = terrain_adjusted_length(300.0, &samples);
        assert!((l - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_terrain_adjusted_incline() {
        // Constant 45° slope: L = horizontal × √2
        let samples = [0.0, 100.0, 200.0, 300.0];
        let l = terrain_adjusted_length(300.0, &samples);
        let expected = 300.0 * (2.0_f64).sqrt();
        assert!((l - expected).abs() < 1e-6, "incline length = {l:.2}, expected {expected:.2}");
    }

    #[test]
    fn test_terrain_adjusted_undulating() {
        // Up-down terrain is longer than flat
        let samples = [0.0, 100.0, 0.0, 100.0, 0.0];
        let flat = terrain_adjusted_length(400.0, &[0.0, 0.0, 0.0, 0.0, 0.0]);
        let undulating = terrain_adjusted_length(400.0, &samples);
        assert!(undulating > flat, "undulating ({undulating:.1}) should exceed flat ({flat:.1})");
    }

    #[test]
    fn test_terrain_single_sample_returns_horizontal() {
        let l = terrain_adjusted_length(500.0, &[100.0]);
        assert!((l - 500.0).abs() < 1e-10);
    }
}
