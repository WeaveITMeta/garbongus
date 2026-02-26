//! # garbongus
//!
//! Fluid mechanics library for pipeline hydraulics, desalination water transport,
//! tunnel flow, suction/vacuum lift, capillary rise, and pipe flow calculations.
//! Supports fresh water and seawater/brine. Pure Rust, zero dependencies.
//!
//! ## Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`fluid`] | Fluid properties: fresh water + seawater/brine (density, viscosity, surface tension, vapor pressure) |
//! | [`flow`] | Flow rate, pipe sizing, unit conversions, Bernoulli, pump power, pressure↔depth |
//! | [`manning`] | Manning equation for tunnel and open-channel gravity flow |
//! | [`pipeline`] | Multi-segment pipeline with elevation profile, friction losses, pump power |
//! | [`pipe`] | Darcy-Weisbach friction loss, Colebrook-White, Reynolds number, flow regime |
//! | [`vacuum`] | Suction/vacuum lift over arbitrary distances |
//! | [`capillary`] | Capillary rise height calculations |
//!
//! ## Quick Example — Seawater Pipeline
//!
//! ```rust
//! use garbongus::fluid::Fluid;
//! use garbongus::flow::{required_diameter, pump_power, FlowRate};
//!
//! // Bay of Bengal feedwater: 32 ppt salinity, 25°C
//! let fluid = Fluid::seawater(25.0, 32.0);
//! println!("Seawater density: {:.1} kg/m³", fluid.density_kg_m3);
//!
//! // Size a pipe for 595 m³/s at 3 m/s
//! let d = required_diameter(595.0, 3.0);
//! println!("Required diameter: {:.1} m", d);
//!
//! // Pump power for 300m elevation lift
//! let pp = pump_power(fluid.density_kg_m3, 595.0, 300.0, 0.85);
//! println!("Pump power: {:.2} GW", pp.shaft_gw());
//!
//! // Flow rate conversions
//! let fr = FlowRate::from_m3s(2.8);
//! println!("Tunnel flow: {:.0} MGD", fr.to_mgd());
//! ```

pub mod capillary;
pub mod fluid;
pub mod flow;
pub mod manning;
pub mod pipe;
pub mod pipeline;
pub mod vacuum;
