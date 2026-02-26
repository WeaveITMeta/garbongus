# garbongus

[![crates.io](https://img.shields.io/crates/v/garbongus.svg)](https://crates.io/crates/garbongus)
[![docs.rs](https://docs.rs/garbongus/badge.svg)](https://docs.rs/garbongus)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Fluid mechanics library for suction/vacuum lift, capillary rise, and pipe flow — pure Rust, zero dependencies, `no_std` compatible.**

---

## What is garbongus?

`garbongus` computes the physics of moving fluids up pipes against gravity. The core question it answers:

> *"I need to pull water up a pipe to height H. What pressure do I need, and is that physically achievable?"*

The key insight: **atmosphere alone can only lift water ~10.3 m** (the "suction lift limit"). Beyond that, you need a pump. `garbongus` calculates the exact pump pressure required for **any distance** — millimeters to kilometers.

---

## Quick Start

```toml
[dependencies]
garbongus = "0.1"
```

```rust
use garbongus::{fluid::Fluid, vacuum::VacuumLift};

fn main() {
    let fluid = Fluid::water(20.0); // water at 20°C

    // How much pump pressure to lift water 50 meters?
    let lift = VacuumLift::new(fluid, 0.05, 50.0); // 50mm radius pipe, 50m height
    let result = lift.calculate();

    println!("Atmospheric max lift:  {:.2} m", result.atmospheric_max_lift_m);
    println!("Achievable by atm:     {}", result.achievable_by_atmosphere);
    println!("Required pump pressure:{:.0} Pa ({:.3} bar)",
        result.required_pump_pressure_pa,
        result.required_pump_pressure_bar()
    );
}
```

Output:
```
Atmospheric max lift:  10.18 m
Achievable by atm:     false
Required pump pressure: 391,650 Pa (3.917 bar)
```

---

## Physics Background

### The Suction Lift Limit

Atmospheric pressure (~101,325 Pa at sea level) acts on the water surface in the source reservoir. This pressure pushes water up a pipe when the top is evacuated. The maximum height this can achieve is:

```
h_max = (P_atm - P_vapor) / (ρ · g)
     = (101325 - 2337) / (998.2 × 9.80665)
     ≈ 10.18 m  (at 20°C)
```

The **vapor pressure** subtracts because if the suction pressure drops below it, the water boils (cavitates) and the column breaks.

### Beyond the Limit: Vacuum Pump

A vacuum pump applies additional pressure differential, allowing lift beyond atmospheric:

```
h_total = h_atm_max + P_pump / (ρ · g)
```

For **any target height H**, the required pump gauge pressure is:

```
P_pump = ρ · g · H + ΔP_friction + P_vapor - P_atm
```

There is **no theoretical upper limit** — with a large enough pump, you can lift water any vertical distance.

### Effect of Temperature

Higher temperature → higher vapor pressure → less available suction:

| Temperature | P_vapor | Max Lift Height |
|-------------|---------|-----------------|
| 0°C  | 611 Pa    | ~10.32 m |
| 20°C | 2,337 Pa  | ~10.18 m |
| 50°C | 12,335 Pa | ~9.28 m  |
| 80°C | 47,390 Pa | ~5.62 m  |
| 100°C| 101,325 Pa| ~0.00 m  |

At 100°C, water is already boiling — no suction lift is possible.

### Altitude Effect

Lower altitude → lower atmospheric pressure → less available suction:

| Altitude | P_atm    | Max Lift (20°C) |
|----------|----------|-----------------|
| 0 m (sea level) | 101,325 Pa | ~10.18 m |
| 1,000 m         | 89,875 Pa  | ~9.01 m  |
| 2,500 m         | 74,682 Pa  | ~7.38 m  |
| 5,000 m (Everest base) | 54,048 Pa | ~5.28 m  |

---

## Modules

### `fluid` — Thermophysical Properties

Water properties as a function of temperature (0–100°C):

```rust
use garbongus::fluid::Fluid;

let f = Fluid::water(20.0);
println!("Density:          {:.2} kg/m³", f.density_kg_m3);
println!("Dynamic viscosity:{:.4e} Pa·s", f.dynamic_viscosity_pa_s);
println!("Surface tension:  {:.5} N/m", f.surface_tension_n_m);
println!("Vapor pressure:   {:.0} Pa", f.vapor_pressure_pa);
println!("Kin. viscosity:   {:.4e} m²/s", f.kinematic_viscosity_m2_s());
```

Custom fluids are also supported:

```rust
let mercury = Fluid::custom(
    20.0,       // temp_c
    13_534.0,   // density kg/m³
    1.526e-3,   // dynamic viscosity Pa·s
    0.4865,     // surface tension N/m
    0.16,       // vapor pressure Pa
);
```

### `capillary` — Capillary Rise (Jurin's Law)

```rust
use garbongus::{fluid::Fluid, capillary::CapillaryRise};

let fluid = Fluid::water(20.0);
let calc = CapillaryRise::new(fluid, 0.001, 0.0); // 1mm radius, 0° contact angle
let result = calc.calculate();

println!("Rise height:         {:.2} mm", result.height_m * 1000.0);
println!("Capillary pressure:  {:.1} Pa", result.capillary_pressure_pa);
```

Formula: `h = 2γ·cos(θ) / (ρ·g·r)`

### `vacuum` — Suction/Vacuum Lift

```rust
use garbongus::{fluid::Fluid, vacuum::VacuumLift};

let fluid = Fluid::water(20.0);
let result = VacuumLift::new(fluid, 0.05, 1000.0) // 50mm pipe, 1 km height
    .flow_velocity(1.0)      // 1 m/s flow with friction
    .roughness(1.5e-6)       // commercial steel
    .ambient_pressure(89_875.0) // 1000m altitude
    .calculate();

println!("Target height:        {:.0} m",   result.target_height_m);
println!("Atm max lift:         {:.2} m",   result.atmospheric_max_lift_m);
println!("Required pump:        {:.0} Pa",  result.required_pump_pressure_pa);
println!("Required pump:        {:.3} bar", result.required_pump_pressure_bar());
println!("Required pump:        {:.2} PSI", result.required_pump_pressure_psi());
println!("Pump head:            {:.2} m",   result.required_pump_head_m(998.2));
```

### `pipe` — Darcy-Weisbach Pipe Flow

```rust
use garbongus::{fluid::Fluid, pipe::DarcyWeisbach};

let fluid = Fluid::water(20.0);
let dw = DarcyWeisbach::new(
    &fluid,
    0.05,    // 50mm diameter
    100.0,   // 100m pipe length
    1.0,     // 1 m/s velocity
    1.5e-6,  // commercial steel roughness
);
let result = dw.calculate();

println!("Reynolds number: {:.0}", result.reynolds_number);
println!("Flow regime:     {:?}", result.flow_regime);
println!("Friction factor: {:.4}", result.friction_factor);
println!("Pressure loss:   {:.1} Pa", result.pressure_loss_pa);
println!("Head loss:       {:.3} m", result.head_loss_m);
```

---

## Benchmarks

All calculations are O(1) — constant time regardless of lift height:

```
vacuum::lift_5m_static      ~45 ns
vacuum::lift_1000m_static   ~45 ns   ← same cost at any distance
vacuum::lift_50m_with_flow  ~380 ns  (includes Colebrook-White iteration)
capillary::jurin_1mm_pipe   ~12 ns
pipe::darcy_weisbach_laminar ~25 ns
pipe::darcy_weisbach_turbulent ~280 ns (Colebrook-White converges in ~4 iter)
fluid::water_properties_20c ~18 ns
```

Run benchmarks:
```bash
cargo bench
# HTML report at: target/criterion/report/index.html
```

---

## Feature Flags

| Flag | Effect | Default |
|------|--------|---------|
| `std` | Use `f64::exp`, `f64::cos`, etc. from std | ✓ |

Without `std`, all math uses built-in Taylor series / iterative implementations (sufficient for engineering accuracy).

---

## License

MIT
