# Garbongus × Eustress WATER — Full Capabilities Audit

## Reference
```
E:\Workspace\EustressEngine\WATER\
├── docs/
│   ├── README.md          — Master physics & 100-year flow model (California example)
│   ├── IGBWP.md           — Indo-Gangetic Basin Water Project (990 km pipeline, 595 m³/s)
│   ├── TUNNEL.md          — Boring Company Tunnel Vision Challenge (U.S. + Prayagraj)
│   ├── WELL_METER.md      — 30M IoT well metering system (ESP32-C3, LoRaWAN, Rust)
│   └── GRIEVANCE_PM.md    — CPGRAMS petition to India PM on IGB aquifer collapse
├── assets/geo/
│   ├── geo.toml           — Eustress terrain config (EPSG:32644, Prayagraj origin)
│   └── vectors/
│       ├── desal_plants.geojson    — 3 nodes: Paradip, Dhamra, Haldia
│       ├── pipeline_route.geojson  — 990 km main trunk + Kanpur/Lucknow/Varanasi branches
│       └── pump_stations.geojson   — Sambalpur, Raipur (main lift), Mirzapur, Prayagraj
└── *.pdf                  — Print versions of above
```

## Objective
Determine exactly what `garbongus v0.1.0` can and cannot do for the WATER project,
what must change, and in what order.

---

## 1. WATER Project Scope (from docs)

The WATER project has **four interlocking workstreams**:

### 1.1 IGBWP — 990 km Desalination Pipeline (IGBWP.md)
- **Problem**: IGB aquifer (#1 most depleted on Earth) losing ~15 km³/year
- **Solution**: 3 desal plants on Bay of Bengal → 990 km pipeline → Prayagraj terminus → aquifer recharge
- **Scale**: 595 m³/s total flow, ~16 m diameter trunk, 9 GW continuous power
- **Segments**: Paradip(0m) → Sambalpur(170m) → Raipur(300m, main lift) → Mirzapur(80m) → Prayagraj(98m)
- **Branches**: Kanpur, Lucknow, Varanasi from Prayagraj hub
- **Terrain multiplier**: 1.30 (mostly flat Gangetic plain, one plateau crossing)
- **Energy**: 3.0 GW desal + 4.1 GW elevation pumping + 1.9 GW friction = 9 GW
- **Cost**: $70–120B over 20–30 years (comparable to China South-North Water Transfer)
- **Feedwater**: Bay of Bengal ~32 g/L salinity (lower than open ocean 35 g/L)

### 1.2 Boring Company Tunnels (TUNNEL.md + IGBWP.md §9)
- **U.S. sites**: Tucson CAP-to-Recharge (#1), Las Vegas Intake (#2), San Diego Desal (#3)
- **India site**: Prayagraj — 800m tunnel under Yamuna, pure alluvium, bypasses barrage allocation
- **Hydraulics**: Manning equation for 3.66m (12 ft) diameter, gravity or pumped flow
- **Scale**: 1 tunnel = 2.8 m³/s (64 MGD); 213 tunnels = full 595 m³/s equilibrium
- **Rust data structures**: `WaterTunnelProposal`, `TunnelHydraulics`, `EvaluationScore` etc.

### 1.3 WELL_METER — 30M IoT Well Sensors (WELL_METER.md)
- **Purpose**: Close the data gap — true pumping rate unknown (5–20× spread in estimates)
- **Hardware**: ESP32-C3 RISC-V, $11.45 BOM, 10-year battery, LoRaWAN, tamper-evident
- **Sensors**: Pressure transducer (d = P/ρg, ±5cm) + paddle-wheel flow (±2% above 0.3 m/s)
- **Physics needed**: `P = ρ·g·d` water column, flow rate Q = k·f·A
- **Backend**: Axum + TimescaleDB + PostGIS, MQTT ingest from ChirpStack
- **Scale**: Phase 0: 10K wells ($230K) → Phase 1: 5M wells ($94M) → Phase 2: 30M wells ($508M)

### 1.4 Policy & Governance (GRIEVANCE_PM.md)
- **Actions**: National tube well metering mandate, independent overdraft audit, minimum environmental flow enforcement
- **Dependency**: No infrastructure commitment beyond Phase 1 is rational until true overdraft known to ±10%

---

## 2. Garbongus v0.1.0 — Current Capabilities

| Module | What It Does | WATER Doc Reference | Coverage |
|--------|-------------|--------------------:|----------|
| `fluid` | Water ρ, μ, σ, P_vapor (0–100°C, 1 atm) | README §Pipe Flow, IGBWP §7, WELL_METER §5 | **Partial** — fresh water only |
| `pipe` | Darcy-Weisbach ΔP, Colebrook-White f, Re, flow regime | README §Pipe Flow (`f·(L/D)·ρv²/2`) | **Core match** — single segment |
| `vacuum` | Suction lift, pump pressure, cavitation check | README §Energy (pumping head) | **Partial** — single-column only |
| `capillary` | Jurin's Law capillary rise | Not referenced | **Irrelevant** to WATER |

### What Works Today for WATER
- ✅ Darcy-Weisbach friction loss — exact formula used in README.md and IGBWP.md
- ✅ Reynolds number, flow regime classification
- ✅ Water density/viscosity at temperature — used in WELL_METER pressure transducer conversion
- ✅ Vapor pressure — cavitation check at pump stations
- ✅ Constants: `G = 9.80665`, `P_ATM = 101_325` — match WELL_METER appendix values

### What Does NOT Work
- ❌ No seawater/brine properties (Bay of Bengal 32 g/L)
- ❌ No multi-segment pipeline (IGBWP has 5 segments with different elevations)
- ❌ No pump station modeling (IGBWP needs 4 stations, pump curves, power calc)
- ❌ No elevation profile / hydraulic grade line
- ❌ No Manning equation (TUNNEL.md uses Manning for open-channel/gravity tunnels)
- ❌ No Bernoulli equation (README §Pipe Flow uses P₁+½ρv₁²+ρgh₁ = P₂+½ρv₂²+ρgh₂)
- ❌ No flow rate / volume conversions (Q=Av, MGD, MLD, acre-feet)
- ❌ No terrain arc-length integration (README §Terrain Curve Extrapolation)
- ❌ No RO desalination energy model (README §Zero-Brine-Discharge 5-stage RO)
- ❌ No well pressure→depth conversion (WELL_METER §5: d = P/(ρg))
- ❌ No pipe diameter sizing from target flow (README: D = 2√(Q/v/π))

---

## 3. Gap Analysis — Mapped to Specific WATER Requirements

### 🔴 P0 — Must Have (blocks IGBWP pipeline sizing)

| # | Gap | WATER Source | Formula/Spec |
|---|-----|-------------|-------------|
| 1 | **Seawater properties** | IGBWP §5, §7 | ρ(T,S), μ(T,S) via UNESCO EOS-80; Bay of Bengal S=32 g/L |
| 2 | **Multi-segment pipeline** | IGBWP §6: 5 segments with elevation | Chain: Paradip→Sambalpur→Raipur→Mirzapur→Prayagraj |
| 3 | **Elevation head** | IGBWP §7: `ρgQH/η` | Total head = elevation + friction; per-station allocation |
| 4 | **Pump power calc** | IGBWP §7: `P = ρgQH/η` | η=0.85, 4.1 GW elevation + 1.9 GW friction = 6 GW pumping |
| 5 | **Manning equation** | TUNNEL.md §3: `Q = (1/n)·A·R^(2/3)·S^(1/2)` | n=0.012 concrete, 3.66m diameter, gravity flow |
| 6 | **Flow rate & volume conversions** | README §100-Year: Q=Av, MGD, AF/yr | `Q_required ≈ 595 m³/s`; unit conversions throughout |
| 7 | **Pipe diameter sizing** | README §Pipe Specs: `D = 2√(A/π)` | Given Q and v, solve for D |

### 🟡 P1 — Should Have (blocks Phase 2+ validation)

| # | Gap | WATER Source | Formula/Spec |
|---|-----|-------------|-------------|
| 8 | **Bernoulli equation** | README §Pipe Flow | `P₁+½ρv₁²+ρgh₁ = P₂+½ρv₂²+ρgh₂` |
| 9 | **Terrain arc-length** | README §Terrain Curve | `L = ∫√(1+(dh/dx)²)dx`, discrete elevation samples |
| 10 | **RO desalination energy** | README §Zero-Brine, IGBWP §7 | 3.0–3.5 kWh/m³ at 32 g/L; osmotic pressure model |
| 11 | **Well pressure→depth** | WELL_METER §5 | `d = P/(ρ·g)` — trivial but needed for firmware validation |
| 12 | **Hydraulic grade line** | IGBWP §7 elevation profile | HGL/EGL plotted over 990 km, detect negative pressure zones |
| 13 | **Pump curves** | IGBWP pump stations | `H(Q) = H_shutoff - a·Q²`, `η(Q)`, system curve intersection |
| 14 | **Water hammer / transient** | 990 km pipeline, valve closures | Joukowsky: `ΔP = ρ·c·Δv`, wave speed in steel/concrete pipe |

### 🟢 P2 — Nice to Have (enables Eustress 3D visualization)

| # | Gap | WATER Source | Formula/Spec |
|---|-----|-------------|-------------|
| 15 | **GeoJSON pipeline ingestion** | `pipeline_route.geojson` | Parse LineString coordinates → pipeline segments |
| 16 | **Bevy integration** | Eustress Engine, `geo.toml` | Pipeline as 3D mesh, pump stations as markers |
| 17 | **Network solver** | IGBWP branches (Kanpur, Lucknow, Varanasi) | Hardy-Cross or gradient method for branched network |
| 18 | **Cost optimization** | IGBWP §11: $70–120B | CAPEX (pipe diameter) vs OPEX (pump energy) tradeoff |
| 19 | **Aquifer recharge model** | IGBWP §8: percolation basins | Recharge rate vs injection rate, aquifer response |

---

## 4. Recommended Implementation Roadmap

### Phase 1: Pipeline Core (v0.2.0) — enables IGBWP §6-§7 calculations

```
New modules:
  src/seawater.rs    — Fluid::seawater(temp_c, salinity_ppt) via UNESCO EOS-80
  src/pipeline.rs    — PipelineSegment chain with elevation, auto-calculates:
                        - Per-segment friction loss (existing Darcy-Weisbach)
                        - Cumulative elevation head
                        - Total pump power requirement (P = ρgQH/η)
  src/manning.rs     — Manning equation for tunnel/open-channel flow
                        Q = (1/n)·A·R^(2/3)·S^(1/2)
  src/flow.rs        — Flow rate helpers:
                        - Q = A·v, A = π·(D/2)², D = 2√(Q/(v·π))
                        - Unit conversions: m³/s ↔ MGD ↔ MLD ↔ AF/yr ↔ L/min
                        - Bernoulli: P₁+½ρv₁²+ρgh₁ = P₂+½ρv₂²+ρgh₂

Changes to existing:
  fluid.rs           — Add Fluid::seawater(temp_c, salinity_ppt)
  pipe.rs            — Add flow_rate_m3s() convenience method
  lib.rs             — Export new modules

Validates against:
  - IGBWP §7 energy budget (3.0 + 4.1 + 1.9 = 9 GW)
  - TUNNEL.md §3 flow scenarios (12-ft tunnel, various head differences)
  - README §100-Year flow rate (595 m³/s for IGB, 1172 m³/s for California)
  - README §Pipe Specs diameter sizing table
```

### Phase 2: Terrain & Pumps (v0.3.0) — enables full IGBWP route analysis

```
New modules:
  src/terrain.rs     — Terrain arc-length integration from elevation samples
                        L = Σ √(Δx² + Δh²), terrain multiplier estimation
  src/pump.rs        — Pump station model:
                        H(Q) curve, η(Q) curve, power = ρgQH/η
                        System curve intersection (operating point)
  src/hgl.rs         — Hydraulic/Energy grade line over pipeline distance
                        Detect negative pressure zones, cavitation risk points
  src/desal.rs       — RO desalination energy model:
                        Osmotic pressure = f(salinity), specific energy (kWh/m³)
                        5-stage cascade per README §Zero-Brine-Discharge

Validates against:
  - IGBWP §6 terrain segments (4 segments, multipliers 1.05–1.20)
  - IGBWP §7 elevation profile (0→170→300→80→98 m)
  - IGBWP §7 energy budget per component
  - pump_stations.geojson (4 stations with elevations)
```

### Phase 3: WELL_METER Physics & Eustress Integration (v0.4.0)

```
New modules:
  src/well.rs        — Well pressure→depth (d = P/ρg), drawdown models
                        Validates WELL_METER §5 sensor physics
  src/transient.rs   — Water hammer: ΔP = ρ·c·Δv, wave speed
  src/geojson.rs     — Parse pipeline_route.geojson → Pipeline segments
                        Parse desal_plants.geojson → node capacities
                        Parse pump_stations.geojson → station elevations

Optional (feature-gated):
  src/bevy.rs        — Bevy Component wrappers per geo.toml spec
  [features]
  bevy = ["dep:bevy"]
```

---

## 5. Specific Validation Targets

These are concrete numbers from the WATER docs that garbongus must reproduce:

### From IGBWP.md

| Calculation | Expected Result | Inputs |
|------------|----------------|--------|
| Required flow rate | ~595 m³/s | 15 km³/yr overdraft, 80% recharge efficiency |
| Pipe diameter | ~15.9 m (or parallel 8m) | Q=595 m³/s, v=3 m/s |
| Elevation pumping power | ~4.1 GW | ρ=1000, g=9.81, Q=595, H=300m, η=0.85 |
| Friction pumping power | ~1.9 GW | f=0.015, L=990km, D=16m, v=3 m/s |
| Desal power | ~3.0 GW | 3.5 kWh/m³ × 595 m³/s × 3600/1000 |
| Total power | ~9 GW | Sum of above |

### From TUNNEL.md

| Calculation | Expected Result | Inputs |
|------------|----------------|--------|
| Tunnel flow (gentle gravity) | 12.6 m³/s, 288 MGD | D=3.66m, n=0.012, S=0.3%, Manning |
| Tunnel flow (pumped high) | 41.0 m³/s, 936 MGD | D=3.66m, S=3.1%, Manning |
| Tucson tunnel velocity | 2.7 m/s | Q=2.8 m³/s, D=3.66m |
| Tunnels for full IGB | 213 | 595 / 2.8 |

### From WELL_METER.md

| Calculation | Expected Result | Inputs |
|------------|----------------|--------|
| Pressure→depth | 5.0 m | P=49050 Pa, ρ=1000, g=9.81 |
| Flow rate DN40 | 2.51 L/s = 151 L/min | A=π×0.02², v=2 m/s |

### From README.md

| Calculation | Expected Result | Inputs |
|------------|----------------|--------|
| California flow rate | ~1,172 m³/s | V₀=3.7e12 m³, R=1%/yr |
| Pipe diameter (California) | ~22.3 m | Q=1172, v=3 m/s |
| Pumping power (California) | ~6.76 GW | Q=1172, H=500m, η=0.85 |

---

## 6. Verdict

| Question | Answer |
|----------|--------|
| Can garbongus aid the WATER project? | **Yes — Darcy-Weisbach, Colebrook-White, fluid properties are the right physics foundation** |
| Is it feasible today at v0.1.0? | **No — 7 critical gaps block pipeline sizing (P0 list above)** |
| What's the #1 missing piece? | **`pipeline.rs` — multi-segment elevation-aware pipeline with pump power calculation** |
| Is the architecture sound for extension? | **Yes — modular, zero-dep, clean separation; new modules plug in directly** |
| Estimated effort to Phase 1? | **1 session — 5 new modules, ~1500 lines** |
| Estimated effort to full WATER coverage? | **3 sessions across Phases 1–3** |
| Blocking issues? | **None — pure computation, no external deps, all formulas documented in WATER docs** |

### Bottom Line

`garbongus` has the **correct physics core** that WATER already uses (Darcy-Weisbach appears
verbatim in README.md §Pipe Flow). The crate needs **7 new modules across 3 phases** to cover
the full WATER scope — from 990 km IGBWP pipeline sizing, to Boring Company tunnel hydraulics
(Manning), to WELL_METER sensor physics (P/ρg). Every formula is explicitly documented in the
WATER `.md` files with expected numerical results — making validation straightforward.

The single highest-value addition is `pipeline.rs` — it unlocks the IGBWP §7 energy budget
calculation (the 9 GW number that determines whether this is a $70B or $120B project).
