# The Passive Baseline: Why Lunar Impact Cratering is Not Charge-and-Release

**A spatial Random Matrix Theory negative control**

Ruqing Chen, GUT Geoservice Inc., Montreal, Canada · June 2026

---

This is the **ninth** study in a unified Random Matrix Theory (RMT) program, and its
essential **negative control**. The program has shown that single-source,
charge-and-release systems show Wigner–Dyson **level repulsion**. A rigorous program
needs the converse test: a system with **no internal memory** should **not** repel.

Lunar impact craters are that system — an externally driven, memoryless bombardment
on a passive surface. If the RMT detector is honest, lunar craters must fail to show
charge-and-release repulsion. **They do.**

## The key correction: the 2D spatial baseline is 0.57, not 0.386

The classic Poisson value **⟨r⟩ = 0.386** is the **1D temporal** baseline (exponential
waiting times). For **2D nearest-neighbor** spacings, complete spatial randomness (CSR)
already gives **⟨r⟩ ≈ 0.57**, because the nearest-neighbor operation excludes a disk
around each point. **Comparing crater spacings to 0.386 would falsely make any random
2D field look "GOE-repulsive."** This pipeline calibrates the correct baseline by Monte
Carlo and is a reusable result for any future spatial RMT analysis.

| Baseline | ⟨r⟩ | meaning |
|---|---|---|
| 1D temporal Poisson | 0.386 | the temporal arm's null |
| **2D spherical CSR** | **0.570** | **correct null for spatial point fields** |
| hard-core exclusion | 0.776 | genuine spatial repulsion |

## Key results

| Subset | n | ⟨r⟩ | Clark-Evans R | KS vs CSR | State |
|---|---|---|---|---|---|
| All craters ≥20 km (global) | 6973 | 0.695 | — | — | mild passive regularity |
| CSR null (same N) | 6973 | 0.572 | ≈1.0 | — | random (reference) |
| **Mare Smythii (old edge) 5–20 km** | 239 | 0.669 | **0.94** | **p=0.12** | **clean CSR ✓** |
| **Mare Australe 8–30 km** | 245 | 0.652 | **0.92** | **p=0.23** | **clean CSR ✓** |
| Oceanus Procellarum (young) 5–20 km | 399 | 0.576 | **0.79 (z=−8)** | p<0.001 | **clustered** |
| Mare Nubium 5–20 km | 156 | 0.599 | 0.86 | p=0.004 | clustered |

### What the numbers say

1. **Global large craters (0.695) sit only mildly above CSR (0.572).** Diagnostics
   reject saturation (no diameter-cut trend) and hard-core exclusion (**52% of nearest
   neighbors physically overlap** — impossible under hard-core). The excess is passive
   geometry (finite crater size + terrain-age density gradients), not memory.

2. **Old limb maria are clean CSR.** Mare Smythii and Mare Australe match complete
   spatial randomness in full distribution shape (Clark-Evans R≈0.93, KS p>0.1).
   Memoryless impacts on a quiet surface → spatial randomness, exactly as predicted.

3. **Young nearside maria are clustered** by **secondary craters** (ejecta chains from
   Copernicus, Kepler, etc.) — Clark-Evans R≈0.79, z=−8. This is correlated triggering,
   not charge-and-release memory. (Note: some clustered subsets have ⟨r⟩≈0.57 by
   coincidence — clustering and finite-size regularity cancel in ⟨r⟩ alone, which is
   why we use Clark-Evans and KS, **independent of ⟨r⟩**.)

### The three-state spectrum

The lunar result completes a spatial spectrum, all separated by one RMT diagnostic:

- **Clustering** (⟨r⟩ depressed, R<1): young maria — secondary-crater correlation, the
  spatial analog of sympathetic solar-flare triggering.
- **Randomness** (⟨r⟩ ≈ CSR): old limb maria — memoryless impacts.
- **Repulsion** (⟨r⟩ > CSR → hard-core): the Earth/orbital charge-and-release systems —
  reservoir depletion / gravitational clearing.

**No charge-and-release mechanism, no level repulsion.** The negative control confirms
the program's core claim from the opposite side.

## Reproduce

```bash
pip install -r requirements.txt
# Global large craters + CSR null + saturation check:
python code/lunar_crater_rmt_pipeline.py data/lunar_craters_ge20km.csv 20
```

## Data

`data/lunar_craters_ge20km.csv` — the 6,973 craters ≥20 km extracted from the full
Robbins (2019) database (1.3 million craters). Full database:
https://astrogeology.usgs.gov/search/map/moon_crater_database_v1_robbins

Columns: lat, lon, diam (km). All distances are great-circle (haversine) on the sphere.

## Repository Structure

```
lunar-crater-rmt-baseline/
├── README.md
├── LICENSE                       # MIT
├── requirements.txt
├── .zenodo.json
├── CITATION.cff
├── paper/
│   ├── paper.tex
│   ├── paper.pdf
│   └── figs/
│       ├── fig1_spectrum.pdf
│       └── fig2_two_anchors.pdf
├── code/
│   └── lunar_crater_rmt_pipeline.py
├── data/
│   └── lunar_craters_ge20km.csv  # 6,973 craters >=20 km (Robbins 2019)
├── figures/
└── results/
    └── lunar_rmt_results.json
```

## The unified program (9 systems)

| # | Domain | State | DOI |
|---|---|---|---|
| 1 | Stratigraphy | GOE | [20774581](https://zenodo.org/records/20774581) |
| 2 | Seismotectonics | scale-dependent | [20768130](https://zenodo.org/records/20768130) |
| 3 | Mantle plumes | GOE | [20768420](https://zenodo.org/records/20768420) |
| 4 | Metallogeny | GOE/GUE | [20768849](https://zenodo.org/records/20768849) |
| 5 | Evolution | GOE | [20783763](https://zenodo.org/records/20783763) |
| 6 | Hydrogeology | GOE/super-GUE | [20780389](https://zenodo.org/records/20780389) |
| 7 | Solar flares | Poisson (depletion limit) | [20784967](https://zenodo.org/records/20784967) |
| 8 | Orbital architecture | GUE (local-frame) | [20785613](https://zenodo.org/records/20785613) |
| **9** | **Lunar craters** | **CSR / clustering (negative control)** | this work |

## License

MIT (code) · Data courtesy Robbins (2019) / USGS Astropedia.
