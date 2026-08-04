# ElectronFlow — Interactive Redox Tower

Interactive companion to Scarampi, Cremer, Soyer et al., *Community Metabolism in Dynamic Redox Landscapes* (under revision).

**Live site:** https://cremerlab.github.io/redox-tower

---

## What the site does

The site shows an **interactive redox tower**: a chart of biologically relevant reduction half-reactions ranked by their reduction potential (mV). Users can adjust pH, gas partial pressures, and metabolite concentrations via sliders. Each half-reaction potential shifts in real time according to the Nernst equation. Pairs of half-reactions (microbial metabolisms) can be selected to display whether they are thermodynamically feasible (ΔG < 0) under the chosen conditions.

The app runs entirely in the browser — no server required. Python code is executed client-side via [Pyodide](https://pyodide.org) and [Panel](https://panel.holoviz.org).

---

## Repository layout

```
app.py                  Main Panel application
app_old.py              Preserved pre-refactor version (for comparison / validation)
environment.yml         Conda environment with all dependencies

data/                   Editable input tables (source of truth for scientific content)
  species.csv           Slider-controlled species: display names, ranges, defaults
  half_rxns.csv         Half-reactions: thermodynamic parameters + stoichiometry
  reactions.csv         Redox pairs (metabolisms): donor/acceptor links + equations
  presets.csv           Named environmental conditions shown in the condition selector

docs/                   Everything served by GitHub Pages
  index.html            Home page
  app.html              Interactive redox tower (generated — do not edit manually)
  app.js                Panel/Pyodide worker bundle (generated — do not edit manually)
  calc.html             Calculation details page
  about.html            About page
  feedback.html         Feedback page
  data/                 Copy of data/ CSVs — fetched by the browser at runtime
    species.csv
    half_rxns.csv
    reactions.csv
    presets.csv
```

---

## Input data

All scientific content lives in CSV files under `data/`. These are the only files you need to edit to add or modify species, half-reactions, metabolisms, or environmental presets.

### `data/species.csv` — all species

One row per species that appears anywhere in a half-reaction. The `slider` column decides whether the user sees a slider for that species or whether it is held at its default value silently.

| Column | Type | Description |
|--------|------|-------------|
| `key` | string | Variable name used in `half_rxns.csv` and throughout the app (e.g. `log_pO2`) |
| `label` | string | Full slider label shown in the UI |
| `disp` | string | Short name used in the conditions summary table |
| `kind` | string | `pH`, `gas` (log₁₀ bar), or `conc` (log₁₀ M) — controls value formatting |
| `start` | float | Slider minimum |
| `end` | float | Slider maximum |
| `step` | float | Slider step size |
| `default` | float | Default value; for fixed species this is the only value ever used |
| `slider` | `true`/`false` | Whether to show a slider in the UI |

The default for all species is 0 (meaning 1 M or 1 bar), which is the biochemical standard state. This ensures that at pH 7 with all sliders at their defaults, each half-reaction gives exactly its tabulated E°′.

To promote a currently fixed species to a slider, change its `slider` value from `false` to `true` — no other changes needed.

Example rows:
```
# Slider species — user can vary this
log_pH2, H₂ partial pressure (log₁₀ bar), H₂, gas, -7.0, 0.0, 0.25, 0.0, true

# Fixed species — always at standard state (1 M), no slider shown
log_cHS, HS⁻ concentration (log₁₀ M), HS⁻, conc, -6.0, 0.0, 0.25, 0.0, false
```

### `data/half_rxns.csv` — half-reaction table

One row per half-reaction displayed on the tower. The Nernst equation is applied generically using the fields in this table.

| Column | Type | Description |
|--------|------|-------------|
| `id` | string | Unique identifier, referenced by reactions.csv |
| `label` | string | Short label shown on the tower |
| `color` | hex string | Bar and label colour |
| `col` | int 0–3 | Tower column (0=Aerobic, 1=Anaerobic resp., 2=Fermentation, 3=Syntrophic/C1) |
| `E_SHE` | float (V) | Standard reduction potential vs SHE at pH 0, all activities = 1 |
| `E_std` | float (mV) | Biochemical standard E°′ at pH 7, unit activities — used for dashed reference lines |
| `ne` | int | Electrons transferred in the reduction direction |
| `nH` | int | Net H⁺ consumed in the reduction direction (positive = H⁺ on left side) |
| `reactants` | JSON dict | `{"species_key": stoichiometric_power}` for species on the left side of the reduction |
| `products` | JSON dict | `{"species_key": stoichiometric_power}` for species on the right side of the reduction |

H⁺ and H₂O are **not** listed here. H⁺ is handled via `nH`; H₂O has activity = 1 by convention and cancels out.

#### The two JSON fields

Both `reactants` and `products` are JSON objects. Each key is a species name from `species.csv`; the value is its stoichiometric coefficient.

- `{}` — no species on that side (other than H⁺ and H₂O)
- `{"log_pO2": 1}` — one O₂ molecule (coefficient 1)
- `{"log_pCO2": 6}` — six CO₂ molecules (coefficient 6)
- `{"log_pCO2": 1, "log_cAc": 1}` — one CO₂ and one acetate

It does not matter whether the species has a slider or not — the notation is the same either way.

#### Three worked examples

**Example 1 — O₂/H₂O** (O₂ has a slider; H₂O not listed)

```
Reaction:   O₂ + 4H⁺ + 4e⁻  →  2H₂O
ne=4, nH=4
reactants = {"log_pO2": 1}   ← O₂ on the left, coefficient 1
products  = {}                ← H₂O not listed (activity = 1 by convention)

log₁₀(Q) = 4·pH − 1·log_pO2
```

**Example 2 — CO₂/Glucose** (both CO₂ and glucose have sliders)

```
Reaction:   6CO₂ + 24H⁺ + 24e⁻  →  Glucose + 6H₂O
ne=24, nH=24
reactants = {"log_pCO2": 6}   ← 6 CO₂ on the left, coefficient 6
products  = {"log_cGlc": 1}   ← Glucose on the right, coefficient 1

log₁₀(Q) = 24·pH  +  1·log_cGlc  −  6·log_pCO2
```

At pH 7 with all species at standard state (log₁₀ = 0): log₁₀(Q) = 168, giving E = −430 mV = E°′ ✓

**Example 3 — SO₄²⁻/HS⁻** (SO₄²⁻ has a slider; HS⁻ is fixed — but the notation is identical)

```
Reaction:   SO₄²⁻ + 9H⁺ + 8e⁻  →  HS⁻ + 4H₂O
ne=8, nH=9
reactants = {"log_cSO4": 1}   ← sulfate on the left
products  = {"log_cHS": 1}    ← HS⁻ on the right

log₁₀(Q) = 9·pH  +  1·log_cHS  −  1·log_cSO4
```

The app looks up `log_cHS` in `species.csv`, sees `slider=false`, and uses its default value (0, i.e. 1 M standard state) silently. If you later want a HS⁻ slider, change `slider` to `true` in `species.csv` — the half-reaction row stays exactly the same.

### `data/reactions.csv` — metabolism (redox pair) table

One row per microbial metabolism shown in the right panel. Each metabolism links two half-reactions by their `id`.

| Column | Type | Description |
|--------|------|-------------|
| `id` | string | Unique identifier |
| `label` | string | Full metabolism name shown in the UI |
| `color` | hex string | Arrow colour |
| `ne` | int | Electrons transferred in the overall reaction |
| `donor` | string | `id` of the donor half-reaction (runs as oxidation) |
| `acceptor` | string | `id` of the acceptor half-reaction (runs as reduction) |
| `eq_don` | string | Donor half-reaction equation (displayed in detail panel) |
| `eq_acc` | string | Acceptor half-reaction equation |
| `eq_net` | string | Overall net reaction |

### `data/presets.csv` — named environmental conditions

One row per preset shown in the **"Selectable conditions"** dropdown. Selecting a preset instantly moves all sliders to the corresponding values.

| Column | Type | Description |
|--------|------|-------------|
| `name` | string | Preset name shown in the dropdown (e.g. `Cow rumen`) |
| `pH` | float | pH value |
| `log_pO2` | float | log₁₀ O₂ partial pressure (bar) |
| `log_pH2` | float | log₁₀ H₂ partial pressure (bar) |
| `log_pCO2` | float | log₁₀ CO₂ partial pressure (bar) |
| `log_pCH4` | float | log₁₀ CH₄ partial pressure (bar) |
| `log_cSO4` | float | log₁₀ SO₄²⁻ concentration (M) |
| `log_cAc` | float | log₁₀ acetate concentration (M) |
| `log_cBut` | float | log₁₀ butyrate concentration (M) |
| `log_cGlc` | float | log₁₀ glucose concentration (M) |
| `log_cProp` | float | log₁₀ propionate concentration (M) |

The column names after `name` match the `key` values of all slider species in `species.csv`. If you add a new slider species, add a corresponding column to `presets.csv` and fill in a value for each row.

**UI behaviour:** selecting a preset sets all sliders to the preset values and shows the preset name in the dropdown. Moving any slider afterwards reverts the dropdown to "— custom —", indicating that conditions no longer match any named preset.

The reference conditions are chosen to illustrate typical scenarios for each environment and should be treated as order-of-magnitude estimates. Actual concentrations and partial pressures can deviate substantially from these representative values depending on the specific site, depth, season, or community composition.

To add a new preset, simply add one row to this file and run the build step. No changes to `app.py` are needed.

---

## How the app works

### Nernst equation (generic)

`app.py` contains a single function `E_nernst(hr, **kw)` that computes the reduction potential for any half-reaction using the data from `half_rxns.csv`:

```
E = E°(SHE) − (RT / nₑF) · ln(Q)
```

where Q = [products] / [reactants]. In log₁₀ form:

```
log₁₀(Q) = nH · pH
          + Σ power × log₁₀(activity)   for all variable and fixed products
          − Σ power × log₁₀(activity)   for all variable and fixed reactants
```

### Free energy

For a complete metabolism pairing donor D with acceptor A:

```
ΔE = E_acceptor − E_donor
ΔG = −nₑ · F · ΔE
```

ΔG < 0 means the metabolism is exergonic (thermodynamically feasible).

### Sliders

Sliders are generated automatically from `species.csv` at startup. Each slider updates all half-reaction potentials and ΔG values in real time.

### Column layout

The tower is divided into four columns defined by `COL_DEFS` in `app.py`. The column assigned to each half-reaction is stored in the `col` field of `half_rxns.csv`.

---

## Development workflow

### Setup

```bash
conda env create -f environment.yml
conda activate redox-stack
```

### Run locally (fast, with autoreload)

```bash
panel serve app.py --show --autoreload
```

This opens the Panel app directly in the browser. Navigation links to the other HTML pages will not work in this mode.

### Run as full website (mirrors GitHub Pages)

```bash
cd docs && python -m http.server 8000
```

Then open http://localhost:8000. The interactive app loads via Pyodide (~30–60 s on first visit while the Python runtime downloads).

### Build the interactive app for deployment

Run this whenever `app.py` or any CSV in `data/` changes:

```bash
cp data/*.csv docs/data/
panel convert app.py --to pyodide-worker --out docs --title "Interactive Redox Tower"
```

This generates `docs/app.html` and `docs/app.js`. Both files must be committed.

### Deploy

Push to `main`. GitHub Pages serves the `docs/` folder automatically.

---

## How to extend the tower

### Add a new species (slider)

1. Add a row to `data/species.csv` with a new `key` (e.g. `log_cLac`).
2. Reference that key in `var_reactants` or `var_products` of any half-reaction in `half_rxns.csv`.
3. Run the build step above.

### Add a new half-reaction

Add one row to `data/half_rxns.csv`. Fields to fill in:

- Look up E°(SHE) from a thermodynamic database (note: this is at pH 0, not the biochemical E°′).
- Count `ne` (electrons) and `nH` (H⁺ consumed) from the balanced reduction half-reaction.
- List all non-H⁺/non-H₂O species in `reactants` and `products` using their `species.csv` keys. If a species is not yet in `species.csv`, add it there first (with `slider=false` initially if you don't want a slider yet).

### Add a new metabolism

Add one row to `data/reactions.csv`, referencing the `id` values of an existing donor and acceptor half-reaction.

### Add a new environmental preset

Add one row to `data/presets.csv` with a `name` and a value for each slider species column. The new preset appears in the dropdown automatically — no changes to `app.py` needed.
