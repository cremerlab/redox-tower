import sys, io, csv, json
import numpy as np
import panel as pn
import matplotlib
matplotlib.use('agg')  # headless backend required for Panel/Pyodide (no display)
import matplotlib.pyplot as plt

pn.extension(sizing_mode='stretch_width')

# ── Constants ─────────────────────────────────────────────────────────────────
R    = 8.314e-3
F    = 96.485
LN10 = np.log(10)

def _rtf(T_K=298.15):
    return R * T_K / F

# ── CSV loader (works locally and in Pyodide/browser) ─────────────────────────
def _read_text(filename):
    if sys.platform == 'emscripten':
        import pyodide.http
        return pyodide.http.open_url(f'data/{filename}').read()
    from pathlib import Path
    return (Path(__file__).parent / 'data' / filename).read_text(encoding='utf-8')

def _load_csv(filename):
    return list(csv.DictReader(io.StringIO(_read_text(filename))))

# ── Species — from data/species.csv ───────────────────────────────────────────
# Columns: key, label, disp, kind, start, end, step, default, slider
# kind  : 'pH' | 'gas' (log₁₀ bar) | 'conc' (log₁₀ M)
# slider: 'true' → FloatSlider shown in UI; 'false' → fixed at default value
SPECIES = {}
for _row in _load_csv('species.csv'):
    SPECIES[_row['key']] = dict(
        label=_row['label'], disp=_row['disp'], kind=_row['kind'],
        start=float(_row['start']), end=float(_row['end']),
        step=float(_row['step']),   default=float(_row['default']),
        slider=(_row['slider'].strip().lower() == 'true'),
    )

# ── Half-reactions — from data/half_rxns.csv ──────────────────────────────────
# Columns: id, label, color, col, E_SHE, ne, nH, E_std, reactants (JSON), products (JSON)
# reactants / products: {"species_key": stoichiometric_power}
# All species (slider or fixed) are referenced the same way.
# H⁺ and H₂O are not listed — H⁺ is handled via nH; H₂O activity = 1 by convention.
HALF_RXNS = []
for _row in _load_csv('half_rxns.csv'):
    for _f in ('reactants', 'products'):    _row[_f] = json.loads(_row[_f])
    for _f in ('col', 'ne', 'nH'):         _row[_f] = int(_row[_f])
    for _f in ('E_SHE', 'E_std'):          _row[_f] = float(_row[_f])
    HALF_RXNS.append(_row)

# ── Generic Nernst equation ───────────────────────────────────────────────────
# E = E°(SHE) − (RT / nₑF) · ln(Q),  Q = [products] / [reactants]
# H⁺ handled via nH: contributes nH·pH to log₁₀(Q) (H⁺ on reactant side).
# kw must contain values for every species key that appears in reactants/products.
def E_nernst(hr, **kw):
    log10_Q = hr['nH'] * kw['pH']
    for sp, p in hr['products'].items():
        log10_Q += p * kw[sp]
    for sp, p in hr['reactants'].items():
        log10_Q -= p * kw[sp]
    return 1000.0 * (hr['E_SHE'] - (_rtf() / hr['ne']) * LN10 * log10_Q)

# x0/x1: left/right x of the bar in axes coordinates; tx: label x; hx: column header x
COL_DEFS = {
    0: dict(x0=0.02, x1=0.11, tx=0.12, hx=0.065,  title='Aerobic'),
    1: dict(x0=0.25, x1=0.34, tx=0.35, hx=0.295,  title='Anaerobic\nrespiration'),
    2: dict(x0=0.50, x1=0.59, tx=0.60, hx=0.545,  title='Fermentation'),
    3: dict(x0=0.76, x1=0.85, tx=0.86, hx=0.805,  title='Syntrophic /\nC1'),
}

HR_COL = {hr['id']: hr['col'] for hr in HALF_RXNS}

# ── Reactions — from data/reactions.csv ───────────────────────────────────────
# Columns: id, label, color, ne, donor, acceptor, eq_don, eq_acc, eq_net
REACTIONS = []
for _row in _load_csv('reactions.csv'):
    _row['ne'] = int(_row['ne'])
    REACTIONS.append(_row)

# ── Defaults and presets ──────────────────────────────────────────────────────
# DEFAULTS covers every species (slider and fixed alike); used as the base kw dict.
# Slider species get their defaults overridden by widget values at runtime.
DEFAULTS = {k: v['default'] for k, v in SPECIES.items()}

# ── Presets — from data/presets.csv ──────────────────────────────────────────
# Columns: name + one column per slider species key
# Each row is a named environmental condition shown in the dropdown selector.
PRESETS = {}
for _row in _load_csv('presets.csv'):
    _name = _row.pop('name')
    PRESETS[_name] = {k: float(v) for k, v in _row.items()}

# ── Sliders (auto-generated for species with slider=true) ─────────────────────
SLIDERS = {
    k: pn.widgets.FloatSlider(
        name=sp['label'], value=sp['default'],
        start=sp['start'], end=sp['end'], step=sp['step'],
        width=250,
    )
    for k, sp in SPECIES.items() if sp['slider']
}

def _build_kw(slider_vals):
    """Merge fixed species defaults with current slider values."""
    kw = dict(DEFAULTS)                          # all species at their defaults
    kw.update(zip(SLIDERS.keys(), slider_vals))  # override with slider values
    return kw

def _set_sliders(preset_values):
    for k, w in SLIDERS.items():
        if k in preset_values:
            w.value = preset_values[k]

_CUSTOM_LABEL = '— select reference condition —'
_PRESET_OPTIONS = [_CUSTOM_LABEL] + list(PRESETS.keys())
_preset_select = pn.widgets.Select(
    name='', options=_PRESET_OPTIONS,
    value=_CUSTOM_LABEL, width=250
)

_applying_preset = [False]  # list so the flag is mutable inside nested functions

def _apply_preset(event):
    name = event.new
    if name == _CUSTOM_LABEL:
        return
    _applying_preset[0] = True
    _set_sliders(PRESETS[name])
    _applying_preset[0] = False

_preset_select.param.watch(_apply_preset, 'value')

def _on_slider_change(event):
    if not _applying_preset[0]:
        _preset_select.value = _CUSTOM_LABEL

for _w in SLIDERS.values():
    _w.param.watch(_on_slider_change, 'value')

rxn_checks = {
    r['id']: pn.widgets.Checkbox(name='', value=(r['id'] == 'aerobic'), width=20)
    for r in REACTIONS
}

def _rxn_feasibility(kw):
    E_now = {hr['id']: E_nernst(hr, **kw) for hr in HALF_RXNS}
    result = {}
    for rxn in REACTIONS:
        dE = E_now[rxn['acceptor']] - E_now[rxn['donor']]
        dG = -rxn['ne'] * F * dE / 1000
        result[rxn['id']] = (dG < 0, dG)
    return result

@pn.depends(*SLIDERS.values())
def rxn_list_pane(*vals):
    kw = _build_kw(vals)
    feasibility = _rxn_feasibility(kw)
    rows = []
    for rxn in REACTIONS:
        feasible, dG = feasibility[rxn['id']]
        sym   = '✓' if feasible else '✗'
        col   = '#2E7D32' if feasible else '#C62828'
        swatch = (f'<span style="display:inline-block;width:10px;height:10px;'
                  f'background:{rxn["color"]};border-radius:2px;'
                  f'vertical-align:middle;margin-right:4px"></span>')
        html  = (swatch
                 + f'<span style="color:{col};font-weight:bold;font-size:14px">{sym}</span>'
                 + f'<span style="font-size:13px;margin-left:3px">{rxn["label"]}</span>')
        rows.append(pn.Row(pn.pane.HTML(html, width=230, height=22),
                           rxn_checks[rxn['id']],
                           margin=(0, 0, 0, 0)))
    return pn.Column(*rows)

# ── Conditions summary (auto-generated from SPECIES) ─────────────────────────
@pn.depends(*SLIDERS.values())
def conditions_display(*vals):
    kw = _build_kw(vals)
    rows = []
    for key, sp in SPECIES.items():
        if not sp['slider']:
            continue
        v = kw[key]
        if sp['kind'] == 'pH':
            display = f'{v:.1f}'
        elif sp['kind'] == 'gas':
            pres = 10 ** v
            display = f'{pres:.3g} bar' if pres >= 0.005 else f'{pres * 1e6:.2g} ppm'
        else:  # conc
            display = f'{10 ** v * 1000:.3g} mM'
        rows.append((sp['disp'], display))
    inner = ''.join(
        f'<tr><td style="color:#555;padding-right:6px">{k}</td>'
        f'<td style="font-weight:bold">{d}</td></tr>'
        for k, d in rows
    )
    return pn.pane.HTML(
        f'<div style="font-size:13px;font-family:monospace;background:#f5f5f5;'
        f'padding:5px 8px;border-radius:4px;border:1px solid #ddd">'
        f'<b style="font-size:13px">Current conditions</b>'
        f'<table style="border-collapse:collapse;margin-top:3px">{inner}</table>'
        f'</div>',
        width=250
    )

# ── Label placement helper ────────────────────────────────────────────────────
def resolve_label_positions(e_vals, min_gap=42):
    if not e_vals:
        return []
    order = sorted(range(len(e_vals)), key=lambda i: e_vals[i])
    pos = [e_vals[i] for i in order]
    for _ in range(300):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i-1] < min_gap:
                mid = (pos[i] + pos[i-1]) / 2
                pos[i-1] = mid - min_gap / 2
                pos[i]   = mid + min_gap / 2
                moved = True
        if not moved:
            break
    result = [0.0] * len(e_vals)
    for rank, orig_idx in enumerate(order):
        result[orig_idx] = pos[rank]
    return result

# ── Plot ──────────────────────────────────────────────────────────────────────
Y_MIN, Y_MAX = -550, 950

def make_figure(kw, active_rxns):
    E_now = {hr['id']: E_nernst(hr, **kw) for hr in HALF_RXNS}

    fig, ax = plt.subplots(figsize=(8.75, 7.5))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(left=0.10, right=0.99, top=0.91, bottom=0.04)

    for cd in COL_DEFS.values():
        ax.text(cd['hx'], Y_MAX - 10, cd['title'],
                ha='center', va='top', fontsize=13, fontweight='bold', color='#444',
                bbox=dict(boxstyle='round,pad=0.2', fc='#f0f0f0', ec='#bbb', lw=0.8))

    for col_idx, cd in COL_DEFS.items():
        hrs = [hr for hr in HALF_RXNS if hr['col'] == col_idx]
        e_vals = [E_now[hr['id']] for hr in hrs]
        label_ys = resolve_label_positions(e_vals)

        for hr, E, label_y in zip(hrs, e_vals, label_ys):
            c = '#555555'
            ax.plot([cd['x0'], cd['x1']], [hr['E_std'], hr['E_std']],
                    color=c, lw=1.2, linestyle='--', alpha=0.30)
            ax.plot([cd['x0'], cd['x1']], [E, E],
                    color=c, lw=1.5, solid_capstyle='butt', alpha=0.90)
            if abs(label_y - E) > 8:
                ax.plot([cd['x1'] + 0.003, cd['tx'] - 0.005],
                        [E, label_y],
                        color=c, lw=0.5, alpha=0.45, linestyle='-',
                        transform=ax.transData)
            ax.text(cd['tx'], label_y, hr['label'],
                    fontsize=12.0, color=c, fontweight='bold',
                    va='center', ha='left')

    same_col_n = {}
    for rxn in REACTIONS:
        if rxn['id'] not in active_rxns:
            continue
        E_don = E_now[rxn['donor']];  E_acc = E_now[rxn['acceptor']]
        col_d = HR_COL[rxn['donor']]; col_a = HR_COL[rxn['acceptor']]
        cd_d  = COL_DEFS[col_d];      cd_a  = COL_DEFS[col_a]
        lw    = 1.5

        dE       = E_acc - E_don
        dG       = -rxn['ne'] * F * dE / 1000
        feasible = dG < 0
        lstyle   = '-'  if feasible else '--'
        alpha    = 0.85 if feasible else 0.50
        col_dG   = '#2E7D32' if feasible else '#C62828'

        if col_d == col_a:
            n  = same_col_n.get(col_d, 0);  same_col_n[col_d] = n + 1
            xb = cd_d['x1'] + 0.008 + n * 0.013
            ymid = (E_don + E_acc) / 2
            ax.plot([xb, xb], [E_don, E_acc],
                    color=rxn['color'], lw=lw, linestyle=lstyle, alpha=alpha)
            ax.annotate('', xy=(xb, E_acc), xytext=(xb, E_don),
                        arrowprops=dict(arrowstyle='->', color=rxn['color'],
                                        lw=lw, mutation_scale=8, linestyle=lstyle))
        else:
            x_don = (cd_d['x0'] + cd_d['x1']) / 2
            x_acc = (cd_a['x0'] + cd_a['x1']) / 2
            xb    = (x_don + x_acc) / 2;  ymid = (E_don + E_acc) / 2
            ax.plot([x_don, x_acc], [E_don, E_acc],
                    color=rxn['color'], lw=lw, linestyle=lstyle, alpha=alpha)
            ax.annotate('', xy=(x_acc, E_acc), xytext=(x_don, E_don),
                        arrowprops=dict(arrowstyle='->', color=rxn['color'],
                                        lw=lw, mutation_scale=8, linestyle=lstyle))

    ax.set_xlim(0, 1)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_ylabel('Reduction Potential E (mV)', fontsize=14)
    ax.tick_params(labelsize=13)
    ax.grid(axis='y', alpha=0.2, linestyle=':')
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    return fig

def make_detail_html(active_rxns, kw):
    if not active_rxns:
        return ''
    E_now = {hr['id']: E_nernst(hr, **kw) for hr in HALF_RXNS}
    html = ''
    for rxn in REACTIONS:
        if rxn['id'] not in active_rxns:
            continue
        E_don = E_now[rxn['donor']];  E_acc = E_now[rxn['acceptor']]
        dE    = E_acc - E_don
        dG    = -rxn['ne'] * F * dE / 1000
        col   = '#2E7D32' if dG < 0 else '#C62828'
        spont = 'exergonic ✓' if dG < 0 else 'endergonic ✗'
        html += (
            f'<div style="font-family:monospace;font-size:13px;padding:6px 10px;'
            f'margin-bottom:6px;border:1px solid {rxn["color"]};border-radius:4px;'
            f'background:#fafafa;line-height:1.6">'
            f'<b style="font-size:15px;color:{rxn["color"]}">{rxn["label"]}</b><br>'
            f'<b>ox:</b>  {rxn["eq_don"]}<br>'
            f'<b>red:</b> {rxn["eq_acc"]}<br>'
            f'<b>net:</b> {rxn["eq_net"]}<br>'
            f'E<sub>donor</sub> = {E_don:+.0f} mV &nbsp;|&nbsp; '
            f'E<sub>acc</sub> = {E_acc:+.0f} mV &nbsp;|&nbsp; '
            f'ΔE = {dE:+.0f} mV &nbsp;|&nbsp; '
            f'<b style="color:{col}">ΔG = {dG:+.1f} kJ/mol ({spont})</b>'
            f'</div>'
        )
    return html

@pn.depends(*SLIDERS.values(), *rxn_checks.values())
def plot_pane(*vals):
    # Panel passes all watched values as a flat tuple: first n_sp are slider values,
    # the rest are checkbox booleans — split accordingly.
    n_sp = len(SLIDERS)
    kw = _build_kw(vals[:n_sp])
    active = [r['id'] for r, v in zip(REACTIONS, vals[n_sp:]) if v]
    fig = make_figure(kw, active)
    pane = pn.pane.Matplotlib(fig, tight=True, format='svg', width=700, height=650)
    plt.close(fig)
    return pane

@pn.depends(*SLIDERS.values(), *rxn_checks.values())
def detail_pane_fn(*vals):
    n_sp = len(SLIDERS)  # same split as plot_pane
    kw = _build_kw(vals[:n_sp])
    active = [r['id'] for r, v in zip(REACTIONS, vals[n_sp:]) if v]
    return pn.pane.HTML(
        make_detail_html(active, kw),
        sizing_mode='stretch_width'
    )

# ── Layout ────────────────────────────────────────────────────────────────────
LEGEND_HTML = """
<div style="font-size:13px;background:#f9f9f9;border:1px solid #ddd;
            border-radius:4px;padding:5px 8px;margin-bottom:4px;line-height:2.2">
  <b style="font-size:13px;color:#333">Legend</b><br>
  <span style="display:inline-block;width:32px;border-bottom:2.5px solid #555;
               vertical-align:middle"></span>
  &nbsp;Potential at chosen conditions<br>
  <span style="display:inline-block;width:32px;border-bottom:2px dashed #aaa;
               vertical-align:middle"></span>
  &nbsp;Standard potential (E°′)<br>
  <span style="color:#2E7D32;font-weight:bold">&#10132; ✓</span>
  &nbsp;Metabolism feasible (ΔG &lt; 0)<br>
  <span style="color:#C62828;font-weight:bold">&#10132; ✗</span>
  &nbsp;Metabolism not feasible (ΔG &gt; 0)
</div>
"""

left_panel = pn.Column(
    pn.pane.HTML('<b style="font-size:15px">Environmental conditions</b>'
                 '<span style="font-size:12px;color:#666;margin-left:4px">'
                 '— use sliders to adjust</span>'),
    *SLIDERS.values(),
    conditions_display,
    pn.pane.HTML('<b style="font-size:14px;margin-top:6px">Selectable conditions</b>'),
    _preset_select,
    width=270, sizing_mode='fixed'
)

right_panel = pn.Column(
    pn.pane.HTML(LEGEND_HTML),
    pn.pane.HTML('<b style="font-size:14px">Select metabolisms to display</b>'),
    rxn_list_pane,
    width=290, sizing_mode='fixed'
)

NAV_HTML = """
<div style="background:#1A3A5C;padding:0 40px;height:56px;display:flex;align-items:center;
            justify-content:space-between;font-family:system-ui,-apple-system,sans-serif;
            margin-bottom:20px">
  <a href="index.html" style="color:#fff;text-decoration:none;font-size:21px;
            font-weight:700;letter-spacing:-0.3px">ElectronFlow</a>
  <div style="display:flex;gap:32px;align-items:center">
    <a href="app.html" style="color:#fff;text-decoration:none;font-size:14px;font-weight:600">
      Interactive Redox Tower</a>
    <a href="calc.html" style="color:rgba(255,255,255,0.65);text-decoration:none;
            font-size:14px;font-weight:500">Calc Details</a>
    <a href="feedback.html" style="color:rgba(255,255,255,0.65);text-decoration:none;
            font-size:14px;font-weight:500">Feedback</a>
    <a href="about.html" style="color:rgba(255,255,255,0.65);text-decoration:none;
            font-size:14px;font-weight:500">About</a>
  </div>
</div>
"""

app = pn.Column(
    pn.pane.HTML(NAV_HTML, sizing_mode='stretch_width', margin=0),
    pn.Column(
        pn.Row(
            left_panel,
            pn.Spacer(width=10),
            pn.Column(plot_pane, detail_pane_fn),
            pn.Spacer(width=40),
            right_panel,
        ),
        sizing_mode='fixed', width=1350, margin=(0, 0, 0, 0)
    ),
    sizing_mode='stretch_width',
    margin=0
)

app.servable()
