import numpy as np
import panel as pn
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

pn.extension(sizing_mode='stretch_width')

# ── Constants ─────────────────────────────────────────────────────────────────
R    = 8.314e-3
F    = 96.485
LN10 = np.log(10)

def _rtf(T_K=298.15):
    return R * T_K / F

# ── Nernst functions ──────────────────────────────────────────────────────────
def E_O2_H2O(pH, log_pO2, **kw):
    return 1000 * (1.229 + _rtf()/4*np.log(10**log_pO2) - _rtf()*LN10*pH)

def E_CO2_Glc(pH, log_pCO2, log_cGlc, **kw):
    return 1000 * (-0.016 + _rtf()/24*np.log((10**log_pCO2)**6/10**log_cGlc) - _rtf()*LN10*pH)

def E_NO3_N2(pH, **kw):
    return 1000 * (1.248 + _rtf()/10*np.log(1e-3**2) - (6/5)*_rtf()*LN10*pH)

def E_NO3_NO2(pH, **kw):
    return 1000 * (0.844 + _rtf()/2*np.log(1.0) - _rtf()*LN10*pH)

def E_NO3_NH4(pH, **kw):
    return 1000 * (0.878 + _rtf()/8*np.log(1e-3/1e-4) - (10/8)*_rtf()*LN10*pH)

def E_Fe(pH, **kw):
    # Fe(OH)3/Fe2+: E°'(pH7) = -100 mV at unit activity; pH-dependent only
    return 1000 * (0.314 - _rtf()*LN10*pH)

def E_SO4_HS(pH, log_cSO4, **kw):
    return 1000 * (0.249 + _rtf()/8*np.log(10**log_cSO4/1e-4) - (9/8)*_rtf()*LN10*pH)

def E_CO2_CH4(pH, log_pCO2, log_pCH4, **kw):
    return 1000 * (0.169 + _rtf()/8*np.log(10**log_pCO2/10**log_pCH4) - _rtf()*LN10*pH)

def E_NAD(pH, **kw):
    return 1000 * (-0.113 + _rtf()/2*np.log(1.0) - (1/2)*_rtf()*LN10*pH)

def E_H2(pH, log_pH2, **kw):
    return 1000 * (-_rtf()*LN10*pH - _rtf()/2*np.log(10**log_pH2))

def E_CO2_Ac(pH, log_pCO2, log_cAc, **kw):
    return 1000 * (0.072 + _rtf()/8*np.log((10**log_pCO2)**2/10**log_cAc) - (7/8)*_rtf()*LN10*pH)

def E_CO2_Prop(pH, log_pCO2, log_cProp, log_cAc, **kw):
    return 1000 * (0.132 + _rtf()/6*np.log(10**log_cAc * 10**log_pCO2 / 10**log_cProp) - _rtf()*LN10*pH)

def E_CO2_But(pH, log_pCO2, log_cBut, **kw):
    return 1000 * (0.113 + _rtf()/20*np.log((10**log_pCO2)**4/10**log_cBut) - (19/20)*_rtf()*LN10*pH)

def E_Crot_But(pH, **kw):
    return 1000 * (0.404 + _rtf()/2*np.log(1.0) - _rtf()*LN10*pH)

def E_Pyr_Lac(pH, **kw):
    return 1000 * (0.229 + _rtf()/2*np.log(1.0) - _rtf()*LN10*pH)

def E_AcAld_EtOH(pH, **kw):
    return 1000 * (0.217 + _rtf()/2*np.log(1e-4/1e-3) - _rtf()*LN10*pH)

def E_Ac_EtOH(pH, log_cAc, **kw):
    # Acetate/Ethanol couple: CH3COO- + 5H+ + 4e- → C2H5OH + H2O
    # E°'(pH0) = 0.130 V  → E°'(pH7) ≈ -388 mV  (from ΔG°'=+10 kJ/mol for syntrophic EtOH ox.)
    # [EtOH] assumed = 1 M; acetate from slider
    return 1000 * (0.130 + _rtf()/4*np.log(10**log_cAc) - (5/4)*_rtf()*LN10*pH)

# ── Half-reactions ────────────────────────────────────────────────────────────
HALF_RXNS = [
    dict(id='o2',      label='O2 / H2O',       color='#E65100', E_std=+816, col=0, func=E_O2_H2O),
    dict(id='no3n2',   label='NO3 / N2',        color='#558B2F', E_std=+751, col=0, func=E_NO3_N2),
    dict(id='no3no2',  label='NO3 / NO2',       color='#827717', E_std=+430, col=0, func=E_NO3_NO2),
    dict(id='glc',     label='CO2 / Glc',       color='#1565C0', E_std=-430, col=0, func=E_CO2_Glc),
    dict(id='no3nh4',  label='NO3 / NH4+',      color='#00838F', E_std=+360, col=1, func=E_NO3_NH4),
    dict(id='fe',      label='Fe(III)/Fe(II)',  color='#8D6E63', E_std=-100, col=1, func=E_Fe),
    dict(id='hs',      label='SO4 / HS-',       color='#6A1B9A', E_std=-217, col=1, func=E_SO4_HS),
    dict(id='ch4',     label='CO2 / CH4',       color='#2E7D32', E_std=-245, col=1, func=E_CO2_CH4),
    dict(id='pyrlac',  label='Pyr / Lac',       color='#AD1457', E_std=-185, col=2, func=E_Pyr_Lac),
    dict(id='etoh',    label='AcAld / EtOH',    color='#00695C', E_std=-197, col=2, func=E_AcAld_EtOH),
    dict(id='co2ac',   label='CO2 / Ac-',       color='#B71C1C', E_std=-290, col=2, func=E_CO2_Ac),
    dict(id='co2but',  label='CO2 / But-',      color='#6200EA', E_std=-280, col=2, func=E_CO2_But),
    dict(id='prop',    label='CO2 / Prop-',     color='#F9A825', E_std=-282, col=3, func=E_CO2_Prop),
    dict(id='nad',     label='NAD+ / NADH',     color='#6D4C41', E_std=-320, col=3, func=E_NAD),
    dict(id='etohac',  label='Ac / EtOH',       color='#00695C', E_std=-388, col=3, func=E_Ac_EtOH),
    dict(id='h2',      label='H+ / H2',         color='#37474F', E_std=-414, col=3, func=E_H2),
]

COL_DEFS = {
    0: dict(x0=0.02, x1=0.11, tx=0.12, hx=0.065,  title='Aerobic'),
    1: dict(x0=0.25, x1=0.34, tx=0.35, hx=0.295,  title='Anaerobic\nrespiration'),
    2: dict(x0=0.50, x1=0.59, tx=0.60, hx=0.545,  title='Fermentation'),
    3: dict(x0=0.76, x1=0.85, tx=0.86, hx=0.805,  title='Syntrophic /\nC1'),
}

HR_COL = {hr['id']: hr['col'] for hr in HALF_RXNS}

# ── Reactions ─────────────────────────────────────────────────────────────────
REACTIONS = [
    # ── Respirations (matching paper Table 1) ─────────────────────────────────
    dict(id='aerobic',   label='Aerobic respiration',          color='#E65100', ne=24,
         donor='glc',    acceptor='o2',
         eq_don='Glucose + 6H2O → 6CO2 + 24H⁺ + 24e⁻',
         eq_acc='6O2 + 24H⁺ + 24e⁻ → 12H2O',
         eq_net='Glucose + 6O2 → 6CO2 + 6H2O'),
    dict(id='denitrif',  label='Denitrification',              color='#558B2F', ne=10,
         donor='h2',     acceptor='no3n2',
         eq_don='5H2 → 10H⁺ + 10e⁻',
         eq_acc='2NO3⁻ + 12H⁺ + 10e⁻ → N2 + 6H2O',
         eq_net='5H2 + 2NO3⁻ + 2H⁺ → N2 + 6H2O'),
    dict(id='fe_red',    label='Iron reduction',               color='#8D6E63', ne=2,
         donor='h2',     acceptor='fe',
         eq_don='H2 → 2H⁺ + 2e⁻',
         eq_acc='2Fe(OH)3 + 4H⁺ + 2e⁻ → 2Fe²⁺ + 6H2O',
         eq_net='H2 + 2Fe(OH)3 + 4H⁺ → 2Fe²⁺ + 6H2O'),
    dict(id='h2sr',      label='H2-driven sulfate reduction',  color='#6A1B9A', ne=8,
         donor='h2',     acceptor='hs',
         eq_don='4H2 → 8H⁺ + 8e⁻',
         eq_acc='SO4²⁻ + 9H⁺ + 8e⁻ → HS⁻ + 4H2O',
         eq_net='4H2 + SO4²⁻ → HS⁻ + H⁺ + 4H2O'),
    dict(id='h2meth',    label='Hydrogenotrophic methanogenesis', color='#2E7D32', ne=8,
         donor='h2',     acceptor='ch4',
         eq_don='4H2 → 8H⁺ + 8e⁻',
         eq_acc='CO2 + 8H⁺ + 8e⁻ → CH4 + 2H2O',
         eq_net='4H2 + CO2 → CH4 + 2H2O'),
    dict(id='ac_meth',   label='Acetoclastic methanogenesis',  color='#33691E', ne=8,
         donor='co2ac',  acceptor='ch4',
         eq_don='CH3COO⁻ + 2H2O → 2CO2 + 7H⁺ + 8e⁻',
         eq_acc='CO2 + 8H⁺ + 8e⁻ → CH4 + 2H2O',
         eq_net='CH3COO⁻ + H2O → CH4 + HCO3⁻'),
    dict(id='aom',       label='Anaerobic oxidation of methane (AOM)', color='#37474F', ne=8,
         donor='ch4',    acceptor='hs',
         eq_don='CH4 + 2H2O → CO2 + 8H⁺ + 8e⁻',
         eq_acc='SO4²⁻ + 9H⁺ + 8e⁻ → HS⁻ + 4H2O',
         eq_net='CH4 + SO4²⁻ → HCO3⁻ + HS⁻ + H2O'),
    # ── Primary fermentations ─────────────────────────────────────────────────
    dict(id='homoacet',  label='Homoacetogenesis (H2 + CO2)',  color='#1B5E20', ne=8,
         donor='h2',     acceptor='co2ac',
         eq_don='4H2 → 8H⁺ + 8e⁻',
         eq_acc='2CO2 + 7H⁺ + 8e⁻ → CH3COO⁻ + 2H2O',
         eq_net='4H2 + 2CO2 → CH3COO⁻ + H⁺ + 2H2O'),
    dict(id='glc_etoh',  label='Ethanol fermentation',         color='#004D40', ne=4,
         donor='glc',    acceptor='etoh',
         eq_don='Glucose → 6CO2 + 24H⁺ + 24e⁻  (approx., ×4/24)',
         eq_acc='2 Acetaldehyde + 4H⁺ + 4e⁻ → 2 Ethanol',
         eq_net='Glucose → 2 Ethanol + 2 CO2  (alcoholic ferm.)'),
    dict(id='glc_but',   label='Butyric acid fermentation',    color='#1565C0', ne=20,
         donor='glc',    acceptor='co2but',
         eq_don='Glucose + 6H2O → 6CO2 + 24H⁺ + 24e⁻  (×20/24)',
         eq_acc='4CO2 + 19H⁺ + 20e⁻ → But⁻ + 6H2O',
         eq_net='Glucose + 2H2O → But⁻ + 2HCO3⁻ + 3H⁺ + 2H2'),
    dict(id='glc_lac',   label='Lactic acid fermentation',     color='#880E4F', ne=4,
         donor='glc',    acceptor='pyrlac',
         eq_don='Glucose → 6CO2 + 24H⁺ + 24e⁻  (approx., ×4/24)',
         eq_acc='2 Pyruvate + 4H⁺ + 4e⁻ → 2 Lactate',
         eq_net='Glucose → 2 Lactate  (glycolysis + LDH)'),
    # ── Secondary fermentations / syntrophic ──────────────────────────────────
    dict(id='prop_synt', label='Syntrophic propionate oxidation', color='#BF360C', ne=6,
         donor='prop',   acceptor='h2',
         eq_don='Prop⁻ + 3H2O → Ac⁻ + HCO3⁻ + 5H⁺ + 6e⁻',
         eq_acc='6H⁺ + 6e⁻ → 3H2',
         eq_net='Prop⁻ + 3H2O → Ac⁻ + HCO3⁻ + H⁺ + 3H2  (requires low H2)'),
    dict(id='but_synt',  label='Syntrophic butyrate oxidation', color='#4A148C', ne=4,
         donor='co2but', acceptor='h2',
         eq_don='But⁻ + 2H2O → 2Ac⁻ + H⁺ + 2H2  (overall)',
         eq_acc='4H⁺ + 4e⁻ → 2H2',
         eq_net='But⁻ + 2H2O → 2Ac⁻ + H⁺ + 2H2  (requires low H2)'),
    dict(id='etoh_synt', label='Syntrophic ethanol oxidation',  color='#00695C', ne=4,
         donor='etohac', acceptor='h2',
         eq_don='C2H5OH + H2O → CH3COO⁻ + 5H⁺ + 4e⁻',
         eq_acc='4H⁺ + 4e⁻ → 2H2',
         eq_net='C2H5OH + H2O → CH3COO⁻ + H⁺ + 2H2  (requires low H2)'),
]

# ── Default values ────────────────────────────────────────────────────────────
DEFAULTS = dict(
    # Biochemical standard state: pH 7, all species at unit activity (1 M / 1 bar)
    # At these values the solid lines coincide with the dashed E°′ reference lines.
    pH=7.0,
    log_pO2=0.0,
    log_pH2=0.0,
    log_pCO2=0.0,
    log_pCH4=0.0,
    log_cSO4=0.0,
    log_cAc=0.0,
    log_cBut=0.0,
    log_cGlc=0.0,
    log_cProp=0.0,
)

# ── Preset conditions ─────────────────────────────────────────────────────────
PRESETS = {
    'Cow rumen': dict(
        pH=6.5,  log_pO2=-7.0,
        log_pH2=-4.5,  log_pCO2=-0.19, log_pCH4=-0.60,
        log_cSO4=-3.0, log_cAc=-1.15,  log_cBut=-1.85, log_cGlc=-4.0, log_cProp=-1.65,
    ),
    'Human colon': dict(
        # High H2: fermentation not fully quenched by methanogenesis
        # Low CH4: many humans have limited methanogenic capacity
        pH=6.8,  log_pO2=-7.0,
        log_pH2=-1.5,  log_pCO2=-0.60, log_pCH4=-3.0,
        log_cSO4=-2.7, log_cAc=-1.5,   log_cBut=-2.0,  log_cGlc=-4.5, log_cProp=-2.0,
    ),
    'Marine sediment': dict(
        pH=7.5,  log_pO2=-7.0,
        log_pH2=-5.0,  log_pCO2=-0.50, log_pCH4=-2.5,
        log_cSO4=-1.55,log_cAc=-2.5,   log_cBut=-3.5,  log_cGlc=-5.0, log_cProp=-3.0,
    ),
    'Freshwater sediment': dict(
        pH=7.0,  log_pO2=-7.0,
        log_pH2=-4.0,  log_pCO2=-0.70, log_pCH4=-1.5,
        log_cSO4=-4.0, log_cAc=-2.5,   log_cBut=-3.0,  log_cGlc=-4.5, log_cProp=-3.0,
    ),
    'Anaerobic digester': dict(
        # Low H2 (~10 ppm) maintained by tight methanogenic coupling
        # → syntrophic butyrate/propionate oxidation becomes feasible
        pH=7.2,  log_pO2=-7.0,
        log_pH2=-5.0,  log_pCO2=-0.40, log_pCH4=-0.20,
        log_cSO4=-3.0, log_cAc=-2.5,   log_cBut=-3.0,  log_cGlc=-5.0, log_cProp=-3.0,
    ),
    'Activated sludge': dict(
        pH=7.2,  log_pO2=-2.0,
        log_pH2=-7.0,  log_pCO2=-1.5,  log_pCH4=-6.0,
        log_cSO4=-2.0, log_cAc=-3.5,   log_cBut=-5.0,  log_cGlc=-4.0, log_cProp=-4.0,
    ),
}

# ── Widgets ───────────────────────────────────────────────────────────────────
w_pH    = pn.widgets.FloatSlider(name='pH',                      value=DEFAULTS['pH'],       start=4.0,  end=10.0, step=0.1,  width=250)
w_lo2   = pn.widgets.FloatSlider(name='O₂ partial pressure (log₁₀ bar)',   value=DEFAULTS['log_pO2'],  start=-10.0,end=0.0,  step=0.5,  width=250)
w_lh2   = pn.widgets.FloatSlider(name='H₂ partial pressure (log₁₀ bar)',   value=DEFAULTS['log_pH2'],  start=-7.0, end=0.0,  step=0.25, width=250)
w_lco2  = pn.widgets.FloatSlider(name='CO₂ partial pressure (log₁₀ bar)',  value=DEFAULTS['log_pCO2'], start=-3.0, end=0.0,  step=0.25, width=250)
w_lch4  = pn.widgets.FloatSlider(name='CH₄ partial pressure (log₁₀ bar)',  value=DEFAULTS['log_pCH4'], start=-4.0, end=0.0,  step=0.25, width=250)
w_lso4  = pn.widgets.FloatSlider(name='SO₄²⁻ concentration (log₁₀ M)',     value=DEFAULTS['log_cSO4'], start=-4.0, end=0.0,  step=0.25, width=250)
w_lac   = pn.widgets.FloatSlider(name='Acetate (log₁₀ M)',                  value=DEFAULTS['log_cAc'],  start=-6.0, end=0.0,  step=0.25, width=250)
w_lbut  = pn.widgets.FloatSlider(name='Butyrate (log₁₀ M)',                 value=DEFAULTS['log_cBut'], start=-6.0, end=0.0,  step=0.25, width=250)
w_lglc  = pn.widgets.FloatSlider(name='Glucose (log₁₀ M)',                  value=DEFAULTS['log_cGlc'], start=-6.0, end=0.0,  step=0.25, width=250)
w_lprop = pn.widgets.FloatSlider(name='Propionate (log₁₀ M)',               value=DEFAULTS['log_cProp'],start=-6.0, end=0.0,  step=0.25, width=250)

def _set_sliders(p):
    w_pH.value    = p['pH']
    w_lo2.value   = p['log_pO2']
    w_lh2.value   = p['log_pH2']
    w_lco2.value  = p['log_pCO2']
    w_lch4.value  = p['log_pCH4']
    w_lso4.value  = p['log_cSO4']
    w_lac.value   = p['log_cAc']
    w_lbut.value  = p['log_cBut']
    w_lglc.value  = p['log_cGlc']
    w_lprop.value = p['log_cProp']

_PRESET_OPTIONS = ['Standard conditions (1M, pH=7)'] + list(PRESETS.keys())
_preset_select = pn.widgets.Select(
    name='', options=_PRESET_OPTIONS,
    value='Standard conditions (1M, pH=7)', width=250
)

def _apply_preset(event):
    name = event.new
    if name == 'Standard conditions (1M, pH=7)':
        _set_sliders(DEFAULTS)
    else:
        _set_sliders(PRESETS[name])

_preset_select.param.watch(_apply_preset, 'value')

rxn_checks = {
    r['id']: pn.widgets.Checkbox(name='', value=(r['id'] == 'aerobic'), width=20)
    for r in REACTIONS
}

def _rxn_feasibility(kw):
    """Return dict of rxn_id -> (feasible, dG)."""
    E_now = {hr['id']: hr['func'](**kw) for hr in HALF_RXNS}
    result = {}
    for rxn in REACTIONS:
        dE = E_now[rxn['acceptor']] - E_now[rxn['donor']]
        dG = -rxn['ne'] * F * dE / 1000
        result[rxn['id']] = (dG < 0, dG)
    return result

@pn.depends(w_pH, w_lo2, w_lh2, w_lco2, w_lch4, w_lso4, w_lac, w_lbut, w_lglc, w_lprop)
def rxn_list_pane(pH, lo2, lh2, lco2, lch4, lso4, lac, lbut, lglc, lprop):
    kw = dict(pH=pH, log_pO2=lo2, log_pH2=lh2, log_pCO2=lco2, log_pCH4=lch4, log_cSO4=lso4,
              log_cAc=lac, log_cBut=lbut, log_cGlc=lglc, log_cProp=lprop)
    feasibility = _rxn_feasibility(kw)
    rows = []
    for rxn in REACTIONS:
        feasible, dG = feasibility[rxn['id']]
        sym   = '✓' if feasible else '✗'
        col   = '#2E7D32' if feasible else '#C62828'
        html  = (f'<span style="color:{col};font-weight:bold;font-size:14px">{sym}</span>'
                 f'<span style="font-size:13px;margin-left:3px">{rxn["label"]}</span>')
        rows.append(pn.Row(pn.pane.HTML(html, width=230, height=22),
                           rxn_checks[rxn['id']],
                           margin=(0, 0, 0, 0)))
    return pn.Column(*rows)

# ── Reactive conditions display (real units) ──────────────────────────────────
@pn.depends(w_pH, w_lo2, w_lh2, w_lco2, w_lch4, w_lso4, w_lac, w_lbut, w_lglc, w_lprop)
def conditions_display(pH, lo2, lh2, lco2, lch4, lso4, lac, lbut, lglc, lprop):
    def fmt_conc(log_val, unit='mM', scale=1000):
        v = 10**log_val * scale
        return f'{v:.3g} {unit}'
    def fmt_pres(log_val):
        v = 10**log_val
        if v >= 0.005:
            return f'{v:.3g} bar'
        return f'{v*1e6:.2g} ppm'
    rows = [
        ('pH', f'{pH:.1f}'),
        ('O₂', fmt_pres(lo2)),
        ('H₂', fmt_pres(lh2)),
        ('CO₂', fmt_pres(lco2)),
        ('CH₄', fmt_pres(lch4)),
        ('SO₄²⁻', fmt_conc(lso4)),
        ('Acetate', fmt_conc(lac)),
        ('Butyrate', fmt_conc(lbut)),
        ('Glucose', fmt_conc(lglc)),
        ('Propionate', fmt_conc(lprop)),
    ]
    inner = ''.join(
        f'<tr><td style="color:#555;padding-right:6px">{k}</td>'
        f'<td style="font-weight:bold">{v}</td></tr>'
        for k, v in rows
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
    """Return y-positions that avoid overlap, preserving relative order."""
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

def make_figure(pH, log_pO2, log_pH2, log_pCO2, log_pCH4, log_cSO4,
                log_cAc, log_cBut, log_cGlc, log_cProp, active_rxns):
    kw = dict(pH=pH, log_pO2=log_pO2, log_pH2=log_pH2, log_pCO2=log_pCO2,
              log_pCH4=log_pCH4, log_cSO4=log_cSO4,
              log_cAc=log_cAc, log_cBut=log_cBut,
              log_cGlc=log_cGlc, log_cProp=log_cProp)
    E_now = {hr['id']: hr['func'](**kw) for hr in HALF_RXNS}

    fig, ax = plt.subplots(figsize=(8.75, 7.5))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(left=0.10, right=0.99, top=0.91, bottom=0.04)

    # ── Column titles ──
    for cd in COL_DEFS.values():
        ax.text(cd['hx'], Y_MAX - 10, cd['title'],
                ha='center', va='top', fontsize=13, fontweight='bold', color='#444',
                bbox=dict(boxstyle='round,pad=0.2', fc='#f0f0f0', ec='#bbb', lw=0.8))

    # ── Half-reaction bars + non-overlapping labels ──
    for col_idx, cd in COL_DEFS.items():
        hrs = [hr for hr in HALF_RXNS if hr['col'] == col_idx]
        e_vals = [E_now[hr['id']] for hr in hrs]
        label_ys = resolve_label_positions(e_vals)

        for hr, E, label_y in zip(hrs, e_vals, label_ys):
            c = hr['color']
            ax.plot([cd['x0'], cd['x1']], [hr['E_std'], hr['E_std']],
                    color=c, lw=1.2, linestyle='--', alpha=0.30)
            ax.plot([cd['x0'], cd['x1']], [E, E],
                    color=c, lw=1.5, solid_capstyle='butt', alpha=0.90)
            # leader line if label is offset from bar
            if abs(label_y - E) > 8:
                ax.plot([cd['x1'] + 0.003, cd['tx'] - 0.005],
                        [E, label_y],
                        color=c, lw=0.5, alpha=0.45, linestyle='-',
                        transform=ax.transData)
            ax.text(cd['tx'], label_y, hr['label'],
                    fontsize=12.0, color=c, fontweight='bold',
                    va='center', ha='left')

    # ── Reaction arrows ──
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

        ax.text(xb + 0.01, ymid, f'{dG:+.0f} kJ', fontsize=12, color=col_dG,
                va='center', ha='left', fontweight='bold')

    ax.set_xlim(0, 1)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_ylabel('Reduction Potential E (mV)', fontsize=14)
    ax.tick_params(labelsize=13)
    ax.grid(axis='y', alpha=0.2, linestyle=':')
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    return fig

def make_detail_html(active_rxns, pH, log_pO2, log_pH2, log_pCO2, log_pCH4, log_cSO4,
                     log_cAc, log_cBut, log_cGlc, log_cProp):
    if not active_rxns:
        return ''
    kw = dict(pH=pH, log_pO2=log_pO2, log_pH2=log_pH2, log_pCO2=log_pCO2,
              log_pCH4=log_pCH4, log_cSO4=log_cSO4,
              log_cAc=log_cAc, log_cBut=log_cBut,
              log_cGlc=log_cGlc, log_cProp=log_cProp)
    E_now = {hr['id']: hr['func'](**kw) for hr in HALF_RXNS}
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

@pn.depends(w_pH, w_lo2, w_lh2, w_lco2, w_lch4, w_lso4,
            w_lac, w_lbut, w_lglc, w_lprop, *rxn_checks.values())
def plot_pane(pH, log_pO2, log_pH2, log_pCO2, log_pCH4, log_cSO4,
              log_cAc, log_cBut, log_cGlc, log_cProp, *check_vals):
    active = [r['id'] for r, v in zip(REACTIONS, check_vals) if v]
    fig = make_figure(pH, log_pO2, log_pH2, log_pCO2, log_pCH4, log_cSO4,
                      log_cAc, log_cBut, log_cGlc, log_cProp, active)
    pane = pn.pane.Matplotlib(fig, tight=True, format='svg', width=700, height=650)
    plt.close(fig)
    return pane

@pn.depends(w_pH, w_lo2, w_lh2, w_lco2, w_lch4, w_lso4,
            w_lac, w_lbut, w_lglc, w_lprop, *rxn_checks.values())
def detail_pane_fn(pH, log_pO2, log_pH2, log_pCO2, log_pCH4, log_cSO4,
                   log_cAc, log_cBut, log_cGlc, log_cProp, *check_vals):
    active = [r['id'] for r, v in zip(REACTIONS, check_vals) if v]
    return pn.pane.HTML(
        make_detail_html(active, pH, log_pO2, log_pH2, log_pCO2, log_pCH4, log_cSO4,
                         log_cAc, log_cBut, log_cGlc, log_cProp),
        sizing_mode='stretch_width'
    )

# ── Layout ────────────────────────────────────────────────────────────────────
LEGEND_HTML = """
<div style="font-size:13px;background:#f9f9f9;border:1px solid #ddd;
            border-radius:4px;padding:5px 8px;margin-bottom:4px;line-height:2.1">
  <b style="font-size:13px;color:#333">Legend</b><br>
  <span style="letter-spacing:2px;color:#888">&#8212; &#8212;</span>
  &nbsp;Dashed line — standard potentials<br>
  <span style="color:#333;font-weight:bold">&#8212;&#8212;&#8212;</span>
  &nbsp;Solid line — potentials at selected conditions<br>
  <span style="color:#2E7D32;font-weight:bold">&#10132; solid green arrow</span>
  &nbsp;— feasible metabolism (ΔG &lt; 0)<br>
  <span style="color:#C62828;font-weight:bold">&#10132; red dashed arrow</span>
  &nbsp;— non feasible metabolism (ΔG &gt; 0)
</div>
"""

left_panel = pn.Column(
    pn.pane.HTML('<b style="font-size:15px">Environmental conditions</b>'
                 '<span style="font-size:12px;color:#666;margin-left:4px">'
                 '— use sliders to adjust</span>'),
    w_pH, w_lo2, w_lh2, w_lco2, w_lch4, w_lso4,
    w_lac, w_lbut, w_lglc, w_lprop,
    conditions_display,
    pn.pane.HTML('<b style="font-size:14px;margin-top:6px">Selectable conditions</b>'),
    _preset_select,
    width=270, sizing_mode='fixed'
)

right_panel = pn.Column(
    pn.pane.HTML('<b style="font-size:15px">Possible microbial metabolisms</b>'
                 '<span style="font-size:12px;color:#666;margin-left:4px">'
                 '— check to display arrow</span>'),
    pn.pane.HTML(LEGEND_HTML),
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


METHODS_HTML = """
<div style="margin-top:18px;padding:12px 16px;background:#f8f9fa;border:1px solid #dee2e6;
            border-radius:6px;font-size:14px;font-family:Georgia,serif;line-height:1.7;
            max-width:1310px">
  <b style="font-size:16px">How condition dependence is calculated</b>
  <p style="margin:6px 0">
    Each redox couple is characterised by its <b>standard reduction potential E°′</b> at pH 7
    (biochemical standard state). To account for actual environmental conditions the
    <b>Nernst equation</b> is applied to each half-reaction:
  </p>
  <div style="background:#fff;border:1px solid #ccc;border-radius:4px;padding:6px 14px;
              font-family:monospace;font-size:14px;margin:6px 0">
    E = E°′ + (RT / n<sub>e</sub>F) · ln(Q)
  </div>
  <p style="margin:6px 0">
    where <i>Q</i> is the reaction quotient for the half reaction, <i>n</i><sub>e</sub> = electrons
    transferred, <i>R</i> = 8.314 J mol⁻¹ K⁻¹, <i>F</i> = 96 485 C mol⁻¹, <i>T</i> = 298 K.
  </p>
  <p style="margin:6px 0">
    For a complete reaction pairing a donor couple (D) with an acceptor couple (A) the free
    energy yield is:
  </p>
  <div style="background:#fff;border:1px solid #ccc;border-radius:4px;padding:6px 14px;
              font-family:monospace;font-size:14px;margin:6px 0">
    ΔG = −n<sub>e</sub> · F · ΔE &nbsp;&nbsp; where &nbsp;&nbsp; ΔE = E<sub>A</sub> − E<sub>D</sub>
  </div>
  <p style="margin:6px 0">
    Note that E<sub>D</sub> appears as negative in the overall reaction because the donor half
    reaction is the reverse of the reduction half reaction, as listed on the redox tower.
  </p>
  <p style="margin:6px 0">
    A reaction is <b style="color:#2E7D32">exergonic (feasible)</b> when ΔG &lt; 0, i.e. when
    the donor half reaction (written in the reduction direction, as it is done on the redox tower)
    sits above the acceptor reaction.
  </p>

  <b style="font-size:15px">Worked example — aerobic respiration (Glucose + O₂)</b>
  <p style="margin:6px 0">
    Two half-reactions are involved, listed in reduction form as on the tower. Q is the
    ratio of reactant to product activities for the reduction half-reaction, with all species
    written out (H₂O activity = 1). We use E°(SHE) as the unambiguous reference here, so that
    Q contains all species explicitly including protons. E°′ is simply E°(SHE) evaluated
    at pH=7 with unit activities for all other species.
  </p>
  <div style="background:#fff;border:1px solid #ccc;border-radius:4px;padding:6px 14px;
              font-family:monospace;font-size:14px;margin:6px 0;line-height:2.1">
    <b>Acceptor</b> &nbsp; O₂ + 4H⁺ + 4e⁻ → 2H₂O &nbsp;&nbsp; (n<sub>e</sub>=4, E°(SHE)=+1229 mV)<br>
    &nbsp;&nbsp;&nbsp; Q = p<sub>O₂</sub> · [H⁺]⁴<br>
    &nbsp;&nbsp;&nbsp; E = 1229 + (RT/4F) · ln( p<sub>O₂</sub> · [H⁺]⁴ )<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = 1229 + (RT/4F)·ln(p<sub>O₂</sub>) + (RT/4F)·4·ln([H⁺])<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = 1229 + (RT/4F)·ln(p<sub>O₂</sub>) − (RT/F)·ln(10)·pH &nbsp; [mV]<br>
    &nbsp; At pH=7, p<sub>O₂</sub>=1: E = 1229 − 7×25.7×ln(10) = 1229 − 414 = <b>+816 mV = E°′</b>
  </div>
  <div style="background:#fff;border:1px solid #ccc;border-radius:4px;padding:6px 14px;
              font-family:monospace;font-size:14px;margin:6px 0;line-height:2.1">
    <b>Donor</b> &nbsp;&nbsp; 6CO₂ + 24H⁺ + 24e⁻ → Glucose + 6H₂O &nbsp;&nbsp; (n<sub>e</sub>=24, E°(SHE)=−16 mV)<br>
    &nbsp;&nbsp;&nbsp; Q = p<sub>CO₂</sub>⁶ · [H⁺]²⁴ / [Glc]<br>
    &nbsp;&nbsp;&nbsp; E = −16 + (RT/24F) · ln( p<sub>CO₂</sub>⁶ · [H⁺]²⁴ / [Glc] )<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = −16 + (RT/4F)·ln(p<sub>CO₂</sub>) + (RT/24F)·24·ln([H⁺]) − (RT/24F)·ln([Glc])<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = −16 + (RT/4F)·ln(p<sub>CO₂</sub>) − (RT/F)·ln(10)·pH − (RT/24F)·ln([Glc]) &nbsp; [mV]<br>
    &nbsp; At pH=7, p<sub>CO₂</sub>=1, [Glc]=1M: E = −16 − 414 = <b>−430 mV = E°′</b>
  </div>
  <p style="margin:6px 0">
    where RT/F ≈ 25.7 mV at 25 °C.
    At standard conditions (pH 7, all species at 1 M or 1 bar):
  </p>
  <div style="background:#fff;border:1px solid #ccc;border-radius:4px;padding:6px 14px;
              font-family:monospace;font-size:14px;margin:6px 0;line-height:2.0">
    E(O₂/H₂O) = +816 mV &nbsp;|&nbsp; E(CO₂/Glc) = −430 mV<br>
    ΔE = (+816) − (−430) = +1246 mV<br>
    ΔG = −24 × 96.485 × 1.246 = <b>−2885 kJ/mol</b>
  </div>
  <p style="margin:6px 0">
    At trace oxygen (p<sub>O₂</sub> = 10⁻⁵ bar = 10 ppm, pH 7, all else standard).
    Atmospheric O₂ is ≈ 0.21 bar, so 10⁻⁵ bar is about 0.005% of atmospheric:
  </p>
  <div style="background:#fff;border:1px solid #ccc;border-radius:4px;padding:6px 14px;
              font-family:monospace;font-size:14px;margin:6px 0;line-height:2.0">
    Q = 10⁻⁵ · (10⁻⁷)⁴ = 10⁻³³ &nbsp;&nbsp; → &nbsp;&nbsp; E = 1229 + (25.7/4)·ln(10⁻³³) = 816 − 74 = <b>+742 mV</b><br>
    ΔE = (+742) − (−430) = +1172 mV &nbsp;|&nbsp; ΔG ≈ <b>−2714 kJ/mol</b> &nbsp;(still strongly exergonic)
  </div>
  <p style="margin:4px 0;font-size:13px;color:#666">
    <i>Note on primary fermentation arrows:</i> glucose fermentation (e.g. → butyrate, → lactate)
    is represented using the overall net reaction ΔG, not a direct electron transfer between two
    couples. The CO₂/glucose couple serves as the donor reference; the ΔG values shown are
    thermodynamically correct for the net reaction but the electron path runs through internal
    carriers (NAD⁺/NADH, ferredoxin) rather than directly between the displayed couples.
  </p>
</div>
"""

INTRO_HTML = """
<div style="max-width:1310px;font-size:14px;font-family:Georgia,serif;line-height:1.7;
            margin-bottom:10px;color:#222">
  This website features an interactive "redox tower", a static version of which is explained in
  Fig. 1 of Scarampi et al. (Citation to follow). In brief, the redox tower lists biologically
  relevant, theoretical redox half reactions on a reduction potential scale going from negative
  values at the top to positive values at the bottom. Each half reaction is listed in the direction
  of "reduction", i.e. acceptance of electrons. To explore possible microbial growth-supporting
  reactions, one should consider pairs of two redox half reactions: On the redox tower each half
  reaction can be reversed and coupled with another half-reaction that is below it on the tower.
  This way, the resulting overall redox reaction would have a positive reduction potential, which
  corresponds to a negative Gibbs free energy, fulfilling the requirement of 2nd Law of
  Thermodynamics. On the static redox tower, it is customary to list all half-reaction reduction
  potentials under biological standard conditions, i.e. 1M concentrations for all metabolites,
  except for hydrogen, listed at pH=7. The interactive version allows varying metabolite
  concentrations, such as hydrogen, thereby exploring shifts in thermodynamic feasibility of
  different microbial growth-supporting metabolisms. In the graphic below dashed lines show
  standard-condition potentials (pH 7); solid lines show actual potentials under the conditions
  selected.
</div>
"""

FOOTER_HTML = """
<div style="margin-top:16px;padding:8px 16px;font-size:13px;color:#666;
            font-family:Georgia,serif;max-width:1310px">
  Source code available on GitHub:
  <a href="https://github.com/cremerlab/redox-tower" target="_blank"
     style="color:#1565C0;text-decoration:none;font-weight:bold">cremerlab/redox-tower</a>.
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
