# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Universality of the power-law index of stellar flare energies
#
# This notebook tests whether the power-law index alpha of the flare energy
# distribution, dN/dE proportional to E^(-alpha), is the same for M, K and G
# dwarfs observed by Kepler. F dwarfs enter as a secondary extension.
# The estimators and the decision rules are described in the text of each
# section and in the accompanying Letter. This file is the source of the
# notebook and the `.ipynb` file is generated from it with jupytext.
#
# The notebook is self-contained. It downloads the catalogs from VizieR when it
# runs, so no data file is stored in the repository. On every run it records the
# catalog identifiers, the access date, the number of rows and a checksum of the
# downloaded data.
#
# Sections:
#
# 1. Setup, parameters and run log
# 2. Data download
# 3. Integrity checks
# 4. Sample selection, dwarf filter and temperature classes
# 5. Estimators: maximum likelihood, xmin by KS, models M1 to M4, goodness of fit
# 6. Validation on synthetic data
# 7. Baseline fits in energy as published
# 8. Luminosity-scaled energy, the primary variable
# 9. Primary fits, bootstrap by star and tests of equality
# 10. Baseline in energy with bootstrap by star
# 11. Sensitivity analysis
# 12. Heterogeneity between stars and reproduction of the published index
# 13. Verdict
# 14. Index in narrow bands of effective temperature
# 15. Additional checks: minimum duration, scale without the radius, bands without
#     the hottest band, permutation test of the whole distribution, Holm correction
# 16. Figures
# 17. Results macros and final checks

# %% [markdown]
# ## 1. Setup, parameters and run log

# %%
import hashlib
import os
import platform
import time
import zlib
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_DIR = Path(".")
CACHE_DIR = Path("data_cache")
USE_CACHE = True
MODE = os.environ.get("FLARE_MODE", "full")
SEED = 20261012
N_JOBS = max(1, min(4, os.cpu_count() or 1))


def seed_for(name):
    """Reproducible seed for a named task, independent of the order of execution."""
    return int(np.random.SeedSequence([SEED, zlib.crc32(name.encode())]).generate_state(1)[0])

# %% [markdown]
# ### Parameters
#
# All tunable choices live in one frozen object so that every sensitivity
# variant is a copy of the baseline with a few fields replaced. Values follow
# the specification of the analysis.

# %%
@dataclass(frozen=True)
class Config:
    # Dwarf filter
    logg_min: float | None = 4.0
    rad_factor: float | None = 1.5
    rad_bin_width: float = 200.0
    rad_bin_min_stars: int = 30
    strict_interval: bool = False
    # Temperature classes: boundaries between M and K, K and G, G and F,
    # and the upper limit of F, in kelvin
    bounds: tuple = (3850.0, 5300.0, 5900.0, 7200.0)
    m_split: float | None = None
    # Exclusions used in sensitivity variants
    drop_single_flare_stars: bool = False
    drop_neighbour_stars: bool = False
    drop_no_rotation: bool = False
    drop_nondsep: bool = False
    # Multiplicative factors on the energies of single classes, as pairs (class, factor)
    class_factor: tuple = ()


BASELINE = Config()
MIN_TAIL = 1000
GRID_DEX = 0.02
MIN_TAIL_SCAN = 1000
VAL_MIN_TAIL = 100
SIZES = {
    "full": {"recovery": 200, "trunc": 100, "gof_size": 100, "gof_synth": 100, "gof_power": 12, "lr": 100,
             "preview_gof": 100, "boot": 1000, "variant": 200, "gof_fixed": 200, "valid_eq": 100},
    "fast": {"recovery": 40, "trunc": 20, "gof_size": 20, "gof_synth": 40, "gof_power": 4, "lr": 20,
             "preview_gof": 20, "boot": 100, "variant": 40, "gof_fixed": 40, "valid_eq": 20},
}[MODE]
# Primary analysis on the luminosity-scaled energy
XMIN_L = 34.5
WINDOW_L = (34.5, 35.5)
WINDOWS_L_SENS = ((34.6, 35.4), (34.8, 35.6), (35.0, 36.0))
ALPHA_TEST = 0.05
DELTA_EQ = 0.10
N_MIN_STAR = 50
SIGMA_E_DEX = 0.2
PRIMARY_CLASSES = ("M", "K", "G")
EXTENDED_CLASSES = ("M", "K", "G", "F")
LN10 = np.log(10.0)
DATA_SOURCES = {
    "flares": "J/ApJS/241/29/table2",
    "stars": "J/ApJS/241/29/table1",
    "params": "J/ApJS/229/30/catalog",
}
EXPECTED_ROWS = {"flares": 162262, "stars": 3420, "params": 197096}

# %%
def package_versions():
    import astropy
    import astroquery
    import matplotlib
    import powerlaw
    import scipy

    mods = (np, pd, scipy, astropy, astroquery, powerlaw, matplotlib)
    out = {"python": platform.python_version()}
    out.update({m.__name__: m.__version__ for m in mods})
    return out


RUN = {
    "started_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
    "mode": MODE,
    "seed": SEED,
    "versions": package_versions(),
}
RUN

# %% [markdown]
# ### Results document
#
# The notebook writes `analysis_results.md` at the end of the run. Each section
# adds its own block of text through this small helper.

# %%
class ResultsDoc:
    def __init__(self):
        self.blocks = []

    def section(self, title):
        self.blocks.append(f"\n## {title}\n")

    def add(self, text):
        self.blocks.append(str(text).rstrip() + "\n")

    def table(self, df, floatfmt=None):
        self.blocks.append(df_to_markdown(df, floatfmt) + "\n")

    def write(self, path):
        Path(path).write_text("\n".join(self.blocks), encoding="utf-8")


def df_to_markdown(df, floatfmt=None):
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        cells = []
        for c in df.columns:
            v = row[c]
            if isinstance(v, (float, np.floating)):
                fmt = (floatfmt or {}).get(c, "{:.2f}")
                cells.append("nan" if np.isnan(v) else fmt.format(v))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


DOC = ResultsDoc()

# %% [markdown]
# ## 2. Data download
#
# Three tables are read from VizieR through astroquery:
#
# * Yang and Liu (2019), `J/ApJS/241/29`: table2 lists 162262 flares with
#   log energy and table1 lists the 3420 flaring stars.
# * Mathur et al. (2017), `J/ApJS/229/30`: stellar parameters of the DR25
#   Kepler catalog. VizieR serves the output parameters in the table named
#   `catalog`.
#
# Failures of the VizieR service stop the notebook with a clear message. There is
# no local copy of the data.

# %%
from astroquery.vizier import Vizier

PARAM_COLUMNS = ["KIC", "Teff", "E_Teff", "e_Teff", "log(g)", "E_log(g)",
                 "e_log(g)", "Rad", "E_Rad", "e_Rad", "Mod"]
COLUMN_REQUEST = {
    "flares": ["**"],
    "stars": ["KIC", "Nf", "Teff", "logg", "Lum", "Prot", "r_Prot", "NS"],
    "params": PARAM_COLUMNS,
}


def fetch_vizier(name, retries=3):
    cat = DATA_SOURCES[name]
    last = None
    for attempt in range(1, retries + 1):
        try:
            viz = Vizier(columns=COLUMN_REQUEST[name], row_limit=-1)
            viz.TIMEOUT = 900
            tables = viz.get_catalogs(cat)
            if cat not in tables.keys():
                raise RuntimeError(f"table {cat} not returned by VizieR")
            return tables[cat].to_pandas()
        except Exception as err:
            last = err
            time.sleep(5 * attempt)
    raise RuntimeError(
        f"Could not download {cat} from VizieR after {retries} attempts. "
        f"Last error: {last}. The notebook needs network access to VizieR."
    )


def canonical_checksum(df, sort_by):
    ordered = df.sort_values(sort_by, kind="mergesort").reset_index(drop=True)
    ordered = ordered[sorted(ordered.columns)]
    payload = ordered.to_csv(index=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_table(name):
    cat = DATA_SOURCES[name]
    cache = CACHE_DIR / (name + ".pkl")
    source = "VizieR"
    df = None
    if USE_CACHE and cache.exists():
        try:
            df = pd.read_pickle(cache)
            source = "local cache of an earlier VizieR download"
        except Exception:
            df = None
    if df is None:
        df = fetch_vizier(name)
        if USE_CACHE:
            CACHE_DIR.mkdir(exist_ok=True)
            df.to_pickle(cache)
    sort_by = ["KIC", "Begin"] if name == "flares" else ["KIC"]
    info = {
        "catalog": cat,
        "rows": len(df),
        "accessed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        "source": source,
        "sha256": canonical_checksum(df, sort_by),
    }
    return df, info


flares_raw, info_flares = load_table("flares")
stars_raw, info_stars = load_table("stars")
params_raw, info_params = load_table("params")
SOURCES_LOG = pd.DataFrame([info_flares, info_stars, info_params])
SOURCES_LOG

# %% [markdown]
# The local cache only speeds up repeated runs and is ignored by git. The
# checksum is computed on the data content, not on the cache file.

# %% [markdown]
# ## 3. Integrity checks
#
# Each check prints a result
# and stores it for the results document. Failed expectations raise a warning in
# the text and do not stop the run, except for structural problems that would
# make later steps meaningless.

# %%
CHECKS = []


def record_check(name, ok, detail):
    CHECKS.append({"check": name, "status": "ok" if ok else "ATTENTION", "detail": detail})


for key, df in (("flares", flares_raw), ("stars", stars_raw), ("params", params_raw)):
    n = len(df)
    record_check(f"row count {key}", n == EXPECTED_ROWS[key],
                 f"{n} rows, expected {EXPECTED_ROWS[key]}")

logE = flares_raw["logE"].to_numpy(dtype=float)
record_check("logE finite", bool(np.isfinite(logE).all()),
             f"{int((~np.isfinite(logE)).sum())} non finite values")
record_check("logE range", 29.3 <= np.nanmin(logE) and np.nanmax(logE) <= 38.6,
             f"min {np.nanmin(logE):.4f}, max {np.nanmax(logE):.4f}, catalog range 29.3 to 38.6")

kic_f = set(flares_raw["KIC"])
kic_s = set(stars_raw["KIC"])
record_check("every flare KIC is in the star table", kic_f <= kic_s,
             f"{len(kic_f - kic_s)} flare KIC values missing from table1")

count = flares_raw.groupby("KIC").size()
nf = stars_raw.set_index("KIC")["Nf"]
mismatch = int((count.reindex(nf.index) != nf).sum())
record_check("Nf equals the number of flares per star", mismatch == 0,
             f"{mismatch} stars with inconsistent counts")

record_check("duplicate stars in table1", not stars_raw["KIC"].duplicated().any(),
             f"{int(stars_raw['KIC'].duplicated().sum())} duplicated KIC values")
record_check("duplicate stars in the Mathur table", not params_raw["KIC"].duplicated().any(),
             f"{int(params_raw['KIC'].duplicated().sum())} duplicated KIC values")

dup_rows = int(flares_raw.duplicated(["KIC", "Begin", "End", "logE"]).sum())
record_check("duplicate flares", True,
             f"{dup_rows} exact duplicate flares (same KIC, begin, end and logE), reported and kept")

ref_star = stars_raw[stars_raw["KIC"] == 8935655]
ref_lum = float(ref_star["Lum"].iloc[0]) if len(ref_star) else float("nan")
record_check("flare to bolometric luminosity ratio reflects the 2020 erratum",
              1e-4 < ref_lum < 3e-4,
              f"KIC 8935655 has Lum {ref_lum:.2e}. The corrected value is about 1.9e-4 and the value before the "
              "erratum is 5.2e-7. Lum is not used in the analysis")

# %%
merged_all = stars_raw.merge(params_raw, on="KIC", how="left", suffixes=("_yang", "_mathur"))
matched = merged_all["Teff_mathur"].notna()
n_unmatched = int((~matched).sum())
record_check("stars matched in the Mathur table", True,
             f"{int(matched.sum())} of {len(merged_all)} matched, {n_unmatched} unmatched and excluded "
             "(unclassified targets are absent from the Mathur table)")

both = merged_all[matched]
dT = (both["Teff_yang"] - both["Teff_mathur"]).to_numpy(dtype=float)
dg = (both["logg"] - both["log(g)"]).to_numpy(dtype=float)
record_check("Teff of table1 against Teff of Mathur", True,
             f"median difference {np.median(dT):.1f} K, standard deviation {np.std(dT):.1f} K, "
             f"{int((np.abs(dT) > 50).sum())} stars differ by more than 50 K")
record_check("log g of table1 against log g of Mathur", True,
             f"median difference {np.median(dg):.3f}, standard deviation {np.std(dg):.4f}")
prov = both["Mod"].value_counts().to_dict()
record_check("provenance of the Mathur parameters among flaring stars", True,
             "counts " + ", ".join(f"{k} {v}" for k, v in prov.items()) +
             ". BTSL means polynomial fits to BT-Settl models for cool dwarfs, and those stars have no distance or extinction")

CHECKS_DF = pd.DataFrame(CHECKS)
CHECKS_DF

# %% [markdown]
# ## 4. Sample selection, dwarf filter and temperature classes
# 5. Estimators: maximum likelihood, xmin by KS, models M1 to M4, goodness of fit
# 6. Validation on synthetic data
# 7. Preview fits of the baseline classes (not final results)
#
# The selection is a function of the configuration, so that every sensitivity
# variant reuses it. The steps and their meaning are:
#
# 1. All flares of table2.
# 2. Flares of stars with a counterpart in the Mathur table.
# 3. Flares of stars that pass the dwarf filter. The filter keeps a star when
#    its surface gravity is at least `logg_min` and its radius is at most
#    `rad_factor` times the median radius of the stars in its 200 K bin of
#    effective temperature. Bins with fewer than 30 stars are merged with the
#    neighbouring bin that has fewer stars. The median is computed on the
#    matched flaring stars before the filter, and it describes the main
#    sequence because most Kepler targets belong to it. The radius is used only
#    for selecting dwarfs.
# 4. Flares of stars with effective temperature inside one of the classes.
# 5. Flares with a finite energy.
#
# Class limits follow Pecaut and Mamajek (2013), Table 5, rounded to 3850, 5300
# and 5900 K. The F class ends at 7200 K, the temperature of an F0V star.

# %%
def radius_reference(teff, rad, width, min_stars):
    edges = np.floor(teff / width).astype(int)
    labels = pd.Series(edges)
    counts = labels.value_counts().sort_index()
    groups = {int(k): [int(k)] for k in counts.index}
    sizes = {int(k): int(v) for k, v in counts.items()}
    while True:
        small = [k for k in groups if sizes[k] < min_stars]
        if not small or len(groups) == 1:
            break
        k = min(small, key=lambda g: sizes[g])
        keys = sorted(groups)
        pos = keys.index(k)
        nbrs = [keys[p] for p in (pos - 1, pos + 1) if 0 <= p < len(keys)]
        tgt = min(nbrs, key=lambda g: sizes[g])
        groups[tgt] = groups[tgt] + groups.pop(k)
        sizes[tgt] += sizes.pop(k)
    lookup = {}
    for g, members in groups.items():
        sel = np.isin(edges, members)
        med = float(np.median(rad[sel]))
        for m in members:
            lookup[m] = med
    return np.array([lookup[e] for e in edges]), groups


def assign_class(teff, cfg):
    mk, kg, gf, ftop = cfg.bounds
    cls = np.full(len(teff), "", dtype=object)
    cls[teff < mk] = "M"
    if cfg.m_split is not None:
        cls[teff < cfg.m_split] = "Mfc"
        cls[(teff >= cfg.m_split) & (teff < mk)] = "Mhi"
    cls[(teff >= mk) & (teff < kg)] = "K"
    cls[(teff >= kg) & (teff < gf)] = "G"
    cls[(teff >= gf) & (teff < ftop)] = "F"
    return cls


def build_star_table(cfg):
    st = merged_all.copy()
    st = st[st["Teff_mathur"].notna()].copy()
    st = st.rename(columns={"Teff_mathur": "Teff", "log(g)": "logg_m", "Rad": "Rad"})
    ref, groups = radius_reference(st["Teff"].to_numpy(float), st["Rad"].to_numpy(float),
                                   cfg.rad_bin_width, cfg.rad_bin_min_stars)
    st["rad_ref"] = ref
    keep = np.ones(len(st), dtype=bool)
    if cfg.logg_min is not None:
        g = st["logg_m"].to_numpy(float)
        if cfg.strict_interval:
            g = g - st["e_log(g)"].fillna(0).to_numpy(float)
        keep &= g >= cfg.logg_min
    if cfg.rad_factor is not None:
        r = st["Rad"].to_numpy(float)
        if cfg.strict_interval:
            r = r + st["E_Rad"].fillna(0).to_numpy(float)
        keep &= r <= cfg.rad_factor * st["rad_ref"].to_numpy(float)
    st["dwarf"] = keep
    st["cls"] = assign_class(st["Teff"].to_numpy(float), cfg)
    excl = np.zeros(len(st), dtype=bool)
    if cfg.drop_single_flare_stars:
        excl |= st["Nf"].to_numpy() == 1
    if cfg.drop_neighbour_stars:
        excl |= st["NS"].astype(str).str.strip().to_numpy() == "p"
    if cfg.drop_no_rotation:
        excl |= st["Prot"].to_numpy(float) < 0
    if cfg.drop_nondsep:
        excl |= st["Mod"].astype(str).str.strip().to_numpy() != "DSEP"
    st["excluded"] = excl
    return st, groups


def select_sample(cfg=BASELINE):
    st, groups = build_star_table(cfg)
    fl = flares_raw[["KIC", "Q", "Begin", "End", "logE"]].merge(
        st[["KIC", "Teff", "logg_m", "Rad", "E_Rad", "e_Rad", "E_Teff", "e_Teff", "Mod", "Nf", "NS",
            "Prot", "dwarf", "cls", "excluded"]],
        on="KIC", how="left")
    in_mathur = fl["Teff"].notna()
    finite = np.isfinite(fl["logE"].to_numpy(float))
    s_dwarf = in_mathur & fl["dwarf"].fillna(False).astype(bool)
    s_excl = s_dwarf & ~fl["excluded"].fillna(True).astype(bool)
    s_class = s_excl & (fl["cls"].fillna("") != "")
    s_fin = s_class & finite
    steps = [
        ("1 all flares in table2", np.ones(len(fl), bool)),
        ("2 star has a counterpart in the Mathur table", in_mathur.to_numpy()),
        ("3 star passes the dwarf filter", s_dwarf.to_numpy()),
        ("3b star not removed by the exclusion options", s_excl.to_numpy()),
        ("4 Teff inside a class", s_class.to_numpy()),
        ("5 finite energy, analysed sample", s_fin.to_numpy()),
    ]
    rows = []
    for label, mask in steps:
        rows.append({"step": label, "stars": int(fl.loc[mask, "KIC"].nunique()),
                     "flares": int(mask.sum())})
    funnel = pd.DataFrame(rows)
    sample = fl[s_fin.to_numpy()].copy()
    for cname, factor in cfg.class_factor:
        sample.loc[sample["cls"] == cname, "logE"] += np.log10(factor)
    sample["E"] = 10.0 ** sample["logE"].to_numpy(float)
    sample = sample.reset_index(drop=True)
    return sample, funnel, st, groups


sample0, funnel0, stars0, rad_groups0 = select_sample(BASELINE)
funnel0

# %% [markdown]
# ### Classes in the baseline sample
#
# The table lists stars, flares and the energy range of each class. The number of
# flares is an upper bound for the number in the power-law tail, which is only
# known after xmin is estimated.

# %%
def class_summary(sample, order=("M", "K", "G", "F")):
    rows = []
    for c in order:
        s = sample[sample["cls"] == c]
        if len(s) == 0:
            continue
        rows.append({
            "class": c,
            "stars": int(s["KIC"].nunique()),
            "flares": len(s),
            "logE min": float(s["logE"].min()),
            "logE median": float(s["logE"].median()),
            "logE max": float(s["logE"].max()),
            "Teff min (K)": float(s["Teff"].min()),
            "Teff max (K)": float(s["Teff"].max()),
        })
    return pd.DataFrame(rows)


summary0 = class_summary(sample0)
summary0

# %%
fmt = {"logE min": "{:.2f}", "logE median": "{:.2f}", "logE max": "{:.2f}",
       "Teff min (K)": "{:.0f}", "Teff max (K)": "{:.0f}"}
radius_bins = pd.DataFrame({
    "bin members (lower edge, K)": [", ".join(str(m * int(BASELINE.rad_bin_width)) for m in v)
                                    for v in rad_groups0.values()],
})

DOC.section("1. Run header")
DOC.add(f"Run started (UTC): {RUN['started_utc']}")
DOC.add(f"Mode: {RUN['mode']}")
DOC.add(f"Seed: {RUN['seed']}")
DOC.add("Package versions: " + ", ".join(f"{k} {v}" for k, v in RUN["versions"].items()))
DOC.add("Data sources:")
DOC.table(SOURCES_LOG[["catalog", "rows", "accessed_utc", "source", "sha256"]])

DOC.section("2. Data integrity")
DOC.table(CHECKS_DF)

DOC.section("3. Selection funnel (baseline configuration)")
DOC.add(f"Configuration: log g at least {BASELINE.logg_min}, radius at most {BASELINE.rad_factor} "
        f"times the median of its {BASELINE.rad_bin_width:.0f} K bin, class limits {BASELINE.bounds} K.")
DOC.table(funnel0, {})
DOC.add("Classes of the analysed sample:")
DOC.table(summary0, fmt)

# %% [markdown]
# ### Dwarf filter and class variants
#
# A quick table of how the sample changes with the filter variants used later as
# sensitivity checks. Nothing is fitted here.

# %%
variants = {
    "baseline": BASELINE,
    "no filter": replace(BASELINE, logg_min=None, rad_factor=None),
    "log g 4.0 only": replace(BASELINE, rad_factor=None),
    "log g 4.2 only": replace(BASELINE, logg_min=4.2, rad_factor=None),
    "radius factor 1.25": replace(BASELINE, rad_factor=1.25),
    "radius factor 2.0": replace(BASELINE, rad_factor=2.0),
    "strict 1 sigma interval": replace(BASELINE, strict_interval=True),
    "Yang class limits": replace(BASELINE, bounds=(4000.0, 5000.0, 6000.0, 7500.0)),
    "limits minus 200 K": replace(BASELINE, bounds=tuple(b - 200.0 for b in BASELINE.bounds)),
    "limits plus 200 K": replace(BASELINE, bounds=tuple(b + 200.0 for b in BASELINE.bounds)),
    "Pecaut limits unrounded": replace(BASELINE, bounds=(3850.0, 5280.0, 5920.0, 7200.0)),
    "drop single flare stars": replace(BASELINE, drop_single_flare_stars=True),
    "drop neighbour stars": replace(BASELINE, drop_neighbour_stars=True),
    "drop stars without rotation": replace(BASELINE, drop_no_rotation=True),
    "drop non DSEP stars": replace(BASELINE, drop_nondsep=True),
}
rows = []
for label, cfg in variants.items():
    smp, _, _, _ = select_sample(cfg)
    row = {"variant": label}
    for c in ("M", "K", "G", "F"):
        s = smp[smp["cls"] == c]
        row[f"{c} stars"] = int(s["KIC"].nunique())
        row[f"{c} flares"] = len(s)
    rows.append(row)
VARIANT_TABLE = pd.DataFrame(rows)
VARIANT_TABLE

# %%
DOC.section("4. Sample size under the filter and class variants")
DOC.add("Stars and flares per class. These counts are upper bounds for the tail used in the fits.")
DOC.table(VARIANT_TABLE, {})

split_cfg = replace(BASELINE, m_split=3200.0)
smp_split, _, _, _ = select_sample(split_cfg)
split_summary = class_summary(smp_split, order=("Mhi", "Mfc"))
DOC.add("Subdivision of the M class at 3200 K (Mhi between 3200 and 3850 K, Mfc below 3200 K):")
DOC.table(split_summary, fmt)

# %% [markdown]
# ## 5. Estimators
#
# All estimators work on the energies in erg, sorted in ascending order.
#
# * **M1**, pure power law above xmin. The maximum likelihood estimator is
#   alpha = 1 + n / sum(ln(x / xmin)) (Clauset et al. 2009, Eq. 3.1).
# * **M2**, power law with a hard upper truncation at xmax, taken as the largest
#   energy of the sample. The density is p(x) = (alpha minus 1) x^(-alpha) /
#   (xmin^(1-alpha) minus xmax^(1-alpha)). The likelihood depends on the data only
#   through the mean of ln(x / xmin), so alpha is the root of a monotonic
#   equation in one variable. alpha can fall below 1 under truncation.
# * **M3**, power law with an exponential cutoff, nested around M1.
# * **M4**, lognormal restricted to x >= xmin.
#
# The value of xmin minimizes the Kolmogorov-Smirnov distance between the data
# above xmin and the fitted model (Clauset et al. 2009, Eq. 3.9). Candidates lie
# on a grid of 0.02 dex and tails must keep at least 1000 points for the real
# data, which is the size from which Clauset et al. (2009, p. 13) report good
# results. Validation on small synthetic samples lowers this limit to 100 points.
# The grid is checked against a full scan of all distinct values in section 6.

# %%
from scipy import integrate, optimize, special, stats


_SHARED = {}


def parallel_map(func, items, n_jobs, shared=None):
    """Map with forked processes. Large arrays travel through the module level dictionary."""
    if shared is not None:
        _SHARED.clear()
        _SHARED.update(shared)
    items = list(items)
    if n_jobs is None or n_jobs <= 1 or len(items) < 2:
        return [func(a) for a in items]
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    ctx = mp.get_context("fork")
    chunk = max(1, len(items) // (4 * n_jobs))
    with ProcessPoolExecutor(max_workers=n_jobs, mp_context=ctx) as ex:
        return list(ex.map(func, items, chunksize=chunk))


def _h(u, L):
    z = u * L
    if abs(z) < 1e-3:
        return L * (0.5 - z / 12.0 + z ** 3 / 720.0)
    return 1.0 / u - L / np.expm1(z)


def solve_u_truncated(s, L):
    """Maximum likelihood u = alpha - 1 for a power law truncated at ln(xmax / xmin) = L.

    s is the mean of ln(x / xmin). The function h is decreasing from L to 0, so the root is unique.
    """
    if not (L > 0):
        return np.nan
    lo, hi = -60.0 / L, 60.0 / L
    if _h(lo, L) - s < 0:
        return lo
    if _h(hi, L) - s > 0:
        return hi
    return optimize.brentq(lambda u: _h(u, L) - s, lo, hi, xtol=1e-12, rtol=1e-12)


def _candidate_indices(xs, grid_dex, min_tail):
    N = len(xs)
    last = N - min_tail
    if last < 1:
        raise ValueError("sample too small for the requested minimum tail")
    if grid_dex is None:
        first = np.r_[True, xs[1:] != xs[:-1]]
        idx = np.flatnonzero(first)
    else:
        lo = np.floor(np.log10(xs[0]) / grid_dex) * grid_dex
        g = 10.0 ** np.arange(lo, np.log10(xs[last]) + grid_dex, grid_dex)
        idx = np.unique(np.searchsorted(xs, g, side="left"))
    return idx[(idx <= last) & (idx < N)]


def _eval_candidate(lx, cs, i, model, AR):
    N = len(lx)
    n = N - i
    t = lx[i:] - lx[i]
    S = cs[i] - n * lx[i]
    if not (S > 0):
        return None
    if model == "M1":
        u = n / S
        F = -np.expm1(-u * t)
    else:
        L = lx[-1] - lx[i]
        if not (L > 0):
            return None
        u = solve_u_truncated(S / n, L)
        if abs(u * L) < 1e-9:
            F = t / L
        else:
            F = np.expm1(-u * t) / np.expm1(-u * L)
    j = AR[:n]
    D = max(np.max((j + 1.0) / n - F), np.max(F - j / n))
    return 1.0 + u, D, n


def scan_xmin(xs, model="M1", grid_dex=GRID_DEX, min_tail=MIN_TAIL_SCAN):
    """Choose xmin by minimizing the KS distance. xs must be sorted in ascending order."""
    lx = np.log(xs)
    cs = np.cumsum(lx[::-1])[::-1]
    AR = np.arange(len(xs), dtype=float)
    best = None
    for i in _candidate_indices(xs, grid_dex, min_tail):
        r = _eval_candidate(lx, cs, i, model, AR)
        if r is None:
            continue
        if best is None or r[1] < best[2]:
            best = (i, r[0], r[1], r[2])
    if best is None:
        raise RuntimeError("no valid xmin candidate")
    i, alpha, D, n = best
    return {"xmin": float(xs[i]), "alpha": float(alpha), "ks": float(D), "n_tail": int(n),
            "xmax": float(xs[-1]) if model == "M2" else np.inf, "model": model}


def fit_fixed_xmin(xs, xmin, model="M1"):
    lx = np.log(xs)
    cs = np.cumsum(lx[::-1])[::-1]
    AR = np.arange(len(xs), dtype=float)
    i = int(np.searchsorted(xs, xmin, side="left"))
    if i >= len(xs) - 1:
        raise RuntimeError("fixed xmin lies above the data")
    r = _eval_candidate(lx, cs, i, model, AR)
    if r is None:
        raise RuntimeError("invalid fixed xmin")
    return {"xmin": float(xs[i]), "alpha": float(r[0]), "ks": float(r[1]), "n_tail": int(r[2]),
            "xmax": float(xs[-1]) if model == "M2" else np.inf, "model": model}


def hill_sigma(alpha, n):
    """Standard error of Clauset et al. 2009, Eq. 3.2, kept only as a reference."""
    return (alpha - 1.0) / np.sqrt(n)


# %% [markdown]
# ### Alternative models and likelihood ratios
#
# Per observation log likelihoods feed the likelihood ratio test of Clauset et
# al. (2009), Appendix C. For models that are not nested the sign of R is
# trusted when the Vuong p value, p = |erfc(R / sqrt(2 n sigma^2))|, is below 0.1.
# For nested models, here M1 inside M3, R follows a chi-squared distribution
# with one degree of freedom for 2R (Appendix C.1). The same convention is used
# by the powerlaw package. Because the cutoff parameter sits on the boundary of
# its range under the null, this p value is conservative.

# %%
def loglik_m1(x, xmin, alpha):
    return np.log(alpha - 1.0) - np.log(xmin) - alpha * np.log(x / xmin)


def loglik_m2(x, xmin, xmax, alpha):
    u = alpha - 1.0
    L = np.log(xmax / xmin)
    norm = 1.0 / L if abs(u * L) < 1e-9 else u / (-np.expm1(-u * L))
    return np.log(norm) - np.log(xmin) - alpha * np.log(x / xmin)


def _logC_m3(alpha, kappa):
    vmax = np.log(800.0 / kappa)
    val, _ = integrate.quad(lambda v: np.exp((1.0 - alpha) * v - kappa * np.exp(v)), 0.0, vmax,
                            epsabs=0.0, epsrel=1e-10, limit=400)
    return np.log(val)


def fit_m3(x, xmin, alpha0):
    """Power law with exponential cutoff above xmin, parameters alpha and ln(kappa) with kappa = lambda xmin."""
    r = x / xmin
    n = len(r)
    sr, sx = np.sum(np.log(r)), np.sum(r)

    def nll(th):
        a, lk = th
        return a * sr + np.exp(lk) * sx + n * _logC_m3(a, np.exp(lk))

    best = None
    for lk0 in (np.log(1.0 / r.max()) - 1.0, np.log(1.0 / r.max()) - 4.0, -20.0):
        res = optimize.minimize(nll, [alpha0, lk0], method="L-BFGS-B",
                                bounds=[(0.2, 6.0), (-40.0, 2.0)])
        if best is None or res.fun < best.fun:
            best = res
    a, lk = best.x
    k = float(np.exp(lk))
    return {"alpha": float(a), "kappa": k, "lambda": k / xmin}


def loglik_m3(x, xmin, alpha, kappa):
    r = x / xmin
    return -alpha * np.log(r) - kappa * r - _logC_m3(alpha, kappa) - np.log(xmin)


def fit_m4(x, xmin):
    """Lognormal restricted to x >= xmin, parameters mu and ln(sigma) of ln x."""
    y = np.log(x)
    n = len(y)
    lxm = np.log(xmin)
    sy, syy = y.sum(), (y ** 2).sum()

    def nll(th):
        mu, ls = th
        sg = np.exp(ls)
        ss = syy - 2 * mu * sy + n * mu ** 2
        return (n * ls + ss / (2 * sg ** 2) + sy + 0.5 * n * np.log(2 * np.pi)
                + n * special.log_ndtr((mu - lxm) / sg))

    best = None
    for mu0, s0 in ((y.mean(), y.std() + 1e-6), (lxm, 2.0), (lxm + 3.0, 3.0)):
        res = optimize.minimize(nll, [mu0, np.log(s0)], method="Nelder-Mead",
                                options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
        if best is None or res.fun < best.fun:
            best = res
    mu, ls = best.x
    return {"mu": float(mu), "sigma": float(np.exp(ls))}


def loglik_m4(x, xmin, mu, sigma):
    y = np.log(x)
    return (-y - np.log(sigma) - 0.5 * np.log(2 * np.pi) - (y - mu) ** 2 / (2 * sigma ** 2)
            - special.log_ndtr((mu - np.log(xmin)) / sigma))


def compare_loglik(ll_a, ll_b, nested=False):
    """Likelihood ratio R of model a against model b and its p value. R positive favours a."""
    n = float(len(ll_a))
    d = ll_a - ll_b
    R = float(d.sum())
    var = float(np.mean((d - d.mean()) ** 2))
    if nested:
        p = float(stats.chi2.sf(abs(2.0 * R), 1))
    else:
        p = float(abs(special.erfc(abs(R) / np.sqrt(2.0 * n * var)))) if var > 0 else 1.0
    return {"R": R, "p": p, "R_norm": float(R / np.sqrt(n * var)) if var > 0 else 0.0}


def likelihood_comparisons(xs, fit):
    """Compare the fitted power law with M3 and M4 above the same xmin."""
    tail = xs[xs >= fit["xmin"]]
    if fit["model"] == "M2":
        ll_pl = loglik_m2(tail, fit["xmin"], fit["xmax"], fit["alpha"])
    else:
        ll_pl = loglik_m1(tail, fit["xmin"], fit["alpha"])
    out = {}
    m4 = fit_m4(tail, fit["xmin"])
    ll4 = loglik_m4(tail, fit["xmin"], m4["mu"], m4["sigma"])
    out["lognormal"] = {**compare_loglik(ll_pl, ll4), **m4}
    a1 = fit["alpha"] if fit["model"] == "M1" else scan_xmin_value_m1(tail, fit["xmin"])
    m3 = fit_m3(tail, fit["xmin"], a1)
    ll3 = loglik_m3(tail, fit["xmin"], m3["alpha"], m3["kappa"])
    ll1 = loglik_m1(tail, fit["xmin"], a1)
    out["cutoff"] = {**compare_loglik(ll3, ll1, nested=True), **m3}
    out["cutoff"]["R"] = -out["cutoff"]["R"]
    return out


def scan_xmin_value_m1(tail, xmin):
    return 1.0 + len(tail) / np.sum(np.log(tail / xmin))


# %% [markdown]
# ### Random numbers and goodness of fit
#
# Power-law numbers come from the transformation method (Clauset et al. 2009,
# Appendix D). The goodness of fit follows the semi-parametric procedure of
# section 4.1 of the same paper. Each synthetic sample keeps the fitted power law
# above xmin and resamples the observed data below xmin, then it is fitted again
# including xmin. The p value is the fraction of synthetic KS distances at least
# as large as the observed one.

# %%
def draw_tail(rng, size, alpha, xmin, xmax=np.inf):
    u = rng.random(size)
    uu = alpha - 1.0
    if np.isfinite(xmax):
        L = np.log(xmax / xmin)
        if abs(uu * L) < 1e-9:
            return xmin * np.exp(u * L)
        top = -np.expm1(-uu * L)
        return xmin * (1.0 - u * top) ** (-1.0 / uu)
    return xmin * (1.0 - u) ** (-1.0 / uu)


def gof_one(args):
    seed, = args
    sh = _SHARED
    rng = np.random.default_rng(seed)
    xs, fit, model = sh["xs"], sh["fit"], sh["model"]
    N = len(xs)
    below = xs[xs < fit["xmin"]]
    m_pl = rng.binomial(N, fit["n_tail"] / N)
    xmax = fit["xmax"] if model == "M2" else np.inf
    pl = draw_tail(rng, m_pl, fit["alpha"], fit["xmin"], xmax)
    body = rng.choice(below, size=N - m_pl, replace=True) if len(below) else np.empty(0)
    syn = np.sort(np.concatenate([pl, body]))
    return scan_xmin(syn, model, sh["grid_dex"], sh["min_tail"])["ks"]


def gof_pvalue(xs, fit, n_synth, seed, grid_dex=GRID_DEX, min_tail=MIN_TAIL_SCAN, n_jobs=1):
    ss = np.random.SeedSequence(seed).spawn(n_synth)
    seeds = [(int(s.generate_state(1)[0]),) for s in ss]
    ks = parallel_map(gof_one, seeds, n_jobs,
                      shared={"xs": xs, "fit": fit, "model": fit["model"],
                              "grid_dex": grid_dex, "min_tail": min_tail})
    ks = np.array(ks)
    return float(np.mean(ks >= fit["ks"])), ks


# %% [markdown]
# ## 6. Validation on synthetic data
#
# Every estimator is tested on data with known truth before it touches the
# catalog. The synthetic samples mimic a flare catalog: a fraction of the points
# follows a power law above xmin and the rest forms an arbitrary body below xmin.
# Tolerances are set before looking at the results and each test reports ok or
# ATTENTION.

# %%
VALIDATION = []


def record_validation(name, ok, detail):
    VALIDATION.append({"test": name, "status": "ok" if ok else "ATTENTION", "detail": detail})


def make_sample(rng, n, alpha, xmin, tail_frac=0.6, xmax=np.inf):
    nt = int(n * tail_frac)
    tail = draw_tail(rng, nt, alpha, xmin, xmax)
    body = xmin * 10.0 ** (-0.8 * np.abs(rng.standard_normal(n - nt)))
    return np.sort(np.concatenate([body, tail]))


# %% [markdown]
# ### V1. Agreement with the powerlaw package and with a full scan of xmin
#
# On synthetic samples, alpha at the same xmin must agree with the powerlaw
# package of Alstott et al. (2014) to numerical precision. The grid of 0.02 dex
# must give the same alpha as a scan over every distinct value, within the
# Hill standard error. The same grid check is repeated on real data in section 7.

# %%
import warnings
import powerlaw

rng = np.random.default_rng(seed_for("v1"))
rows = []
for k in range(5):
    x = make_sample(rng, 3000, 2.0 + 0.1 * k, 1e33)
    g = scan_xmin(x, "M1", min_tail=VAL_MIN_TAIL)
    f = scan_xmin(x, "M1", grid_dex=None, min_tail=VAL_MIN_TAIL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pl = powerlaw.Fit(x[x >= g["xmin"]], xmin=g["xmin"], verbose=False)
    rows.append({"alpha grid": g["alpha"], "alpha full scan": f["alpha"],
                 "alpha powerlaw at grid xmin": float(pl.power_law.alpha),
                 "log10 xmin grid": np.log10(g["xmin"]), "log10 xmin full": np.log10(f["xmin"]),
                 "Hill sigma": hill_sigma(g["alpha"], g["n_tail"])})
V1 = pd.DataFrame(rows)
ok_pl = bool((np.abs(V1["alpha grid"] - V1["alpha powerlaw at grid xmin"]) < 1e-6).all())
ok_grid = bool((np.abs(V1["alpha grid"] - V1["alpha full scan"]) < V1["Hill sigma"]).all())
record_validation("alpha agrees with powerlaw at the same xmin", ok_pl,
                  f"largest difference {np.abs(V1['alpha grid'] - V1['alpha powerlaw at grid xmin']).max():.2e}")
record_validation("grid scan agrees with the full scan", ok_grid,
                  f"largest difference {np.abs(V1['alpha grid'] - V1['alpha full scan']).max():.4f}, "
                  f"smallest Hill sigma {V1['Hill sigma'].min():.4f}")
V1

# %% [markdown]
# ### V2. Recovery of alpha and xmin
#
# True alpha values of 1.7, 2.0 and 2.4, xmin of 1e33 erg, 5000 points, 60 per
# cent of them in the tail. The table gives the mean bias, the scatter, the mean
# Hill error and the fraction of fits within one Hill error of the truth. The
# Hill error ignores the uncertainty on xmin, so a coverage below 68 per cent is
# expected and is the reason for using the bootstrap in the final analysis.

# %%
rows = []
for alpha in (1.7, 2.0, 2.4):
    rng = np.random.default_rng(seed_for(f"v2-{alpha}"))
    est, sig, xm = [], [], []
    for _ in range(SIZES["recovery"]):
        x = make_sample(rng, 5000, alpha, 1e33)
        f = scan_xmin(x, "M1", min_tail=VAL_MIN_TAIL)
        est.append(f["alpha"])
        sig.append(hill_sigma(f["alpha"], f["n_tail"]))
        xm.append(np.log10(f["xmin"]))
    est, sig, xm = np.array(est), np.array(sig), np.array(xm)
    rows.append({"true alpha": alpha, "mean bias": est.mean() - alpha, "scatter": est.std(),
                 "mean Hill error": sig.mean(), "coverage of 1 Hill error": np.mean(np.abs(est - alpha) < sig),
                 "mean log10 xmin minus 33": xm.mean() - 33.0})
V2 = pd.DataFrame(rows)
record_validation("recovery of alpha without bias", bool((np.abs(V2["mean bias"]) < 0.03).all()),
                  f"largest mean bias {np.abs(V2['mean bias']).max():.4f} over {SIZES['recovery']} samples per alpha")
V2

# %% [markdown]
# ### V3. Truncation
#
# A power law with alpha 1.8 truncated at 1e36 erg, 20000 points. A pure power
# law fit (M1) must be biased towards steeper values and the truncated fit (M2)
# must recover the truth.

# %%
rng = np.random.default_rng(seed_for("v3"))
e1, e2 = [], []
for _ in range(SIZES["trunc"]):
    x = make_sample(rng, 20000, 1.8, 1e33, 0.7, xmax=1e36)
    e1.append(scan_xmin(x, "M1", min_tail=VAL_MIN_TAIL)["alpha"])
    e2.append(scan_xmin(x, "M2", min_tail=VAL_MIN_TAIL)["alpha"])
V3 = pd.DataFrame({"model": ["M1 pure", "M2 truncated"],
                   "mean alpha": [np.mean(e1), np.mean(e2)],
                   "scatter": [np.std(e1), np.std(e2)],
                   "bias": [np.mean(e1) - 1.8, np.mean(e2) - 1.8]})
record_validation("truncated fit recovers alpha and the pure fit is biased",
                  abs(V3["bias"][1]) < 0.03 and V3["bias"][0] > abs(V3["bias"][1]),
                  f"bias M1 {V3['bias'][0]:.4f}, bias M2 {V3['bias'][1]:.4f}")
V3

# %% [markdown]
# ### V4. Goodness of fit: size and power
#
# Under a true power law the rejection rate at p below 0.1 should stay close to
# its nominal level. Lognormal data should be rejected, but the power of any such
# test depends on how many points reach the tail and is limited when xmin comes
# from the KS distance, because that choice picks the part of the tail that looks
# most like a power law. Power is reported for two sample sizes. The requirement
# is that power exceeds the size of the test and grows with the sample size.
# Clauset et al. (2009, p. 20 and 21) show the same dependence on sample size
# for the likelihood ratio. Each p value uses synthetic samples drawn with the semi-parametric
# procedure.

# %%
def gof_rate(make, n_rep, seed_name, **kw):
    rng = np.random.default_rng(seed_for(seed_name))
    out = []
    for k in range(n_rep):
        x = make(rng)
        f = scan_xmin(x, "M1", min_tail=VAL_MIN_TAIL)
        p, _ = gof_pvalue(x, f, SIZES["gof_synth"], seed_for(f"{seed_name}-gof-{k}"),
                          min_tail=VAL_MIN_TAIL, n_jobs=N_JOBS)
        out.append((p, f["n_tail"]))
    return np.array(out)


def lognormal_sample(n):
    return lambda rng: np.sort(np.exp(rng.normal(np.log(1e33), 1.6, n)))


size = gof_rate(lambda rng: make_sample(rng, 2000, 2.0, 1e33), SIZES["gof_size"], "v4-size")
rej01, rej05 = float(np.mean(size[:, 0] < 0.1)), float(np.mean(size[:, 0] < 0.05))
record_validation("goodness of fit keeps its size under a true power law", 0.02 <= rej01 <= 0.25,
                  f"{len(size)} samples, rejection rate {rej01:.2f} at 0.1 and {rej05:.2f} at 0.05, mean p {size[:, 0].mean():.2f}")
power_rows = []
for n in (30000, 200000):
    r = gof_rate(lognormal_sample(n), SIZES["gof_power"], f"v4-power-{n}")
    power_rows.append({"lognormal sample size": n, "samples": len(r),
                       "mean points in the tail": float(r[:, 1].mean()),
                       "rejection rate at 0.1": float(np.mean(r[:, 0] < 0.1))})
V4 = pd.DataFrame(power_rows)
record_validation("goodness of fit power against the lognormal is above its size and grows with the sample size",
                  bool(V4["rejection rate at 0.1"].iloc[-1] > rej01 and
                       V4["rejection rate at 0.1"].iloc[-1] >= V4["rejection rate at 0.1"].iloc[0]),
                  f"rejection rate {V4['rejection rate at 0.1'].iloc[0]:.2f} at {V4['lognormal sample size'].iloc[0]} points "
                  f"and {V4['rejection rate at 0.1'].iloc[-1]:.2f} at {V4['lognormal sample size'].iloc[-1]} points, "
                  f"size {rej01:.2f}. Limited power is a documented limitation of the test at the KS xmin")
V4

# %% [markdown]
# ### V5. Likelihood ratio tests
#
# Power-law data must not be separated from the lognormal alternative with a
# reliable sign in favour of the lognormal, lognormal data must favour the
# lognormal with growing power as the sample size increases, and data with an exponential cutoff must be
# separated from the pure power law by the nested test. The lognormal fit can
# degenerate on power-law data, with a very large sigma, because a wide lognormal
# looks like a power law over a limited range. That case is expected and shows as
# p near 1.

# %%
rng = np.random.default_rng(seed_for("v5"))
n_rep = SIZES["lr"]
wrong = 0
for _ in range(n_rep):
    x = make_sample(rng, 6000, 2.0, 1e33)
    f = scan_xmin(x, "M1", min_tail=VAL_MIN_TAIL)
    c = likelihood_comparisons(x, f)["lognormal"]
    wrong += int(c["R"] < 0 and c["p"] < 0.1)
ln_rows = []
for n in (30000, 200000):
    reps = n_rep if n == 30000 else max(10, n_rep // 4)
    ok = 0
    for _ in range(reps):
        y = np.sort(np.exp(rng.normal(np.log(1e33), 1.6, n)))
        f = scan_xmin(y, "M1", min_tail=VAL_MIN_TAIL)
        c = likelihood_comparisons(y, f)["lognormal"]
        ok += int(c["R"] < 0 and c["p"] < 0.1)
    ln_rows.append({"lognormal sample size": n, "samples": reps, "lognormal favoured with p below 0.1": ok / reps})
V5 = pd.DataFrame(ln_rows)
cut_ok = 0
n_cut = max(5, n_rep // 5)
for _ in range(n_cut):
    tail = draw_tail(rng, 20000, 1.5, 1e33)
    tail = tail[rng.random(len(tail)) < np.exp(-tail / 3e35)]
    z = np.sort(np.concatenate([1e33 * 10.0 ** (-0.8 * np.abs(rng.standard_normal(10000))), tail]))
    f = scan_xmin(z, "M1", min_tail=VAL_MIN_TAIL)
    c = likelihood_comparisons(z, f)["cutoff"]
    cut_ok += int(c["R"] < 0 and c["p"] < 0.1)
record_validation("power-law data are not wrongly assigned to the lognormal", wrong / n_rep <= 0.2,
                  f"{wrong} of {n_rep} samples with R negative and p below 0.1 (expected near 0.1 or lower)")
record_validation("likelihood ratio power against the lognormal grows with the sample size",
                  bool(V5["lognormal favoured with p below 0.1"].iloc[-1] > 0.1 and
                       V5["lognormal favoured with p below 0.1"].iloc[-1] >= V5["lognormal favoured with p below 0.1"].iloc[0]),
                  f"fraction {V5['lognormal favoured with p below 0.1'].iloc[0]:.2f} at {V5['lognormal sample size'].iloc[0]} points "
                  f"and {V5['lognormal favoured with p below 0.1'].iloc[-1]:.2f} at {V5['lognormal sample size'].iloc[-1]} points")
record_validation("cutoff data are separated from the pure power law", cut_ok / n_cut >= 0.9,
                  f"{cut_ok} of {n_cut} samples with the nested test p below 0.1")
V5

# %% [markdown]
# ### V6. Timing
#
# Time of one full fit on a sample of the size of the largest class, used to
# plan the bootstrap and the calibration of the equality test.

# %%
rng = np.random.default_rng(seed_for("v6"))
big = make_sample(rng, 60000, 2.0, 1e33, 0.5)
TIMES = {}
for m in ("M1", "M2"):
    t0 = time.perf_counter()
    for _ in range(3):
        scan_xmin(big, m)
    TIMES[m] = (time.perf_counter() - t0) / 3.0
pd.DataFrame({"model": list(TIMES), "seconds per fit, 60000 points": list(TIMES.values())})

# %%
VAL_DF = pd.DataFrame(VALIDATION)
DOC.section("5. Validation on synthetic data")
DOC.add("Each estimator is tested on samples with known truth. Sample sizes depend on the run mode.")
DOC.table(VAL_DF)
DOC.add("Recovery of alpha and xmin (V2):")
DOC.table(V2, {"mean bias": "{:.4f}", "scatter": "{:.4f}", "mean Hill error": "{:.4f}",
               "coverage of 1 Hill error": "{:.2f}", "mean log10 xmin minus 33": "{:.3f}"})
DOC.add("Power of the goodness of fit test for lognormal data (V4):")
DOC.table(V4, {"lognormal sample size": "{:.0f}", "samples": "{:.0f}", "mean points in the tail": "{:.0f}",
               "rejection rate at 0.1": "{:.2f}"})
DOC.add("Likelihood ratio against the lognormal for lognormal data (V5):")
DOC.table(V5, {"lognormal sample size": "{:.0f}", "samples": "{:.0f}", "lognormal favoured with p below 0.1": "{:.2f}"})
DOC.add("Truncated data (V3):")
DOC.table(V3, {"mean alpha": "{:.4f}", "scatter": "{:.4f}", "bias": "{:.4f}"})
VAL_DF

# %% [markdown]
# ## 7. Baseline fits in energy as published
#
# This section applies the estimators to the baseline sample in energy as
# published, with xmin from the KS distance. These are the baseline fits used to
# compare with the literature. The errors from the bootstrap by star are added in
# section 10. The primary analysis of the universality test uses the
# luminosity-scaled energy of section 8.

# %%
def fit_class(x, model, n_gof, seed, n_jobs):
    xs = np.sort(np.asarray(x, dtype=float))
    f = scan_xmin(xs, model)
    out = dict(f)
    out["n_total"] = len(xs)
    out["hill_sigma"] = hill_sigma(f["alpha"], f["n_tail"])
    if n_gof:
        out["gof_p"], _ = gof_pvalue(xs, f, n_gof, seed, n_jobs=n_jobs)
    out["comparisons"] = likelihood_comparisons(xs, f)
    return xs, out


PREVIEW = {}
timing = []
for c in ("M", "K", "G", "F"):
    xc = sample0.loc[sample0["cls"] == c, "E"].to_numpy()
    for model in ("M1", "M2"):
        t0 = time.perf_counter()
        xs, res = fit_class(xc, model, SIZES["preview_gof"], seed_for(f"prev-{c}-{model}"), N_JOBS)
        res["seconds"] = time.perf_counter() - t0
        PREVIEW[(c, model)] = res
        timing.append({"class": c, "model": model, "points": len(xs), "seconds incl. GoF": res["seconds"]})
pd.DataFrame(timing)

# %%
rows = []
for (c, model), r in PREVIEW.items():
    comp = r["comparisons"]
    rows.append({
        "class": c, "model": model, "points": r["n_total"], "tail": r["n_tail"],
        "log10 xmin": np.log10(r["xmin"]), "alpha": r["alpha"], "Hill error": r["hill_sigma"],
        "KS D": r["ks"], "GoF p": r["gof_p"],
        "R vs lognormal": comp["lognormal"]["R"], "p (Vuong)": comp["lognormal"]["p"],
        "R vs cutoff": comp["cutoff"]["R"], "p (nested)": comp["cutoff"]["p"],
    })
PREVIEW_DF = pd.DataFrame(rows)
PREVIEW_DF

# %% [markdown]
# ### Checks on the real data
#
# The grid of 0.02 dex is compared with a full scan of distinct values on a random
# subsample of each class, and alpha at the same xmin is compared with the
# powerlaw package. The subsample keeps the full scan affordable.

# %%
rows = []
for c in ("M", "K", "G", "F"):
    xc = sample0.loc[sample0["cls"] == c, "E"].to_numpy()
    rng = np.random.default_rng(seed_for(f"grid-check-{c}"))
    sub = np.sort(rng.choice(xc, size=min(8000, len(xc)), replace=False))
    g = scan_xmin(sub, "M1")
    t0 = time.perf_counter()
    f = scan_xmin(sub, "M1", grid_dex=None)
    t_full = time.perf_counter() - t0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pl = powerlaw.Fit(sub[sub >= g["xmin"]], xmin=g["xmin"], verbose=False)
    rows.append({"class": c, "points": len(sub), "alpha grid": g["alpha"], "alpha full scan": f["alpha"],
                 "log10 xmin grid": np.log10(g["xmin"]), "log10 xmin full": np.log10(f["xmin"]),
                 "alpha powerlaw": float(pl.power_law.alpha), "Hill error": hill_sigma(g["alpha"], g["n_tail"]),
                 "full scan seconds": t_full})
GRID_CHECK = pd.DataFrame(rows)
record_validation("grid and full scan agree on real data",
                  bool((np.abs(GRID_CHECK["alpha grid"] - GRID_CHECK["alpha full scan"]) < GRID_CHECK["Hill error"]).all()),
                  "largest difference in alpha " + f"{np.abs(GRID_CHECK['alpha grid'] - GRID_CHECK['alpha full scan']).max():.4f}")
GRID_CHECK

# %% [markdown]
# ### Alpha as a function of xmin
#
# Fixing xmin on a grid shows how alpha depends on the lower limit of the fit
# range. A true power law gives a plateau. The table also gives the number of
# flares above each xmin.

# %%
rows = []
for c in ("M", "K", "G", "F"):
    xc = np.sort(sample0.loc[sample0["cls"] == c, "E"].to_numpy())
    for lxm in np.arange(33.0, 35.51, 0.5):
        try:
            r = fit_fixed_xmin(xc, 10.0 ** lxm, "M1")
        except RuntimeError:
            continue
        if r["n_tail"] >= 300:
            rows.append({"class": c, "log10 xmin": lxm, "alpha": r["alpha"], "tail": r["n_tail"],
                         "Hill error": hill_sigma(r["alpha"], r["n_tail"]), "KS D": r["ks"]})
PROFILE = pd.DataFrame(rows)
PROFILE_WIDE = PROFILE.pivot(index="log10 xmin", columns="class", values="alpha").reindex(columns=["M", "K", "G", "F"])
PROFILE_WIDE

# %%
DOC.section("6. Baseline fits in energy as published")
DOC.add("Maximum likelihood fits with xmin from the KS distance. M1 is the pure power law and M2 has an upper "
        "truncation at the largest energy of the class. The Hill error ignores the uncertainty on xmin and is "
        "only a reference. The goodness of fit p uses " + str(SIZES["preview_gof"]) + " synthetic samples, "
        "so its precision is limited. Positive R favours the power law over the alternative.")
DOC.table(PREVIEW_DF, {"log10 xmin": "{:.2f}", "alpha": "{:.3f}", "Hill error": "{:.3f}", "KS D": "{:.4f}",
                       "GoF p": "{:.2f}", "R vs lognormal": "{:.1f}", "p (Vuong)": "{:.3f}",
                       "R vs cutoff": "{:.1f}", "p (nested)": "{:.3g}"})
DOC.add("Alpha of the pure power law (M1) for fixed values of xmin, with at least 300 flares above xmin. A plateau would "
        "indicate a power law over that range:")
DOC.table(PROFILE, {"log10 xmin": "{:.1f}", "alpha": "{:.3f}", "Hill error": "{:.3f}", "KS D": "{:.4f}"})
DOC.add("Grid and full scan of xmin on random subsamples of 8000 flares, and alpha from the powerlaw package at the same xmin:")
DOC.table(GRID_CHECK, {"alpha grid": "{:.3f}", "alpha full scan": "{:.3f}", "log10 xmin grid": "{:.2f}",
                       "log10 xmin full": "{:.2f}", "alpha powerlaw": "{:.3f}", "Hill error": "{:.3f}",
                       "full scan seconds": "{:.1f}"})

# %% [markdown]
# ### Run time planning
#
# One fit of the largest class takes the time given below, so the bootstrap of
# the next section needs B fits per class and the calibration of the equality
# test needs B_NULL fits for every class.

# %%
per_fit = {k: v for k, v in TIMES.items()}
plan = pd.DataFrame({
    "task": ["bootstrap, one class, B = 1000", "all four classes, B = 1000", "equality test calibration, B_NULL = 500, four classes"],
    "estimated minutes (serial)": [1000 * per_fit["M2"] / 60, 4000 * per_fit["M2"] / 60, 2000 * per_fit["M2"] / 60],
})
plan["estimated minutes (parallel)"] = plan["estimated minutes (serial)"] / N_JOBS
DOC.add("Run time planning from the fit time on 60000 points (" + ", ".join(f"{k} {v:.2f} s" for k, v in TIMES.items()) + "):")
DOC.table(plan, {"estimated minutes (serial)": "{:.1f}", "estimated minutes (parallel)": "{:.1f}"})
plan


# %% [markdown]
# ## 8. Luminosity-scaled energy, the primary variable
#
# The rollover of the flare energy distribution depends on the stellar class when
# energy is measured in erg, and it lines up when each energy is divided by the
# luminosity of its star. A detection threshold is a relative amplitude, and the
# energy of a flare is the equivalent duration times the stellar luminosity in the
# Kepler band (Davenport 2016, p. 5), so the completeness limit in erg grows with
# the luminosity. The scaled energy is
#
# E_L = E / (L / L_sun), with L / L_sun = (R / R_sun)^2 (Teff / 5772 K)^4.
#
# The radius and the temperature come from the Mathur catalog. The bolometric
# luminosity is not the Kepler band luminosity, so a factor that depends on Teff
# through the instrument response remains. The radius has a typical uncertainty
# of 27 per cent, which enters as noise in E_L and is tested in section 11.
# E_L is the primary variable of the test of universality.

# %%
TEFF_SUN = 5772.0


def add_scaled_energy(sample):
    out = sample.copy()
    lum = out["Rad"].to_numpy(float) ** 2 * (out["Teff"].to_numpy(float) / TEFF_SUN) ** 4
    out["logL"] = np.log10(lum)
    out["logEL"] = out["logE"] - out["logL"]
    return out


def histogram_peak(values, width):
    edges = np.arange(values.min(), values.max() + width, width)
    h, e = np.histogram(values, bins=edges)
    return float(e[int(h.argmax())] + width / 2.0)


sampleL = add_scaled_energy(sample0)
rows = []
for c in EXTENDED_CLASSES:
    d = sampleL[sampleL["cls"] == c]
    rows.append({"class": c, "median log L": d["logL"].median(),
                 "peak log E, 0.1 dex": histogram_peak(d["logE"], 0.1),
                 "peak log E, 0.3 dex": histogram_peak(d["logE"], 0.3),
                 "peak log E_L, 0.1 dex": histogram_peak(d["logEL"], 0.1),
                 "peak log E_L, 0.3 dex": histogram_peak(d["logEL"], 0.3)})
PEAKS = pd.DataFrame(rows)
SPREAD = {k: float(PEAKS[k].max() - PEAKS[k].min()) for k in PEAKS.columns[2:]}
PEAKS

# %%
DOC.section("7. Luminosity-scaled energy")
DOC.add("Peak of the differential flare energy distribution in energy E and in the luminosity-scaled energy E_L. "
        "The spread between classes is " + ", ".join(f"{k}: {v:.2f} dex" for k, v in SPREAD.items()) + ".")
DOC.table(PEAKS, {c: "{:.2f}" for c in PEAKS.columns[1:]})

# %% [markdown]
# ## 9. Primary fits, bootstrap by star and tests of equality
#
# **Data model.** For every class the flares above the common lower limit
# `XMIN_L` of E_L are summarized by three sums per star: the number of flares,
# the sum of ln(E_L / xmin) and the sum of E_L / xmin. Every estimator below is a
# function of these three totals, so a bootstrap that resamples whole stars only
# needs to resample the per star sums. Stars without flares above the limit stay
# in the resampling units.
#
# **Estimators.**
#
# * Power law with an exponential cutoff above xmin (M3), the primary model. It
#   has an index alpha_c and a cutoff scale lambda_c for every class.
# * Truncated power law in the window `WINDOW_L` (M2 with fixed limits), a
#   secondary estimate.
# * Pure power law above xmin (M1), for reference.
#
# **Test of equality.** The classes are independent samples. The statistic is
# the weighted sum of squared deviations of the indices from their weighted mean,
# with weights from the bootstrap variances. Its null distribution comes from the
# bootstrap by star with the observed indices subtracted, which respects the
# correlation between flares of the same star. The p value is
# (1 + number of bootstrap values at least as large as the observed one) /
# (B + 1). The chi-squared p value with the same statistic is also given as a
# reference. Pairwise differences use percentile intervals and the correction of
# Holm. The nested likelihood ratio is reported as well, but
# it treats flares as independent and is therefore anticonservative.

# %%
def class_units(sample, c, col, lo, hi=None):
    """Per star sums (number, sum of ln r, sum of r) with r = x / xmin, for flares in [lo, hi)."""
    d = sample.loc[sample["cls"] == c, ["KIC", col]]
    kic = d["KIC"].to_numpy()
    v = d[col].to_numpy(float)
    units, inv = np.unique(kic, return_inverse=True)
    sel = v >= lo
    if hi is not None:
        sel &= v < hi
    ln_r = (v[sel] - lo) * LN10
    k = len(units)
    n = np.bincount(inv[sel], minlength=k).astype(float)
    s1 = np.bincount(inv[sel], weights=ln_r, minlength=k)
    s2 = np.bincount(inv[sel], weights=np.exp(ln_r), minlength=k)
    return np.column_stack([n, s1, s2])


def fit_m3_sums(n, s1, s2, starts=None):
    """Power law with exponential cutoff from the sums. Returns alpha, kappa = lambda xmin and the log likelihood."""
    def nll(th):
        a, lk = th
        return a * s1 + np.exp(lk) * s2 + n * _logC_m3(a, np.exp(lk))

    if starts is None:
        a1 = 1.0 + n / s1
        mean_r = s2 / n
        starts = [(a1, np.log(1.0 / (20.0 * mean_r))), (a1, np.log(1.0 / (200.0 * mean_r))), (a1, -20.0)]
    best = None
    for a0, lk0 in starts:
        res = optimize.minimize(nll, [a0, lk0], method="L-BFGS-B", bounds=[(0.2, 6.0), (-40.0, 2.0)])
        if best is None or res.fun < best.fun:
            best = res
    a, lk = best.x
    return {"alpha": float(a), "kappa": float(np.exp(lk)), "ll": -float(best.fun)}


def m3_boot_one(args):
    n, s1, s2, a0, lk0 = args
    if n < 30:
        return (np.nan, np.nan)
    try:
        r = fit_m3_sums(n, s1, s2, starts=[(a0, lk0)])
        return (r["alpha"], r["kappa"])
    except Exception:
        return (np.nan, np.nan)


def window_alpha(n, s1, L):
    if n < 30 or not (s1 > 0):
        return np.nan
    return 1.0 + solve_u_truncated(s1 / n, L)


def primary_class_fit(sample, c, lo, hi, B, seed, n_jobs=N_JOBS, var="logEL"):
    """Point estimates and bootstrap by star of the M3, window M2 and M1 indices for one class."""
    Uf = class_units(sample, c, var, lo)
    Uw = class_units(sample, c, var, lo, hi)
    K = len(Uf)
    tot, totw = Uf.sum(axis=0), Uw.sum(axis=0)
    L = (hi - lo) * LN10
    m3 = fit_m3_sums(*tot)
    est = {"cls": c, "stars": K, "tail": int(tot[0]), "window": int(totw[0]),
           "alpha_c": m3["alpha"], "kappa": m3["kappa"], "alpha_1": 1.0 + tot[0] / tot[1],
           "alpha_w": window_alpha(totw[0], totw[1], L)}
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(K, np.full(K, 1.0 / K), size=B)
    Sb, Sw = counts @ Uf, counts @ Uw
    lk0 = np.log(m3["kappa"])
    args = [(Sb[i, 0], Sb[i, 1], Sb[i, 2], m3["alpha"], lk0) for i in range(B)]
    boot_m3 = np.array(parallel_map(m3_boot_one, args, n_jobs))
    est["boot_alpha_c"] = boot_m3[:, 0]
    est["boot_kappa"] = boot_m3[:, 1]
    est["boot_alpha_1"] = 1.0 + Sb[:, 0] / Sb[:, 1]
    est["boot_alpha_w"] = np.array([window_alpha(Sw[i, 0], Sw[i, 1], L) for i in range(B)])
    est["hill_se"] = hill_sigma(est["alpha_1"], est["tail"])
    est["units_tot"] = tot
    return est


def holm_adjust(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def equality_test(alpha_hat, alpha_boot, names, margin=DELTA_EQ, pairs=True):
    """Bootstrap test of equal indices for independent classes, with pairwise differences."""
    alpha_hat = np.asarray(alpha_hat, float)
    ab = np.asarray(alpha_boot, float)
    ab = ab[np.isfinite(ab).all(axis=1)]
    B, C = ab.shape
    w = 1.0 / ab.var(axis=0, ddof=1)

    def stat(dev):
        m = (dev * w).sum(axis=-1, keepdims=True) / w.sum()
        return (w * (dev - m) ** 2).sum(axis=-1)

    t_obs = float(stat(alpha_hat))
    t_boot = stat(ab - alpha_hat)
    out = {"B": B, "stat": t_obs, "p_boot": float((1 + np.sum(t_boot >= t_obs)) / (B + 1)),
           "p_chi2": float(stats.chi2.sf(t_obs, C - 1)),
           "alpha_common": float((alpha_hat * w).sum() / w.sum())}
    rows = []
    for i in range(C if pairs else 0):
        for j in range(i + 1, C):
            d_hat = alpha_hat[i] - alpha_hat[j]
            d_b = ab[:, i] - ab[:, j]
            p = (1 + np.sum(np.abs(d_b - d_hat) >= abs(d_hat))) / (B + 1)
            ci95 = np.percentile(d_b, [2.5, 97.5])
            ci90 = np.percentile(d_b, [5.0, 95.0])
            rows.append({"pair": f"{names[i]} minus {names[j]}", "difference": d_hat, "se": float(d_b.std(ddof=1)),
                         "ci95 low": ci95[0], "ci95 high": ci95[1], "p": float(p),
                         "equivalent within margin": bool(ci90[0] > -margin and ci90[1] < margin)})
    if pairs:
        pair_table = pd.DataFrame(rows)
        pair_table["p Holm"] = holm_adjust(pair_table["p"].to_numpy())
        out["pairs"] = pair_table
        out["max_abs_diff"] = float(pair_table["difference"].abs().max())
    else:
        out["max_abs_diff"] = float(np.ptp(alpha_hat))
    return out


def verdict_label(test, level=ALPHA_TEST, margin=DELTA_EQ):
    if test["p_boot"] >= level:
        return "H0 not rejected"
    if test["max_abs_diff"] < margin:
        return "H0 rejected, differences below the margin"
    return "H0 rejected"


# %%
def table_fit(res_by_class, key, boot_key, classes, label):
    rows = []
    for c in classes:
        r = res_by_class[c]
        b = r[boot_key][np.isfinite(r[boot_key])]
        lo16, hi84 = np.percentile(b, [16, 84])
        lo99, hi99 = np.percentile(b, [0.5, 99.5])
        rows.append({"class": c, "stars": r["stars"], "flares": r["tail"] if key != "alpha_w" else r["window"],
                     label: r[key], "se": float(b.std(ddof=1)), "16 to 84 low": lo16, "16 to 84 high": hi84,
                     "99 low": lo99, "99 high": hi99})
    return pd.DataFrame(rows)


PRIMARY = {}
for c in EXTENDED_CLASSES:
    PRIMARY[c] = primary_class_fit(sampleL, c, WINDOW_L[0], WINDOW_L[1], SIZES["boot"], seed_for(f"primary-{c}"))

fmt_fit = {"alpha": "{:.3f}", "se": "{:.3f}", "16 to 84 low": "{:.3f}", "16 to 84 high": "{:.3f}",
           "99 low": "{:.3f}", "99 high": "{:.3f}"}
TAB_C = table_fit(PRIMARY, "alpha_c", "boot_alpha_c", EXTENDED_CLASSES, "alpha")
TAB_W = table_fit(PRIMARY, "alpha_w", "boot_alpha_w", EXTENDED_CLASSES, "alpha")
TAB_1 = table_fit(PRIMARY, "alpha_1", "boot_alpha_1", EXTENDED_CLASSES, "alpha")
TAB_C

# %%
def run_equality(res, key, classes_sets):
    out = {}
    for label, cl in classes_sets.items():
        ah = [res[c][key] for c in cl]
        ab = np.column_stack([res[c]["boot_" + key] for c in cl])
        out[label] = equality_test(ah, ab, list(cl))
    return out


SETS = {"M, K, G": PRIMARY_CLASSES, "M, K, G, F": EXTENDED_CLASSES}
EQ_C = run_equality(PRIMARY, "alpha_c", SETS)
EQ_W = run_equality(PRIMARY, "alpha_w", SETS)
EQ_1 = run_equality(PRIMARY, "alpha_1", SETS)
summary_rows = []
for name, eq in (("cutoff index (primary)", EQ_C), ("window truncated index", EQ_W), ("pure power law index", EQ_1)):
    for label, t in eq.items():
        summary_rows.append({"estimator": name, "classes": label, "statistic": t["stat"], "p bootstrap": t["p_boot"],
                             "p chi-squared": t["p_chi2"], "common index": t["alpha_common"],
                             "largest difference": t["max_abs_diff"], "verdict": verdict_label(t)})
EQ_TABLE = pd.DataFrame(summary_rows)
EQ_TABLE

# %%
EQ_C["M, K, G"]["pairs"]

# %% [markdown]
# ### Comparison with the value of the fractal-diffusive model
#
# H0c states that the index equals 3/2 (Aschwanden 2014, p. 8, Eqs. 25 to 27).
# The hypothesis is rejected for a class when 3/2 falls outside the 99 per cent
# bootstrap interval. No statement about a universality class follows from it.

# %%
rows = []
for c in EXTENDED_CLASSES:
    for name, key, bkey in (("cutoff", "alpha_c", "boot_alpha_c"), ("window", "alpha_w", "boot_alpha_w")):
        b = PRIMARY[c][bkey]
        b = b[np.isfinite(b)]
        lo99, hi99 = np.percentile(b, [0.5, 99.5])
        rows.append({"class": c, "estimator": name, "alpha": PRIMARY[c][key], "alpha minus 3/2": PRIMARY[c][key] - 1.5,
                     "99 low": lo99, "99 high": hi99, "3/2 rejected": bool(not (lo99 <= 1.5 <= hi99))})
H0C = pd.DataFrame(rows)
H0C

# %% [markdown]
# ### Nested likelihood ratio and design effect
#
# The nested likelihood ratio compares separate indices with one common index
# for the cutoff model, each class keeping its own cutoff. Because the flares of
# one star are correlated, the chi-squared reference is anticonservative, and the
# design effect below shows by how much. It is the ratio of the bootstrap standard
# error by star to the Hill error that assumes independent flares.

# %%
def fit_m3_common(sums, alpha0, lk0):
    C = len(sums)

    def nll(th):
        a = th[0]
        tot = 0.0
        for k, (n, s1, s2) in enumerate(sums):
            lk = th[1 + k]
            tot += a * s1 + np.exp(lk) * s2 + n * _logC_m3(a, np.exp(lk))
        return tot

    res = optimize.minimize(nll, [alpha0] + list(lk0), method="L-BFGS-B",
                            bounds=[(0.2, 6.0)] + [(-40.0, 2.0)] * C)
    return {"alpha": float(res.x[0]), "ll": -float(res.fun)}


def lr_common(res, classes):
    sums = [tuple(res[c]["units_tot"]) for c in classes]
    ll_sep = sum(fit_m3_sums(*s)["ll"] for s in sums)
    a0 = float(np.mean([res[c]["alpha_c"] for c in classes]))
    lk0 = [np.log(res[c]["kappa"]) for c in classes]
    com = fit_m3_common(sums, a0, lk0)
    lam = 2.0 * (ll_sep - com["ll"])
    return {"Lambda": lam, "p naive": float(stats.chi2.sf(max(lam, 0.0), len(classes) - 1)), "alpha common": com["alpha"]}


LR = {k: lr_common(PRIMARY, v) for k, v in SETS.items()}
DEFF = pd.DataFrame([{"class": c, "Hill error, independent flares": PRIMARY[c]["hill_se"],
                      "bootstrap error by star, pure power law":
                          float(np.nanstd(PRIMARY[c]["boot_alpha_1"], ddof=1))} for c in EXTENDED_CLASSES])
DEFF["design effect, standard error ratio"] = DEFF.iloc[:, 2] / DEFF.iloc[:, 1]
pd.DataFrame([{"classes": k, **v} for k, v in LR.items()]), DEFF

# %% [markdown]
# ### Goodness of fit of the cutoff model at a fixed lower limit
#
# With xmin fixed there is nothing to refit except the two parameters. The
# synthetic samples are drawn from the fitted cutoff model with the observed
# number of flares, then refitted, and the p value is the fraction of synthetic
# KS distances at least as large as the observed one. The procedure treats flares
# as independent.

# %%
def m3_table(alpha, kappa, npts=6000):
    vmax = min(np.log(800.0 / kappa), 60.0)
    v = np.linspace(0.0, vmax, npts)
    f = np.exp((1.0 - alpha) * v - kappa * np.exp(v))
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(v))])
    return v, cdf / cdf[-1]


def ks_against_table(v_obs, v, cdf):
    x = np.sort(v_obs)
    n = len(x)
    F = np.interp(x, v, cdf)
    j = np.arange(n, dtype=float)
    return float(max(np.max((j + 1.0) / n - F), np.max(F - j / n)))


def gof_m3_one(args):
    seed, = args
    sh = _SHARED
    rng = np.random.default_rng(seed)
    v, cdf = m3_table(sh["alpha"], sh["kappa"])
    sim = np.interp(rng.random(sh["n"]), cdf, v)
    s1, s2 = sim.sum(), np.exp(sim).sum()
    fit = fit_m3_sums(sh["n"], s1, s2, starts=[(sh["alpha"], np.log(sh["kappa"]))])
    v2, cdf2 = m3_table(fit["alpha"], fit["kappa"])
    return ks_against_table(sim, v2, cdf2)


def gof_m3(sample, c, lo, n_synth, seed, var="logEL"):
    d = sample.loc[(sample["cls"] == c) & (sample[var] >= lo), var].to_numpy(float)
    v_obs = (d - lo) * LN10
    fit = fit_m3_sums(len(v_obs), v_obs.sum(), np.exp(v_obs).sum())
    v, cdf = m3_table(fit["alpha"], fit["kappa"])
    d_obs = ks_against_table(v_obs, v, cdf)
    ss = np.random.SeedSequence(seed).spawn(n_synth)
    seeds = [(int(s.generate_state(1)[0]),) for s in ss]
    ks = parallel_map(gof_m3_one, seeds, N_JOBS,
                      shared={"alpha": fit["alpha"], "kappa": fit["kappa"], "n": len(v_obs)})
    return {"KS D": d_obs, "p": float(np.mean(np.array(ks) >= d_obs)), "n": len(v_obs)}


GOF_C = pd.DataFrame([{"class": c, **gof_m3(sampleL, c, XMIN_L, SIZES["gof_fixed"], seed_for(f"gofc-{c}"))}
                      for c in EXTENDED_CLASSES])
GOF_C

# %%
DOC.section("8. Primary analysis in the luminosity-scaled energy")
DOC.add(f"Lower limit of E_L: log10 {XMIN_L} (all classes), window of the truncated fit: log10 {WINDOW_L[0]} to "
        f"{WINDOW_L[1]}. Bootstrap by star with B = {SIZES['boot']}. Errors are bootstrap standard errors, and the "
        "intervals are percentile intervals.")
DOC.add("Index of the power law with exponential cutoff (primary model):")
DOC.table(TAB_C, fmt_fit)
DOC.add("Index of the truncated power law in the window:")
DOC.table(TAB_W, fmt_fit)
DOC.add("Index of the pure power law above the lower limit:")
DOC.table(TAB_1, fmt_fit)
DOC.add("Tests of equality of the index between classes:")
DOC.table(EQ_TABLE, {"statistic": "{:.2f}", "p bootstrap": "{:.4f}", "p chi-squared": "{:.3g}",
                     "common index": "{:.3f}", "largest difference": "{:.3f}"})
DOC.add("Pairwise differences of the cutoff index, classes M, K and G:")
DOC.table(EQ_C["M, K, G"]["pairs"], {"difference": "{:.3f}", "se": "{:.3f}", "ci95 low": "{:.3f}",
                                      "ci95 high": "{:.3f}", "p": "{:.4f}", "p Holm": "{:.4f}"})
DOC.add("Comparison with 3/2:")
DOC.table(H0C, {"alpha": "{:.3f}", "alpha minus 3/2": "{:.3f}", "99 low": "{:.3f}", "99 high": "{:.3f}"})
DOC.add("Nested likelihood ratio, separate against common index for the cutoff model, and the design effect of the "
        "clustering of flares in stars:")
DOC.table(pd.DataFrame([{"classes": k, **v} for k, v in LR.items()]),
          {"Lambda": "{:.2f}", "p naive": "{:.3g}", "alpha common": "{:.3f}"})
DOC.table(DEFF, {c: "{:.3f}" for c in DEFF.columns[1:]})
DOC.add("Goodness of fit of the cutoff model at the fixed lower limit, " + str(SIZES["gof_fixed"]) +
        " synthetic samples (flares treated as independent):")
DOC.table(GOF_C, {"KS D": "{:.4f}", "p": "{:.3f}"})

# %% [markdown]
# ### Validation of the cutoff fit and of the bootstrap test
#
# Two checks on synthetic data before the tests are used on the catalog.
#
# * The cutoff model fitted from per star sums must recover the generating
#   index. The samples come from the numerical cutoff distribution.
# * The test of equality must keep its size when the index is the same in all
#   classes, also when every star has its own index around the common value, and
#   it must detect a difference. The stars have different numbers of flares, as in
#   the catalog. The same data are also tested with a Wald statistic that treats
#   the flares as independent, to show what the bootstrap by star corrects.
#   The expected size lies between 0.02 and 0.10 at the level 0.05 and the
#   expected power for a difference of 0.2 in the index is at least 0.8.

# %%
rng = np.random.default_rng(seed_for("v7"))
bias_a, est_a = [], []
for _ in range(max(10, SIZES["recovery"] // 5)):
    v, cdf = m3_table(1.9, 0.02)
    sim = np.interp(rng.random(20000), cdf, v)
    fit = fit_m3_sums(len(sim), sim.sum(), np.exp(sim).sum())
    est_a.append((fit["alpha"], fit["kappa"]))
est_a = np.array(est_a)
V7 = pd.DataFrame({"parameter": ["alpha", "kappa"], "truth": [1.9, 0.02],
                   "mean estimate": est_a.mean(axis=0), "scatter": est_a.std(axis=0)})
record_validation("cutoff fit from per star sums recovers the index", abs(est_a[:, 0].mean() - 1.9) < 0.05,
                  f"mean alpha {est_a[:, 0].mean():.3f} against 1.9, scatter {est_a[:, 0].std():.3f}, "
                  f"{len(est_a)} samples of 20000 flares")


def simulate_classes(rng, alphas, n_stars=300, mean_flares=60.0, tau=0.15, L=np.log(10.0)):
    """Window power-law flares for stars with their own index around the class value."""
    sums = []
    for a in alphas:
        n_i = rng.poisson(rng.lognormal(np.log(mean_flares) - 0.5, 1.0, n_stars))
        a_i = a + tau * rng.standard_normal(n_stars)
        a_rep = np.repeat(a_i, n_i)
        uu = a_rep - 1.0
        u = rng.random(len(a_rep))
        r = np.exp(-np.log1p(-u * (-np.expm1(-uu * L))) / uu)
        star = np.repeat(np.arange(n_stars), n_i)
        s1 = np.bincount(star, weights=np.log(r), minlength=n_stars)
        sums.append(np.column_stack([n_i.astype(float), s1]))
    return sums, L


def test_replicates(alphas, reps, B=300):
    rng = np.random.default_rng(seed_for("v8-" + "-".join(map(str, alphas))))
    rej_boot, rej_naive = 0, 0
    for _ in range(reps):
        sums, L = simulate_classes(rng, alphas)
        ah, ab, naive_w = [], [], []
        for U in sums:
            K = len(U)
            tot = U.sum(axis=0)
            ah.append(window_alpha(tot[0], tot[1], L))
            counts = rng.multinomial(K, np.full(K, 1.0 / K), size=B)
            Sb = counts @ U
            ab.append([window_alpha(Sb[i, 0], Sb[i, 1], L) for i in range(B)])
            eps = 1e-4
            info = (window_alpha(tot[0], tot[1] + eps * tot[0], L) - ah[-1]) / eps
            naive_w.append(None)
        ab = np.array(ab).T
        eq = equality_test(ah, ab, ["a", "b", "c"])
        rej_boot += int(eq["p_boot"] < ALPHA_TEST)
        v_naive = np.array([1.0 / (tot_n * (1.0 / (a - 1.0) ** 2 - L ** 2 * np.exp((a - 1.0) * L)
                                              / np.expm1((a - 1.0) * L) ** 2))
                            for a, tot_n in zip(ah, [U[:, 0].sum() for U in sums])])
        w = 1.0 / v_naive
        am = float((np.array(ah) * w).sum() / w.sum())
        chi = float((w * (np.array(ah) - am) ** 2).sum())
        rej_naive += int(stats.chi2.sf(chi, 2) < ALPHA_TEST)
    return rej_boot / reps, rej_naive / reps


size_b, size_n = test_replicates((2.0, 2.0, 2.0), SIZES["valid_eq"])
power_b, power_n = test_replicates((2.0, 2.0, 2.2), SIZES["valid_eq"])
V8 = pd.DataFrame({"scenario": ["same index, size at 0.05", "index 0.2 higher in one class, power"],
                   "bootstrap by star": [size_b, power_b], "Wald with independent flares": [size_n, power_n]})
record_validation("equality test keeps its size with correlated flares", 0.02 <= size_b <= 0.10,
                  f"rejection rate {size_b:.2f} at 0.05 with the bootstrap by star against {size_n:.2f} for the independent "
                  f"flares Wald test, {SIZES['valid_eq']} replicates")
record_validation("equality test detects a difference of 0.2", power_b >= 0.8,
                  f"rejection rate {power_b:.2f} against {power_n:.2f} for the independent flares Wald test")
V8

# %% [markdown]
# ## 10. Baseline in energy with bootstrap by star
#
# The baseline of the analysis fits the truncated power law (M2)
# above the xmin chosen by the KS distance, in energy as published. Every
# bootstrap replicate resamples whole stars and chooses xmin again. The test of
# equality uses the same procedure as in section 9, here comparing indices that
# belong to different energy ranges in each class. The goodness of fit and the
# likelihood ratios come from section 7.

# %%
def baseline_arrays(sample, c):
    d = sample.loc[sample["cls"] == c, ["KIC", "E"]].sort_values("KIC", kind="mergesort")
    kic = d["KIC"].to_numpy()
    arr = d["E"].to_numpy(float)
    _, start = np.unique(kic, return_index=True)
    return arr, np.r_[start, len(arr)]


def baseline_boot_one(args):
    seed, = args
    sh = _SHARED
    rng = np.random.default_rng(seed)
    arr, offs = sh["arr"], sh["offs"]
    K = len(offs) - 1
    idx = rng.integers(0, K, K)
    x = np.concatenate([arr[offs[i]:offs[i + 1]] for i in idx])
    if sh["noise"] > 0:
        x = x * 10.0 ** rng.normal(0.0, sh["noise"], len(x))
    x.sort()
    f = scan_xmin(x, sh["model"], GRID_DEX, MIN_TAIL_SCAN)
    return f["alpha"], f["xmin"], f["n_tail"]


def baseline_class(sample, c, model, B, seed, noise=0.0):
    arr, offs = baseline_arrays(sample, c)
    f = scan_xmin(np.sort(arr), model, GRID_DEX, MIN_TAIL_SCAN)
    ss = np.random.SeedSequence(seed).spawn(B)
    seeds = [(int(s.generate_state(1)[0]),) for s in ss]
    out = np.array(parallel_map(baseline_boot_one, seeds, N_JOBS,
                                shared={"arr": arr, "offs": offs, "model": model, "noise": noise}))
    return {"cls": c, "fit": f, "stars": len(offs) - 1, "alpha": f["alpha"], "boot_alpha": out[:, 0],
            "boot_xmin": out[:, 1], "boot_tail": out[:, 2]}


BASE = {c: baseline_class(sampleL, c, "M2", SIZES["boot"], seed_for(f"base-{c}")) for c in EXTENDED_CLASSES}
rows = []
for c in EXTENDED_CLASSES:
    b = BASE[c]
    lo16, hi84 = np.percentile(b["boot_alpha"], [16, 84])
    rows.append({"class": c, "stars": b["stars"], "tail": b["fit"]["n_tail"], "log10 xmin": np.log10(b["fit"]["xmin"]),
                 "se of log10 xmin": float(np.std(np.log10(b["boot_xmin"]), ddof=1)),
                 "alpha": b["alpha"], "se": float(b["boot_alpha"].std(ddof=1)),
                 "16 to 84 low": lo16, "16 to 84 high": hi84,
                 "GoF p": PREVIEW[(c, "M2")]["gof_p"]})
BASE_TAB = pd.DataFrame(rows)
BASE_EQ = {label: equality_test([BASE[c]["alpha"] for c in cl], np.column_stack([BASE[c]["boot_alpha"] for c in cl]), list(cl))
           for label, cl in SETS.items()}
BASE_EQ_TAB = pd.DataFrame([{"classes": k, "statistic": v["stat"], "p bootstrap": v["p_boot"], "p chi-squared": v["p_chi2"],
                             "common index": v["alpha_common"], "largest difference": v["max_abs_diff"],
                             "verdict": verdict_label(v)} for k, v in BASE_EQ.items()])
BASE_TAB, BASE_EQ_TAB

# %%
DOC.section("9. Baseline in energy as published, with bootstrap by star")
DOC.add(f"Truncated power law (M2) above the xmin of the KS distance, B = {SIZES['boot']} resamplings of stars with xmin "
        "estimated again in each. The goodness of fit comes from the fit of section 6.")
DOC.table(BASE_TAB, {"log10 xmin": "{:.2f}", "se of log10 xmin": "{:.2f}", "alpha": "{:.3f}", "se": "{:.3f}",
                     "16 to 84 low": "{:.3f}", "16 to 84 high": "{:.3f}", "GoF p": "{:.2f}"})
DOC.table(BASE_EQ_TAB, {"statistic": "{:.2f}", "p bootstrap": "{:.4f}", "p chi-squared": "{:.3g}",
                        "common index": "{:.3f}", "largest difference": "{:.3f}"})

# %% [markdown]
# ## 11. Sensitivity analysis
#
# Every variant repeats the primary analysis of section 9 with a smaller number of
# bootstrap resamplings. The verdict of a variant is the one for the cutoff index
# of classes M, K and G. The tables report the index of each class, the p value
# of the test of equality, and whether all three classes keep at least 1000
# flares above the lower limit.
#
# * S1: fixed lower limits in energy as published, which show what happens
#   without the luminosity scaling.
# * S2: class limits of Yang, shifts of 200 K, unrounded Pecaut limits and the
#   split of class M at 3200 K.
# * S3: variants of the dwarf filter.
# * S4: multiplicative factors of 1.25 on the energy of one class, and noise on
#   energy and luminosity in the bootstrap.
# * S5: other windows and lower limits of E_L.
# * S7: exclusions of stars.
# * S10: class F limited to Teff below 6500 K.

# %%
def variant_fit(smp, classes, lo, hi, B, seed_name, var="logEL"):
    res = {}
    for c in classes:
        if (smp["cls"] == c).sum() < 200:
            continue
        res[c] = primary_class_fit(smp, c, lo, hi, B, seed_for(f"{seed_name}-{c}"), var=var)
    return res


def variant_row(label, res, main=PRIMARY_CLASSES):
    row = {"variant": label}
    for c in res:
        row[f"alpha {c}"] = res[c]["alpha_c"]
    have = [c for c in main if c in res]
    if len(have) == len(main):
        ec = equality_test([res[c]["alpha_c"] for c in have], np.column_stack([res[c]["boot_alpha_c"] for c in have]), have)
        ew = equality_test([res[c]["alpha_w"] for c in have], np.column_stack([res[c]["boot_alpha_w"] for c in have]), have)
        row.update({"p cutoff": ec["p_boot"], "p window": ew["p_boot"], "largest difference": ec["max_abs_diff"],
                    "min tail": int(min(res[c]["tail"] for c in have)), "verdict": verdict_label(ec)})
    return row


VARIANTS = []
base_cfg = BASELINE
for lo, hi in ((34.6, 35.4), (34.8, 35.6), (35.0, 36.0), (34.3, 35.3), (34.7, 35.7)):
    VARIANTS.append((f"S5 window {lo} to {hi}", base_cfg, lo, hi, PRIMARY_CLASSES + ("F",)))
yang = (4000.0, 5000.0, 6000.0, 7500.0)
for label, cfg in (
    ("S2 Yang class limits", replace(base_cfg, bounds=yang)),
    ("S2 limits minus 200 K", replace(base_cfg, bounds=tuple(b - 200.0 for b in base_cfg.bounds))),
    ("S2 limits plus 200 K", replace(base_cfg, bounds=tuple(b + 200.0 for b in base_cfg.bounds))),
    ("S2 Pecaut limits unrounded", replace(base_cfg, bounds=(3850.0, 5280.0, 5920.0, 7200.0))),
    ("S3 no dwarf filter", replace(base_cfg, logg_min=None, rad_factor=None)),
    ("S3 log g 4.0 only", replace(base_cfg, rad_factor=None)),
    ("S3 log g 4.2 only", replace(base_cfg, logg_min=4.2, rad_factor=None)),
    ("S3 radius factor 1.25", replace(base_cfg, rad_factor=1.25)),
    ("S3 radius factor 2.0", replace(base_cfg, rad_factor=2.0)),
    ("S3 strict one sigma interval", replace(base_cfg, strict_interval=True)),
    ("S4 energy of M times 1.25", replace(base_cfg, class_factor=(("M", 1.25),))),
    ("S4 energy of M divided by 1.25", replace(base_cfg, class_factor=(("M", 0.8),))),
    ("S4 energy of K times 1.25", replace(base_cfg, class_factor=(("K", 1.25),))),
    ("S4 energy of K divided by 1.25", replace(base_cfg, class_factor=(("K", 0.8),))),
    ("S4 energy of G times 1.25", replace(base_cfg, class_factor=(("G", 1.25),))),
    ("S4 energy of G divided by 1.25", replace(base_cfg, class_factor=(("G", 0.8),))),
    ("S7 drop single flare stars", replace(base_cfg, drop_single_flare_stars=True)),
    ("S7 drop neighbour stars", replace(base_cfg, drop_neighbour_stars=True)),
    ("S7 drop stars without rotation", replace(base_cfg, drop_no_rotation=True)),
    ("S7 drop non DSEP stars", replace(base_cfg, drop_nondsep=True)),
    ("S10 class F below 6500 K", replace(base_cfg, bounds=(3850.0, 5300.0, 5900.0, 6500.0))),
):
    VARIANTS.append((label, cfg, XMIN_L, WINDOW_L[1], PRIMARY_CLASSES + ("F",)))

VAR_ROWS = [variant_row("baseline", {c: PRIMARY[c] for c in EXTENDED_CLASSES})]
for label, cfg, lo, hi, classes in VARIANTS:
    smp = add_scaled_energy(select_sample(cfg)[0])
    res = variant_fit(smp, classes, lo, hi, SIZES["variant"], "var-" + label)
    VAR_ROWS.append(variant_row(label, res))
VAR_TAB = pd.DataFrame(VAR_ROWS)
VAR_TAB

# %%
# Split of class M at 3200 K, tested with the upper part of M against K and G
split_smp = add_scaled_energy(select_sample(replace(base_cfg, m_split=3200.0))[0])
split_res = variant_fit(split_smp, ("Mhi", "Mfc", "K", "G", "F"), XMIN_L, WINDOW_L[1], SIZES["variant"], "var-split")
split_row = variant_row("S2 M split at 3200 K, Mhi against K and G", split_res, main=("Mhi", "K", "G"))
split_extra = {f"alpha {c}": split_res[c]["alpha_c"] for c in split_res}
pd.DataFrame([split_row])

# %%
# S1: fixed lower limits in energy as published, no luminosity scaling
s1_rows = []
for lo in (33.0, 33.5, 34.0):
    res = variant_fit(sampleL, PRIMARY_CLASSES + ("F",), lo, lo + 1.0, SIZES["variant"], f"s1-{lo}", var="logE")
    row = variant_row(f"S1 energy as published, lower limit log10 {lo}", res)
    row["alpha pure power law M, K, G, F"] = ", ".join(f"{res[c]['alpha_1']:.2f}" for c in PRIMARY_CLASSES + ("F",) if c in res)
    s1_rows.append(row)
S1_TAB = pd.DataFrame(s1_rows)
S1_TAB

# %%
# S4: noise on the energy and on the luminosity inside the bootstrap
def primary_class_noise(sample, c, lo, hi, B, seed, sig_e, with_L):
    d = sample.loc[sample["cls"] == c]
    units, inv = np.unique(d["KIC"].to_numpy(), return_inverse=True)
    K = len(units)
    logE_ = d["logE"].to_numpy(float)
    logL_ = d["logL"].to_numpy(float)
    rad = d["Rad"].to_numpy(float)
    teff = d["Teff"].to_numpy(float)
    d_r = 0.5 * (d["E_Rad"].to_numpy(float) + d["e_Rad"].to_numpy(float)) / rad
    d_t = 0.5 * (d["E_Teff"].to_numpy(float) + d["e_Teff"].to_numpy(float)) / teff
    sig_star = np.zeros(K)
    sig_flare = np.sqrt((2.0 * d_r) ** 2 + (4.0 * d_t) ** 2) / LN10
    sig_star[inv] = sig_flare
    rng = np.random.default_rng(seed)
    L = (hi - lo) * LN10
    full, win = [], []
    for _ in range(B):
        counts = rng.multinomial(K, np.full(K, 1.0 / K))
        w = counts[inv].astype(float)
        e = logE_ + (rng.normal(0.0, sig_e, len(logE_)) if sig_e > 0 else 0.0)
        l = logL_ + ((rng.standard_normal(K) * sig_star)[inv] if with_L else 0.0)
        v = e - l
        sel = (v >= lo) & (w > 0)
        r = (v[sel] - lo) * LN10
        ws = w[sel]
        full.append((ws.sum(), (ws * r).sum(), (ws * np.exp(r)).sum()))
        selw = sel & (v < hi)
        rw = (v[selw] - lo) * LN10
        wsw = w[selw]
        win.append((wsw.sum(), (wsw * rw).sum()))
    point = primary_class_fit(sample, c, lo, hi, 10, seed, n_jobs=1)
    lk0 = np.log(point["kappa"])
    boot = np.array(parallel_map(m3_boot_one, [(a, b, cc, point["alpha_c"], lk0) for a, b, cc in full], N_JOBS))
    point["boot_alpha_c"] = boot[:, 0]
    point["boot_alpha_w"] = np.array([window_alpha(n, s, L) for n, s in win])
    return point


def equality_centered(res, key, classes):
    ab = np.column_stack([res[c]["boot_" + key] for c in classes])
    ab = ab[np.isfinite(ab).all(axis=1)]
    shifted = ab - ab.mean(axis=0) + np.array([res[c][key] for c in classes])
    return equality_test([res[c][key] for c in classes], shifted, list(classes))


noise_rows = []
for label, sig_e, with_L in (("S4 noise of 0.2 dex on the energy", SIGMA_E_DEX, False),
                             ("S4 noise on the luminosity from the Mathur intervals", 0.0, True),
                             ("S4 both noises", SIGMA_E_DEX, True)):
    res = {c: primary_class_noise(sampleL, c, XMIN_L, WINDOW_L[1], SIZES["variant"], seed_for(f"noise-{label}-{c}"),
                                  sig_e, with_L) for c in PRIMARY_CLASSES + ("F",)}
    row = {"variant": label}
    for c in res:
        row[f"alpha {c}"] = res[c]["alpha_c"]
        row[f"mean bootstrap alpha {c}"] = float(np.nanmean(res[c]["boot_alpha_c"]))
    ec = equality_centered(res, "alpha_c", PRIMARY_CLASSES)
    ew = equality_centered(res, "alpha_w", PRIMARY_CLASSES)
    row.update({"p cutoff": ec["p_boot"], "p window": ew["p_boot"], "largest difference": ec["max_abs_diff"],
                "min tail": int(min(res[c]["tail"] for c in PRIMARY_CLASSES)), "verdict": verdict_label(ec)})
    noise_rows.append(row)
NOISE_TAB = pd.DataFrame(noise_rows)
NOISE_TAB

# %%
# S4 in energy as published: noise of 0.2 dex in the bootstrap of the baseline
BASE_NOISE = {c: baseline_class(sampleL, c, "M2", max(SIZES["variant"], 100), seed_for(f"basenoise-{c}"), noise=SIGMA_E_DEX)
              for c in PRIMARY_CLASSES}
BASE_NOISE_TAB = pd.DataFrame([{"class": c, "alpha": BASE_NOISE[c]["alpha"],
                                "mean bootstrap alpha with noise": float(BASE_NOISE[c]["boot_alpha"].mean()),
                                "se with noise": float(BASE_NOISE[c]["boot_alpha"].std(ddof=1)),
                                "se without noise": float(BASE[c]["boot_alpha"].std(ddof=1))} for c in PRIMARY_CLASSES])
BASE_NOISE_TAB

# %%
fmt_var = {"p cutoff": "{:.3f}", "p window": "{:.3f}", "largest difference": "{:.3f}"}
for col in list(VAR_TAB.columns):
    if col.startswith("alpha"):
        fmt_var[col] = "{:.3f}"
DOC.section("10. Sensitivity analysis")
DOC.add(f"Bootstrap by star with B = {SIZES['variant']} resamplings in every variant. The verdict refers to the cutoff index of "
        "classes M, K and G.")
DOC.table(VAR_TAB, fmt_var)
DOC.add("Split of class M at 3200 K, upper part of M against K and G, with the index of the fully convective part:")
DOC.table(pd.DataFrame([{**split_row, **{k: v for k, v in split_extra.items() if k not in split_row}}]), fmt_var)
DOC.add("Lower limits in energy as published (no luminosity scaling):")
fmt_s1 = dict(fmt_var)
DOC.table(S1_TAB, fmt_s1)
DOC.add("Noise in the bootstrap, centred on the bootstrap mean:")
fmt_n = {c: "{:.3f}" for c in NOISE_TAB.columns if c.startswith(("alpha", "mean", "p ", "largest"))}
DOC.table(NOISE_TAB, fmt_n)
DOC.add("Noise of 0.2 dex on the energy in the bootstrap of the baseline in energy as published:")
DOC.table(BASE_NOISE_TAB, {c: "{:.3f}" for c in BASE_NOISE_TAB.columns[1:]})

# %% [markdown]
# ## 12. Heterogeneity between stars and reproduction of the published index
#
# **Heterogeneity.** The main result pools the flares of each class. As a check,
# every star with at least `N_MIN_STAR` flares in the window of E_L gets its own
# truncated index. The weighted mean, the heterogeneity statistic Q of Cochran
# and the between star scatter of DerSimonian and Laird describe whether the
# pooled index hides stars with different indices. The variances come from the
# observed information of each star and assume independent flares within a star.
#
# **Reproduction of Yang and Liu.** The published indices come from the
# bisector of the ordinary least squares lines (Isobe et al. 1990, p. 106,
# Table 1) fitted to the differential distribution in fixed energy ranges, with
# class limits of 4000, 5000, 6000 and 7500 K and no dwarf filter. This check
# shows how much of the published result is due to the estimator and the ranges.

# %%
def window_information(n, u, L):
    z = u * L
    if abs(z) < 1e-3:
        return n * L ** 2 / 12.0
    return n * (1.0 / u ** 2 - L ** 2 * np.exp(z) / np.expm1(z) ** 2)


def star_indices(sample, c, lo, hi, nmin):
    d = sample.loc[(sample["cls"] == c) & (sample["logEL"] >= lo) & (sample["logEL"] < hi), ["KIC", "logEL"]]
    L = (hi - lo) * LN10
    rows = []
    for kic, g in d.groupby("KIC"):
        if len(g) < nmin:
            continue
        s = float(((g["logEL"].to_numpy() - lo) * LN10).mean())
        u = solve_u_truncated(s, L)
        info = window_information(len(g), u, L)
        if info > 0:
            rows.append((kic, 1.0 + u, 1.0 / info, len(g)))
    return pd.DataFrame(rows, columns=["KIC", "alpha", "var", "flares"])


def heterogeneity(df):
    w = 1.0 / df["var"].to_numpy()
    a = df["alpha"].to_numpy()
    am = (w * a).sum() / w.sum()
    q = float((w * (a - am) ** 2).sum())
    k = len(df)
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (q - (k - 1)) / c)
    return {"stars": k, "mean index": float(am), "Q": q, "p Q": float(stats.chi2.sf(q, k - 1)),
            "between star scatter": float(np.sqrt(tau2)), "median star index": float(np.median(a))}


het_rows = []
for c in EXTENDED_CLASSES:
    df = star_indices(sampleL, c, WINDOW_L[0], WINDOW_L[1], N_MIN_STAR)
    if len(df) >= 5:
        het_rows.append({"class": c, **heterogeneity(df)})
HET = pd.DataFrame(het_rows)
HET

# %%
def ols_bisector(x, y):
    xm, ym = x.mean(), y.mean()
    sxx = np.sum((x - xm) ** 2)
    sxy = np.sum((x - xm) * (y - ym))
    syy = np.sum((y - ym) ** 2)
    b1, b2 = sxy / sxx, syy / sxy
    return (b1 * b2 - 1.0 + np.sqrt((1.0 + b1 ** 2) * (1.0 + b2 ** 2))) / (b1 + b2)


def yang_index(logE, lo, hi, width, min_counts=1, poisson_weights=False):
    edges = np.arange(np.floor(logE.min() / width) * width, logE.max() + width, width)
    h, e = np.histogram(logE, bins=edges)
    centre = 0.5 * (e[1:] + e[:-1])
    de = 10.0 ** e[1:] - 10.0 ** e[:-1]
    ok = (h >= max(min_counts, 1)) & (e[:-1] >= lo - 1e-9) & (e[1:] <= hi + 1e-9)
    x, y = centre[ok], np.log10(h[ok] / de[ok])
    if poisson_weights:
        X = np.column_stack([np.ones(len(x)), x])
        beta = np.linalg.solve(X.T @ (h[ok][:, None] * X), X.T @ (h[ok] * y))
        return -float(beta[1]), int(h[ok].sum())
    return -float(ols_bisector(x, y)), int(h[ok].sum())


yang_cfg = replace(BASELINE, logg_min=None, rad_factor=None, bounds=yang)
yang_smp, _, _, _ = select_sample(yang_cfg)
yang_all = flares_raw[["KIC", "logE"]].merge(stars_raw[["KIC", "Teff"]], on="KIC")
yang_all["cls"] = assign_class(yang_all["Teff"].to_numpy(float), yang_cfg)
yang_counts = {"F": 8869, "G": 55259, "K": 47112, "M": 50439}
yang_ranges = {"F": (34.5, 36.0), "G": (34.5, 36.0), "K": (33.5, 35.0), "M": (33.5, 35.0)}
yang_published = {"F": 2.11, "G": 1.96, "K": 1.78, "M": 2.13}
rows = []
for c in ("M", "K", "G", "F"):
    lg = yang_smp.loc[yang_smp["cls"] == c, "logE"].to_numpy(float)
    row = {"class": c, "published index": yang_published[c], "fit range": f"{yang_ranges[c][0]} to {yang_ranges[c][1]}"}
    for wd in (0.2, 0.25, 0.3):
        a, n_in = yang_index(lg, *yang_ranges[c], wd)
        row[f"bin {wd} dex"] = a
    row["bin 0.25 dex, at least 20 flares per bin"] = yang_index(lg, *yang_ranges[c], 0.25, min_counts=20)[0]
    row["bin 0.25 dex, Poisson weighted line"] = yang_index(lg, *yang_ranges[c], 0.25, poisson_weights=True)[0]
    lg_all = yang_all.loc[yang_all["cls"] == c, "logE"].to_numpy(float)
    row["bin 0.25 dex, all stars of table 1"] = yang_index(lg_all, *yang_ranges[c], 0.25)[0]
    row["flares, Mathur matched"] = len(lg)
    row["flares, all of table 1"] = len(lg_all)
    row["flares, Fig. 3 of Yang"] = yang_counts[c]
    row["flares in range, 0.25 dex bins"] = yang_index(lg, *yang_ranges[c], 0.25)[1]
    rows.append(row)
YANG = pd.DataFrame(rows)
YANG

# %%
DOC.section("11. Heterogeneity between stars and reproduction of the published index")
DOC.add(f"Truncated index per star in the E_L window, stars with at least {N_MIN_STAR} flares in it.")
DOC.table(HET, {"mean index": "{:.3f}", "Q": "{:.1f}", "p Q": "{:.3g}", "between star scatter": "{:.3f}",
                "median star index": "{:.3f}"})
DOC.add("Bisector of the least squares lines (Isobe et al. 1990) on the differential distribution in energy, with the class "
        "limits of Yang and Liu (4000, 5000, 6000 and 7500 K), no dwarf filter and the energy ranges of their Fig. 3. "
        "The index of the table is the negative slope. The last two variants use only bins with at least 20 flares and a line "
        "weighted by the Poisson counts:")
DOC.add("Largest difference to the published index over all variants of each class: " + ", ".join(
    f"{r['class']} {max(abs(r[k] - r['published index']) for k in YANG.columns if k.startswith('bin')):.2f}"
    for _, r in YANG.iterrows()) + ". The last bin column uses all stars of table 1 of Yang and Liu, including those without Mathur "
    "parameters. The flare counts of the classes are compared with Fig. 3 of Yang and Liu in the columns that start with flares.")
DOC.table(YANG, {c: "{:.2f}" for c in YANG.columns if c.startswith("bin") or c == "published index"})

# %% [markdown]
# ## 13. Verdict
#
# The verdict on H0 follows the decision rules fixed before the analysis,
# adapted to the primary analysis in E_L.
#
# * H0 is the equality of the cutoff index of classes M, K and G. It is rejected
#   when the bootstrap p value of the test by star is below 0.05.
# * When H0 is rejected but every pairwise difference is below the margin of
#   0.10 in the index, the difference is called detectable and small.
# * The verdict is robust when it is the same in all variants of the sensitivity
#   analysis in which the three classes keep at least 1000 flares above the lower
#   limit. Variants that change the verdict are listed.
# * Index and 3/2 are compared through the 99 per cent bootstrap interval.

# %%
main_verdict = verdict_label(EQ_C["M, K, G"])
usable = VAR_TAB.dropna(subset=["verdict"])
usable = usable[(usable["min tail"] >= MIN_TAIL) & (usable["variant"] != "baseline")]  # the baseline is not a variant of itself
changed = usable[usable["verdict"] != main_verdict]
ROBUST = len(changed) == 0
n_rej_3_2 = int(H0C[(H0C["estimator"] == "cutoff") & (H0C["class"].isin(PRIMARY_CLASSES))]["3/2 rejected"].sum())
VERDICT = {
    "verdict": main_verdict,
    "p": EQ_C["M, K, G"]["p_boot"],
    "common index": EQ_C["M, K, G"]["alpha_common"],
    "largest difference": EQ_C["M, K, G"]["max_abs_diff"],
    "variants used": int(len(usable)),
    "variants with another verdict": int(len(changed)),
    "robust": ROBUST,
    "classes with 3/2 rejected": n_rej_3_2,
}
pd.DataFrame([VERDICT]), changed[["variant", "p cutoff", "verdict"]] if len(changed) else "no variant changes the verdict"

# %%
VAL_ALL = pd.DataFrame(VALIDATION)
DOC.section("12. Verdict")
DOC.add(f"Primary analysis: {VERDICT['verdict']}. Bootstrap p value {VERDICT['p']:.3f}, common cutoff index "
        f"{VERDICT['common index']:.3f}, largest pairwise difference {VERDICT['largest difference']:.3f}.")
DOC.add(f"Variants with at least {MIN_TAIL} flares above the lower limit in each class: {VERDICT['variants used']}. "
        f"Variants with another verdict: {VERDICT['variants with another verdict']}. "
        + ("The verdict is robust." if ROBUST else "The verdict depends on the choices of the analysis:"))
if len(changed):
    DOC.table(changed[["variant", "p cutoff", "verdict"]], {"p cutoff": "{:.3f}"})
DOC.add(f"The value 3/2 is outside the 99 per cent interval of the cutoff index in {n_rej_3_2} of 3 classes.")

# %% [markdown]
# ## 14. Index in narrow bands of effective temperature
#
# The equality of the index between the classes M, K and G was not robust, and
# the class limits of Yang and Liu and the lower limit of E_L changed the verdict.
# This section asks whether the index varies with Teff on a scale finer than the
# classes.
#
# **Bands.** Twelve bands with the same number of stars, set by the quantiles of
# the Teff of all analysed stars. The band limits are fixed in the bootstrap
# replicates. Variants use 8 and 16 bands and bands of 250 K. A band with fewer
# than 1000 flares above the lower limit is shown and flagged, and stays out of
# the tests (Clauset et al. 2009, p. 13).
#
# **Estimators.** The cutoff model (primary), the truncated power law in the
# window and the pure power law of section 9, in every band, with the common lower
# limit `XMIN_L` and the window `WINDOW_L`. The bootstrap resamples stars inside
# each band.
#
# **Tests.**
#
# * Equality of the index in all valid bands, with the statistic and the
#   bootstrap of section 9.
# * Trend: weighted linear regression of the index on the median Teff of the
#   band, with weights from the bootstrap variances. The regression is repeated
#   in every bootstrap replicate to obtain the error of the slope.
# * Shape: weighted chi-squared of three models for the index as a function of the
#   band: a constant, a straight line and one value for each class. The chi-squared
#   p values use the bootstrap variances and a normal approximation, so they are
#   only a reference.
#
# **Completeness diagnostic.** The peak of the differential distribution of E_L
# in each band, and its slope with Teff. A slope of the peak would mean that the
# fixed lower limit sits at a different depth of completeness in each band.

# %%
BAND_COUNTS = (12,)
CLASS_EDGES = (3850.0, 5300.0, 5900.0)


def star_table(sample):
    return sample.groupby("KIC", as_index=False).agg(Teff=("Teff", "first"))


def quantile_edges(sample, nb):
    teff = star_table(sample)["Teff"].to_numpy(float)
    edges = np.quantile(teff, np.linspace(0.0, 1.0, nb + 1))
    edges[-1] += 1.0
    return edges


def fixed_edges(sample, width, lo, min_tail):
    """Bands of fixed width, merged with a neighbour until each has enough flares above lo."""
    teff = sample["Teff"].to_numpy(float)
    edges = list(np.arange(np.floor(teff.min() / width) * width, teff.max() + width, width))
    edges[-1] = max(edges[-1], teff.max() + 1.0)
    while len(edges) > 2:
        idx = np.clip(np.searchsorted(edges, teff, side="right") - 1, 0, len(edges) - 2)
        tail = np.bincount(idx[sample["logEL"].to_numpy(float) >= lo], minlength=len(edges) - 1)
        k = int(np.argmin(tail))
        if tail[k] >= min_tail:
            break
        nbr = k - 1 if k == len(tail) - 1 or (k > 0 and tail[k - 1] <= tail[k + 1]) else k + 1
        edges.pop(max(k, nbr))
    return np.array(edges)


def band_index(teff, edges):
    return np.clip(np.searchsorted(edges, teff, side="right") - 1, 0, len(edges) - 2)


def peak_per_band(sample, edges, width=0.2):
    band = band_index(sample["Teff"].to_numpy(float), edges)
    return np.array([histogram_peak(sample.loc[band == k, "logEL"], width) if (band == k).sum() > 50 else np.nan
                     for k in range(len(edges) - 1)])


def band_fits(sample, edges, lo, hi, B, seed_name, los=None):
    """Fits of every estimator in the bands defined by edges. los gives a lower limit for each band."""
    band = band_index(sample["Teff"].to_numpy(float), edges)
    out = []
    for k in range(len(edges) - 1):
        sel = sample.loc[band == k].assign(cls="X")
        if sel["KIC"].nunique() < 10:
            continue
        lo_k = lo if los is None else float(los[k])
        res = primary_class_fit(sel, "X", lo_k, lo_k + (hi - lo), B, seed_for(f"{seed_name}-{k}"))
        star_teff = star_table(sel)["Teff"]
        res.update({"band": k, "edge low": float(edges[k]), "edge high": float(edges[k + 1]),
                    "Teff median": float(star_teff.median()), "valid": bool(res["tail"] >= MIN_TAIL),
                    "lower limit": lo_k})
        out.append(res)
    return out


def weighted_trend(x, y, w):
    X = np.column_stack([np.ones(len(x)), x - x.mean()])
    A = X.T @ (w[:, None] * X)
    beta = np.linalg.solve(A, X.T @ (w * y))
    return float(beta[1]), float(beta[0])


def trend_test(alpha_hat, alpha_boot, x):
    ab = np.asarray(alpha_boot, float)
    ab = ab[np.isfinite(ab).all(axis=1)]
    w = 1.0 / ab.var(axis=0, ddof=1)
    s_hat, _ = weighted_trend(x, np.asarray(alpha_hat, float), w)
    s_b = np.array([weighted_trend(x, row, w)[0] for row in ab])
    p = float((1 + np.sum(np.abs(s_b - s_hat) >= abs(s_hat))) / (len(s_b) + 1))
    lo, hi = np.percentile(s_b, [2.5, 97.5])
    return {"slope": s_hat, "se": float(s_b.std(ddof=1)), "low": float(lo), "high": float(hi), "p": p}


def shape_models(alpha_hat, se, x, teff_med):
    a = np.asarray(alpha_hat, float)
    w = 1.0 / np.asarray(se, float) ** 2
    n = len(a)
    chi_c = float((w * (a - (w * a).sum() / w.sum()) ** 2).sum())
    s, c0 = weighted_trend(np.asarray(x, float), a, w)
    chi_l = float((w * (a - (c0 + s * (np.asarray(x, float) - np.mean(x)))) ** 2).sum())
    grp = np.searchsorted(np.array(CLASS_EDGES), np.asarray(teff_med, float), side="right")
    chi_s, used = 0.0, 0
    for g in np.unique(grp):
        m = grp == g
        chi_s += float((w[m] * (a[m] - (w[m] * a[m]).sum() / w[m].sum()) ** 2).sum())
        used += 1
    rows = [("constant", chi_c, n - 1), ("straight line", chi_l, n - 2), ("one value per class", chi_s, n - used)]
    return pd.DataFrame([{"model": m, "chi-squared": c, "dof": d, "reference p": float(stats.chi2.sf(c, d)) if d > 0 else np.nan}
                         for m, c, d in rows])


def band_summary(bands, estimator, label):
    """Equality, trend and shape of one estimator over the valid bands."""
    key = {"cutoff": "alpha_c", "window": "alpha_w", "pure": "alpha_1"}[estimator]
    use = [b for b in bands if b["valid"]]
    ah = np.array([b[key] for b in use])
    ab = np.column_stack([b["boot_" + key] for b in use])
    x = np.array([b["Teff median"] for b in use]) / 1000.0
    eq = equality_test(ah, ab, [f"band {b['band']}" for b in use], pairs=False)
    tr = trend_test(ah, ab, x)
    se = np.nanstd(ab, axis=0, ddof=1)
    return {"label": label, "estimator": estimator, "valid bands": len(use), "equality": eq, "trend": tr,
            "shape": shape_models(ah, se, x, x * 1000.0)}


def band_table(bands):
    rows = []
    for b in bands:
        rows.append({"band": b["band"], "Teff from": b["edge low"], "Teff to": b["edge high"], "Teff median": b["Teff median"],
                     "stars": b["stars"], "tail": b["tail"], "valid": "yes" if b["valid"] else "no",
                     "alpha cutoff": b["alpha_c"], "se cutoff": float(np.nanstd(b["boot_alpha_c"], ddof=1)),
                     "alpha window": b["alpha_w"], "se window": float(np.nanstd(b["boot_alpha_w"], ddof=1)),
                     "alpha pure": b["alpha_1"], "se pure": float(np.nanstd(b["boot_alpha_1"], ddof=1))})
    return pd.DataFrame(rows)


EDGES12 = quantile_edges(sampleL, 12)
BANDS12 = band_fits(sampleL, EDGES12, XMIN_L, WINDOW_L[1], SIZES["boot"], "bands12")
BAND_TAB = band_table(BANDS12)
BAND_SUM = {e: band_summary(BANDS12, e, "12 quantile bands") for e in ("cutoff", "window", "pure")}
BAND_TAB

# %%
SUM_TAB = pd.DataFrame([{"estimator": e, "valid bands": s["valid bands"], "equality statistic": s["equality"]["stat"],
                         "p bootstrap": s["equality"]["p_boot"], "p chi-squared": s["equality"]["p_chi2"],
                         "slope per 1000 K": s["trend"]["slope"], "se": s["trend"]["se"], "95 low": s["trend"]["low"],
                         "95 high": s["trend"]["high"], "p slope": s["trend"]["p"]} for e, s in BAND_SUM.items()])
SHAPE_TAB = pd.concat([s["shape"].assign(estimator=e) for e, s in BAND_SUM.items()], ignore_index=True)
SUM_TAB, SHAPE_TAB

# %%
# Completeness diagnostic: peak of the E_L distribution in each band and its slope with Teff
def peak_boot(sample, edges, B, seed, width=0.2):
    band = band_index(sample["Teff"].to_numpy(float), edges)
    rng = np.random.default_rng(seed)
    nb = len(edges) - 1
    peaks = np.full((B, nb), np.nan)
    prep = []
    for k in range(nb):
        d = sample.loc[band == k, ["KIC", "logEL"]]
        units, inv = np.unique(d["KIC"].to_numpy(), return_inverse=True)
        v = d["logEL"].to_numpy(float)
        bins = np.arange(v.min(), v.max() + width, width)
        prep.append((len(units), inv, v, bins))
    for b in range(B):
        for k, (K, inv, v, bins) in enumerate(prep):
            w = rng.multinomial(K, np.full(K, 1.0 / K))[inv]
            h, e = np.histogram(v, bins=bins, weights=w)
            peaks[b, k] = e[int(h.argmax())] + width / 2.0
    return peaks


PEAKS_B = peak_boot(sampleL, EDGES12, SIZES["variant"], seed_for("peaks"))
PEAK_POINT = peak_per_band(sampleL, EDGES12)
x_med = np.array([b["Teff median"] for b in BANDS12]) / 1000.0
pk_slopes = np.array([np.polyfit(x_med, row, 1)[0] for row in PEAKS_B])
PEAK_DIAG = {"slope": float(np.polyfit(x_med, PEAK_POINT, 1)[0]), "se": float(pk_slopes.std(ddof=1)),
             "min": float(np.nanmin(PEAK_POINT)), "max": float(np.nanmax(PEAK_POINT))}
PEAK_TAB = pd.DataFrame({"band": [b["band"] for b in BANDS12], "Teff median": x_med * 1000.0, "peak log E_L": PEAK_POINT,
                         "se of the peak": PEAKS_B.std(axis=0, ddof=1)})
rot = sampleL.groupby("KIC").agg(Teff=("Teff", "first"), Prot=("Prot", "first")).reset_index()
rot["band"] = band_index(rot["Teff"].to_numpy(float), EDGES12)
PEAK_TAB["median Prot (d)"] = [float(rot.loc[(rot["band"] == k) & (rot["Prot"] > 0), "Prot"].median()) for k in range(len(EDGES12) - 1)]
PEAK_DIAG, PEAK_TAB

# %% [markdown]
# ### Variants of the bands
#
# Variants of the number of bands, of the lower limit, of the sample and of the
# measurement of Teff. Each repeats the fits with fewer bootstrap resamplings.
# The verdict of a variant is the one for the equality of the cutoff index in the
# valid bands, and the slope is the one of the cutoff index.

# %%
def variant_summary(label, bands):
    out = {"variant": label, "valid bands": sum(b["valid"] for b in bands), "bands": len(bands)}
    if out["valid bands"] >= 4:
        sc = band_summary(bands, "cutoff", label)
        sw = band_summary(bands, "window", label)
        out.update({"p equality cutoff": sc["equality"]["p_boot"], "slope cutoff": sc["trend"]["slope"],
                    "p slope cutoff": sc["trend"]["p"], "p equality window": sw["equality"]["p_boot"],
                    "slope window": sw["trend"]["slope"], "p slope window": sw["trend"]["p"],
                    "verdict": "H0T rejected" if sc["equality"]["p_boot"] < ALPHA_TEST else "H0T not rejected"})
    return out


BV = SIZES["variant"]
VBAND = [variant_summary("baseline, 12 bands", BANDS12)]
for nb in (8, 16):
    VBAND.append(variant_summary(f"{nb} quantile bands",
                                 band_fits(sampleL, quantile_edges(sampleL, nb), XMIN_L, WINDOW_L[1], BV, f"bands{nb}")))
VBAND.append(variant_summary("bands of 250 K, merged", band_fits(sampleL, fixed_edges(sampleL, 250.0, XMIN_L, MIN_TAIL), XMIN_L,
                                                                  WINDOW_L[1], BV, "bands250")))
BAND_LOWLIM = {}
for lo in (34.6, 34.8):
    BAND_LOWLIM[f"lo{lo}"] = band_fits(sampleL, EDGES12, lo, lo + 1.0, BV, f"lo{lo}")
    VBAND.append(variant_summary(f"lower limit {lo}", BAND_LOWLIM[f"lo{lo}"]))
for dl in (0.3, 0.5):
    los = np.where(np.isfinite(PEAK_POINT), PEAK_POINT + dl, XMIN_L)
    BAND_LOWLIM[f"peak{dl}"] = band_fits(sampleL, EDGES12, XMIN_L, XMIN_L + 1.0, BV, f"peak{dl}", los=los)
    VBAND.append(variant_summary(f"lower limit at the peak of each band plus {dl}", BAND_LOWLIM[f"peak{dl}"]))
for label, cfg in (("no dwarf filter", replace(BASELINE, logg_min=None, rad_factor=None)),
                   ("without neighbour stars", replace(BASELINE, drop_neighbour_stars=True)),
                   ("without stars lacking rotation", replace(BASELINE, drop_no_rotation=True)),
                   ("without non DSEP stars", replace(BASELINE, drop_nondsep=True))):
    smp = add_scaled_energy(select_sample(cfg)[0])
    VBAND.append(variant_summary(label, band_fits(smp, quantile_edges(smp, 12), XMIN_L, WINDOW_L[1], BV, "bv-" + label)))
VBAND_TAB = pd.DataFrame(VBAND)
VBAND_TAB

# %%
# Measurement noise in Teff: stars are resampled, their Teff is drawn inside the Mathur intervals and they are
# assigned to the fixed bands again in every replicate
def band_fits_teff_noise(sample, edges, lo, hi, B, seed_name):
    d = sample.assign(cls="X")
    Uf = class_units(d, "X", "logEL", lo)
    Uw = class_units(d, "X", "logEL", lo, hi)
    stars = d.groupby("KIC").agg(Teff=("Teff", "first"), up=("E_Teff", "first"), down=("e_Teff", "first")).sort_index()
    teff, up, down = (stars[c].to_numpy(float) for c in ("Teff", "up", "down"))
    K, nb = len(stars), len(edges) - 1
    L = (hi - lo) * LN10
    rng = np.random.default_rng(seed_for(seed_name))
    tasks, win = [], []
    for b in range(B):
        counts = rng.multinomial(K, np.full(K, 1.0 / K)).astype(float)
        u = rng.standard_normal(K)
        tf = teff + np.where(u >= 0, u * up, u * down)
        band = band_index(tf, edges)
        for k in range(nb):
            m = band == k
            tasks.append((b, k, (counts[m, None] * Uf[m]).sum(axis=0), (counts[m, None] * Uw[m]).sum(axis=0)))
    return tasks, nb, L


def noise_m3(args):
    n, s1, s2, a0, lk0 = args
    return m3_boot_one((n, s1, s2, a0, lk0))


NOISE_B = SIZES["variant"]
tasks, nb_, L_ = band_fits_teff_noise(sampleL, EDGES12, XMIN_L, WINDOW_L[1], NOISE_B, "teffnoise")
args = [(tk[2][0], tk[2][1], tk[2][2], BANDS12[tk[1]]["alpha_c"], np.log(BANDS12[tk[1]]["kappa"])) for tk in tasks]
fits = np.array(parallel_map(noise_m3, args, N_JOBS))
noise_cut = fits[:, 0].reshape(NOISE_B, nb_)
noise_win = np.array([window_alpha(tk[3][0], tk[3][1], L_) for tk in tasks]).reshape(NOISE_B, nb_)
x_all = np.array([b["Teff median"] for b in BANDS12]) / 1000.0
valid_all = np.array([b["valid"] for b in BANDS12])
rows = []
for name, arr, key in (("cutoff", noise_cut, "alpha_c"), ("window", noise_win, "alpha_w")):
    ah = np.array([b[key] for b in BANDS12])
    ab = arr[:, valid_all]
    ab = ab[np.isfinite(ab).all(axis=1)]
    centred = ab - ab.mean(axis=0) + ah[valid_all]
    eq = equality_test(ah[valid_all], centred, [str(i) for i in range(valid_all.sum())], pairs=False)
    tr = trend_test(ah[valid_all], centred, x_all[valid_all])
    rows.append({"estimator": name, "mean change of the index in a band": float(np.nanmean(arr[:, valid_all].mean(axis=0) - ah[valid_all])),
                 "p equality": eq["p_boot"], "slope per 1000 K": tr["slope"], "se": tr["se"], "p slope": tr["p"]})
TEFF_NOISE = pd.DataFrame(rows)
TEFF_NOISE

# %%
# Validation on synthetic data: size and power of the band tests, and recovery of a slope
def band_test_replicates(slope, reps, nb=12, B=200, n_stars=200):
    rng = np.random.default_rng(seed_for(f"v9-{slope}"))
    xs = np.linspace(3.3, 6.3, nb)
    alphas = 2.0 + slope * (xs - xs.mean())
    rej, rej_slope, est = 0, 0, []
    for _ in range(reps):
        sums, L = simulate_classes(rng, tuple(alphas), n_stars=n_stars)
        ah, ab = [], []
        for U in sums:
            K = len(U)
            tot = U.sum(axis=0)
            ah.append(window_alpha(tot[0], tot[1], L))
            counts = rng.multinomial(K, np.full(K, 1.0 / K), size=B)
            Sb = counts @ U
            ab.append([window_alpha(Sb[i, 0], Sb[i, 1], L) for i in range(B)])
        ab = np.array(ab).T
        eq = equality_test(ah, ab, list(range(nb)), pairs=False)
        tr = trend_test(ah, ab, xs)
        rej += int(eq["p_boot"] < ALPHA_TEST)
        rej_slope += int(tr["p"] < ALPHA_TEST)
        est.append(tr["slope"])
    return rej / reps, rej_slope / reps, float(np.mean(est) - slope)


size_g, size_s, bias0 = band_test_replicates(0.0, SIZES["valid_eq"])
power_g, power_s, bias1 = band_test_replicates(0.1, SIZES["valid_eq"])
record_validation("band equality test keeps its size without a trend", 0.02 <= size_g <= 0.10,
                  f"rejection rate {size_g:.2f} at 0.05, {SIZES['valid_eq']} replicates of 12 bands")
record_validation("band trend test keeps its size without a trend", 0.02 <= size_s <= 0.10,
                  f"rejection rate {size_s:.2f} at 0.05")
record_validation("band equality test detects a trend of 0.1 per 1000 K", power_g >= 0.8,
                  f"rejection rate {power_g:.2f}, and {power_s:.2f} for the trend test")
record_validation("slope of the index in Teff is recovered", max(abs(bias0), abs(bias1)) <= 0.02,
                  f"mean bias {bias0:.3f} without trend and {bias1:.3f} with a trend of 0.1 per 1000 K")
pd.DataFrame({"scenario": ["no trend", "trend of 0.1 per 1000 K"], "equality test": [size_g, power_g],
              "trend test": [size_s, power_s], "mean bias of the slope": [bias0, bias1]})

# %%
band_usable = VBAND_TAB[VBAND_TAB["valid bands"] >= 5].dropna(subset=["verdict"])
band_usable = band_usable[~band_usable["variant"].str.startswith("baseline")]  # the baseline is not a variant of itself
BAND_VERDICT = {"verdict": "H0T rejected" if BAND_SUM["cutoff"]["equality"]["p_boot"] < ALPHA_TEST else "H0T not rejected",
                "p": BAND_SUM["cutoff"]["equality"]["p_boot"], "slope": BAND_SUM["cutoff"]["trend"]["slope"],
                "slope se": BAND_SUM["cutoff"]["trend"]["se"], "slope p": BAND_SUM["cutoff"]["trend"]["p"],
                "variants used": int(len(band_usable))}
BAND_VERDICT["variants changed"] = int((band_usable["verdict"] != BAND_VERDICT["verdict"]).sum())
BAND_VERDICT["robust"] = BAND_VERDICT["variants changed"] == 0
BAND_VERDICT["slope sign changes"] = int((np.sign(band_usable["slope cutoff"]) != np.sign(BAND_VERDICT["slope"])).sum())
BAND_VERDICT["slope variants significant"] = int((band_usable["p slope cutoff"] < ALPHA_TEST).sum())
BAND_VERDICT["slope robust"] = bool(BAND_VERDICT["slope variants significant"] == BAND_VERDICT["variants used"]
                                    and BAND_VERDICT["slope sign changes"] == 0)
pd.DataFrame([BAND_VERDICT])

# %%
fmt_b = {"Teff from": "{:.0f}", "Teff to": "{:.0f}", "Teff median": "{:.0f}", "alpha cutoff": "{:.3f}", "se cutoff": "{:.3f}",
         "alpha window": "{:.3f}", "se window": "{:.3f}", "alpha pure": "{:.3f}", "se pure": "{:.3f}"}
DOC.section("14. Index in narrow bands of effective temperature")
DOC.add(f"Twelve bands with the same number of stars. Lower limit of E_L log10 {XMIN_L}, window {WINDOW_L[0]} to {WINDOW_L[1]}, "
        f"bootstrap by star with B = {SIZES['boot']}. A band with fewer than {MIN_TAIL} flares in the tail is flagged and stays "
        "out of the tests.")
DOC.table(BAND_TAB, fmt_b)
DOC.add("Equality of the index in the valid bands, and the slope with Teff:")
DOC.table(SUM_TAB, {"equality statistic": "{:.2f}", "p bootstrap": "{:.4f}", "p chi-squared": "{:.3g}",
                    "slope per 1000 K": "{:.3f}", "se": "{:.3f}", "95 low": "{:.3f}", "95 high": "{:.3f}", "p slope": "{:.4f}"})
DOC.add("Shape of the index as a function of the band, weighted chi-squared with a normal reference p:")
DOC.table(SHAPE_TAB[["estimator", "model", "chi-squared", "dof", "reference p"]], {"chi-squared": "{:.1f}", "reference p": "{:.3g}"})
DOC.add(f"Completeness diagnostic: peak of the E_L distribution in each band. Slope of the peak {PEAK_DIAG['slope']:.3f} "
        f"(error {PEAK_DIAG['se']:.3f}) dex per 1000 K, range {PEAK_DIAG['min']:.2f} to {PEAK_DIAG['max']:.2f}. The median rotation "
        "period is given as context only.")
DOC.table(PEAK_TAB, {"band": "{:.0f}", "Teff median": "{:.0f}", "peak log E_L": "{:.2f}", "se of the peak": "{:.2f}",
                     "median Prot (d)": "{:.1f}"})
DOC.add(f"Variants of the bands, B = {BV} resamplings each:")
DOC.table(VBAND_TAB, {c: ("{:.3f}" if c.startswith(("p ", "slope")) else "{:.0f}") for c in VBAND_TAB.columns
                      if c.startswith(("p ", "slope", "valid", "bands"))})
DOC.add(f"Measurement noise in Teff, {NOISE_B} replicates, centred on the bootstrap mean:")
DOC.table(TEFF_NOISE, {"mean change of the index in a band": "{:.3f}", "p equality": "{:.3f}", "slope per 1000 K": "{:.3f}",
                       "se": "{:.3f}", "p slope": "{:.3f}"})
DOC.add(f"Verdict on the equality of the cutoff index across the bands: {BAND_VERDICT['verdict']} (p {BAND_VERDICT['p']:.3f}). "
        f"Slope {BAND_VERDICT['slope']:.3f} per 1000 K with error {BAND_VERDICT['slope se']:.3f} (p {BAND_VERDICT['slope p']:.3f}). "
        f"{BAND_VERDICT['variants changed']} of {BAND_VERDICT['variants used']} variants with at least 5 valid bands change "
        "the verdict, and the sign of the slope changes in " + str(BAND_VERDICT["slope sign changes"]) + ". The slope of the cutoff index is "
        f"significant in {BAND_VERDICT['slope variants significant']} of {BAND_VERDICT['variants used']} variants, so the linear trend is "
        + ("robust." if BAND_VERDICT["slope robust"] else "not robust."))

# %% [markdown]
# ## 15. Additional checks
#
# Additional checks on the selection of the flares and on the test.
#
# * Common minimum duration. The catalog needs at least three consecutive points
#   above the threshold, so the shortest detectable flare lasts the same in every
#   star. Restricting the flares to a common minimum duration above that floor
#   tests whether the selection by duration produces the dependence of the index
#   on Teff.
# * Scale without the radius. The luminosity of Eq. (1) is replaced by one that
#   depends on Teff alone, with the median radius of the dwarfs in the 200 K bin
#   of the star. This removes the radius noise from the scale.
# * Bands without the hottest valid band.
# * A test of equality of the whole distribution, which does not depend on the
#   cutoff model. The statistic compares the empirical distribution of E_L above
#   the lower limit among groups of stars. The null distribution comes from
#   permuting the group labels among stars, so the flares of one star always stay
#   together.
# * Holm correction for the main tests.
# * Range of the index over the windows and over all variants.

# %%
CADENCE_DAYS = 29.4 / 1440.0
DUR_VARIANTS = ((4, "at least 4 points"), (5, "at least 5 points"))


def with_points(smp):
    out = smp.copy()
    out["points"] = np.rint((out["End"].to_numpy(float) - out["Begin"].to_numpy(float)) / CADENCE_DAYS).astype(int) + 1
    return out


def with_teff_scale(smp):
    ref = stars0[["KIC", "rad_ref"]].drop_duplicates("KIC")
    out = smp.merge(ref, on="KIC", how="left")
    lum = out["rad_ref"].to_numpy(float) ** 2 * (out["Teff"].to_numpy(float) / TEFF_SUN) ** 4
    out["logL"] = np.log10(lum)
    out["logEL"] = out["logE"] - out["logL"]
    return out


def band_alpha_range(bands):
    use = [b for b in bands if b["valid"]]
    a = np.array([b["alpha_c"] for b in use])
    k_min, k_max = int(np.argmin(a)), int(np.argmax(a))
    return {"alpha min": float(a[k_min]), "Teff at min": use[k_min]["Teff median"],
            "alpha max": float(a[k_max]), "Teff at max": use[k_max]["Teff median"]}


def band_row(label, smp, B, lo=XMIN_L, nb=12):
    bands = band_fits(smp, quantile_edges(smp, nb), lo, lo + 1.0, B, "x14b-" + label)
    row = {"variant": label, **variant_summary(label, bands)}
    row.update(band_alpha_range(bands))
    s = band_summary(bands, "cutoff", label)
    row["statistic"] = s["equality"]["stat"]
    row["dof"] = s["valid bands"] - 1
    return row, bands


def class_row(label, smp, B):
    res = variant_fit(smp, EXTENDED_CLASSES, XMIN_L, WINDOW_L[1], B, "x14c-" + label)
    return variant_row(label, res)


# %% [markdown]
# ### Common minimum duration and scale without the radius

# %%
X_CLASS, X_BAND = [], []
BX = SIZES["variant"]
sample_pts = with_points(sampleL)
cases = [("baseline", sampleL)]
for npts, lab in DUR_VARIANTS:
    cases.append((f"duration {lab}", sample_pts[sample_pts["points"] >= npts].reset_index(drop=True)))
cases.append(("scale by Teff only", with_teff_scale(sampleL)))
cases.append(("scale by Teff only, duration at least 4 points",
              with_teff_scale(sample_pts[sample_pts["points"] >= 4].reset_index(drop=True))))
BAND_SETS = {}
for lab, smp in cases:
    X_CLASS.append(class_row(lab, smp, BX))
    row, bands = band_row(lab, smp, BX)
    X_BAND.append(row)
    BAND_SETS[lab] = bands
X_CLASS_TAB = pd.DataFrame(X_CLASS)
X_BAND_TAB = pd.DataFrame(X_BAND)
DUR_SHARE = pd.DataFrame({"points": sorted(sample_pts["points"].unique())[:6]})
DUR_SHARE["share of flares"] = [float((sample_pts["points"] == p).mean()) for p in DUR_SHARE["points"]]
DUR_BY_CLASS = pd.DataFrame([{"class": c, "share with at least 4 points": float((sample_pts.loc[(sample_pts["cls"] == c) & (sample_pts["logEL"] >= XMIN_L), "points"] >= 4).mean())}
                             for c in EXTENDED_CLASSES])
X_CLASS_TAB, X_BAND_TAB

# %% [markdown]
# ### Bands without the hottest valid band

# %%
def drop_hottest(bands):
    use = [b for b in bands if b["valid"]]
    hot = max(use, key=lambda b: b["Teff median"])
    return [b for b in bands if b["band"] != hot["band"]], hot


NOHOT = []
for lab, nb in (("12 bands", 12), ("8 bands", 8), ("16 bands", 16)):
    bands = BANDS12 if nb == 12 else band_fits(sampleL, quantile_edges(sampleL, nb), XMIN_L, WINDOW_L[1], BX, f"x14d-{nb}")
    kept, hot = drop_hottest(bands)
    for est in ("cutoff", "window", "pure"):
        sm = band_summary(kept, est, lab)
        NOHOT.append({"bands": lab, "estimator": est, "dropped band median Teff": hot["Teff median"],
                      "dropped band tail": hot["tail"], "valid bands": sm["valid bands"],
                      "statistic": sm["equality"]["stat"], "dof": sm["valid bands"] - 1,
                      "p bootstrap": sm["equality"]["p_boot"], "p chi-squared": sm["equality"]["p_chi2"]})
NOHOT_TAB = pd.DataFrame(NOHOT)
NOHOT_TAB

# %% [markdown]
# ### Permutation test of the whole distribution
#
# Cramer-von Mises type statistic over groups g,
# T = sum_g n_g sum_j p_j (F_gj - F_j)^2, with F_g the empirical distribution of
# E_L above the lower limit in group g, F the pooled one and p_j the pooled
# probability of bin j, in 40 bins of equal pooled probability. Group labels are
# permuted among stars. A star with no flare above the limit does not enter. The
# test compares the whole distribution (index, cutoff and scale together), so a
# rejection does not by itself point to the index.

# %%
def star_bin_counts(smp, lo, nbins=40):
    t = smp.loc[smp["logEL"] >= lo, ["KIC", "logEL", "Teff", "cls"]]
    edges = np.quantile(t["logEL"].to_numpy(float), np.linspace(0, 1, nbins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    b = np.clip(np.searchsorted(edges, t["logEL"].to_numpy(float), side="right") - 1, 0, nbins - 1)
    kic, inv = np.unique(t["KIC"].to_numpy(), return_inverse=True)
    C = np.bincount(inv * nbins + b, minlength=len(kic) * nbins).reshape(len(kic), nbins).astype(float)
    info = t.groupby("KIC").agg(Teff=("Teff", "first"), cls=("cls", "first")).reindex(kic)
    return C, info


def cvm_stat(C, labels, groups):
    tot = C.sum(axis=0)
    N = tot.sum()
    p = tot / N
    F = np.cumsum(tot) / N
    s = 0.0
    for g in groups:
        Gc = C[labels == g].sum(axis=0)
        n = Gc.sum()
        if n == 0:
            continue
        s += n * float(np.sum(p * (np.cumsum(Gc) / n - F) ** 2))
    return s


def perm_test(C, labels, groups, P, seed):
    keep = np.isin(labels, groups)
    C, lab = C[keep], np.asarray(labels)[keep]
    t_obs = cvm_stat(C, lab, groups)
    rng = np.random.default_rng(seed)
    t_p = np.empty(P)
    for i in range(P):
        t_p[i] = cvm_stat(C, rng.permutation(lab), groups)
    return {"stat": t_obs, "p": float((1 + np.sum(t_p >= t_obs)) / (P + 1)), "P": P,
            "stars": int(keep.sum())}


PERM_P = SIZES["boot"] * 2
PERM = []
for lo in (XMIN_L, 34.8):
    C, info = star_bin_counts(sampleL, lo)
    lab_c = info["cls"].to_numpy(object)
    r = perm_test(C, lab_c, list(PRIMARY_CLASSES), PERM_P, seed_for(f"perm-cls-{lo}"))
    PERM.append({"test": "classes M, K, G", "lower limit": lo, **r})
    pair_p = []
    for a, b in (("M", "K"), ("M", "G"), ("K", "G")):
        rr = perm_test(C, lab_c, [a, b], PERM_P, seed_for(f"perm-{a}{b}-{lo}"))
        pair_p.append(rr["p"])
        PERM.append({"test": f"pair {a} and {b}", "lower limit": lo, **rr, "p Holm": np.nan})
    adj = holm_adjust(pair_p)
    for k in range(3):
        PERM[-3 + k]["p Holm"] = float(adj[k])
    edges = EDGES12
    lab_b = band_index(info["Teff"].to_numpy(float), edges)
    valid_b = [b["band"] for b in BANDS12 if b["valid"]]
    r = perm_test(C, lab_b, valid_b, PERM_P, seed_for(f"perm-band-{lo}"))
    PERM.append({"test": "valid bands of Teff", "lower limit": lo, **r})
    hot = max((b for b in BANDS12 if b["valid"]), key=lambda b: b["Teff median"])["band"]
    r = perm_test(C, lab_b, [k for k in valid_b if k != hot], PERM_P, seed_for(f"perm-bandnohot-{lo}"))
    PERM.append({"test": "valid bands without the hottest", "lower limit": lo, **r})
PERM_TAB = pd.DataFrame(PERM)
PERM_TAB

# %% [markdown]
# ### Holm correction of the main tests and ranges of the index

# %%
def p_against(boot, hat, value):
    b = np.asarray(boot, float)
    b = b[np.isfinite(b)]
    return float((1 + np.sum(np.abs(b - hat) >= abs(hat - value))) / (len(b) + 1))


MAIN_TESTS = [("equality of the cutoff index in M, K and G", EQ_C["M, K, G"]["p_boot"]),
              ("equality of the cutoff index across the valid bands", BAND_SUM["cutoff"]["equality"]["p_boot"])]
for c in PRIMARY_CLASSES:
    MAIN_TESTS.append((f"index of class {c} equal to 3/2", p_against(PRIMARY[c]["boot_alpha_c"], PRIMARY[c]["alpha_c"], 1.5)))
_adj = holm_adjust([p for _, p in MAIN_TESTS])
HOLM_TAB = pd.DataFrame({"test": [t for t, _ in MAIN_TESTS], "p": [p for _, p in MAIN_TESTS], "p Holm": _adj})

_win = VAR_TAB[VAR_TAB["variant"].isin(["baseline"]) | VAR_TAB["variant"].str.startswith("S5")]
RANGE_TAB = pd.DataFrame([
    {"class": c, "baseline": float(VAR_TAB.loc[VAR_TAB["variant"] == "baseline", f"alpha {c}"].iloc[0]),
     "window min": float(_win[f"alpha {c}"].min()), "window max": float(_win[f"alpha {c}"].max()),
     "all variants min": float(VAR_TAB[f"alpha {c}"].min()), "all variants max": float(VAR_TAB[f"alpha {c}"].max())}
    for c in PRIMARY_CLASSES])
HOLM_TAB, RANGE_TAB

# %%
DOC.section("15. Additional checks")
DOC.add("Common minimum duration. The shortest flare of the catalog has three points above the threshold. Share of the flares "
        "by number of points in the analysed sample:")
DOC.table(DUR_SHARE, {"points": "{:.0f}", "share of flares": "{:.3f}"})
DOC.add("Share of the flares above the lower limit with at least 4 points, by class:")
DOC.table(DUR_BY_CLASS, {"share with at least 4 points": "{:.3f}"})
DOC.add(f"Classes with a common minimum duration and with the scale that uses Teff alone, B = {BX} resamplings:")
DOC.table(X_CLASS_TAB, {c: ("{:.3f}" if c.startswith(("alpha", "p ", "largest")) else "{:.0f}") for c in X_CLASS_TAB.columns
                        if c.startswith(("alpha", "p ", "largest", "min"))})
DOC.add(f"Twelve bands with the same variants, B = {BX} resamplings:")
DOC.table(X_BAND_TAB, {c: ("{:.3f}" if c.startswith(("p ", "slope", "alpha")) else "{:.0f}" if c.startswith(("valid", "bands", "Teff", "dof")) else "{:.1f}")
                       for c in X_BAND_TAB.columns if c not in ("variant", "verdict")})
DOC.add("Equality of the index across the bands after dropping the hottest valid band:")
DOC.table(NOHOT_TAB, {"dropped band median Teff": "{:.0f}", "dropped band tail": "{:.0f}", "valid bands": "{:.0f}",
                      "statistic": "{:.1f}", "dof": "{:.0f}", "p bootstrap": "{:.4f}", "p chi-squared": "{:.3g}"})
DOC.add(f"Permutation test of the whole distribution of E_L above the lower limit, {PERM_P} permutations of the labels among stars. "
        "It tests the equality of the entire distribution and not only of the index. Pairwise p values are corrected with Holm.")
DOC.table(PERM_TAB, {"lower limit": "{:.1f}", "stat": "{:.4f}", "p": "{:.4f}", "P": "{:.0f}", "stars": "{:.0f}", "p Holm": "{:.4f}"})
DOC.add("Holm correction over the main tests (equality by classes, equality by bands, and the comparison with 3/2 in each class):")
DOC.table(HOLM_TAB, {"p": "{:.4f}", "p Holm": "{:.4f}"})
DOC.add("Range of the cutoff index over the five windows and over all variants of the sensitivity analysis:")
DOC.table(RANGE_TAB, {c: "{:.2f}" for c in RANGE_TAB.columns if c != "class"})

# %% [markdown]
# ## 16. Figures
#
# The figures use a layout in the style of Super Mongo: tick marks pointing inwards
# on all four sides, no grid, black text and axes, and colour only for the data. All text uses Courier Prime,
# which is released under the SIL Open Font License. The font files are in the
# `fonts` folder next to the notebook. The notebook stops with a clear message when
# the font is not found and checks that the exported files contain no other font.
#
# Requirements of the journal followed here: files in EPS with a PDF copy, width
# of 80 mm for one column or 165 mm for two columns, lines of at least 0.3 pt,
# curves told apart by colour (colour-blind safe palette) and by line style and symbol, never by grey, no figure number or
# caption inside the file, names `figure1.eps`, `figure2.eps`, `figure3.eps`.
# Greek letters are written out because Courier Prime has no Greek glyphs. Subscripts use mathtext set in Courier Prime.

# %%
import re

FIG_DIR = OUTPUT_DIR / "figures"
MM = 1.0 / 25.4
LINESTYLES = {"M": "-", "K": "--", "G": ":", "F": "-."}
MARKERS = {"M": "o", "K": "s", "G": "^", "F": "D"}
# Colour-blind safe palette (Okabe and Ito). Colour is added on top of line style and symbol, so the
# curves can still be told apart in a black and white print. Text, axes and reference lines stay black.
CLASS_COLOR = {"M": "#0072B2", "K": "#009E73", "G": "#E69F00", "F": "#CC79A7"}
EST_COLOR = {"cutoff": "#0072B2", "window": "#D55E00", "baseline": "#009E73", "published": "#CC79A7"}


def find_font_dir():
    candidates = []
    if os.environ.get("FLARE_FONT_DIR"):
        candidates.append(Path(os.environ["FLARE_FONT_DIR"]))
    candidates += [Path("fonts"), Path("..") / "fonts"]
    for c in candidates:
        if (c / "CourierPrime-Regular.ttf").exists():
            return c
    raise RuntimeError("Courier Prime was not found. Put the CourierPrime TTF files in a folder named fonts next to the "
                       "notebook or set the environment variable FLARE_FONT_DIR. The figures cannot be made without it.")


def setup_figure_style():
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm

    font_dir = find_font_dir()
    for f in sorted(font_dir.glob("CourierPrime-*.ttf")):
        fm.fontManager.addfont(str(f))
    if "Courier Prime" not in {f.name for f in fm.fontManager.ttflist}:
        raise RuntimeError("Courier Prime could not be registered from " + str(font_dir))
    plt.rcParams.update({
        "font.family": "Courier Prime", "font.size": 7.0, "axes.linewidth": 0.6, "lines.linewidth": 0.8,
        "lines.markersize": 3.0, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
        "ytick.right": True, "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.minor.size": 1.6,
        "ytick.minor.size": 1.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5, "xtick.minor.visible": True, "ytick.minor.visible": True,
        "axes.grid": False, "axes.unicode_minus": False, "ps.fonttype": 42, "pdf.fonttype": 42,
        "mathtext.fontset": "custom", "mathtext.rm": "Courier Prime", "mathtext.it": "Courier Prime:italic",
        "mathtext.bf": "Courier Prime:bold", "mathtext.default": "it",
        "legend.frameon": False, "legend.fontsize": 6.0, "figure.dpi": 120, "savefig.dpi": 300,
    })
    return plt


def save_figure(fig, name):
    FIG_DIR.mkdir(exist_ok=True)
    eps = FIG_DIR / f"{name}.eps"
    fig.savefig(eps)
    fig.savefig(FIG_DIR / f"{name}.pdf")
    text = eps.read_bytes().decode("latin-1")
    fonts = set(re.findall(r"/FontName /(\S+)", text))
    bad = sorted(f for f in fonts if not f.startswith("CourierPrime"))
    if bad:
        raise RuntimeError(f"{name}.eps contains fonts other than Courier Prime: {bad}")
    return {"file": f"{name}.eps", "fonts": ", ".join(sorted(fonts)),
            "width mm": round(fig.get_size_inches()[0] / MM, 1)}


plt = setup_figure_style()
FIG_INFO = []

# %%
# Figure 1: complementary distribution of E_L above the common limit, with the cutoff fits
fig = plt.figure(figsize=(80 * MM, 66 * MM))
ax = fig.add_axes([0.2, 0.15, 0.77, 0.82])
for c in ("M", "K", "G", "F"):
    v = sampleL.loc[(sampleL["cls"] == c) & (sampleL["logEL"] >= XMIN_L), "logEL"].to_numpy(float)
    v = np.sort(v)
    n = len(v)
    ccdf = (n - np.arange(n)) / n
    above = np.unique(np.round(np.logspace(0, np.log10(n), 45)).astype(int))
    idx = n - above
    ax.plot(v[idx], ccdf[idx], linestyle="none", marker=MARKERS[c], markerfacecolor="white", markeredgecolor=CLASS_COLOR[c],
            markeredgewidth=0.6, markersize=2.6)
    vt, cdf = m3_table(PRIMARY[c]["alpha_c"], PRIMARY[c]["kappa"])
    ax.plot(XMIN_L + vt / LN10, 1.0 - cdf, color=CLASS_COLOR[c], linestyle=LINESTYLES[c], linewidth=0.8)
ax.set_yscale("log")
ax.set_xlim(XMIN_L - 0.05, 36.9)
ax.set_ylim(1.5e-5, 1.5)
ax.set_xlabel("log $E_L$ (erg)")
ax.set_ylabel("fraction of flares above $E_L$")
handles = [plt.Line2D([], [], color=CLASS_COLOR[c], linestyle=LINESTYLES[c], marker=MARKERS[c], markerfacecolor="white",
                      markeredgewidth=0.6, markersize=2.6, linewidth=0.8) for c in ("M", "K", "G", "F")]
ax.legend(handles, ["M", "K", "G", "F"], loc="upper right", handlelength=2.6)
FIG_INFO.append(save_figure(fig, "figure1"))
plt.close(fig)

# %%
# Figure 2: indices by class for the primary analysis, the secondary and baseline estimators, and the published values
published = {"M": (2.13, 0.05), "K": (1.78, 0.02), "G": (1.96, 0.04), "F": (2.11, 0.09)}
usable_variants = VAR_TAB[(VAR_TAB["min tail"] >= MIN_TAIL)] if "min tail" in VAR_TAB else VAR_TAB
fig = plt.figure(figsize=(80 * MM, 72 * MM))
ax = fig.add_axes([0.17, 0.12, 0.8, 0.85])
for k, c in enumerate(("M", "K", "G", "F")):
    r = PRIMARY[c]
    b = r["boot_alpha_c"][np.isfinite(r["boot_alpha_c"])]
    lo16, hi84 = np.percentile(b, [16, 84])
    lo99, hi99 = np.percentile(b, [0.5, 99.5])
    x0 = k - 0.24
    ax.plot([x0, x0], [lo99, hi99], color=EST_COLOR["cutoff"], linewidth=0.5)
    ax.plot([x0, x0], [lo16, hi84], color=EST_COLOR["cutoff"], linewidth=1.4)
    ax.plot(x0, r["alpha_c"], marker="o", color=EST_COLOR["cutoff"], markersize=3.2, linestyle="none")
    col = f"alpha {c}"
    if col in usable_variants:
        vals = usable_variants[col].dropna()
        ax.plot([k - 0.36, k - 0.36], [vals.min(), vals.max()], color="black", linestyle=":", linewidth=0.8)
    x1 = k - 0.08
    ax.errorbar(x1, r["alpha_w"], yerr=np.nanstd(r["boot_alpha_w"], ddof=1), fmt="s", color=EST_COLOR["window"], markerfacecolor="white",
                markersize=3.0, elinewidth=0.6, capsize=1.5, capthick=0.5, markeredgewidth=0.6)
    bb = BASE[c]["boot_alpha"]
    lo16b, hi84b = np.percentile(bb, [16, 84])
    x2 = k + 0.08
    ax.errorbar(x2, BASE[c]["alpha"], yerr=[[BASE[c]["alpha"] - lo16b], [hi84b - BASE[c]["alpha"]]], fmt="^", color=EST_COLOR["baseline"],
                markerfacecolor="white", markersize=3.2, elinewidth=0.6, capsize=1.5, capthick=0.5, markeredgewidth=0.6)
    x3 = k + 0.24
    ax.errorbar(x3, published[c][0], yerr=published[c][1], fmt="x", color=EST_COLOR["published"], markersize=3.2, elinewidth=0.6,
                capsize=1.5, capthick=0.5, markeredgewidth=0.7)
ax.axhline(1.5, color="black", linestyle="--", linewidth=0.6)
ax.text(3.46, 1.53, "3/2", ha="right", va="bottom", fontsize=6)
ax.set_xticks(range(4))
ax.set_xticklabels(["M", "K", "G", "F"])
ax.set_xlim(-0.55, 3.55)
ax.set_ylim(1.2, 3.75)
ax.set_xlabel("class")
ax.set_ylabel("index alpha")
ax.xaxis.set_minor_locator(plt.NullLocator())
hs = [plt.Line2D([], [], color=EST_COLOR["cutoff"], marker="o", markersize=3.2, linestyle="none"),
      plt.Line2D([], [], color=EST_COLOR["window"], marker="s", markerfacecolor="white", markersize=3.0, linestyle="none", markeredgewidth=0.6),
      plt.Line2D([], [], color=EST_COLOR["baseline"], marker="^", markerfacecolor="white", markersize=3.2, linestyle="none", markeredgewidth=0.6),
      plt.Line2D([], [], color=EST_COLOR["published"], marker="x", markersize=3.2, linestyle="none", markeredgewidth=0.7),
      plt.Line2D([], [], color="black", linestyle=":", linewidth=0.8)]
ax.legend(hs, ["cutoff, $E_L$", "window, $E_L$", "baseline, $E$", "Yang and Liu", "range over variants"], loc="upper left",
          ncol=2, columnspacing=0.8, handletextpad=0.3, borderaxespad=0.2)
FIG_INFO.append(save_figure(fig, "figure2"))
plt.close(fig)

# %%
# Figure 3: bootstrap distributions of the differences and of the test statistic
def equality_null(res, key, classes):
    ah = np.array([res[c][key] for c in classes])
    ab = np.column_stack([res[c]["boot_" + key] for c in classes])
    ab = ab[np.isfinite(ab).all(axis=1)]
    w = 1.0 / ab.var(axis=0, ddof=1)

    def stat(dev):
        m = (dev * w).sum(axis=-1, keepdims=True) / w.sum()
        return (w * (dev - m) ** 2).sum(axis=-1)

    return float(stat(ah)), stat(ab - ah)


fig = plt.figure(figsize=(165 * MM, 66 * MM))
axa = fig.add_axes([0.085, 0.16, 0.395, 0.8])
axb = fig.add_axes([0.585, 0.16, 0.395, 0.8])
pair_style = {("M", "K"): "-", ("M", "G"): "--", ("K", "G"): ":"}
pair_color = {("M", "K"): "#0072B2", ("M", "G"): "#D55E00", ("K", "G"): "#009E73"}
for (ci, cj), ls in pair_style.items():
    bi, bj = PRIMARY[ci]["boot_alpha_c"], PRIMARY[cj]["boot_alpha_c"]
    n = min(len(bi), len(bj))
    d = bi[:n] - bj[:n]
    d = d[np.isfinite(d)]
    h, e = np.histogram(d, bins=np.linspace(-0.6, 0.6, 49), density=True)
    axa.stairs(h, e, color=pair_color[(ci, cj)], linestyle=ls, linewidth=0.9)
axa.axvline(0.0, color="black", linewidth=0.5)
for s in (-DELTA_EQ, DELTA_EQ):
    axa.axvline(s, color="black", linewidth=0.4, linestyle=(0, (1, 3)))
axa.set_xlabel("difference of the cutoff index between classes")
axa.set_ylabel("bootstrap density")
axa.set_xlim(-0.6, 0.6)
ha = [plt.Line2D([], [], color=pair_color[k], linestyle=ls, linewidth=0.9) for k, ls in pair_style.items()]
axa.legend(ha, ["M minus K", "M minus G", "K minus G"], loc="upper right", handlelength=2.6)
axa.text(0.03, 0.97, "A", transform=axa.transAxes, ha="left", va="top")
for key, ls, lab in (("alpha_c", "-", "cutoff"), ("alpha_w", "--", "window")):
    t_obs, t_boot = equality_null(PRIMARY, key, PRIMARY_CLASSES)
    h, e = np.histogram(np.log10(np.maximum(t_boot, 1e-3)), bins=np.linspace(-3, 2, 51), density=True)
    axb.stairs(h, e, color=EST_COLOR["cutoff" if key == "alpha_c" else "window"], linestyle=ls, linewidth=0.9)
    axb.axvline(np.log10(max(t_obs, 1e-3)), color=EST_COLOR["cutoff" if key == "alpha_c" else "window"], linestyle=ls,
                linewidth=0.6)
axb.set_xlabel("log $T$")
axb.set_ylabel("bootstrap density")
axb.set_xlim(-2.0, 2.0)
hb = [plt.Line2D([], [], color=EST_COLOR[k], linestyle=ls, linewidth=0.9) for k, ls in (("cutoff", "-"), ("window", "--"))]
axb.legend(hb, ["cutoff", "window"], loc="upper right", handlelength=2.6)
axb.text(0.03, 0.97, "B", transform=axb.transAxes, ha="left", va="top")
FIG_INFO.append(save_figure(fig, "figure3"))
plt.close(fig)

# %%
# Figure 4: index against the effective temperature in narrow bands, with the median rotation period of each band
def band_teff_axes(fig, rect, xlabel=True):
    ax_ = fig.add_axes(rect)
    ax_.set_xlim(3000, 6700)
    for edge in CLASS_EDGES:
        ax_.axvline(edge, color="black", linestyle=":", linewidth=0.5)
    if xlabel:
        ax_.set_xlabel("$T_{\\mathrm{eff}}$ (K)")
    return ax_


fig = plt.figure(figsize=(80 * MM, 104 * MM))
ax = band_teff_axes(fig, [0.17, 0.37, 0.8, 0.60], xlabel=False)
ax.tick_params(labelbottom=False)
for b in BANDS12:
    flagged = not b["valid"]
    se_c = float(np.nanstd(b["boot_alpha_c"], ddof=1))
    se_w = float(np.nanstd(b["boot_alpha_w"], ddof=1))
    eb = ax.errorbar(b["Teff median"] - 30.0, b["alpha_c"], yerr=se_c, fmt="o", color=EST_COLOR["cutoff"],
                     markerfacecolor="white" if flagged else EST_COLOR["cutoff"], markersize=3.0, elinewidth=0.6, capsize=1.3,
                     capthick=0.5, markeredgewidth=0.6)
    ew = ax.errorbar(b["Teff median"] + 30.0, b["alpha_w"], yerr=se_w, fmt="s", color=EST_COLOR["window"], markerfacecolor="white",
                     markersize=3.0, elinewidth=0.6, capsize=1.3, capthick=0.5, markeredgewidth=0.6)
    if flagged:
        for e in (eb, ew):
            for line in e[2]:
                line.set_linestyle(":")
ax.axhline(1.5, color="black", linestyle="--", linewidth=0.6)
ax.text(6680, 1.53, "3/2", ha="right", va="bottom", fontsize=6)
for xc, lab in ((3450, "M"), (4600, "K"), (5600, "G"), (6300, "F")):
    ax.text(xc, 3.22, lab, ha="center", va="top", fontsize=6)
ax.set_ylim(1.2, 3.3)
ax.set_ylabel("index alpha")
ax.text(0.02, 0.97, "A", transform=ax.transAxes, ha="left", va="top")
hs = [plt.Line2D([], [], color=EST_COLOR["cutoff"], marker="o", markersize=3.0, linestyle="none"),
      plt.Line2D([], [], color=EST_COLOR["window"], marker="s", markerfacecolor="white", markersize=3.0, linestyle="none",
                 markeredgewidth=0.6)]
ax.legend(hs, ["cutoff, $E_L$", "window, $E_L$"], loc="upper left", bbox_to_anchor=(0.0, 0.93), borderaxespad=0.3)
axr = band_teff_axes(fig, [0.17, 0.09, 0.8, 0.24])
for b, prot in zip(BANDS12, PEAK_TAB["median Prot (d)"]):
    axr.plot(b["Teff median"], prot, marker="o", color="black", markersize=3.0, linestyle="none",
             markerfacecolor="white" if not b["valid"] else "black", markeredgewidth=0.6)
axr.set_ylim(0, 13)
axr.set_ylabel("$P_{\\mathrm{rot}}$ (d)")
axr.text(0.02, 0.95, "B", transform=axr.transAxes, ha="left", va="top")
FIG_INFO.append(save_figure(fig, "figure4"))
plt.close(fig)

# %%
# Figure 5: cutoff index in the bands for three lower limits
fig = plt.figure(figsize=(80 * MM, 66 * MM))
ax = band_teff_axes(fig, [0.17, 0.15, 0.8, 0.82])
SERIES = ((BANDS12, "log $E_L$ from 34.5", "#0072B2", "o", -40.0), (BAND_LOWLIM["lo34.8"], "log $E_L$ from 34.8", "#009E73", "^", 0.0),
          (BAND_LOWLIM["peak0.5"], "peak of each band plus 0.5 dex", "#CC79A7", "D", 40.0))
for bands, lab, col, mk, dx in SERIES:
    for b in bands:
        se = float(np.nanstd(b["boot_alpha_c"], ddof=1))
        e = ax.errorbar(b["Teff median"] + dx, b["alpha_c"], yerr=se, fmt=mk, color=col,
                        markerfacecolor="white" if (not b["valid"] or mk != "o") else col, markersize=2.8, elinewidth=0.6,
                        capsize=1.2, capthick=0.5, markeredgewidth=0.6)
        if not b["valid"]:
            for line in e[2]:
                line.set_linestyle(":")
ax.axhline(1.5, color="black", linestyle="--", linewidth=0.6)
ax.text(6680, 1.53, "3/2", ha="right", va="bottom", fontsize=6)
ax.set_ylim(1.0, 3.6)
ax.set_ylabel("cutoff index alpha")
hs = [plt.Line2D([], [], color=col, marker=mk, markersize=2.8, linestyle="none", markeredgewidth=0.6,
                 markerfacecolor=col if mk == "o" else "white") for _, _, col, mk, _ in SERIES]
ax.legend(hs, [s[1] for s in SERIES], loc="upper left", borderaxespad=0.3)
FIG_INFO.append(save_figure(fig, "figure5"))
plt.close(fig)
pd.DataFrame(FIG_INFO)

# %% [markdown]
# ## 17. Results macros and final checks
#
# `results_macros.tex` holds one macro for each number used in the manuscript, so
# that text and results cannot drift apart. Names contain letters only, in the
# pattern quantity followed by the class letter. Numbers use a decimal point and
# no thousands separator. The file is rebuilt from scratch on every run. The header
# states the run mode, and a run in fast mode is marked as not valid for
# publication.

# %%
def tex_num(x, d=2):
    if x is None or not np.isfinite(x):
        return "nan"
    s = f"{abs(x):.{d}f}"
    return "\\ensuremath{-}" + s if x < 0 and float(s) != 0.0 else s


def tex_p(p, B):
    """p value from B resamplings. A value at the floor 1 / (B + 1) is written as a bound."""
    if p is None or not np.isfinite(p):
        return "nan"
    floor = np.ceil(1000.0 / (B + 1)) / 1000.0
    if p <= 1.0 / (B + 1) + 1e-12:
        return "\\ensuremath{<}" + f"{floor:.3f}"
    return "\\ensuremath{<}0.001" if p < 0.001 else f"{p:.3f}"


def tex_int(x):
    return f"{int(round(x))}"


MACROS = {}


def macro(name, value):
    assert re.fullmatch(r"[A-Za-z]+", name), name
    MACROS[name] = value


B_MAIN = SIZES["boot"]
macro("RunMode", MODE)
macro("RunValid", "yes" if MODE == "full" else "no, fast mode")
macro("BootB", tex_int(B_MAIN))
macro("XminL", tex_num(XMIN_L, 1))
macro("WindowLow", tex_num(WINDOW_L[0], 1))
macro("WindowHigh", tex_num(WINDOW_L[1], 1))
macro("MarginEq", tex_num(DELTA_EQ, 2))
funnel = funnel0.set_index("step")
macro("NflaresAll", tex_int(funnel.iloc[0]["flares"]))
macro("NstarsAll", tex_int(funnel.iloc[0]["stars"]))
macro("NflaresMatched", tex_int(funnel.iloc[1]["flares"]))
macro("NstarsMatched", tex_int(funnel.iloc[1]["stars"]))
macro("NstarsUnmatched", tex_int(funnel.iloc[0]["stars"] - funnel.iloc[1]["stars"]))
macro("NflaresDwarf", tex_int(funnel.iloc[2]["flares"]))
macro("NstarsDwarf", tex_int(funnel.iloc[2]["stars"]))
macro("NflaresSample", tex_int(funnel.iloc[-1]["flares"]))
macro("NstarsSample", tex_int(funnel.iloc[-1]["stars"]))
for c in EXTENDED_CLASSES:
    r = PRIMARY[c]
    sc = sampleL[sampleL["cls"] == c]
    macro(f"Nstars{c}", tex_int(sc["KIC"].nunique()))
    macro(f"Nflares{c}", tex_int(len(sc)))
    macro(f"Ntail{c}", tex_int(r["tail"]))
    macro(f"Nwindow{c}", tex_int(r["window"]))
    macro(f"AlphaCut{c}", tex_num(r["alpha_c"]))
    macro(f"AlphaCutErr{c}", tex_num(float(np.nanstd(r["boot_alpha_c"], ddof=1))))
    macro(f"LogCutScale{c}", tex_num(XMIN_L - np.log10(r["kappa"]), 1))
    macro(f"AlphaWin{c}", tex_num(r["alpha_w"]))
    macro(f"AlphaWinErr{c}", tex_num(float(np.nanstd(r["boot_alpha_w"], ddof=1))))
    macro(f"AlphaOne{c}", tex_num(r["alpha_1"]))
    macro(f"AlphaOneErr{c}", tex_num(float(np.nanstd(r["boot_alpha_1"], ddof=1))))
    macro(f"AlphaBase{c}", tex_num(BASE[c]["alpha"]))
    macro(f"AlphaBaseErr{c}", tex_num(float(BASE[c]["boot_alpha"].std(ddof=1))))
    macro(f"LogXminBase{c}", tex_num(np.log10(BASE[c]["fit"]["xmin"])))
    macro(f"NtailBase{c}", tex_int(BASE[c]["fit"]["n_tail"]))
    g = GOF_C.set_index("class").loc[c]
    macro(f"GofCut{c}", tex_p(float(g["p"]), SIZES["gof_fixed"]))
    macro(f"KsCut{c}", tex_num(float(g["KS D"]), 3))
    macro(f"GofBase{c}", tex_p(PREVIEW[(c, "M2")]["gof_p"], SIZES["preview_gof"]))
    macro(f"DesignEff{c}", tex_num(float(DEFF.set_index("class").loc[c].iloc[-1]), 1))
    row = H0C[(H0C["class"] == c) & (H0C["estimator"] == "cutoff")].iloc[0]
    macro(f"CutMinusHalf{c}", tex_num(float(row["alpha minus 3/2"])))
for c in PRIMARY_CLASSES:
    hrow = HET[HET["class"] == c]
    if len(hrow):
        macro(f"HetStars{c}", tex_int(hrow.iloc[0]["stars"]))
        macro(f"HetScatter{c}", tex_num(float(hrow.iloc[0]["between star scatter"])))
for key, name in (("M, K, G", ""), ("M, K, G, F", "F")):
    macro(f"PEqualCut{name}", tex_p(EQ_C[key]["p_boot"], B_MAIN))
    macro(f"PEqualWin{name}", tex_p(EQ_W[key]["p_boot"], B_MAIN))
    macro(f"PEqualOne{name}", tex_p(EQ_1[key]["p_boot"], B_MAIN))
    macro(f"PEqualBase{name}", tex_p(BASE_EQ[key]["p_boot"], B_MAIN))
macro("AlphaCommonCut", tex_num(EQ_C["M, K, G"]["alpha_common"]))
macro("DeltaMaxCut", tex_num(EQ_C["M, K, G"]["max_abs_diff"]))
for _, pr in EQ_C["M, K, G"]["pairs"].iterrows():
    a, _, b = pr["pair"].split()
    macro(f"DeltaCut{a}{b}", tex_num(pr["difference"]))
    macro(f"DeltaCutLow{a}{b}", tex_num(pr["ci95 low"]))
    macro(f"DeltaCutHigh{a}{b}", tex_num(pr["ci95 high"]))
    macro(f"PHolmCut{a}{b}", tex_p(pr["p Holm"], B_MAIN))
macro("LambdaNaive", tex_num(LR["M, K, G"]["Lambda"], 1))
macro("PNaiveLR", tex_p(LR["M, K, G"]["p naive"], 10 ** 9))
macro("VerdictText", VERDICT["verdict"])
macro("NvariantsTotal", tex_int(VERDICT["variants used"]))
macro("NvariantsChanged", tex_int(VERDICT["variants with another verdict"]))
macro("NvariantsSameVerdict", tex_int(VERDICT["variants used"] - VERDICT["variants with another verdict"]))
macro("RejectThreeHalves", tex_int(VERDICT["classes with 3/2 rejected"]))
macro("VerdictRobust", "yes" if VERDICT["robust"] else "no")
pk = PEAKS.set_index("class")
for c in EXTENDED_CLASSES:
    macro(f"PeakE{c}", tex_num(float(pk.loc[c, "peak log E, 0.1 dex"]), 2))
    macro(f"PeakEL{c}", tex_num(float(pk.loc[c, "peak log E_L, 0.1 dex"]), 2))
macro("PeakSpreadE", tex_num(SPREAD["peak log E, 0.1 dex"], 2))
macro("PeakSpreadEL", tex_num(SPREAD["peak log E_L, 0.1 dex"], 2))
for _, yr in YANG.iterrows():
    macro(f"YangRepro{yr['class']}", tex_num(float(yr["bin 0.25 dex"]), 2))
macro("NbandsTeff", tex_int(len(BANDS12)))
macro("NbandsValid", tex_int(sum(b["valid"] for b in BANDS12)))
for est, nm in (("cutoff", "Cut"), ("window", "Win"), ("pure", "One")):
    sm = BAND_SUM[est]
    macro(f"PEqualBands{nm}", tex_p(sm["equality"]["p_boot"], B_MAIN))
    macro(f"SlopeBands{nm}", tex_num(sm["trend"]["slope"], 3))
    macro(f"SlopeBandsErr{nm}", tex_num(sm["trend"]["se"], 3))
    macro(f"SlopeBandsLow{nm}", tex_num(sm["trend"]["low"], 3))
    macro(f"SlopeBandsHigh{nm}", tex_num(sm["trend"]["high"], 3))
    macro(f"PSlopeBands{nm}", tex_p(sm["trend"]["p"], B_MAIN))
    shp = sm["shape"].set_index("model")
    for key, mn in (("constant", "Const"), ("straight line", "Line"), ("one value per class", "Step")):
        macro(f"ChiBands{nm}{mn}", tex_num(float(shp.loc[key, "chi-squared"]), 1))
        macro(f"DofBands{mn}", tex_int(shp.loc[key, "dof"]))
macro("PeakSlope", tex_num(PEAK_DIAG["slope"], 3))
macro("PeakSlopeErr", tex_num(PEAK_DIAG["se"], 3))
macro("PeakMin", tex_num(PEAK_DIAG["min"], 2))
macro("PeakMax", tex_num(PEAK_DIAG["max"], 2))
macro("BandsVariantsUsed", tex_int(BAND_VERDICT["variants used"]))
macro("BandsVariantsChanged", tex_int(BAND_VERDICT["variants changed"]))
macro("BandsRobust", "yes" if BAND_VERDICT["robust"] else "no")
macro("BandsSlopeVariantsSignificant", tex_int(BAND_VERDICT["slope variants significant"]))
macro("BandsSlopeRobust", "yes" if BAND_VERDICT["slope robust"] else "no")


# Extremes of the cutoff index among the valid bands and the tail of the band left out of the tests
_valid = [b for b in BANDS12 if b["valid"]]
_lo = min(_valid, key=lambda b: b["alpha_c"])
_hi = max(_valid, key=lambda b: b["alpha_c"])
for _tag, _b in (("Min", _lo), ("Max", _hi)):
    macro(f"Band{_tag}Alpha", tex_num(_b["alpha_c"], 2))
    macro(f"Band{_tag}Err", tex_num(float(np.nanstd(_b["boot_alpha_c"], ddof=1)), 2))
    macro(f"Band{_tag}Teff", tex_int(round(_b["Teff median"] / 10.0) * 10))
macro("BandDroppedTail", tex_int(min(b["tail"] for b in BANDS12 if not b["valid"])))
_prot = PEAK_TAB["median Prot (d)"].to_numpy(float)
macro("BandHotValidTail", tex_int(BANDS12[max(i for i, b in enumerate(BANDS12) if b["valid"])]["tail"]))
macro("RotBandOne", tex_num(_prot[0], 1))
macro("RotBandTwo", tex_num(_prot[1], 1))
macro("RotBandHot", tex_num(float(_prot[max(i for i, b in enumerate(BANDS12) if b["valid"])]), 1))


def clean_label(s):
    s = re.sub(r"^S\d+ ", "", s)
    s = re.sub(r"^lower limit at the peak of each band plus ([0-9.]+)", r"limit at band peak plus \1 dex", s)
    s = s.replace("lower limit ", "lower limit log E_L ").replace("strict one sigma interval", "strict one-sigma interval")
    s = s.replace("drop ", "without ").replace("without single flare stars", "without single-flare stars")
    s = s.replace("log E_L", "$\\log E_L$")
    return s[0].upper() + s[1:]


def write_variant_tables(folder):
    """Bodies of the tables of the sensitivity variants (classes and bands), read by the manuscript."""
    Bv = SIZES["variant"]
    keep = VAR_TAB.dropna(subset=["verdict"])
    keep = keep[(keep["min tail"] >= MIN_TAIL) | (keep["variant"] == "baseline")]
    items = []
    for _, r in keep.iterrows():
        name = "Baseline" if r["variant"] == "baseline" else clean_label(r["variant"])
        mark = "" if "not rejected" in r["verdict"] else "$^{\\dagger}$"
        pv = tex_p(float(r["p cutoff"]), B_MAIN if r["variant"] == "baseline" else Bv)
        items.append(f"{name}{mark} & {tex_num(r['alpha M'], 2)} & {tex_num(r['alpha K'], 2)} & {tex_num(r['alpha G'], 2)} & {pv} \\\\")
    head = ["\\begin{tabular}{@{}lrrrr@{}}", "\\toprule", "Variant & M & K & G & $p$ \\\\", "\\midrule"]
    foot = ["\\bottomrule", "\\end{tabular}"]
    body = head + items + foot
    Path(folder / "table_variants.tex").write_text("% Rebuilt on every run of the notebook. Do not edit by hand.\n" + "\n".join(body) + "\n",
                                                   encoding="utf-8")
    rows = ["\\begin{tabular}{@{}lrrr@{}}", "\\toprule", "Variant & Bands & $p$ & Slope \\\\", "\\midrule"]
    for _, r in VBAND_TAB.dropna(subset=["verdict"]).iterrows():
        if r["valid bands"] < 5:
            continue
        pv = tex_p(float(r["p equality cutoff"]), B_MAIN if r["variant"].startswith("baseline") else Bv)
        rows.append(f"{clean_label(r['variant'])} & {int(r['valid bands'])} & {pv} & {tex_num(r['slope cutoff'], 3)} \\\\")
    rows += ["\\bottomrule", "\\end{tabular}"]
    Path(folder / "table_bandvariants.tex").write_text("% Rebuilt on every run of the notebook. Do not edit by hand.\n" + "\n".join(rows) + "\n",
                                                       encoding="utf-8")


def write_band_table(path):
    """Body of the table of the temperature bands, written by the notebook and read by the manuscript."""
    rows = ["\\begin{tabular}{@{}rrrrr@{}}", "\\toprule", "$T_{\\mathrm{eff}}$ (K) & Stars & $N$ & Cutoff & Window \\\\", "\\midrule"]
    for b in BANDS12:
        mark = "" if b["valid"] else "$^{a}$"
        se_c = float(np.nanstd(b["boot_alpha_c"], ddof=1))
        se_w = float(np.nanstd(b["boot_alpha_w"], ddof=1))
        rows.append(f"{int(round(b['Teff median']))}{mark} & {int(b['stars'])} & {int(b['tail'])} & "
                    f"{tex_num(b['alpha_c'], 2)}$\\pm${tex_num(se_c, 2)} & {tex_num(b['alpha_w'], 2)}$\\pm${tex_num(se_w, 2)} \\\\")
    rows += ["\\bottomrule", "\\end{tabular}"]
    Path(path).write_text("% Rebuilt on every run of the notebook. Do not edit by hand.\n" + "\n".join(rows) + "\n", encoding="utf-8")


# Additional checks (section 15)
_xk = {"duration at least 4 points": "DurFour", "duration at least 5 points": "DurFive",
       "scale by Teff only": "TeffOnly", "scale by Teff only, duration at least 4 points": "TeffOnlyDurFour"}
_xc = X_CLASS_TAB.set_index("variant")
_xb = X_BAND_TAB.set_index("variant")
for lab, nm in _xk.items():
    for c in PRIMARY_CLASSES:
        macro(f"Alpha{nm}{c}", tex_num(float(_xc.loc[lab, f"alpha {c}"]), 2))
    macro(f"PEqual{nm}", tex_p(float(_xc.loc[lab, "p cutoff"]), BX))
    macro(f"BandMinAlpha{nm}", tex_num(float(_xb.loc[lab, "alpha min"]), 2))
    macro(f"BandMaxAlpha{nm}", tex_num(float(_xb.loc[lab, "alpha max"]), 2))
    macro(f"PEqualBands{nm}", tex_p(float(_xb.loc[lab, "p equality cutoff"]), BX))
    macro(f"ChiBands{nm}", tex_num(float(_xb.loc[lab, "statistic"]), 1))
    macro(f"SlopeBands{nm}", tex_num(float(_xb.loc[lab, "slope cutoff"]), 3))
_db = DUR_BY_CLASS.set_index("class")
for c in EXTENDED_CLASSES:
    macro(f"ShareFourPoints{c}", tex_num(float(_db.loc[c, "share with at least 4 points"]), 2))
macro("ShareThreePoints", tex_num(float(DUR_SHARE.loc[DUR_SHARE["points"] == 3, "share of flares"].iloc[0]), 2))
_nh = NOHOT_TAB[NOHOT_TAB["bands"] == "12 bands"].set_index("estimator")
for est, nm in (("cutoff", "Cut"), ("window", "Win"), ("pure", "One")):
    macro(f"ChiBandsNoHot{nm}", tex_num(float(_nh.loc[est, "statistic"]), 1))
    macro(f"PEqualBandsNoHot{nm}", tex_p(float(_nh.loc[est, "p bootstrap"]), BX))
macro("DofBandsNoHot", tex_int(_nh.iloc[0]["dof"]))
macro("NbandsNoHot", tex_int(_nh.iloc[0]["valid bands"]))
macro("NoHotAllReject", "yes" if bool((NOHOT_TAB["p bootstrap"] < ALPHA_TEST).all()) else "no")
_pm = PERM_TAB[PERM_TAB["lower limit"] == XMIN_L].set_index("test")
macro("PermClasses", tex_p(float(_pm.loc["classes M, K, G", "p"]), PERM_P))
macro("PermBands", tex_p(float(_pm.loc["valid bands of Teff", "p"]), PERM_P))
macro("PermBandsNoHot", tex_p(float(_pm.loc["valid bands without the hottest", "p"]), PERM_P))
macro("PermPerm", tex_int(PERM_P))
for pair, nm in (("pair M and K", "MK"), ("pair M and G", "MG"), ("pair K and G", "KG")):
    macro(f"PermHolm{nm}", tex_p(float(_pm.loc[pair, "p Holm"]), PERM_P))
macro("PermAllReject", "yes" if bool((PERM_TAB["p"] < ALPHA_TEST).all()) else "no")
_hl = HOLM_TAB.set_index("test")
macro("HolmClasses", tex_p(float(_hl.loc["equality of the cutoff index in M, K and G", "p Holm"]), B_MAIN))
macro("HolmBands", tex_p(float(_hl.loc["equality of the cutoff index across the valid bands", "p Holm"]), B_MAIN))
macro("HolmThreeHalvesMax", tex_p(float(max(_hl.loc[f"index of class {c} equal to 3/2", "p Holm"] for c in PRIMARY_CLASSES)), B_MAIN))
_rg = RANGE_TAB.set_index("class")
for c in PRIMARY_CLASSES:
    macro(f"WindowMin{c}", tex_num(float(_rg.loc[c, "window min"]), 2))
    macro(f"WindowMax{c}", tex_num(float(_rg.loc[c, "window max"]), 2))
    macro(f"VariantsMin{c}", tex_num(float(_rg.loc[c, "all variants min"]), 2))
    macro(f"VariantsMax{c}", tex_num(float(_rg.loc[c, "all variants max"]), 2))


def write_macros(path):
    lines = ["% Numbers of the analysis, rebuilt on every run of the notebook. Do not edit by hand.",
             f"% Run started (UTC): {RUN['started_utc']}, mode: {MODE}, seed: {SEED}",
             "% Data checksums (sha256):"]
    for _, s in SOURCES_LOG.iterrows():
        lines.append(f"% {s['catalog']} {s['sha256']}")
    if MODE != "full":
        lines.append("% WARNING: fast mode. These numbers are for testing and not for publication.")
    for k in sorted(MACROS):
        lines.append("\\newcommand{\\" + k + "}{" + str(MACROS[k]) + "}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


write_macros(OUTPUT_DIR / "results_macros.tex")
write_band_table(OUTPUT_DIR / "table_bands.tex")
write_variant_tables(OUTPUT_DIR)
print(len(MACROS), "macros written")

# %%
DOC.section("17. Figures")
DOC.add("Files in the folder figures, in EPS and PDF. Figure 1 shows the complementary distribution of E_L above the "
        "common lower limit with the cutoff fits. Figure 2 shows the indices of each class for the primary analysis, "
        "the window estimate, the baseline in energy, the values of Yang and Liu and the range over the sensitivity "
        "variants. Figure 3 shows the bootstrap distributions of the differences between classes and of the test statistic. "
        "Figure 4 shows the index in the narrow bands of effective temperature for the cutoff and window estimators, with the class limits, and the median rotation period of each band. "
        "Figure 5 shows the cutoff index in the bands for three lower limits. "
        "The check of the fonts found only Courier Prime in the files.")
DOC.table(pd.DataFrame(FIG_INFO))

# Limitations found automatically
lim = []
for c in EXTENDED_CLASSES:
    if PRIMARY[c]["tail"] < MIN_TAIL:
        lim.append(f"Class {c} has {PRIMARY[c]['tail']} flares above the lower limit, below {MIN_TAIL}.")
for _, g in GOF_C.iterrows():
    if g["p"] < 0.1:
        p_txt = f"below {1.0 / (SIZES['gof_fixed'] + 1):.3f}" if g["p"] <= 1.0 / (SIZES["gof_fixed"] + 1) + 1e-12 else f"{g['p']:.3f}"
        lim.append(f"The cutoff model is rejected by the goodness of fit test for class {g['class']} (p {p_txt}, "
                   f"KS D {g['KS D']:.3f}). The test treats flares as independent and rejects small departures when the sample is large.")
for c in PRIMARY_CLASSES:
    hrow = HET[HET["class"] == c]
    if len(hrow) and hrow.iloc[0]["p Q"] < 0.05:
        lim.append(f"The index varies between stars of class {c} (between star scatter {hrow.iloc[0]['between star scatter']:.2f}).")
if not VERDICT["robust"]:
    lim.append("The verdict changes in " + str(VERDICT["variants with another verdict"]) + " variants of the sensitivity analysis.")
flagged_bands = [b["band"] for b in BANDS12 if not b["valid"]]
if flagged_bands:
    lim.append(f"Bands {flagged_bands} of the twelve have fewer than {MIN_TAIL} flares in the tail and stay out of the band tests.")
if not BAND_VERDICT["robust"]:
    lim.append(f"The verdict on the equality of the index across the bands changes in {BAND_VERDICT['variants changed']} of "
               f"{BAND_VERDICT['variants used']} variants.")
if abs(PEAK_DIAG["slope"]) > 2.0 * PEAK_DIAG["se"]:
    lim.append(f"The peak of the E_L distribution varies with Teff (slope {PEAK_DIAG['slope']:.3f} dex per 1000 K), so the fixed lower "
               "limit sits at different depths of completeness in the bands.")
if not BAND_VERDICT["slope robust"]:
    lim.append("The linear trend of the cutoff index with Teff is not robust: the slope is significant in "
               f"{BAND_VERDICT['slope variants significant']} of {BAND_VERDICT['variants used']} variants and depends on the "
               "lower limit of E_L. The window estimate has a slope compatible with zero.")
lim.append("The cutoff model trades the index against the cutoff scale, which raises the error of the index in some bands. "
           "The window estimate is more precise but is not the primary model.")
lim.append("The dependence of the index on rotation is not analysed, and the rotation period distribution changes with Teff.")
lim.append("The luminosity used in E_L is bolometric, from the radius and temperature of the Mathur catalog, and not the "
           "Kepler band luminosity. The typical radius uncertainty is 27 per cent.")
lim.append("The catalog has no injection and recovery test, so completeness is handled through the lower limit of E_L.")
DOC.section("18. Limitations found by the run")
for item in lim:
    DOC.add("* " + item)

DOC.section("19. Validation summary")
DOC.table(pd.DataFrame(VALIDATION))
DOC.write(OUTPUT_DIR / "analysis_results.md")
print("analysis_results.md written with", len(DOC.blocks), "blocks")
