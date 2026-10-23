# Flare energy distributions of Kepler M, K and G dwarfs

Code and results for the Letter "Is the power-law index of flare energies universal across
M, K and G dwarfs?". The notebook tests whether the power-law index of the flare energy
distribution is the same for M, K and G dwarfs in the Kepler flare catalogue of Yang and Liu
(2019), and how it changes with effective temperature. Every number, table and figure of the
Letter is written by the notebook.

## Contents

| File or folder | What it holds |
|---|---|
| `flare_universality.py` | Source of the notebook in jupytext percent format. |
| `flare_universality.ipynb` | The same notebook, generated from the `.py` file. It is stored without outputs. |
| `analysis_results.md` | Every table of the analysis, the validation summary and the limitations found by the run. Result of the full run. |
| `results_macros.tex` | One LaTeX macro for each number used in the Letter. |
| `table_bands.tex`, `table_variants.tex`, `table_bandvariants.tex` | Tables of the Letter. |
| `figures/` | Figures 1 to 5 in EPS and PDF. |
| `fonts/` | Courier Prime, the only font of the figures (SIL Open Font License). |
| `requirements.txt` | Python packages. |

## Running

```
pip install -r requirements.txt
jupytext --to ipynb flare_universality.py
jupyter nbconvert --to notebook --execute flare_universality.ipynb
```

Run the commands from this folder, so that the notebook finds the `fonts` folder. A full run
takes about 45 minutes on four cores. Setting the environment variable `FLARE_MODE=fast` runs the
same steps with fewer replicates, in a few minutes. A fast run marks its outputs as not valid for
publication, and it overwrites `analysis_results.md`, `results_macros.tex`, the tables and the
figures, so run it in a copy of the folder.

## Data

The notebook downloads the catalogues from VizieR, so it needs network access:

* Yang and Liu (2019), J/ApJS/241/29, tables `table1` (stars) and `table2` (flares).
* Mathur et al. (2017), J/ApJS/229/30, stellar properties of the Kepler DR25 targets.

No data file is stored in this repository. Downloads are cached in `data_cache/`, which git
ignores. On every run the notebook records the catalogue identifiers, the access date, the
number of rows and a checksum of the data, and they appear in `analysis_results.md` and at the
top of `results_macros.tex`. The random seed is fixed.

## Licence and citation

The code is released under the MIT licence (see `LICENSE`). The fonts keep their own licence.
Please cite the Letter and this repository (see `CITATION.cff`). Archived version: see the DOI
in the Letter.
