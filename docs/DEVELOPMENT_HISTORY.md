# TrickleBed-Dynamics v2.9 — Two-tracer identifiability, with explicit synthetic KCl (27 October 2026)

**Evidence status:** Sales et al. (2005) Figure 5 supplies a published, vector-digitized **sodium-gluconate** inlet/output record. The source mentions auxiliary **KCl** in its moment analysis, but **does not publish a KCl time series** in the user-supplied PDF. No experimental KCl observations were found and none are supplied or inferred here. **Synthetic KCl curves never enter the fit to real sodium-gluconate data.**

## New: what v2.9 actually tests

- `dual_tracer_v29.py` performs a two-parameter exploratory in-sample fit to the **real graph-extracted NaGl** Figure 5, and separately creates an explicitly **synthetic** controlled NaGl/KCl experiment for Fisher information and parameter identifiability.
- Controlled synthetic assumptions: KCl has **zero adsorption** and an **assumed** pore diffusivity `1.3e-9 m²/s`; NaGl adsorption and pore accessibility are joint unknowns. Their true synthetic values are known because they are generated from the same forward model. This is a demonstration with an **inverse crime**, not experimental validation.
- `results/v29_dual_tracer/REPORT.md`, `metrics.json`, `fisher_summary.csv`, `local_fisher_robustness.csv`, `fig5_real_vs_fitted.csv`, and `hypothetical_two_tracers.csv` provide numerical outputs and provenance.
- `results/v29_dual_tracer/NaGl_experimental_fit.png`, `fisher_information.png`, `synthetic_dual_tracer.png` visualize results. The latter two plots contain hypothetical information, not KCl measurements.
- `tests/test_v29_dual_tracer.py` adds structural checks; 128 tests passed across the full project at packaging time.
- `KCL_DATA_REQUEST.md` is an **unsent** draft for requesting authentic KCl and sodium-gluconate measurements from article authors.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python dual_tracer_v29.py
python -m pytest -q
```

*Reference:* Sales et al. (2005), **Experimental Dynamic Evaluation of Three-Phase Reactors**, BJCE, DOI 10.1590/S0104-66322005000300012. This is a nonreactive RTD numerical study, not an experimentally validated reactive CFD model or a literal implementation of the paper's analytical Laplace transfer function.

---

# TrickleBed-Dynamics v2.8 — Diagnosing model discrepancy using published RTD figures

**26 October 2026.** This release is a scientific discrepancy analysis of the v2.7 heterogeneous non-reactive tracer model using **Sales et al. (2005), Figure 5**. It is not a validated industrial reactor model or a reproduction of the full reactive CFD of Alopaeus et al. (2012).

## Main finding

The uncalibrated pore-exchange model over-retains tracer (graphical RMSE ~0.04938), whereas the simpler axial-only model is ~0.01487. Reducing assumed *accessible wetted pore volume* to 25% gives ~0.01343; a limited one-dimensional grid gives ~0.01233 at 10%. A two-parameter exploration of accessibility and axial dispersion reached ~0.00949 (in-sample). **These are not measured pore fractions or independent validation scores.** The more flexible empirical Gamma curve of v2.6 still had a lower in-sample graphical RMSE (~0.00704) with three fitted parameters.

The baseline pore + adsorption capacity implied by our uncertain assumptions is roughly 56% of the mobile liquid storage, providing a physically interpretable candidate for the excess delayed response. This is a hypothesis about the **model closure**, not proof that the real catalyst has that effective accessible pore volume. The digitized Figure-5 data remain correlated graphical samples, not raw lab readings.

## Run

```bash
python -m pip install -r requirements.txt
python v28_diagnosis.py
python -m pytest -q
```

## New artifacts

- `v28_diagnosis.py`: discrete physical-process ablations, pore-availability grid, and two-parameter sensitivity.
- `results/v28_root_cause/REPORT.md`: scientific report and limits, including Figure-5-only within-curve holdout.
- `results/v28_root_cause/ablation.csv`: fixed-choice comparisons.
- `results/v28_root_cause/pore_accessibility_profile.csv`: seven assumed accessible volume fractions.
- `results/v28_root_cause/parameter_tradeoff.csv`: accessibility–dispersion grid.
- `results/v28_root_cause/ablation_comparison.png`, `accessibility_profile.png`, `parameter_tradeoff.png`: graphs.
- `tests/test_v28_diagnosis.py`: structural tests of modified fraction, conservation and time response.

## Source and permissions

Sales, F. G., Maranhão, L. C. A., Pereira, J. A. F. R. & Abreu, C. A. M. (2005), *Experimental Dynamic Evaluation of Three-Phase Reactors*, Brazilian Journal of Chemical Engineering, 22(3), 443–452, DOI 10.1590/S0104-66322005000300012. The user supplied `download.pdf`; it is **not included** in the package. The v2.6 **vector-digitized** Figure-5 curves with provenance are retained, but are not instrument-level raw data. No independent second-flow validation is claimed because Figure-6 inlet waveform was not published.

---

# TrickleBed-Dynamics v2.7 — published-figure physical RTD model comparison

**26 October 2026.** We compare the v2.6 graph-digitized Sales (2005) Figure 5 inlet/output with a conservative one-dimensional finite-volume *nonreactive* tracer model inspired by the heterogeneous equations (7)–(9) of Sales et al., including axial dispersion, partially wetted particle exchange, radial pore diffusion, and reversible linear adsorption. The implementation is an approximation, **not** a literal reproduction of Sales Eqs. (13)–(15) or Alopaeus's 2012 reactive CFD.

**Run**: `python heterogeneous_rtd_v27.py && python -m pytest -q`

**Outputs**: `results/v27_heterogeneous_rtd/REPORT.md`, `figure5_heterogeneous_comparison.png`, `metrics.json`, `parameter_sensitivity.csv`.

**Experimental-source status**: inputs and outputs are vector-digitized *published figures*, not the authors' original instrument data. The 10 L/h Figure 6 experimental inlet signal was not published; independent cross-flow predictive validation is **not** claimed. Paper-reported kinetic and diffusion estimates conflict in units and an assumed adsorption rate and bed voidage remain explicit.

---

# TrickleBed-Dynamics v2.6 — PDF-vector extraction of published RTD curves

**Current upgrade (26 October 2026):** The user supplied the original 10-page PDF of Sales et al. (2005). Version 2.6 extracts embedded *vector glyph geometry* directly from Figures 5 and 6, calibrates figure axes against vector tick marks, and fits an exploratory gamma vs exponential inlet-to-outlet residence-time model for the 5 L/h test. **These are digitized published figures, not raw instrument data. The 10 L/h inlet time series is not published and cannot support independent transfer-function validation.**

**Key files:** `vector_digitize_v26.py`, `fit_vector_rtd_v26.py`, `results/v26_vector_extraction/REPORT.md`, `PROVENANCE.json`, `figure5_*_vector.csv`, `figure6_*_vector.csv`.

To regenerate from the uploaded PDF (PDF not redistributed in archive):

```bash
python -m pip install -r requirements.txt
python vector_digitize_v26.py --pdf /path/to/download.pdf --out results/v26_vector_extraction
python fit_vector_rtd_v26.py
python -m pytest -q
```

The retained v0.1–v2.5 code and outputs remain below. All earlier synthetic/demonstration results keep their original limitations.

---

# TrickleBed-Dynamics v2.5 — published-figure provenance and calibrated-digitization tooling

**Status:** data audit and a manual pixel-coordinate calibration tool, NOT new experimentally validated samples. Publication PDF was viewed through a browser, but no original plot bitmap could be downloaded to the runtime.

Run: `python figure_provenance_v25.py` and `python -m pytest -q`. For user-supplied plot images run `python pixel_digitize_v25.py --help`.

See `results/v25_provenance/REPORT.md` for exact scientific limitations and calibration steps. The previous modules (v0.1–v2.4) are retained.

Source: Sales et al. (2005) DOI 10.1590/S0104-66322005000300012.
