# RTD calibration & uncertainty: synthetic demonstration (v1.3)

**Not experimental calibration, not a validated reactor model, and not reproduced 2012 Alopaeus CFD data.** No experimental RTD dataset was available in the repository. Results must never be presented as measured or as demonstrated generalization to real equipment.

## Protocol
- Generate tracer impulse profiles from a gamma tanks-in-series model, with known synthetic truth `tau=85 s`, `N=5.2`, `delay=7 s` and deterministic random Gaussian measurement noise (seed 20261026).
- Fit `tau`, effective number of tanks `N` and delay on synthetic flow `Q/Qref=1.0` only by constrained nonlinear least squares.
- Predict the **independently generated**, previously unseen dataset at `Q/Qref=1.35` with **no refitting**. Assume tau scales with inverse flow while N and delay stay constant; this extrapolation law is an assumption.
- Compare to a single-CSTR (`N=1`) constrained baseline, fitted only on calibration data.
- Residual bootstrap (`n=100`, deterministic seed) gives approximate confidence bands for parameters. Plot shaded bands are *parameter uncertainty*, not observation intervals and not calibrated predictive uncertainty.

## Results
- RTD model: tau=86.666 s, N=5.444, delay=4.741 s.
- Calibration RMSE=0.0001684 s^-1.
- Held-out RMSE=0.0002491 s^-1.
- Held-out single-CSTR RMSE=0.0023809 s^-1.
- 95% residual bootstrap parameter bands (in order: tau, N, delay): [(np.float64(83.832), np.float64(88.694)), (np.float64(5.044), np.float64(5.748)), (np.float64(2.597), np.float64(7.95))].

## Interpretation and limitations
This is a **method/pipeline check** with ground truth intentionally in the fitted model family. Strong held-out performance is expected and is not scientific proof of physical generalization. The gamma RTD is a reduced-order hydraulic surrogate with no reaction, wetting mechanics, explicit catalyst, heat, or gas dynamics. It is NOT yet coupled to the v1.2 four-compartment simulation. Data must be replaced by independent published or experimental tracer observations for external validation, with sample-time and excitation protocol recorded. Parameter identifiability requires further checking on real measurements. No claim of experimental validation is supported.

## Reproduction
`python rtd_experiment_v13.py` then `python -m pytest -q`.
