# Sales (2005) PDF-vector RTD extraction and transfer fit — v2.6

## Source and provenance

- Sales, Maranhão, Pereira & Abreu (2005), *Experimental Dynamic Evaluation of Three-Phase Reactors*, Brazilian Journal of Chemical Engineering 22(3), pp. 443–452, printed page **450**, Fig. **5** and **6**; DOI 10.1590/S0104-66322005000300012.
- User provided publisher-style **vector PDF**, not an image screenshot. Script `vector_digitize_v26.py` reads the original vector-triangle and circle marker coordinates with PyMuPDF; no OCR or invented observations. Source PDF SHA-256 is in `PROVENANCE.json`.
- Figure 5: **59 resampled points per displayed series** reconstructed from high-density vector marker geometry. *Do not interpret them as 59 statistically independent laboratory measurements.* Figure 6: **61 printed triangle centers**.
- PDF point axes are calibrated using literal plot tick marks: `PROVENANCE.json` stores full geometry and a conservative ±1pt/±2.2pt positioning scenario (not experimental uncertainties). The raw plotted markers may extend slightly beyond [0,1] due to glyph offset; model fitting clips signals to that interval, while CSV retains direct unmodified coordinates.

## Figure 5: QL=5 L/h, QG=300 L/h — input-output model

We use the actually plotted *negative step inlet* together with a non-negative normalized residence-time transfer kernel (gamma distribution). The output is calculated by convolving the **change in inlet** with the delayed RTD kernel. Fit uses SciPy nonlinear least squares.

| Model | Fitted RTD mean (s) | Delay (s) | Graphical curve RMSE |
|---|---:|---:|---:|
| Gamma | 74.349 | 73.000 | 0.007039 |
| Exponential | 89.624 | 56.307 | 0.008531 |
| No residence | N/A | N/A | 0.136170 |

Gamma shape estimate = 0.600. Figure 5 half-concentration crossing inlet = 325.7s; outlet = 476.7s. RMSE is a **graphical-fit diagnostic** not measurement validation; the dense trace points are correlated.

## Figure 6: QL=10 L/h

- **61 observed plotted triangles** extracted from the original PDF. Output half-concentration crossing = 233.5s.
- The **authors' own fitted-model curve** is separately reconstructed from the printed vector polyline. Its RMSE compared to the plotted triangle centers is **0.02763** normalized signal units. This is a comparison of the AUTHORS' in-sample fit to their own figure, **not a validation of our gamma kernel**.
- **No numerical inlet waveform** for this independent run is published in Figure 6, so we do **not** assert a valid cross-flow input-output prediction or calculate a misleading independent transfer-model RMSE.
- An experimentally defensible follow-up requires the actual inlet recording or a documented input forcing protocol. A prompt step cannot be equated to the measured inlet without evidence.

## Scientific interpretation and limitations

1. We have **published graph-derived experimental observations**, not true *raw instrument observations*. This is a substantial provenance improvement over v2.5 hand estimates.
2. This validates the extraction method and enables exploratory calibration of a phenomenological transfer kernel at **one flow condition**; it does **not** validate a fully coupled reactive/thermal industrial reactor model.
3. The principal bottleneck for next-stage publication-quality validation is an independent paired inlet/outlet RTD dataset for a second condition, plus known line width/calibration/instrument errors.
4. The figure in Sales (2005) is a step-down tracer response, **not** a pulse-response probability density.

## Reproduce

```bash
python vector_digitize_v26.py --pdf /path/to/download.pdf --out results/v26_vector_extraction
python fit_vector_rtd_v26.py
pytest -q
```

The source PDF is **not bundled** for copyright reasons; it is user-supplied and must be supplied explicitly. No rounding or inferred pseudo-experiments are added to the published figures.
