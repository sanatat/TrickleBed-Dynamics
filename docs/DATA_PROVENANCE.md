# Data provenance and reuse

## Experimental publication, not instrument data

* **Source:** Sales et al. (2005), *Experimental Dynamic Evaluation of Three-Phase Reactors*, DOI `10.1590/S0104-66322005000300012`.
* **Figure 5:** sodium-gluconate inlet and outlet at liquid flow 5 L/h. Traces extracted from PDF vector geometry and resampled; samples are **correlated points from a published drawing**, not raw, independent instrument measurements.
* **Figure 6:** experimental sodium-gluconate outlet markers at liquid flow 10 L/h; an independently measured inlet waveform was not published in the supplied PDF.
* **Table 3:** author-fitted hydrodynamic quantities (`f_e`, `h_L`, `Pe`) derived from experiments, **not** raw RTD curves.
* **KCl:** the source mentions an auxiliary tracer in its methodology but does not give a released time series. **All repository KCl responses are purely hypothetical generated signals.**

Published-figure derivation details, original PDF SHA-256, coordinate calibration, and uncertainty notes are in [`results/v26_vector_extraction/PROVENANCE.json`](../results/v26_vector_extraction/PROVENANCE.json). All data and graph-derived files should retain attribution to the article. No original PDF is distributed.

## File interpretation

| File or directory | Data category |
|---|---|
| `results/v26_vector_extraction/fig5_*_vector.csv` | graph-derived published inlet/outlet signal (resampling) |
| `results/v26_vector_extraction/fig6_experimental_exit_vector.csv` | published figure experimental marker coordinates |
| `results/v19_sales_experimental/sales2005_table3.csv` | transcribed *author-optimized* parameters |
| `results/v29_dual_tracer/fig5_real_vs_fitted.csv` | graph-derived signal and fitted numerical output |
| `results/v29_dual_tracer/hypothetical_two_tracers.csv` | synthetic NaGl and KCl simulation |
| `results/v29_dual_tracer/local_fisher_robustness.csv` | local hypothetical identifiability sensitivity |

## Reproducibility warnings

- Do **not** label vector-digitized points as original raw lab data.
- Do **not** call plotted/resampled points independent laboratory replicates.
- Do **not** describe the KCl numerical demonstration as external validation.
- Do **not** describe the code as validated CFD or approved industrial control logic.
- For independent verification of vector extraction, download the publication yourself and use `SALES2005_PDF=/path/to/file.pdf` when running the tests.
