# Experimental-literature evidence module — v1.9

## What is actually real?
The six records in `sales2005_table3.csv` are **published, experimentally derived OPTIMIZED HYDRODYNAMIC PARAMETERS** (not raw measured tracer time/concentration). Transcribed from Sales, Maranhão, Pereira & Abreu (2005), *Experimental dynamic evaluation of three-phase reactors*, *Brazilian Journal of Chemical Engineering* 22(3), Table 3 (printed page 450 / PDF page 8), DOI 10.1590/S0104-66322005000300012.

Variables: liquid superficial velocity ×10^-3 m/s, gas superficial velocity ×10^-2 m/s, fitted catalyst wetting efficiency, liquid holdup, and liquid Peclet number. The paper reports measured RTD step-response curves in Figures 5 and 6, but **numerical raw points were not accessible**. No curves were manufactured or presented as digitized measurements.

## What was computed here?
A reproducible six-observation leave-one-out comparison: training-only mean baseline, liquid-flow-only ridge linear model, and liquid-plus-gas-flow ridge linear model. All units and centering are explicit in the script `sales2005_experimental_v19.py`. CSV gives every held-out prediction. With only 6 values, results are exploratory and may be unstable. This does NOT calibrate the v1.8 mechanistic RTD solver or test it against measured RTD points.

## What is not claimed?
- No time-domain RTD calibration to measured profiles; Figures 5/6 are publication figures only.
- No demonstration of experimental reactor temperature or kinetics accuracy.
- No ability to identify all wetting and dispersion parameters from six fitted results.
- No ability to transfer estimated coefficients to the 2012 Alopaeus/Neste octene reactor without physical mismatch checks.

## Next step
Obtain raw paired inlet/outlet tracer observations from authors or digitize *legible* printed curves with x/y mapping and uncertainty; fit the 1-D tracer model against the **measured output conditional on the measured input**, rather than an idealized delta pulse. Verify consistent experiment flow rate and parameter metadata. Keep published theoretical-model curve separate from experimental markers.
