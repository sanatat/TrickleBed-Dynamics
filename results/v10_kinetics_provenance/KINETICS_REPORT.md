# v1.0 Literature-supported temperature dependence

**NOT a reproduction of Gorshkova et al. (2012).**

- Lylykangas et al. (2003), AIChE Journal 49:1508-1515: publisher abstract explicitly reports activation energies of **34 kJ/mol for TMP-1** and **49 kJ/mol for TMP-2**, on Ni/Al2O3.
- Computed Arrhenius *relative-rate* curves are physically meaningful conditional on Arrhenius temperature dependence, but they do **not** determine the absolute reaction rates.
- Legacy v0.9 model represents isooctenes as one lumped concentration, so the two Ea scenarios are **proxies** applied separately to the same lumped system, NOT a chemical simulation of both isomers.
- The adsorption coefficients, pre-exponential factors, mass transfer, heat transfer and gas-liquid equilibrium closures remain **uncalibrated illustrations**.
- The model may be especially misleading about absolute thermal rise and reaction conversion. This release does not claim a quantitative fit to Fig. 2 (2012).

## Next scientifically necessary work

1. Obtain the **complete original Lylykangas 2003 paper with parameter table** and confirm units, rate law, and both-isomer constants.
2. Implement separate TMP-1, TMP-2, isooctane states and any measured isomerization pathways.
3. Recover the full wetting closure of Lappalainen (2008) and physical thermal/phase-equilibrium properties.
4. Compare digitized 2012 Fig. 2 values only after closures are verified.
