# Joint reactive–RTD simulation (v1.5)

Numerically integrate a passive tracer simultaneously with the v6 four-compartment reactive model. Axial dispersion now affects the mobile-liquid reactant/product fields as well as the tracer, providing a shared hydrodynamic transport hypothesis. A step change in liquid flow at 180 s changes residence time and prescribed wetting.

**All inputs and comparison curves are synthetic.** This is *numerical coupling*, not experimental validation or a fully consistent multifluid hydrodynamic closure. The tracer cannot physically change chemical kinetics (it is passive); the shared transport operator links the responses. The calibrated dispersion value from v1.4 is weakly identifiable and not validated. See metrics.json and code.
