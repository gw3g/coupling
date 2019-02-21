# coupling
`L`-loop running (L=1,...,5), for QCD with any number of massless quarks.

## Instructions

To use this module, `Qcd(nf)` should be run to setup the proper group factors
and beta function coefficients. (By default it is the quenched limit.)

The function `Solver(t,L)` will return α(t) evaluated at t=ln(μ²/Λ²), and
determined by integrating the `L` order renormalisation group (RG) equations from t=10³.
As a boundary condition, we use the UV-approximation also available in `A_asymp(t, L)'.

As an example, one can compare the result of the RG evolution in the 2-loop case
with the exact (implicit) solution:
```python
from alphas import *
Qcd(3) # nf=3
Solver(.5,2)
A_2loop(.5)
```

We include also an effective coupling for space- and time-like momenta. `A_analy(t)`.
It is a 1-loop level approximation based on [hep-ph/9512336](https://arxiv.org/abs/hep-ph/9512336).

## Plots

Included in `plot.py`

## Interpolation

To avoid re-running the RG solver every time the coupling is evaluted at a scale μ², it is faster
to use an interpolated function. For that purpose we provide the function `ReadTable(nf,l)`
which can read files like "table_Coupling_{nf0,3-loop}.dat". The output is a table
with the first column being μ/Λ, and the second column contains the associated α(...).
The inerpolation function is `interpolate_Alpha(tbl,mu)`.
