# coupling
L-loop running (L=1,...,5)

## Instructions

To use this module, `qcd(nf)` should be used. (By default it is the quenched limit.)

For example, one can compare the result of the RG evolution in the 2-loop case
with the exact (implicit) solution:
```python
from alphas import *
qcd(3) # nf=3
solver(.5,2)
A_2loop(.5)
```
## Plots

Included in `plot.py`
