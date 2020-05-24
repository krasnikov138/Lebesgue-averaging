# Lebesgue-averaging
The attempt of reproducing Shilkov's article "Lebesgue moment method for solving the neutron transport equation".

### Implemented features
-  Unifed reader of endf files for following __MT__ sections:
   -  `MT=3` for cross section data
   -  `MT=6` for energy-angle distribution of reaction products
-  Tool for splitting of energy interval on carriers of resonances (using the algorithm from the page 12 in the article)
-  Calculation of Legesgue's variable for each carrier (see section _Legesgue's averaging_)
-  Auxiliary formulas for interpolating on two dimensional energy-angle distribution in case of angular part given as Legengre polinomials. Calculated coefs of transitions between carriers of resonances.
### To do
-  Obtain optimal grid for Legesgue's variable using system of exponential or basis of rational functions.
-  Solve transport equation on this grid.
-  Backmapping neutron current in E-space using Legesgue's moments.
### Building
Type in root project's directory:
```bash
make release
```
The standard directory for object files is `run/bin`. Executable file name is `run/main`.
To see info about supported command variables, type:
```bash
run/main --help
```
