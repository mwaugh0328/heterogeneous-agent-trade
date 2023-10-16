Here is a basic description of files.

- [coumpute_baseline.jl](compute_baseline.jl) This file first generates moments from the EK data set, takes in calibrated parameter values computes the world equilibrium which is then compared to model counterparts, and then generates output to compare trade in the model vs. the data (Figure 2 in the paper), the trade elasticities (Figure 3) in the paper, and micro moments are reported. The plotting file is [plot-calibration.ipynb](../../notebooks/plot-calibration.ipynb)

- [coumpute_baseline-log.jl](compute_baseline-log.jl) Same thing as above, but now for the log preference case.

