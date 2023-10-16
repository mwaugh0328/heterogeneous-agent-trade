Here is a basic description of files that I use results from in the paper. A more detailed description of everything will be populated soon.

- [two-country.ipynb](../../notebooks/two-country.ipynb) This is a jupyter notebook (julia) used to compute a symmetric two country model and illustrate some properties as to how everything works. Figure 1a and Figure 1b are then created from the output and this is plotted in a jupyter notebook (python) in [plot-micro-elasticity.ipynb](../../notebooks/plot-micro-elasticity.ipynb)

- [calibrate-gravity-as-guide.jl](calibrate-gravity-as-guide.jl) This is one of the main driver files to calibrate the model as described in the paper. If one were to run it, this takes time and resources. 

- [calibrate-all.jl](calibration-all.jl) This is an alternative and much faster approach to calibrating the model. It looks for a parameter vector and a vector of wages and an interest rate that satisfies (i) the moment conditions and (ii) that markets clear, all in one step. If is, however, far more sensitive to the initial guess.

- [coumpute-baseline.jl](compute-baseline.jl) This file first generates moments from the EK data set, takes in calibrated parameter values and computes a world equilibrium. Outcomes from this are then compared to the data, and output to plot trade in the model vs. the data (Figure 2 in the paper), the trade elasticities (Figure 3) in the paper, and micro moments are reported. The plotting file is [plot-calibration.ipynb](../../notebooks/plot-calibration.ipynb)

- [coumpute-baseline-log.jl](compute-baseline-log.jl) Same thing as above, but now for the log preference case.

