import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from hat_support import *
from hat_solution import *

# Parameters
β = 0.96        # Discount Rate
σ = 3           # EIS
α = 0.36        # Output Elasticity of Capital
δ = 0.08        # Depreciation Rate
φ = -3          # Borrowing Constraint
σ_ε = 0.25      # Standard deviation of logit shock
N = 5           # Number of countries

# Set up grids
# ρ = 0.6 and  σ_z = 0.3
log_zgrid, P, estar = tauchen(5, 0.2, 0.3919)
zgrid = np.exp(log_zgrid)
agrid = np.linspace(φ, 8, 100) 

# Prices
p = 1.5 * np.ones(N)
p[0] = 1

# Labor suppy
L = np.ones(N)

# Create Model Instance
m = HATmodel(β, σ, α, φ, δ, σ_ε, N, P, L)

# Find equilibrium
eqm = m.find_equilibrium(agrid, zgrid, 1.029, 1, p)






