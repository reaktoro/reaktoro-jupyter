# ---
# jupyter:
#   jupytext:
#     formats: notebooks//ipynb,scripts//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Evaluating standard thermodynamic properties of substances and reactions
#
# This tutorial demonstrates how to use Reaktoro to evaluate standard thermodynamic properties of substances and reactions.
#
# > *If your main interest is on computing thermodynamic properties, rather than chemical equilibrium and kinetics modeling, 
# > you may want to check [ThermoFun](https://thermohub.org/thermofun/thermofun/), which is an excellent project dedicated for this task.*

# %% [markdown]
# First, we import everything from the `reaktoro` package by

# %%
from reaktoro import *  

# %% [markdown]
# We need a thermodynamic database that enables us to compute thermodynamic properties of species and reactions. The 
# object `bd` is an instance of the class `Database`:

# %%
db = Database('supcrt98.xml')  

# %% [markdown]
# To access the thermodynamic property evaluations, we create a `Thermo` object:

# %%
thermo = Thermo(db)

# %% [markdown]
# For computing the standarGIT d Gibbs energy of $\mathrm{Na+}$ at given temperature and pressure, we use function
# `standardPartialMolarGibbsEnergy()` 

# %%
T = 300.0  # temperature in K
P = 10.0e5  # pressure in Pa (equivalent to 10 bar)

G0 = thermo.standardPartialMolarGibbsEnergy(T, P, 'Na+')  

# %% [markdown]
# Similarly, we can compute the $\log(K)$ of given reaction at given $(T, P)$ by using function `logEquilibriumConstant()`: 

# %%
logK = thermo.logEquilibriumConstant(T, P, 'H2O(l) = H+ + OH-')  

# %% [markdown]
# To print computed values, `.val` must be used:

# %%
print(f'G0(Na+) = {G0.val} J/mol')  # use G0.ddT or G0.ddP to get temperature and pressure derivatives
print(f'logK(H2O(l) = H+ + OH-) = {logK.val}')  # use logK.ddT or logK.ddP to get temperature and pressure derivatives

# %% [markdown]
# **Note!** To access partial derivatives of `G0` with respect to temperature or pressure, use `G0.ddT` or `G0.ddP`.
