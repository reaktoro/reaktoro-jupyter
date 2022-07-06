# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Evaluating standard thermodynamic properties of substances and reactions
#
# This tutorial demonstrates how to use Reaktoro to evaluate the standard thermodynamic properties of substances and
# reactions.
#
# > **Note:** If your main interest is on computing thermodynamic properties, rather than chemical equilibrium and
# kinetics modeling, you may want to check [ThermoFun](https://thermohub.org/thermofun/thermofun/), which is an
# excellent project dedicated for this task.
#
# First, we import the `reaktoro` package:

from reaktoro import *

# Then we create an object of class [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html) to have
# access to a thermodynamic database that contains the necessary data and allows us to compute thermodynamic properties
# at a given temperature and pressure condition:

db = Database('supcrt98.xml')

# To evaluate the thermodynamic properties, we create a
# [Thermo](https://reaktoro.org/cpp/classReaktoro_1_1Thermo.html) object:

thermo = Thermo(db)

# Below, we show how the standard Gibbs energy of Na<sup>+</sup> is computed at 360 K and 10 bar:

T = 360.0  # temperature in K
P = 10.0e5  # pressure in Pa (equivalent to 10 bar)
G0 = thermo.standardPartialMolarGibbsEnergy(T, P, 'Na+')
print(f'G0(Na+) = {G0.val} J/mol')

# > **Note:** Use `G0.ddT` or `G0.ddP` to get temperature or pressure derivatives.
#
# We can also compute the log(*K*) of a reaction at given *T* and *P* as follows:

logK = thermo.logEquilibriumConstant(T, P, 'Ca++ + 2*Cl- = CaCl2(aq)')
print(f'logK(Ca++ + 2*Cl- = CaCl2(aq)) = {logK.val}')

# > **Note:** Use `logK.ddT` or `logK.ddP` to get temperature or pressure derivatives.
#
# The rest of the properties, e.g., the standard partial molar enthalpy, entropy, or internal energy of say
# Ca<sup>2+</sup> at 250 K and 1 bar can be similarly accessed by the corresponding functions:

# +
T = 250.0
P = 1.0e5
H = thermo.standardPartialMolarEnthalpy(T, P, 'Ca++')
S = thermo.standardPartialMolarEntropy(T, P, 'Ca++')
U = thermo.standardPartialMolarInternalEnergy(T, P, 'Ca++')

print(f'H(Ca++) = {H.val} J/mol')
print(f'S(Ca++) = {S.val} J/K')
print(f'U(Ca++) = {U.val} J/mol')
