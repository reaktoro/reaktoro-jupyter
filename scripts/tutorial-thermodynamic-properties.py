# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,../notebooks//ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
# ---

# # Evaluating standard thermodynamic properties of substances and reactions

# This tutorial demonstrates how to use Reaktoro to evaluate standard thermodynamic properties of substances and reactions.

# > **Note:** If your main interest is on computing thermodynamic properties, rather than chemical equilibrium and kinetics modeling, you may want to check [ThermoFun], which is an excellent project dedicated for this task.

# First, we import the `reaktoro` package:

from reaktoro import *

# Then we create an object of class [Database] to have access to a thermodynamic database that contains the necesary data to allows us to compute thermodynamic properties at a given temperature and pressure condition:

db = Database('supcrt98.xml')

# To evaluate the thermodynamic properties, we create a [Thermo] object:

thermo = Thermo(db)

# Below we show how the standard Gibbs energy of `Na+` is computed at 360 K and 10 bar:

T = 360.0  # temperature in K
P = 10.0e5  # pressure in Pa (equivalent to 10 bar)
G0 = thermo.standardPartialMolarGibbsEnergy(T, P, 'Na+')
print(f'G0(Na+) = {G0.val} J/mol')

# > **Note:** Use G0.ddT or G0.ddP to get temperature or pressure derivatives.

# We can also compute the log(*K*) of a reaction at given *T* and *P* as follows:

logK = thermo.logEquilibriumConstant(T, P, 'Ca++ + 2*Cl- = CaCl2(aq)')
print(f'logK(Ca++ + 2*Cl- = CaCl2(aq)) = {logK.val}')

# > **Note:** Use logK.ddT or logK.ddP to get temperature or pressure derivatives.

# [Database]: https://reaktoro.org/cpp/classReaktoro_1_1Database.html
# [Thermo]: https://reaktoro.org/cpp/classReaktoro_1_1Thermo.html
