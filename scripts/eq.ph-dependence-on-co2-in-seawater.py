# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
# ---

# # Dependence of the pH on the CO<sub>2</sub>(g) amount in seawater
#
# This tutorial demonstrates how pH is dependent on the added CO<sub>2</sub>(g) amount in seawater.

# We start by importing the **reaktoro** package:

from reaktoro import *

# We initialize a thermodynamic database by:

db = Database("supcrt98.xml")

# The chemical system includes aqueous and gaseous phases:

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O C Ca Mg K Cl Na S N")
editor.addGaseousPhase(["CO2(g)"])
system = ChemicalSystem(editor)

# Initializing chemical problem corresponding to the seawater content:

# +
T = 25 + 273.15
P = 1e5
water_kg = 1.0

problem = EquilibriumProblem(system)
problem.setTemperature(T, "celsius")
problem.setPressure(P, "bar")
problem.add("H2O", water_kg, "kg")
problem.add("Ca++", 412.3 * water_kg, "mg")
problem.add("Mg++", 1290 * water_kg, "mg")
problem.add("Na+", 10768.0 * water_kg, "mg")
problem.add("K+", 399.1 * water_kg, "mg")
problem.add("Cl-", 19353.0 * water_kg, "mg")
problem.add("HCO3-", 141.682 * water_kg, "mg")
problem.add("SO4--", 2712.0 * water_kg, "mg")
# -

# Next, we define equilibrium solver to be used for range of equilibrium problems:

solver = EquilibriumSolver(system)
state = ChemicalState(system)
solver.solve(state, T, P, problem.elementAmounts())

# The function for the pH evaluation is defined by:

evaluate_pH = ChemicalProperty.pH(system)

# Finally, we define the auxiliary lists with amounts of CO<sub>2</sub> in the chemical state and corresponding to
# that state pH values:

# +
co2_initial = 0.0
co2_delta = 0.1
nsteps = 50
co2_amounts = [co2_initial]
phs = [evaluate_pH(state.properties()).val]

for i in range(nsteps):

    # Add more CO2 to the problem
    problem.add("CO2", co2_delta, "mmol")

    # Equilibrate state with updated problem
    solver.solve(state, T, P, problem.elementAmounts())

    # Append new ph
    phs.append(evaluate_pH(state.properties()).val)

    # Append new CO2 amount
    co2_amounts.append(co2_amounts[-1] + co2_delta)
# -

# Plot pH as a function of the CO<sub>2</sub> amount added into the seawater:

# +
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

ax.plot(co2_amounts, phs, label=f"pH", color='C2')
ax.legend(loc="upper right")
ax.grid(True)
ax.set(xlabel=r'CO$_2$ amount [mol]')
ax.set(ylabel='pH [-]')
ax.set(title=r'Dependence of pH on the CO$_2$ amount in seawater')
fig.savefig('ph-dependence-on-co2-amount-in-seawater.png', bbox_inches='tight')
