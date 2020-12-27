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

# # Solubility of the table salt with water
#
# Importing the **reaktoro** package:

from reaktoro import *

# Initialize a thermodynamic database:

db = Database("supcrt98.xml")

# Initializing chemical editor:

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Na Cl")
editor.addMineralPhase("Halite") # sodium chloride (solid)

# Initializing chemical system:

system = ChemicalSystem(editor)

# To evaluate the content of the chemical system:

print("system = ", system)

# Initializing chemical problem:
T = 25 + 273.15
P = 1e5

problem = EquilibriumProblem(system)
problem.setTemperature(T, "kelvin")
problem.setPressure(P, "pascal")
problem.add("H2O", 1.0, "kg") # water

solver = EquilibriumSolver(system)
state = ChemicalState(system)
solver.solve(state, T, P, problem.elementAmounts())

evaluate_pH = ChemicalProperty.pH(system)

hcl_initial = 0.0

hcl_amounts = [hcl_initial]
phs = [evaluate_pH(state.properties()).val]

hcl_delta = 0.1
nsteps = 50

for i in range(nsteps):

    # Add more halite to the problem
    problem.add("HCl", hcl_delta, "mmol")

    # Equilibrate state with updated problem
    solver.solve(state, T, P, problem.elementAmounts())

    # Append new ph
    phs.append(evaluate_pH(state.properties()).val)

    # Append new halite mass
    hcl_amounts.append(hcl_amounts[-1] + hcl_delta)

print(hcl_amounts)
print(phs)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()

ax.plot(hcl_amounts, phs, label=f"ph")
ax.legend(loc="upper right")
ax.grid(True)
ax.set(xlabel='HCl amount [mol]')
ax.set(ylabel='pH [-]')
ax.set(title='Dependence of pH on th HCl amount')
fig.savefig('ph-dependence-on-hcl-amount.png')

