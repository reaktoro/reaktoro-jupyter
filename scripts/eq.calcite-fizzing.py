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

# # Calcite fizzing
#
#
#
# Import necessary packages:

from reaktoro import *
import numpy as np
import matplotlib.pyplot as plt

# Function that simulates chemical equilibrium of 10 mol calcite with given solute (defined by the instance `problem`):

def calcite_fizzing(mHCl, system, T, P):

    n0Calcite = 10.0

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")
    problem.add("HCl", mHCl, "mmol")
    problem.add("Calcite", n0Calcite, "mol")

    solver = EquilibriumSolver(system)
    state = ChemicalState(system)

    # Equilibrate pure water
    solver.solve(state, problem)

    # Fetch the amount of final calcite in the equilibrium state
    nCO2g = state.speciesAmount("CO2(g)")

    return nCO2g

# Initialize chemical system with aqueous, gaseous, and calcite phases:

db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O C Cl")
editor.addGaseousPhase(["CO2(g)"]).setChemicalModelSpycherPruessEnnis()
editor.addMineralPhase("Calcite")
system = ChemicalSystem(editor)

# Initialize the array of molar amount of NaCl to be mixed with water and calculate the amount of released CO2(g)
# when limestone is reacting with acidic brine:

amounts_hcl = np.arange(0.0, 5.0, 0.1)
P = 1.0
T = 25.0
mCO2 = [calcite_fizzing(mHCl, system, T, P) for mHCl in amounts_hcl]
print(mCO2)

# Ask Allan: why CO2(g) is not releasing?
