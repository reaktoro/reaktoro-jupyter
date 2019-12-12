# ---
# jupyter:
#   jupytext:
#     formats: py:light,../notebooks//ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Performing a chemical equilibrium calculation with customized activity model
#
# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation with customized activity
# models. First, we import everything from the `reaktoro` package by

from reaktoro import *
import numpy as np
import matplotlib.pyplot as plt

# To indicate phases and corresponding to them species, as well as models used to evaluate activities of the aqueous
# species, we use `ChemicalEditor`. Here, the default database SUPCRT98 is used.

editor = ChemicalEditor()

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setActivityModelDrummondCO2()
editor.addGaseousPhase(["H2O(g)", "CO2(g)"]) \
    .setChemicalModelSpycherPruessEnnis()
editor.addMineralPhase("Halite")

system = ChemicalSystem(editor)

problem = EquilibriumProblem(system)
problem.setTemperature(60, "celsius")
problem.setPressure(300, "bar")
problem.add("H2O", 1, "kg")
problem.add("CO2", 100, "g")
problem.add("NaCl", 1, "mol")

state = equilibrate(problem)

print(state)
