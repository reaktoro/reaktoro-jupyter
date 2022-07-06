# -*- coding: utf-8 -*-
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

# # Solubility of sodium-chloride in water
#
# This tutorial demonstrates how to use Reaktoro to study the solubility of sodium chloride salt in water. In particular,
# we simulate the precipitation of halite for very concentrated NaCl-brines by preventing excessive complexation of
# Na<sup>+</sup> and Cl<sup>-</sup> into NaCl(aq).
#
# First, we import everything from the `reaktoro` package and load the SUPCRT98 database

# +
from reaktoro import *

db = Database("supcrt98.xml")
# -

# To specify aqueous and mineral phases, we use the code below.  We set the Setschenow activity model to observe Halite
# precipitation. Otherwise, all the NaCl gets split into Na<sup>+</sup>, Cl<sup>-</sup> and NaCl(aq) species. Unlike
# the default HKF activity model, the Setschenow model uses the ionic strength to compute the activity coefficient
# of neutral species (e.g., NaCl(aq)). This is similar to the approach used in PHREEQC package.

# +
editor = ChemicalEditor(db)
aqueousphase = editor.addAqueousPhaseWithElements("H O Na Cl")
aqueousphase.setActivityModelSetschenow("NaCl(aq)", 1.0)
editor.addMineralPhase("Halite")

system = ChemicalSystem(editor)
# -

# The temperature and pressure are chosen to be 30 &deg;C and 300 bar (which is realistic for the oil and gas sites).
# Moreover, we consider a very saline NaCl-brine, the mixture of 1 kg of water and 14.0 mol of NaCl salt.

# +
problem = EquilibriumProblem(system)
problem.setTemperature(30.0, "celsius")
problem.setPressure(300.0, "bar")
problem.add("H2O", 1.0, "kg")
problem.add("NaCl", 14.0, "mol")

state = equilibrate(problem)
# -

# The full result of the calculation is output below. The obtained ionic strength is about 6.15 molal, and the amount of
# halite has 7.85 mol.

print(state)

# To access the halite molality only, we use a function below, where molality is calculate as amount of halite (in mol)
# divided by the mass of water (set to 1 kg the definition of
# [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html)):

print("Halite molality = {:4.2f} molal".format(state.speciesAmount("Halite")))

# To calculate the ionic strength directly, we need first to create function `evaluate_I` that accepts the chemical
# state properties. For that purpose, the code below can be used:

# Access chemical state properties
props = state.properties()
# Create function for evaluating the ionic strength
evaluate_I = ChemicalProperty.ionicStrength(system)
# Calculate the ionic strength
I = evaluate_I(props)
print("I = {:4.2f} molal".format(I.val))
