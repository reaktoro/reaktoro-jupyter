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
#       jupytext_version: 1.3.2
# ---

# # Carbon dioxide gas solubility in the NaCl-brine
#
# This tutorial demonstrates how to simulate the solubility of CO2 gas in the NaCl-brine and its dependence of
# salinity of the brine (also called as a salting-out effect).
#
# First, we import all the necessary packages for further simulations:

from reaktoro import *
import numpy as np
import matplotlib.pyplot as plt

# Function `solubility_co2()` returns the concentration of  CO<sup>2</sup>(g) that was dissolved in the brined:

def solubility_co2(system, T, P, n0CO2g, mNaCl):

    # Define equilibrium problem as a mixture of NaCl-brine with given salinity and
    # CO2 of a given initial concentration at fixed T and P
    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")
    problem.add("CO2", n0CO2g, "mol")
    problem.add("NaCl", mNaCl, "mol")

    # Equilibrate chemical problem
    state = equilibrate(problem)

    # Return the difference of initial amount of the CO2(g) and remaining one after equilibration
    return (n0CO2g - state.speciesAmount("CO2(g)"))

# Define database
database = Database("supcrt98.xml")

# Initialize phases with chemical editor
editor = ChemicalEditor(database)
editor.addAqueousPhaseWithElements("H O C Na Cl")
editor.addGaseousPhase(["CO2(g)"])
editor.addMineralPhase("Halite")

# Create chemical system
system = ChemicalSystem(editor)

# Initialize temperature (in celsius)
T = np.arange(20.0, 150.0, 5.0)
# Initialize pressure (in bar)
P = 1.0  # 100.0
# Initial amount of CO2(g) (in mol)
n0CO2g = 10.0

# Generate the lists of CO2(g) mols that got dissolved in NaCl-brine of different salinity
deltaCO2_nacl1 = [solubility_co2(system, x, P, n0CO2g, mNaCl=1.0) for x in T]
deltaCO2_nacl2 = [solubility_co2(system, x, P, n0CO2g, mNaCl=2.0) for x in T]
deltaCO2_nacl4 = [solubility_co2(system, x, P, n0CO2g, mNaCl=4.0) for x in T]
deltaCO2_nacl6 = [solubility_co2(system, x, P, n0CO2g, mNaCl=6.0) for x in T]

# Plot solubility of CO2(g) as function of temperature for different salinities of NaCl-brine:

fig, ax = plt.subplots()
ax.plot(T, deltaCO2_nacl1, label=f"1 NaCl molal")
ax.plot(T, deltaCO2_nacl2, label=f"2 NaCl molal")
ax.plot(T, deltaCO2_nacl4, label=f"4 NaCl molal")
ax.plot(T, deltaCO2_nacl6, label=f"6 NaCl molal")
ax.legend(loc="upper right")
ax.grid(True)
ax.set(xlabel='Temperature [°C]')
ax.set(ylabel='Solubility [mol/kgw]')
ax.set(title='Solubility of CO2 in NaCl brine, P = ' + str(P) + ' bar')

# We see the illustration of the so-called salting-out effect. It indicates lower solubility of the CO2(g) for more
# saline NaCl-brines.