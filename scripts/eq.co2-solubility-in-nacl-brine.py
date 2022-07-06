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
# ---

# # Carbon dioxide gas solubility in the NaCl-brine
#
# This tutorial demonstrates how to simulate the solubility of CO<sub>2</sub> gas in the NaCl-brine and its dependence on
# the salinity of the brine (also called a salting-out effect).
#
# First, we import all the necessary packages for further simulations:

from reaktoro import *
import numpy as np
import matplotlib.pyplot as plt

# Function `solubility_co2()` returns the concentration of  CO<sub>2</sub>(g) that was dissolved in the brined:

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

# Below, we set up the chemical system for running the solubility calculations:

# +
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
P = 1.0
# Initial amount of CO2(g) (in mol)
n0CO2g = 10.0

# Generate the lists of CO2(g) mols that got dissolved in NaCl-brine of different salinity
deltaCO2_nacl1 = [solubility_co2(system, t, P, n0CO2g, mNaCl=1.0) for t in T]
deltaCO2_nacl2 = [solubility_co2(system, t, P, n0CO2g, mNaCl=2.0) for t in T]
deltaCO2_nacl4 = [solubility_co2(system, t, P, n0CO2g, mNaCl=4.0) for t in T]
deltaCO2_nacl6 = [solubility_co2(system, t, P, n0CO2g, mNaCl=6.0) for t in T]
# -

# Plot solubility of CO<sub>2</sub>(g) as a function of temperature for different salinities of NaCl-brine:

fig, ax = plt.subplots()
ax.plot(T, deltaCO2_nacl1, label=f"1 NaCl molal")
ax.plot(T, deltaCO2_nacl2, label=f"2 NaCl molal")
ax.plot(T, deltaCO2_nacl4, label=f"4 NaCl molal")
ax.plot(T, deltaCO2_nacl6, label=f"6 NaCl molal")
ax.legend(loc="best")
ax.grid(True)
ax.set(xlabel='Temperature [°C]')
ax.set(ylabel='Solubility [mol/kgw]')
ax.set(title='Solubility of CO2 in NaCl brine, P = ' + str(P) + ' bar')
plt.savefig('co2-solubility-nacl-h2o-vs-temperature-1bar.png', bbox_inches='tight')

# We see the illustration of the so-called salting-out effect. It indicates lower solubility of the CO<sub>2</sub>(g)
# for more saline NaCl-brines. Moreover, we see that the solubility of carbon dioxide decreases with the growth of the
# temperature.
#
# Alternatively, we can study the dependence of the CO<sub>2</sub>(g) solubility on the pressure increase:

# Initialize pressure range (in bar)
P = np.arange(1.0, 300.0, 1.0)
# Initialize pressure (in celsius)
T = 300.0
# Generate the lists of CO2(g) mols that got dissolved in NaCl-brine of different salinity
deltaCO2_nacl1 = [solubility_co2(system, T, p, n0CO2g, mNaCl=1.0) for p in P]
deltaCO2_nacl2 = [solubility_co2(system, T, p, n0CO2g, mNaCl=2.0) for p in P]
deltaCO2_nacl4 = [solubility_co2(system, T, p, n0CO2g, mNaCl=4.0) for p in P]
deltaCO2_nacl6 = [solubility_co2(system, T, p, n0CO2g, mNaCl=6.0) for p in P]

# Below, we plot solubility of CO<sub>2</sub>(g) as a function of pressure for different salinities of NaCl-brine:

fig, ax = plt.subplots()
ax.plot(P, deltaCO2_nacl1, label=f"1 NaCl molal")
ax.plot(P, deltaCO2_nacl2, label=f"2 NaCl molal")
ax.plot(P, deltaCO2_nacl4, label=f"4 NaCl molal")
ax.plot(P, deltaCO2_nacl6, label=f"6 NaCl molal")
ax.legend(loc="best")
ax.grid(True)
ax.set(xlabel='Pressure [bar]')
ax.set(ylabel='Solubility [mol/kgw]')
ax.set(title='Solubility of CO2 in NaCl brine, T = ' + str(T) + ' celsius')
plt.savefig('co2-solubility-nacl-h2o-vs-pressure-300celsius.png', bbox_inches='tight')
