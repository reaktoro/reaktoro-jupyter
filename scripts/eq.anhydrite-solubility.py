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

# # Gypsum/anhydrite solubility in water
#
# In this tutorial, we investigate the dependence of the sulfate mineral anhydrite (CaSO<sub>4</sub>) solubility in
# water for different temperatures and pressures.
#
# > **Note**: In the databases available in Reaktoro, CaSO<sub>4</sub> is always referred to as
# > anhydrite. When exposed to water, anhydrite transforms by the absorption of water to the more commonly
# > known/occurring form gypsum (CaSO<sub>4</sub> · 2H<sub>2</sub>O).
#
# First, we import necessary packages for the presented simulations:

from reaktoro import *
import numpy as np
import math as math
import matplotlib.pyplot as plt

# Below, function `water_problem()` defines the chemical problem corresponding to pure water (closed system):

def water_problem(system, T, P):

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")

    return problem

# Function `solubility_of_anhydrite()` simulates chemical equilibrium of 5 mol anhydrite with given solute (defined by the
# instance `problem`):

def solubility_of_anhydrite(problem, system, reaction):

    n0Anhydrite = 5.0

    solver = EquilibriumSolver(system)
    state = ChemicalState(system)

    # Equilibrate pure water
    solver.solve(state, problem)

    # Add 10mol of anhydrite
    state.setSpeciesAmount("Anhydrite", n0Anhydrite, "mol")

    # Equilibrate pure water with anhydrite
    solver.solve(state)

    # Calculate ph of the current state
    evaluate_I = ChemicalProperty.ionicStrength(system)
    I = evaluate_I(state.properties()).val

    # Calculate ph of the current state
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(state.properties()).val

    # Fetch chemical properties
    props = state.properties()
    # Calculate equilibrium constant
    lnK = reaction.lnEquilibriumConstant(props)
    # Calculate reaction quotient
    lnQ = reaction.lnReactionQuotient(props)

    # Calculate saturation ratio
    lnSR = lnQ.val - lnK.val
    # Calculate saturation index as log10(SR)
    SI = lnSR / math.log(10)

    print(f"P = {problem.pressure() * 1e-5:.1f} bar, T = {problem.temperature() - 273.15} C: "
          f"ph = {pH:.2f}, I = {I * 1e3:.2f} mmolal, lnK = {lnK.val:.4f}, SI = {SI:e}")

    # Fetch the amount of final anhydrite in the equilibrium state
    nAnhydrite = state.speciesAmount("Anhydrite")

    return n0Anhydrite - nAnhydrite

# Next, we initialize chemical system with aqueous, gaseous, and anhydrite phases:

db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Ca S")
editor.addMineralPhase("Anhydrite")
system = ChemicalSystem(editor)

# The anhydrite chemical reaction is defined based on the reaction equation and chemical system:

reaction = Reaction(ReactionEquation("Anhydrite = SO4-- + Ca++"), system)

# Temperatures from 0 till 90 &deg;C are initialized with the list:

temperatures = np.arange(0.0, 91.0, 5.0)

# Below, we calculate solubilities of anhydrite in water for the pressure P = 1 bar and save it in txt-file:

P = 1.0 # in bar
print(f"Solubility in water, P = 1 bar:")
delta_anhydrite_water_P1 = [solubility_of_anhydrite(water_problem(system, T, P), system, reaction) for T in temperatures]
np.savetxt('reaktoro-water-delta-anhydrite-p-' + str(P) + '.txt', delta_anhydrite_water_P1)

# Inside the function `solubility_of_anhydrite()`, we also evaluate ph, ionic strength I, equilibrium constant lnK,
# and saturation index SI. If the SI > 0, the solution is supersaturated with anhydrite, whereas if the SI < 0,
# the solution is undersaturated with it. Finally, SI = 0 indicates equilibrium. The obtained solubility in
# water is 2.308049 g/L (25 &deg;C), which is considerably higher than for calcite. We use 172.17 g/mol as a
# molar mass of anhydrite:

print(f"Solubility of anhydrite in water = {delta_anhydrite_water_P1[1]:.6f} mol/kgw = ... = "
      f"{delta_anhydrite_water_P1[1] * 0.17217 * 1e3:.6f} g/L")

# According to the calculated values, the lnK at the 25 &deg;C is lnK = -9.9159, which corresponds to the database
# values of `log_k = -4.36` (see, for example, `phreeqc.dat`)

lnK = -9.9159
log10K = lnK * math.log10(math.exp(1))
print("Anhydrite/gypsum logK = ", log10K)

# Calculate solubilities of anhydrite in water for the pressures 100 and 1000 bar and save it in txt-file:

# +
P = 100.0 # in bar
print(f"Solubility in water, P = 100 bar:")
delta_anhydrite_water_P100 = [solubility_of_anhydrite(water_problem(system, T, P), system, reaction) for T in temperatures]
np.savetxt('reaktoro-water-delta-anhydrite-p-' + str(P) + '.txt', delta_anhydrite_water_P100)

P = 1000.0 # in bar
print(f"Solubility in water, P = 1000 bar:")
delta_anhydrite_water_P1000 = [solubility_of_anhydrite(water_problem(system, T, P), system, reaction) for T in
                               temperatures]
np.savetxt('reaktoro-water-delta-anhydrite-p-' + str(P) + '.txt', delta_anhydrite_water_P1000)
# -


# Let us now plot solubilities of anhydrite as function of different temperatures for different pressure:

# +
fig, ax = plt.subplots()

ax.plot(temperatures, delta_anhydrite_water_P1, label="P = 1", color='C3')
ax.plot(temperatures, delta_anhydrite_water_P100, label="P = 100", color='C4')
ax.plot(temperatures, delta_anhydrite_water_P1000, label="P = 1000", color='C5')
ax.legend(loc="best")
ax.set_title(r'Anhydrite solubility in water')
ax.grid(True)
ax.set_ylabel('Solubility [mol/kgw]')
ax.set_xlabel(r'Temperature [$^{\circ}$C]')
fig.savefig('anhydrite-solubility.png', bbox_inches='tight')
# -

# We see that anhydrite solubility decreases with increasing temperature.
# Increasing pressure also increases the solubility of calcium sulfate.
