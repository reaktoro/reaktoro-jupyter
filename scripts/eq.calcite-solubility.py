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

# # Calcite solubility in water and CO2-saturated rainwater
#
# In this tutorial, we investigate the dependence of calcite solubility in water (closed system) and carbon-dioxide
# saturated rainwater (open system) on temperature and pressure change.
#
# We begin by importing necessary packages:

from reaktoro import *
import numpy as np
import matplotlib.pyplot as plt

# Function below defines the chemical problem corresponding to the pure water (also referred as a closed system):

def water_problem(system, T, P):

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")

    return problem

# Function `rainwater_problem()` defines the chemical problem corresponding to rainwater saturated with carbon-dioxide:

def rainwater_problem(system, T, P):

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    # Rainwater composition
    problem.add("H2O", 1.0, "kg")
    problem.add("Na", 2.05, "mg") # Sodium, 2.05 ppm = 2.05 mg/L ~ 2.05 mg/kgw
    problem.add("K", 0.35, "mg") # Potassium
    problem.add("Ca", 1.42, "mg") # Calcium
    problem.add("Mg", 0.39, "mg") # Magnesium
    problem.add("Cl", 3.47, "mg") # Chloride
    problem.add("SO4", 2.19, "mg")
    problem.add("NO3", 0.27, "mg")
    problem.add("NH4", 0.41, "mg")
    #  Saturated with carbon dioxide
    problem.add("CO2", 0.36, "mol")  # amount of carbon dioxide to saturate water

    return problem

# Finally, function `water_co2_problem` defines the chemical problem corresponding to an open system, i.e.,
# the carbon-dioxide saturated water with the partial pressure of the atmosphere pCO2 = 3.408. To convert the amount
# of CO2 gas from millimeters of mercury to parts per million (ppm) we use the instruction on the
# [following website](https://sciencing.com/calculate-ppm-vapor-pressure-6457861.html). In particular, pCO2 = 3.408
# corresponds to (3.408 / 760) * 106 = 0.475 ppm, which converts to mol/L by formula 0.475 / 44.01 = 0.010793 mol/L,
# where 44.01 g/mol is the CO2 molar amount.

def water_co2_problem(system, T, P):

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")
    problem.add("CO2", 0.010793, "mol")

    return problem

# Function that numerically models chemical equilibrium of 10 mol calcite with a given solute (defined by the input
# instance `problem`):

def solubility_of_calcite(problem, system):

    n0Calcite = 10.0

    # Define solver to calculate equilibrium and the initial chemical state
    solver = EquilibriumSolver(system)
    state = ChemicalState(system)

    # Equilibrate the solution given by the chemical problem `problem`
    solver.solve(state, problem)

    # Add `n0Calcite` amount of calcite
    state.setSpeciesAmount("Calcite", n0Calcite, "mol")

    # Equilibrate solution with added calcite
    solver.solve(state)

    # Calculate ph of the current state
    evaluate_pH = ChemicalProperty.pH(system)
    print(f"P = {problem.pressure() * 1e-5:.1f} bar, "
          f"T = {problem.temperature() - 273.15} C: "
          f"ph = {evaluate_pH(state.properties()).val:.2f}")

    # Fetch the amount of final calcite in the equilibrium state
    nCalcite = state.speciesAmount("Calcite")

    # Return the difference between the initially added and remaining calcite
    return n0Calcite - nCalcite

# Initialize chemical system with aqueous, gaseous, and calcite phases:

db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O C Ca Mg K Cl Na S N")
editor.addGaseousPhase(["CO2(g)"])
editor.addMineralPhase("Calcite")
system = ChemicalSystem(editor)

# Initialize the array of temperatures from 20 &deg;C till 90 &deg;C:

temperatures = np.arange(20.0, 91.0, 5.0)

# Calculate solubilities of calcite in water and CO2-saturated rainwater for pressure P = 1 bar and save it in txt-file:

P = 1.0 # in bar
print(f"Solubility in water (closed system):")
delta_calcite_water_P1 = [solubility_of_calcite(water_problem(system, T, P), system) for T in temperatures]
print(f"Solubility in rainwater:")
delta_calcite_rainwater_P1 = [solubility_of_calcite(rainwater_problem(system, T, P), system) for T in temperatures]
np.savetxt('reaktoro-water-delta-calcite-p-' + str(P) + '.txt', delta_calcite_water_P1)
np.savetxt('reaktoro-rainwater-delta-calcite-p-' + str(P) + '.txt', delta_calcite_rainwater_P1)

# Check if the obtained value for 25 &deg;C (second element of the list `delta_calcite_water_P1`) corresponds to the
# values of Wikipedia, i.e., the solubility in water equals to 0.013 g/L (25 &deg;C). Below, we use 100.0869 g/mol as
# the calcite molar mass:

print(f"Solubility of calcite in water (closed system) equals to {delta_calcite_water_P1[1]:.6f} mol/kgw = ... = "
      f"{delta_calcite_water_P1[1] * 0.1000869 * 1e3:.6f} g/L")

# Thus, in the closed system (with the pure water), approximately 0.124 mmol of calcite dissolve.
# The amount of calcite that dissolves is independent of the initial value (provided that it exceeds the solubility
# limit).

# Calculate solubilities of calcite in water and CO2-saturated rainwater for pressure P = 100 bar and save it in
# txt-file:

P = 100.0 # in bar
print(f"Solubility in water (closed system):")
delta_calcite_water_P100 = [solubility_of_calcite(water_problem(system, T, P), system) for T in temperatures]
print(f"Solubility in rainwater (open system):")
delta_calcite_rainwater_P100 = [solubility_of_calcite(rainwater_problem(system, T, P), system) for T in temperatures]
np.savetxt('reaktoro-water-delta-calcite-p-' + str(P) + '.txt', delta_calcite_water_P100)
np.savetxt('reaktoro-rainwater-delta-calcite-p-' + str(P) + '.txt', delta_calcite_rainwater_P100)

# Plot solubilities of calcite in water and CO2-saturated rainwater for pressure P = 1:

fig, ax = plt.subplots()
ax.plot(temperatures, delta_calcite_water_P1, label=f"Calcite in water", color='C1')
ax.plot(temperatures, delta_calcite_rainwater_P1, label=f"Calcite in rainwater", color='C3')
ax.legend(loc="upper right")
ax.set(xlabel=r'Temperature [$^\circ$C]',
       ylabel='Solubility [mol/kgw]',
       title='Comparison of calcite solubility')
ax.grid()
fig.savefig('calcite-solubility-P-1bar.png')

# The plot illustrates that calcium carbonate has very low solubility in pure water, but in rainwater saturated with
# carbon dioxide, its solubility increases due to the formation of more soluble calcium bicarbonate. Calcium carbonate
# is unusual in that its solubility increases as the temperature of the rainwater decreases.
#
# Let us now plot solubilities on the different scales and for different pressure:

# +
fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.plot(temperatures, delta_calcite_water_P1, label="P = 1", color='C1')
ax1.plot(temperatures, delta_calcite_water_P100, label="P = 100", color='C2')
ax1.legend(loc="best")
ax1.set_title('Calcite in water')
ax1.grid()
ax1.set_ylabel('Solubility [mol/kgw]')

ax2.plot(temperatures, delta_calcite_rainwater_P1, label=f"P = 1", color='C3')
ax2.plot(temperatures, delta_calcite_rainwater_P100, label=f"P = 100", color='C4')
ax2.set_title('Calcite in rainwater')
ax2.legend(loc="best")
ax2.grid()
ax2.set_ylabel('Solubility [mol/kgw]')
ax2.set_xlabel(r'Temperature [$^\circ$C]')

fig.tight_layout()
fig.savefig('calcite-solubility.png', bbox_inches='tight')
# -

# We see that increasing pressure also increases the solubility of calcium carbonate.

# Finally, the solubility of calcite in the open system is simulated using the function `water_co2_problem()`,
# representing pure water and CO<sup>2</sup>. The pressure corresponding to the partial pressure of the atmosphere
# pCO2 = 3.408 is equal to 39 = $\mathsf{10^{-3.408}}$ Pa.

P = 39 * 1e-5 # in bar
print(f"Solubility in water and CO2 (open system, P = 39 Pascal):")
delta_calcite_water_co2_P1 = [solubility_of_calcite(water_co2_problem(system, T, P), system) for T in temperatures]
print(f"Solubility of calcite in open system equals to {delta_calcite_water_co2_P1[1]:.6f} mol/kgw = ... = "
   f"{delta_calcite_water_co2_P1[1] * 0.1000869 * 1e3:.6f} g/L")

# We see that the solubility of calcite in the open system is about four times higher than in the closed system (
# 0.000124 mol/kgw or 0.012455 g/L).
