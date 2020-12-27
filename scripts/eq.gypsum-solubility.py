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

# # Gypsum solubility in water
#
# In this tutorial, we investigate the gypsum (CaSO<sup>4<\sup>) solubility in water. In the databases available in
# Reaktoro, CaSO<sup>4<\sup> is always referred to as anhydrite. When exposed to water, anhydrite transforms
# to the more commonly known/occurring form gypsum, (CaSO<sup>4<\sup> ·2H<sup>2<\sup>O), by the absorption of water.
#
# Import necessary packages:

from reaktoro import *
import numpy as np
import math as math
import matplotlib.pyplot as plt

# Function that defines the chemical problem corresponding to pure water (closed system):

def water_problem(system, T, P):

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")

    return problem

# Function that simulates chemical equilibrium of 5 mol gypsum with given solute (defined by the instance `problem`):

def solubility_of_gypsum(problem, system, reaction):

    n0Gypsum = 5.0

    solver = EquilibriumSolver(system)
    state = ChemicalState(system)

    # Equilibrate pure water
    solver.solve(state, problem)

    # Add 10mol of gypsum
    state.setSpeciesAmount("Anhydrite", n0Gypsum, "mol")

    # Equilibrate pure water with gypsum
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
    SR = lnQ.val - lnK.val
    # Calculate saturation index
    SI = math.exp(SR)

    print(f"P = {problem.pressure() * 1e-5:.1f} bar, T = {problem.temperature() - 273.15} C: "
          f"ph = {pH:.2f}, I = {I * 1e3:.2f} mmolal, lnK = {lnK.val:.4f}, SI = {SI:.2f}")

    # Fetch the amount of final gypsum in the equilibrium state
    nGypsum = state.speciesAmount("Anhydrite")

    return n0Gypsum - nGypsum

# Initialize chemical system with aqueous, gaseous, and gypsum phases:

db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Ca S")
editor.addMineralPhase("Anhydrite")
system = ChemicalSystem(editor)

# Define gypsum chemical reaction based on the reaction equation and chemical system:

reaction = Reaction(ReactionEquation("Anhydrite = SO4-- + Ca++"), system)

# Initialize the array of temperatures from 20 &deg;C till 90 &deg;C:

temperatures = np.arange(20.0, 91.0, 5.0)

# Calculate solubilities of gypsum in water and CO2-saturated rainwater for pressure P = 1 bar and save it in txt
# file:

P = 1.0 # in bar
print(f"Solubility in water, P = 1 bar:")
delta_gypsum_water_P1 = [solubility_of_gypsum(water_problem(system, T, P), system, reaction) for T in temperatures]
np.savetxt('reaktoro-water-delta-gypsum-p-' + str(P) + '.txt', delta_gypsum_water_P1)

# Inside the function `solubility_of_gypsum()`, we also evaluate ph, ionic strength I, equilibrium constant lnK,
# adn saturation index SI. Since the SI > 0, the solution is supersaturated with gypsum.
# The obtained solubility in water is 2.308049 g/L (25 &deg;C), which is considerably higher than for calcite.
# Below, we use 172.17 g/mol as a molar mass of gypsum:

print(f"Solubility of gypsum in water = {delta_gypsum_water_P1[1]:.6f} mol/kgw = ... = "
      f"{delta_gypsum_water_P1[1] * 0.17217 * 1e3:.6f} g/L")

# According to the calculated values, the lnK at the 25 &deg;C is lnK = -9.9159, which corresponds to the database
# values of `log_k = -4.36` (see, for example, `phreeqc.dat`)

lnK = -9.9159
log10K = lnK * math.log10(math.exp(1))
print("Anhydrite/gypsum logK = ", log10K)

# Calculate solubilities of gypsum in water and CO2-saturated rainwater for pressure P = 100 bar and save it in txt
# file:

P = 100.0 # in bar
print(f"Solubility in water, P = 100 bar:")
delta_gypsum_water_P100 = [solubility_of_gypsum(water_problem(system, T, P), system, reaction) for T in temperatures]
np.savetxt('reaktoro-water-delta-gypsum-p-' + str(P) + '.txt', delta_gypsum_water_P1)

# Let us now plot solubilities on the different scales and for different pressure:

# +
fig, ax = plt.subplots()

ax.plot(temperatures, delta_gypsum_water_P1, label="P=1", color='C1')
ax.plot(temperatures, delta_gypsum_water_P100, label="P=100", color='C2')
ax.legend(loc="best")
ax.set_title('gypsum in water')
ax.grid()
ax.set_ylabel('Solubility [mol/kgw]')

fig.tight_layout()
fig.savefig('gypsum-solubility.png', bbox_inches='tight')

# We see again the solubility increase with decreasing temperature in CO2-saturated rainwater.
# Increasing pressure also increases the solubility of calcium carbonate.

# The plot illustrates that calcium carbonate has very low solubility in pure water, but in rainwater saturated with
# carbon dioxide, its solubility increases due to the formation of more soluble calcium bicarbonate. Calcium
# carbonate is unusual in that its solubility increases as the temperature of the water decreases.
#
