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

# # Mass balance and mass action equations and related chemical properties
#
# In this tutorial, we clarify how to access certain basic properties of the chemical system and chemical equilibrium
# state, such as mass balance and mass action equations:

# We start by defining H<sub>2</sub>O-CO<sub>2</sub> chemical system defined as a mixture of 100 mol of H<sub>2</sub>0
# and 2 mols of CO<sub>2</sub> at T = 100 &deg;C and P = 50 bar:

# +
from reaktoro import *
db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O C")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])

system = ChemicalSystem(editor)

T = 100 # in celsius
P = 50 # in bar

problem = EquilibriumProblem(system)
problem.setPressure(P, "bar")
problem.setTemperature(T, "celsius")
problem.add("H2O", 100, "mol")
problem.add("CO2", 2, "mol")

# Equilibrate chemical problem
state = equilibrate(problem)
# -

# Obtain chemical species, chemical amounts, and formula matrix:

# +
b = state.elementAmounts()
n = state.speciesAmounts()
A = system.formulaMatrix()

print("b = ", b)
print("n = ", n)
print("A = ", A)
# -

# To evaluate the satisfaction of the mass balance equation

# +
# Import numpy package to work with arrays
import numpy as np

# Calculate the residual of the mass balance equation
r = b - np.dot(A, n)
# Calculate the norm of the residual
r_norm = np.linalg.norm(r)
print("||r|| = ", r_norm)
# -

# How much of the CO<sub>2</sub>(g) is dissolved as CO<sub>2</sub>(aq)?

print(f"CO2(aq) amount is {state.speciesAmount('CO2(aq)'):6.4e} mol")

# How much of the H<sub>2</sub>O(l) has evaporated as H<sub>2</sub>O(g)?

print(f"H2O(g) amount is {state.speciesAmount('H2O(g)'):6.4e} mol")

# What is the amount of H<sup>+</sup> species?

print(f"H+ amount is {state.speciesAmount('H+'):6.4e} mol")

# A nicer output of the formula matrix (where one can control the spacing and format):

rows, cols = A.shape
for i in range(rows):
    for j in range(cols):
        print(f"{A[i][j]:4.0f}", end="")
    print("\n")

# Rank is the maximal number of linearly independent columns of A,
# and it is equal to the dimension of the vector space spanned by its rows.

rank = np.linalg.matrix_rank(A)
print("Rank of A is", rank)

# Which phases exist in the equilibrium state?

# Collect phases names in the list
phases_names = [phase.name() for phase in system.phases()]
# Fetch stability indices from the state
stability_indices = state.phaseStabilityIndices()
print("Phase   : Phase amounts : Stability indices")
for name, si in zip(phases_names, stability_indices):
    print(f"{name:>7} : {state.phaseAmount(name):13.4f} : {si:6.4e}")

# > **Note**: The stability index is
# > * *zero* (very close to zero) if the phase is stable,
# > * *negative* if the phase is under-saturated, and
# > * *positive* if the phase is over-saturated.

# To access the molar masses of the elements in the system and evaluate their mass:

# Collect elements names and molar masses
element_names, molar_masses = zip(*[(element.name(), element.molarMass()) for element in system.elements()])
print("\nElement : Molar mass (g/mol) : Mass (g)")
for name, molar_mass, amount in zip(element_names, molar_masses, b):
    print(f"{name:>7} : {molar_mass * 1e3:18.2e} : {molar_mass * 1e3 * amount:8.2e}")

# To evaluate the mass of the species:

# +
# Fetch species names
species_names = [speices.name() for speices in system.species()]

# Species counter
i = 0
print("\nSpecies  : Amount (mol) :   Mass (g)")
for name, amount in zip(species_names, n):
    # Calculate species molar mass as the multiplication of the formula matrix column and element molar masses (in g)
    species_molar_mass = np.dot(A[:, i], molar_masses) * 1e3
    # Calculate species mass
    mass = amount * species_molar_mass

    print(f"{name:>8} : {amount:12.4f} : {mass:9.4f}")

    # Increase the species counter
    i += 1
# -

# To evaluate the properties (i.e., chemical potentials, logarithms of activities) of the system:

# +
TKelvin = T + 273.15 # in kelvin
PPascal = P * 1e5
props = system.properties(TKelvin, PPascal, n)
print("\nChemical potentials of the species:")
for mu, species, index in zip(props.chemicalPotentials().val,
                              system.species(),
                              list(range(1, system.numSpecies()+1))):
    print(f"\u03BC_{index} ({species.name():>8}) = {mu:12.4f} (J/mol)")

print("\nLogarithms of activities of the species:")
for lna, species, index in zip(props.lnActivities().val,
                              system.species(),
                              list(range(1, system.numSpecies()+1))):
    print(f"ln(a_{index} ({species.name():>8}) = {lna:8.4f}")
# -

# To evaluate equilibrium constants for the reactions:

# Initialize reaction equations
equations = ["H2O(l) = H+ + OH-",
             "HCO3- + H+ = CO2(aq) + H2O(l)",
             "H2O(l) + CO2(aq) = CO3-- + 2*H+",
             "CO2(aq) = CO2(g)"]
# Initialize reactions
reactions = [Reaction(ReactionEquation(equation), system) for equation in equations]
# Fetch equilibrium constants for each reaction
lnKs = [reaction.lnEquilibriumConstant(props) for reaction in reactions]
print("\nEquilibrium constants of reactions:")
for equation, reaction, lnK in zip(equations,
                              reactions,
                              lnKs):
    print(f"lnK ({equation:>32}) = {lnK.val:6.4f}")

# To control whether these constants correspond to the definition via the standard chemical potential,
# let us consider the equation `H2O(l) = H+ + OH-`:

# +
# Standard chemical potentials
mu0_H = 0.0
mu0_H2O = -242.992 * 1e3
mu0_OH = -155.559 * 1e3

R = 8.314 # J / (mol * K)
lnK = - 1 / R / TKelvin * (mu0_OH + mu0_H- mu0_H2O)

print("\nEquilibrium constants via standard chemical potentials:")
print("lnK (H2O(l) = H+ + OH-) = ", lnK)
