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

# # Calculation of ionic strength and activity coefficients of aqueous species
#
# In this tutorial, we explain how to access the ionic strength of the equilibrium state as well as clarify how it
# can be calculated manually by accessing specific properties of the chemical state. We also calculate the activity
# coefficients of aqueous species and solvent water.
#
# First, we import all the required python packages:

from reaktoro import *
import numpy as np
from math import *

# Define database, initialize phases with chemical editor, and create chemical system:

database = Database("supcrt98.xml")
editor = ChemicalEditor(database)
editor.addAqueousPhaseWithElements("H O C Na Cl Ca").setChemicalModelDebyeHuckel()
editor.addGaseousPhase(["H2O(g)"])
system = ChemicalSystem(editor)
print(system)

# Set thermodynamic conditions:

T = 25 # in celsius
P = 1 # in bar

# First problem simulates mixing sodium chlorite with water:

problem1 = EquilibriumProblem(system)
problem1.setTemperature(T, "celsius")
problem1.setPressure(P, "bar")
problem1.add("H2O", 1.0, "kg")
problem1.add("NaCl", 1.0, "mol")

# The second one demonstrates experiment of mixing water with CaCl<sub>2</sub>:

problem2 = EquilibriumProblem(system)
problem2.setTemperature(T, "celsius")
problem2.setPressure(P, "bar")
problem2.add("H2O", 1.0, "kg")
problem2.add("CaCl2", 1, "mol")

# Equilibration of the water mixed 1 mol of sodium chlorite and water mixed with 1 mol of CaCl<sub>2</sub> results into
# the following two states:

state1 = equilibrate(problem1)
state2 = equilibrate(problem2)

# To evaluate the ionic strength, we need to define a corresponding function `ionic_strength()`

# +
ionic_strength = ChemicalProperty.ionicStrength(system)
I1 = ionic_strength(state1.properties()).val
I2 = ionic_strength(state2.properties()).val

print(f"Ionic strength of state1 is {I1:f} molal")
print(f"Ionic strength of state2 is {I2:f} molal")
# -

# We see that the ionic strength of the second mix is higher, which can be explained by the fact that CaCl<sub>2</sub>
# contains two ions of Cl<sup>-</sup>.
#
# ## Calculating the ionic strength
#
# Below, we explain, which information one needs to fetch from chemical state to be able to calculate ionic strength. \
# First, let us fix the molality of 1 kg of solvent water (where 18.0154 * 1e-3 kg/mol is a molar mass of water):

mw_h2o = 1 / 18.0154 / 1e-3
print(f"mw_h2o = {mw_h2o:f} molal")

# Next we collect the list of species and their amounts in the chemical system:

species = system.species()
n2 = state2.speciesAmounts()

# Since out of two phases (aqueous and gaseous) we need to consider only species from the aqueous phase,
# we fetch indices of the aqueous species.

indx_aqueous_phase = system.indexPhase("Aqueous")
indx_gaseous_phase = system.indexPhase("Gaseous")
indx_aqueous_species = system.indicesSpeciesInPhases([indx_aqueous_phase])
indx_all_species = system.indicesSpeciesInPhases([indx_aqueous_phase, indx_gaseous_phase])
print(f"Indices of aq. species:", indx_aqueous_species)
print(f"Indices of all species:", indx_all_species)

# We see that the difference between these two lists is only in the last index.
#
# Out the list `species` (with all species), we collect lists of only aqueous species, their amounts, names, and
# corresponding charges:

species_aq = [species[i] for i in indx_aqueous_species]
n2_aq = [n2[i] for i in indx_aqueous_species]
names_aq = [species.name() for species in species_aq]
z_aq = [species.charge() for species in species_aq]

# Amount of the water is obtained by:

n_h2o = state2.speciesAmount("H2O(l)")

# Next, we calculate the molalities of aqueous species and print their names, charges, and molalities:

m_aq = mw_h2o * np.divide(n2_aq, n_h2o)
print(f"   Species : charge, molalities")
for name, Z, m in zip(names_aq, z_aq, m_aq):
    print(f"{name:>10} : {Z:6.0f}, {m:6.2e} molal")

# The ionic strength can be calculated by:

I2 = 1/2 * sum([m * Z**2 for m, Z in zip(m_aq, z_aq)])
print(f"Ionic strength of state2 is {I2:f} molal (calculated manually)")

# ## Calculating the activity coefficients for aqueous ionic species (Davis model)
#
# Calculating and outputting the activity coefficients for aqueous ionic species is done by:

A_gamma = 0.5095
gammas = [10**(-A_gamma * z**2 * (sqrt(I2) / (1 + sqrt(I2)) - 0.3 * I2)) for z in z_aq]
print(f"   Species : Activity coefficients")
for name, gamma in zip(names_aq, gammas):
    print(f"{name:>10} : {gamma:2.4f}")

# We see that many of the activity coefficients are away from $\gamma_i$ = 1 (which corresponds to an ideal solution).
#
# ## Calculating activity of the water solvent
#
# To calculate the activity of the water solvent, we need fractions of the species, which are stored in the class
# [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) accessed from chemical
# state by the method `properties()`.

properties = state2.properties()
fractions = properties.moleFractions().val

# Let us output only those species that have fractions bigger than machine precision (set here to 10<sup>-16</sup>):

machine_precision = 1e-16
print(f"   Species : Mole fractions")
for name, x in zip(names_aq, fractions):
    if x > machine_precision:
        print(f"{name:>10} : {x:6.4e}")

# We see that solvent water possesses the biggest fraction as well as CaCl<sub>2</sub>(aq) and ions Ca<sup>2+</sup>,
# CaCl<sup>+</sup>, and Cl<sup>-</sup> (caused by addition of CaCl<sub>2</sub> to the water). The fraction of solvent
# water
# can be accessed via index of this species:

indx_h2o = system.indexSpecies("H2O(l)")
x_h2o = fractions[indx_h2o]
print(f"Index of the water solvent is {indx_h2o}")
print(f"Fraction of the water solvent is {x_h2o:6.4f}")

# Finally, we calculate activity coefficient of solvent water by:

ln10 = numpy.log(10.0)
sqrtI2 = numpy.sqrt(I2)
gamma_h2o_davis = exp(ln10/55.5084*A_gamma*(2*(I2 + 2*sqrtI2)/(1 + sqrtI2) - 4 * numpy.log(1 + sqrtI2) - 0.3 * I2**2)
                      - (1 - x_h2o)/x_h2o)
gamma_h2o_ideal = exp(- (1 - x_h2o)/x_h2o)
print(f"Activity coefficient of water solvent (Davis model) is {gamma_h2o_davis:6.4f}")
print(f"Activity coefficient of water solvent (ideal model) is {gamma_h2o_ideal:6.4f}")

# ## Demonstration of Coulomb’s law
#
# According to Coulomb’s law, the activity coefficient decreases as the concentration increases because the
# electrostatic forces become stronger as the ions approach. Thus, for more concentrated solutions, the repulsion
# effect seems to dominate. Let us demonstrate how it can be seen in Reaktoro simulations. First, we access the
# activity coefficients of the `state2` via its properties obtained earlier:

gamma_1_mol = np.exp(properties.lnActivityCoefficients().val)

# Next, we increase the concentration of CaCl<sub>2</sub> in the mixture and recalculate activity coefficients:

problem2.add("CaCl2", 2, "mol")
state2 = equilibrate(problem2)
properties = state2.properties()
gamma_2_mol = np.exp(properties.lnActivityCoefficients().val)

print(f"Species with decreased activity coeffs. after adding more CaCl2 to the water:")
for name, gamma_1_mol, gamma_2_mol  in zip(names_aq, gamma_1_mol, gamma_2_mol):
    if gamma_1_mol > gamma_2_mol:
        print(f"{name:>10} : {gamma_1_mol:6.4e} -> {gamma_2_mol:6.4e}")
