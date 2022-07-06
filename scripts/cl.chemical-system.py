# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
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

# # Functionality of ChemicalSystem class
#
# In this tutorial, we provide an explanation on the functionality of the class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html)
# that represents a chemical system and its attributes and properties. Below,
# we provide a tutorial of the methods that can be used to study the characteristics of the considered chemical system.
#
# A chemical system is a description of the phases of interest in the modeling problem and the chemical species that
# compose those phases. For example, when modeling the chemistry of an aqueous solution, only one phase should be
# enough, an *aqueous phase*. If one is interested in modeling the solubility of gases in an aqueous solution,
# then it makes sense to also define a *gaseous phase* with one or more gases. When modeling aqueous solutions and
# minerals, under a variety of temperature, pressure, and chemical conditions, it might make sense to define the
# chemical system with many *mineral phases*.
# Assume that we have defined `system`, an instance of
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) object, by the code below:

# +
# Import the reaktoro Python package
from reaktoro import *

# Initialize a thermodynamic database
db = Database("supcrt98.xml")

# Define the chemical system
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements(
    "H O Na Cl C"
)  # add aqueous phase by the all possible combination of provided elements
editor.addGaseousPhase(["CO2(g)"])  # add one gaseous species

# Construct the chemical system
system = ChemicalSystem(editor)
# -

# > Check previous tutorials to learn the  steps above! For example,
# > [Basics of equilibrium calculation](eq.equilibrium-basics.ipynb),
# > [Equilibrium calculation of carbonate species](eq.equilibrium-carbonates.ipynb),
# > [Equilibrium calculations using equilibrium solver](eq.co2-brine-using-equilibrium-solver.ipynb), or
# > [Custom activity model for equilibrium calculations](eq.custom-activity-models.ipynb).
#
# The most general print-out of the chemical system can be done by

print(system)

# However, to access some limited information, such as the list of species, phases, or elements separately, the following
# code can be used

import numpy as np

print(f"List with {system.numSpecies()} species:")
for species in system.species():
    print(species.name())

print(f"List with {system.numPhases()} phases:")
for phase in system.phases():
    print(phase.name())

print(f"List with {(system.numElements())} elements:")
for element in system.elements():
    print(element.name())

# To output additional information about phases, for instance, one can use

print("List of phases with number of species in each phase:")
for phase in system.phases():
    print(f" * Phase {phase.name()} contains {phase.numSpecies()} species")

# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html)
# provides the formula matrix (whose entry *(j, i)* is given by the number of atoms of its *j*th element in
# its *i*th species). To access it, one need to use

matrix = system.formulaMatrix()
print(f"Formula matrix of the size {matrix.shape[0]} x {matrix.shape[1]}:\n", matrix)

# For instance, the matrix printed above is the matrix of the size *6 x 23*, i.e.,
#
# \begin{bmatrix}
# 1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
# 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
# 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 2 & 2 & 2 & 1 & 1 & 1 & 1 & 1 & 0 & 0 & 1 & 0 & 1 & 0 \\
# 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 0 \\
# 1 & 2 & 3 & 0 & 1 & 2 & 3 & 4 & 0 & 0 & 1 & 2 & 3 & 0 & 1 & 2 & 2 & 0 & 0 & 1 & 2 & 1 & 2 \\
# 0 & 0 & -2 & -1 & -1 & -1 & -1 & -1 & 1 & 0 & 0 & 0 & -1 & 0 & 0 & 0 & -1 & 1 & 0 & 0 & 0 & -1 & 0 \\
# \end{bmatrix}.
#
# Assume the enumeration starts from 0 (which is the case in Python language). Here, the row with index 1 corresponds to
# the element Cl, which is only present in species
# Cl<sup>-</sup>,
# ClO<sup>-</sup>,
# ClO<sup>2-</sup>,
# ClO<sup>3-</sup>,
# ClO<sup>4-</sup>
# (with indices 3, 4, 5, 6, and 7),
# HCl(aq),
# HClO(aq),
# HClO_2(aq)
# (with indices 13, 14, and 15),
# and
# NaCl(aq) (with index 18).
# Since only one atom of Cl contributes to each species, the second row contains only 1 in non-zero values.
#
# To get the index of certain element, phase, or species, functions `index__()` or `index__WithError()` (the latter
# results in system throwing an exception if the element does not exist), where instead of `__` one can used `Element`,
# `Phase`, or `Spieces`:

print("Index of the element H: ", system.indexElement("H"))
print("Index of the phase Aqueous: ", system.indexPhase("Aqueous"))
print("Index of the species Cl-: ", system.indexSpecies("Cl-"))

# When working with a  set of species, one can request a set of corresponding indices. Let us collect all the species
# with chlorine and retrieve indices of corresponding species:

species = [
    "Cl-",
    "ClO-",
    "ClO2-",
    "ClO3-",
    "ClO4-",
    "HCl(aq)",
    "HClO(aq)",
    "HClO2(aq)",
    "NaCl(aq)",
]
print("Indices of species with Cl: ", system.indicesSpecies(species))

# They must correspond with positions of the non-zero elements in the row 1 of formula matrix discussed-above.
# Having the instance of chemical system, we can calculate the molar amounts of the elements (in units of mol), when the
# molar amount of species is provided:

n = np.ones(system.numSpecies())
elements_amount = system.elementAmounts(n)
hydrogen_amount = system.elementAmount(2, n)
print(
    "Element amounts (in mol) provided 1 molal for all species: ",
    elements_amount,
)
print(
    "Hydrogen amounts (in mol) provided 1 molal for all species: ",
    hydrogen_amount,
)

# To study the chemical system even further, one can access the class
# [ThermoProperties](https://reaktoro.org/cpp/classReaktoro_1_1ThermoProperties.html) by providing temperature and
# pressure, i.e.,

T = 60
P = 100
thermo_properties = system.properties(T, P)

# Object `thermo_properties` contains information about
# * standard partial molar Gibbs energies of the species (in units of J/mol),
# * standard partial molar enthalpies of the species (in units of J/mol),
# * standard partial molar volumes of the species (in units of m3/mol),
# * standard partial molar entropies of the species (in units of J/(mol*K)),
# * standard partial molar internal energies of the species (in units of J/mol),
# * standard partial molar Helmholtz energies of the species (in units of J/mol),
# * standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)),
# * standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
#
# For instance, partial molar Gibbs energies or enthalpies can be accessed as follows:

print("List of standard partial molar Gibbs energies of the species:")
for energies, species in zip(
    thermo_properties.standardPartialMolarGibbsEnergies().val,
    system.species()
):
    print(f"\u03B4G ({species.name():>10}) = {energies}")

print("List of standard partial molar enthalpies of the species:")
for enthalpies, species in zip(
    thermo_properties.standardPartialMolarEnthalpies().val,
    system.species()
):
    print(f"\u03B4H ({species.name():>10}) = {enthalpies}")

# Alternatively, by providing the vector with molar amounts of the species (in units of mol) class
# [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) can be accessed:

chemical_properties = system.properties(T, P, n)

# This object contains various chemical properties, such as mole fractions, the logarithm of activities, chemical
# potentials of the species, in addition to thermodynamic properties listed above:

print("Chemical potentials of the species:")
for potential, species, index in zip(
    chemical_properties.chemicalPotentials().val,
    system.species(),
    list(range(1, system.numSpecies()+1))
):
    print(f"\u03BC_{index} ({species.name():>10}) = {potential}")

print("Logarithms of activities of the species:")
for activity, species, index in zip(
    chemical_properties.lnActivities().val,
    system.species(),
    list(range(1, system.numSpecies()+1))
):
    print(f"ln (a_{index}) ({species.name():>10}) = {activity}")


# ### Definition of chemical system using GEMs
#
# Chemical systems can be initialize using alternative backends, such as GEMs and PHREEQC, using corresponding classes
# [Gems](https://reaktoro.org/cpp/classReaktoro_1_1Gems.html) and
# [Phreeqc](https://reaktoro.org/cpp/classReaktoro_1_1Phreeqc.html)
# `gems = Gems("some-gems-project-file.lst")` or `phreeqc = Phreeqc("phreeqc.dat")`,
# which allows to define instance of [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html)
# class by respective `system = ChemicalSystem(gems)` or `system = ChemicalSystem(phreeqc)`.
