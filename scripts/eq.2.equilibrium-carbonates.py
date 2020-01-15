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
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---



# # Performing a chemical equilibrium calculation of carbonate species
#
# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation with carbon species.
# We start by importing the `reaktoro` package:
#

from reaktoro import *

# ### Chemical system definition
#
# The default thermodynamic databases embedded into Reaktoro is SUPCRT92, so you do not have to initialize the
# database `db = Database('supcrt98.xml')`, unless you use an alternative one.

# Class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html)
# provides convenient operations to initialize
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) instances.

editor = ChemicalEditor()

# Alternatively, the editor can initialize from the file:

# Initialize a thermodynamic database with supcrt98.xml
db = Database("supcrt98.xml")
#  Define the editor of the chemical system
editor = ChemicalEditor(db)

# To define a chemical system, aqueous, gaseous, and mineral phases must be added. For aqueous phase, it can be done
# from a list of chemical element names. The database will be searched for all species that could be formed out of those
# elements.

editor.addAqueousPhaseWithElements("H O C Ca Cl Mg")

# To add gaseous phase, we provide a specific list of gaseous species that must be used:

editor.addGaseousPhase(["H2O(g)", "CO2(g)", "H2(g)", "O2(g)", "CH4(g)"])

# For the [MineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1MineralPhase.html) object, we add two pure mineral
# phases:

editor.addMineralPhase("Calcite")
editor.addMineralPhase("Dolomite")

# To create an object of class [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html)
# using the chemical system definition details stored in the object editor, we use

system = ChemicalSystem(editor)

# To output the details of the chemical system to the console, i.e., the list of elements, species, and different
# phases, we can execute

print(system)

# For the specific information about the chemical system, one can use the following methods:

print("Number of chemical species: ", len(system.species()))
print("Number of phases: ", len(system.phases()))
print("Number of elements: ", len(system.elements()))

print("List of species in chemical system: \n")
for species in system.species():
    print(species.name())

# ### Defining the chemical equilibrium problem
#
# Next, we create an equilibrium problem with our prescribed equilibrium conditions for
# amounts of elements that are consistent with our intention of calculating reaction of calcite with
# injected 0.002 molal MaCl<sub>2</sub> brine. Both temperature and pressure are assumed to be default values, i.e.,
# 25 &deg;C and 1 bar, respectively.

problem = EquilibriumProblem(system)
problem.add("H2O", 1, "kg")
problem.add("MgCl2", 0.002, "mol")
problem.add("CaCO3", 1, "mol")

# ### Calculating the chemical equilibrium state
#
# In this step, we use the `equilibrate()` function to calculate the chemical
# equilibrium state of the system with the given equilibrium conditions stored in the object problem.

state = equilibrate(problem)

# Reaktoro uses an efficient **Gibbs energy minimization** computation to determine the species amounts that correspond
# to a state of minimum Gibbs energy in the system, while satisfying the prescribed amount conditions for the
# temperature, pressure, and element amounts. The result is stored in the object `state`, of class
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html).

# To output the result of equilibration to the console, we use

print(state)

# We will obtain table describing the chemical state of the system. For example, the molar amounts, molar fractions,
# activities, activity coefficients, and chemical potentials of the species. The molar volumes of the phases,
# the amounts of each element in the phases, and also the total phase molar amounts.

# Alternatively, to save equilibrated state into a file, one can use method `output()`:

state.output("result.txt")  # Output the equilibrium state into a file result.txt

# To print the amounts of some specific aqueous species, one can use

# Print the amounts of some aqueous species
print("Amount of CO2(aq):", state.speciesAmount("CO2(aq)"))
print("Amount of HCO3-:", state.speciesAmount("HCO3-"))
print("Amount of CO3--:", state.speciesAmount("CO3--"))

# Similarly, one can also print the amounts of certain element, say carbon, in both aqueous and gaseous phases

# Print the amounts of element C in both aqueous and gaseous phases
print("Amount of C in aqueous phase:", state.elementAmountInPhase("C", "Aqueous"))
print("Amount of C in gaseous phase:", state.elementAmountInPhase("C", "Gaseous"))
print("Amount of C in calcite phase:", state.elementAmountInPhase("C", "Calcite"))
