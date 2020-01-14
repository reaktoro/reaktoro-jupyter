# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: py:light,../notebooks//ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Performing a chemical equilibrium calculation

# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation. We start by importing the `reaktoro` package:

from reaktoro import *

# Next, we need a thermodynamic database that enable us to compute thermodynamic properties of species and reactions. For this, we create an object of class [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html):

db = Database("supcrt98.xml")

# To indicate the phases of interest (as well as their species) that may potentially exist at equilibrium, we create an object of class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html):

editor = ChemicalEditor(db)

# We consider an aqueous phase composed of all aqueous species in the database that can be formed with the given chemical elements below:

editor.addAqueousPhaseWithElements("H O Na Cl C Ca Mg Si")

# > **Note:** This automatic selection of chemical species for a phase can result in a large number of them. This potentially increases the computing cost of the chemical reaction calculations. If you are using Reaktoro in demanding applications, you might want to manually specify the chemical species of each phase in your chemical system. This can be achieved by providing an explicit list of species names, e.g., `editor.addAqueousPhase("H2O(l) H+ OH- CO2(aq)")`. Note, however, that care is required here to ensure relevant species are not missing. The just given example is a bad one in fact, with important species such as `HCO3-` and `CO3--` missing in the list!.

# We are interested in a gaseous phase containing exactly the following gases (which may or not exist in positive amounts at the end of our equilibrium calculation later!):

editor.addGaseousPhase("H2O(g) CO2(g)")

# Finally, we consider some pure minerals that could exist in positive amounts in our equilibrium calculations:

editor.addMineralPhase("Halite")
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Magnesite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Quartz")

# Next follows an important step, where we create the chemical system with the information so far collected in the [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object `editor`:

system = ChemicalSystem(editor)

# The [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) class is one of the most
# important classes in Reaktoro. It is the class used to computationally represent a chemical system, with all its
# phases, species, and elements. It is also the class used for computation of thermodynamic properties of phases and
# species, such as activities, chemical potentials, standard Gibbs energies, enthalpies, phase molar volumes, densities,
# and many others. Many classes in Reaktoro require an instance of
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) for their initialization,
# since any chemical calculation needs to know the definition of the chemical system and the thermodynamic models
# describing the non-ideal behavior of the phases.

problem = EquilibriumProblem(system)

# Reaktoro provides the class [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html)
# for convenient description of equilibrium conditions. Using this class allows one to set the temperature and
# pressure at equilibrium, and a recipe that describes a mixture of substances and their amounts, which can be seen
# as initial conditions for the equilibrium calculation.

problem.setTemperature(70, "celsius")
problem.setPressure(100, "bar")
problem.add("H2O", 1.0, "kg")
problem.add("CO2", 2.0, "mol")
problem.add("NaCl", 1.0, "mol")
problem.add("CaCO3", 10.0, "g")
problem.add("MgCO3", 5.0, "g")
problem.add("Quartz", 1.0, "mol")

# > **Note:** The substance names above can either be chemical formulas, such as CaCO3 and CaCl2, as well as names of species that can be found in the database, such as Quartz. Reaktoro will break down the chemical formulas of the substances and calculate the amount of each chemical element in the system. These element amounts are inputs to the equilibrium calculation. In the future, we will only allow species names to be provided, since this is a safer way of preventing unfeasible elemental mass conditions to be imposed (e.g., there are *x* moles of C and *y* moles of O, and distributing these among the species always produce an excess of either C or O).

# We now perform a fast Gibbs energy minimization calculation to compute the chemical equilibrium state of the system at given conditions stored in `problem`. For this, we use the convenient function `equilibrate`:

state = equilibrate(problem)

# ## Using class ChemicalState to inspect species amounts

# The result of the `equilibrate` call before, `state`, is an object of class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html). This object contains the temperature, pressure, and amounts of the species at the computed chemical equilibrium state. We can access these properties as follows:

T = state.temperature()
P = state.pressure()
n = state.speciesAmounts()

print(f"T = {T} K")
print(f"P = {P} Pa")
print(f"n = (in mol)\n{n}")

# To print the name of each species and its amount (in mol), we execute the following loop:

for species in system.species():
    name = species.name()
    amount = state.speciesAmount(name)
    print(f"{name:>15} = {amount}")

# You can also output the chemical state to a file

state.output("results/state.txt")

# ## Using the class ChemicalProperties to obtain species activities

# If you require chemical properties of a system that depend on temperature (*T*), pressure (*P*), and composition (*n*), then [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) class is what you need

properties = ChemicalProperties(system)

# We can compute the chemical properties of the system at the state of equilibrium we found before:

properties.update(T, P, n)

# Alternatively, we could also have done:

properties = state.properties()

# > **Note:** The call above creates a new object of [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) each time. If you are using Reaktoro in a simulator that needs the chemical properties of the system at millions/billions of states each time step, prefer to the `update` method of an existing [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) object.

# Once we have computed the chemical properties, we can query for some of them. Below we get the natural log of the species activities:

lna = properties.lnActivities().val

# > **Note:** The use of `.ddT`, `.ddP`, and `.ddn`, instead of `.val`, extracts the derivatives of the activities (or any other chemical property) with respect to *T*, *P*, and *n*, respectively.

# To compute the actual activities (not their natural log), and print them one by one, we do

a = numpy.exp(lna)
for i, species in enumerate(system.species()):
    print(f"{species.name():>15} = {a[i]}")

# ## Calculating the pH of the aqueous solution

# Let's create a pH function that computes the pH of the aqueous solution given the chemical properties of the system (this will be soon simplified!)

evaluate_pH = ChemicalProperty.pH(system)
pH = evaluate_pH(properties)
print(f"The pH of the aqueous phase is {pH.val}.")
print(f"Its sensitivity with respect to speciation, ∂(pH)/∂n, is:")
for i, species in enumerate(system.species()):
    print(f"{species.name():>15} = {pH.ddn[i]}")