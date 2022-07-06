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
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Chemical equilibrium calculations

# This tutorial demonstrates how to use Reaktoro to perform chemical equilibrium calculations. We start by importing the `reaktoro` package:

from reaktoro import *

# ## Initializing thermodynamic database

# Next, we need a thermodynamic database that enables us to compute the thermodynamic properties of species and
# reactions. For this, we create an object of class [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html):

db = Database("supcrt98.xml")

# > For more detailed overview on the functionality of this class, please check the tutorial
# [**Database**](cl.database.ipynb).

# ## Initializing chemical system

# To indicate the phases of interest (as well as their species) that may potentially exist at equilibrium,
# we create an object of class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html):

editor = ChemicalEditor(db)

# We consider an aqueous phase composed of all aqueous species in the database that can be formed with the given
# chemical elements below:

editor.addAqueousPhaseWithElements("H O Na Cl C Ca Mg Si")

# > **Note:** This automatic selection of chemical species for a phase can result in a large number of them. This
# potentially increases the computing cost of the chemical reaction calculations. If you are using Reaktoro in
# demanding applications, you might want to manually specify the chemical species of each phase in your chemical
# system. This can be achieved by providing an explicit list of species names, e.g.,
# `editor.addAqueousPhase("H2O(l) H+ OH- CO2(aq)")`. Note, however, that care is required here to ensure relevant
# species are not missing. The just
# given example is a bad one in fact, with important species such as `HCO3-` and `CO3--` missing in the list.

# We are interested in a gaseous phase containing exactly the following gases (which may not exist in positive
# amounts at the end of our equilibrium calculation later):

editor.addGaseousPhase("H2O(g) CO2(g)")

# Finally, we consider some pure minerals that could exist in positive amounts in our equilibrium calculations:

editor.addMineralPhase("Halite")
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Magnesite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Quartz")

# > See the tutorial [**ChemicalEditor**](cl.chemical-editor.ipynb) for studying further capabilities of
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class.

# ## Initializing chemical system

# Next, follows an important step with the creation of the chemical system with the information so far collected in the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object `editor`:

system = ChemicalSystem(editor)

# The [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) class is one of the most
# important classes in Reaktoro. It is the class used to computationally represent a chemical system, with all its
# phases, species, and elements. It is also the class used for computation of thermodynamic properties of phases and
# species, such as activities, chemical potentials, standard Gibbs energies, enthalpies, phase molar volumes, densities,
# and many others. Many classes in Reaktoro require an instance of
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) for their initialization,
# since any chemical calculation needs to know the definition of the chemical system and the thermodynamic models
# describing the non-ideal behavior of the phases.

# > See [**ChemicalSystem**](cl.chemical-system.ipynb) for the explanation on functionality of the class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html).

# ## Initializing equilibrium problem

# We use class [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html) to specify the
# conditions at which our system should be in equilibrium.

problem = EquilibriumProblem(system)

# Reaktoro provides the class [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html)
# for the convenient description of equilibrium conditions. Using this class allows one to set the temperature and
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

# > **Note:** The substance names above can either be chemical formulas, such as `CaCO3` and `CaCl2`, as well as names of
# species that can be found in the database, such as `Quartz`. Reaktoro will break down the chemical formulas of the
# substances and calculate the amount of each chemical element in the system. These element amounts are inputs to the
# equilibrium calculation. In the future, we will only allow species names to be provided since this is a safer way
# of preventing unfeasible elemental mass conditions to be imposed (e.g., there are *x* moles of C and *y* moles of
# O, and distributing these among the species always produce an excess of either C or O).

# ## Equilibration of chemical problem

# We now perform a fast Gibbs energy minimization calculation to compute the chemical equilibrium state of the system
# at given conditions stored in `problem`. For this, we use the convenient function
# [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e):

state = equilibrate(problem)

# ## Analyzing species amounts

# The result of the `equilibrate` call before, `state`, is an object of class [ChemicalState](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html). This object contains the temperature, pressure,
# and amounts of the species at the computed chemical equilibrium state. We can access these properties as follows:

T = state.temperature()
P = state.pressure()
n = state.speciesAmounts()

print(f"T = {T} K")
print(f"P = {P} Pa")
print(f"n (in mol) = \n{n}")

# To print the name of each species and its amount (in mol), we execute the following loop:

print("Species names : n (in mol)")
for species in system.species():
    name = species.name()
    amount = state.speciesAmount(name)
    print(f"{name:>13} = {amount}")

# To examine the complete information about the chemical state, you can also output it to a file

state.output("state.txt")

# ## Analyzing chemical properties

# If you require chemical properties of a system that depend on temperature (*T*), pressure (*P*), and composition (*n*),
# then [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) class is what you need

properties = ChemicalProperties(system)

# We can compute the chemical properties of the system at the state of equilibrium we found before:

properties.update(T, P, n)

# Alternatively, one can also do:

properties = state.properties()

# > **Note:** The call above creates a new object of [
# > ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html) each time. If you are using
# > Reaktoro in a simulator that needs the chemical properties of the system at millions/billions of states each time
# > step, instead of populating many instances of ChemicalProperties, prefer to use
# > [ChemicalProperties::update](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html#af923d85484865039fa56889c1a2f36c9)
# > method of an existing [ChemicalProperties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalProperties.html)
# > object.

# ### Analyzing activities

# Once we have computed the chemical properties, we can query for some of them. Below we get the natural log of
# species activities:

lna = properties.lnActivities().val

# > **Note:** The use of `.ddT`, `.ddP`, and `.ddn`, instead of `.val`, extracts the derivatives of the activities
# (or any other chemical property) with respect to *T*, *P*, and *n*, respectively.

# To compute the actual activities (not their natural log), and print them one by one, we do

a = numpy.exp(lna)
print("Species names : activities")
for i, species in enumerate(system.species()):
    print(f"{species.name():>13} : {a[i]}")

# ### Analyzing chemical potentials

# Similarly, we can inspect chemical potentials:

mu = properties.chemicalPotentials().val
print("Species names : potentials (in kJ/mol)")
for i, species in enumerate(system.species()):
    print(f"{species.name():>13} : {mu[i]}")

# ### Calculating the pH of the aqueous solution

# Let's create a pH function that computes the pH of the aqueous solution given the chemical properties of the system.

evaluate_pH = ChemicalProperty.pH(system)
pH = evaluate_pH(properties)

# > **Note:** This will be soon simplified!

# Besides its value, which can be obtained by  `pH.val`, `ph` also contains derivatives with respect to the species
# amounts. It can be accessed by `pH.ddn`:

print(f"The pH of the aqueous phase is {pH.val}.")
print(f"Its sensitivity with respect to speciation is:")
print("Species names : ∂(pH)/∂n (in 1/mol)")
for i, species in enumerate(system.species()):
    print(f"{species.name():>13} : {pH.ddn[i]}")

