# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ../notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
# ---

# # Performing a chemical equilibrium calculation with full evaporation of the water

# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation when the whole mineral
# phase is evaporating as the result of reacting with brine. We start by importing the `reaktoro` package:

from reaktoro import *

# Define the phases and corresponding to these phases species, which the chemical system should have.
# It is done using [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object.

editor = ChemicalEditor()
editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--", "Ca++"]) \
    .setChemicalModelDebyeHuckel() \
    .setActivityModelDrummondCO2()
editor.addGaseousPhase(["CO2(g)", "H2O(g)"]) \
    .setChemicalModelSpycherPruessEnnis()
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Halite")

# Here, [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html) is created by specifying the
# list of exact names of the species. These names must coincide with those used in the database that was specified
# during the initialization of the [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html)
# object, otherwise, an exception will be thrown. In this case, the default database SUPCRT92 is used.

# To initialize the chemical model of the [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html)
# with the Debye-Huckel equation of state, we use method [setChemicalModelDebyeHuckel](
# https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html#aa3f53d5cb5ae7adfb50e563c7a198ce6). However,
# we specify that the Drummond (1981) activity model must be used to model $\mathrm{CO2(aq)} using
# [setActivityModelDrummondCO2](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html
# #a8d98d8294d81b26043e3a8d43e386c21).

# Then, [GaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1GaseousPhase.html) is composed from the names of
# the provided gaseous species $\mathrm{H_2O(g)}$ and $\mathrm{CO_2(g)}$. These names must conform to those
# used in the database that was specified  during the initialization of the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object, otherwise, an exception
# will be thrown.
# Here, method `setChemicalModelSpycherPruessEnnis()` sets Spycher et al. (2003) equation of state. This model only
# supports the gaseous species $\mathrm{H_2O(g)}$ and $\mathrm{CO_2(g)}$.

# The [MineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1MineralPhase.html) object is created by specifying the
# names of the mineral species one by one. Analogously to the gaseous species, provided names must
# coincide with those used in the database (specified during the initialization of
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object) ,
# otherwise, an exception will be thrown. In this case, method
# [addMineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#a05b263aa9d797a105feb9b83e05e1b86)
# is used to create two pure mineral phases with calcite and halite.

# ### Chemical and reaction system definition

# To initialize the chemical system, we use class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html), which requires the
# instance of [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) defined earlier.

system = ChemicalSystem(editor)

# Class [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) serves as a system of the
# chemical reaction by a collection of [Reaction](https://reaktoro.org/cpp/classReaktoro_1_1Reaction.html) class
# instances. It provides convenient methods that calculate the equilibrium constants, reaction quotients,
# and rates of the reactions.

reactions = ReactionSystem(editor)

# ### Chemical problem definition
#
# We define an equilibrium problem providing amounts of compounds. In particular, we mix 1 kg of water with 1 mol of
# sodium-chloride and 200 kg of carbon-dioxide. The amount of calcite in the system is set to 10 mol.

problem = EquilibriumProblem(system)
problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CaCO3", 10, "mol")
problem.add("CO2", 200, "kg")

# To customize options for the equilibrium calculations, class
# [EquilibriumOptions](https://reaktoro.org/cpp/structReaktoro_1_1EquilibriumOptions.html) can be used.
# For instance, in [NonlinearOptions](https://reaktoro.org/cpp/structReaktoro_1_1NonlinearOptions.html), it contains
# information about the nonlinear solver used.

options = EquilibriumOptions()
options.optimum.output.active = True
options.epsilon = 1e-20

# Here, we set the field of [OutputterOptions](https://reaktoro.org/cpp/structReaktoro_1_1OutputterOptions.html)
# class to be `True` in order to determine whether the intermediate values of equilibrium simulations must be output
# to the console.
# Then, we set the parameter $\varepsilon$ (used for the numerical representation of a zero molar amount) to be equal
# 1e-50. The molar amount of the `i`-th species is considered zero if $n[i] < \varepsilon \cdot \min b$, where `b` is
# the vector of element molar amounts.

# Finally, we use function [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#a908245bfa7d236d8c556241dc87f489e)
# proving not only the instance of equilibrium problem but also the specified `options`.

state = equilibrate(problem, options)

# **Note:** Function `equilibrate` is intended for convenience only. For performance critical applications, use class
# [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html).

# Run of the `equilibrium` method generates the following error *Error: Could not calculate the equilibrium state of the system.*
# The reason of this error follows after the the error:
# *Reason: Convergence could not be established with given equilibrium conditions, initial guess, and(or) numerical
# parameters*. If we review the amount of species, which were output to the console, we will see that `n[H2O(l)] =
# 0`. This means that all the water used in the mixing got evaporated, therefore, equilibration cannot be proceeded.
