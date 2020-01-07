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

# # Performing a chemical equilibrium calculation with fixed properties

# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation when some of the
# properties are held fixed. We start by importing the `reaktoro` package:

from reaktoro import *

# ## Initializing chemical editor
#
# Class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html)
# provides convenient operations to initialize
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) instances.
# Below, we define the editor of the chemical system from the default database SUPCRT92:

editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")

# Here, [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html) is created by specifying the
# list of compound or substance names that might not necessarily represent names of species in the database. The list
# of compounds will be broken into a list of element names, and the database will be similarly searched for all species
# that could be formed out of those elements.

# Then, [GaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1GaseousPhase.html) is composed from the names of
# the provided gaseous species. These names must conform to those used in the database that was specified  during the
# initialization of the `ChemicalEditor` object, otherwise, an exception will be thrown.

# The [MineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1MineralPhase.html) object is created by specifying the
# names of the mineral species one by one. Analogously to the gaseous species, provided names must
# coincide with those used in the database (specified during the initialization of
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object) ,
# otherwise, an exception will be thrown. In this case, method
# [addMineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#a05b263aa9d797a105feb9b83e05e1b86)
# is used to create one pure mineral phases with calcite.

# ### Chemical system definition
#
# Construction of chemical system is done by calling

system = ChemicalSystem(editor)

# ### Inverse equilibrium problems
#
# Generally, the inverse equilibrium problem is handled by the class [EquilibriumInverseProblem](
# https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html).
# In an inverse equilibrium problem, not all elements have known molar amounts. Their amount constraints are
# replaced by other equilibrium constraints such as fixed species amount or activity, or the volume or total amount
# of a phase. Since the amounts of elements are not known a priori, an inverse equilibrium calculation tries to
# determine amounts of titrants that can control the specified equilibrium constraints. The amount of the titrants
# are unknown, and its addition or removal is done over the calculation so that the equilibrium state is driven
# towards a state where all given equilibrium constraints are satisfied.

# First problem, which we consider, is the problem with fixed mass of mineral (in this case, calcite) and
# molar amount of species $\mathrm{CO_2(g)}$ in equilibrium.
# Besides, we initialized the amount of $\mathrm{H_2O}$ by 1 kg and sodium-chloride $\mathrm{NaCl}$
# by 0.1 mol.
#
# Using function [equilibrate()](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e),
# we calculate the chemical equilibrium state of the system with the given equilibrium conditions stored in the
# object `problem1`.

problem1 = EquilibriumInverseProblem(system)
problem1.add("H2O", 1, "kg")
problem1.add("NaCl", 0.1, "mol")
problem1.fixSpeciesMass("Calcite", 100, "g")
problem1.fixSpeciesAmount("CO2(g)", 1.0, "mol")

state1 = equilibrate(problem1)

# To output the properties of the [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html)
# obtained by the equilibration, we use method `output()`:

state1.output('state1.txt')

# Indeed, in the file with the resulting information about `state1`, in the column *Species* we see 1 mol for the
# amount of $\mathrm{CO_2(g)}$, whereas in the column *Phases* for calcite the mass is set to 0.1 kg. Alkalinity
# in this case becomes equal to 0.0214083 $\mathrm{[eq/L]}$.

# In the second equilibrium inverse problem, which we consider similar conditions used for the first problem. However,
# in addition, we fix the total alkalinity of the aqueous solution by the method
# [alkalinity](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html#a870a693edb5134d8c8d7aff325e98cde),
# where we provide the actual value of the total alkalinity of the aqueous solution,
# the units of the total alkalinity (must be convertible to eq/L), and the name of titrant that control the solution
# alkalinity. Again, the problem is equilibrated by the standard method `equilibrate` and the obtained chemical state
# is then output for the comparison with the earlier obtained chemical state.

problem2 = EquilibriumInverseProblem(system)
problem2.add("H2O", 1, "kg")
problem2.add("NaCl", 0.1, "mol")
problem2.fixSpeciesMass("Calcite", 100, "g")
problem2.fixSpeciesAmount("CO2(g)", 1.0, "mol")
problem2.alkalinity(25.0, "meq/L", "Cl")

state2 = equilibrate(problem2)

state2.output('state2.txt')

# Due to the fixed alkalinity, we see 0.0250008 $\mathrm{[eq/L]}$ in the resulting file. In this case, ph is slightly
# higher, i.e., pH = 6.12591. Moreover, unlike the earlier simulations, here we obtain negative reduction potential,
# i. e., pE = -6.51608.

# Third inverse problem is initialized analogously with 1 kg of $\mathrm{H_2O}$ and 0.1 mol of sodium-chloride
# $\mathrm{NaCl}$. In the case, however, only the amount of calcite is fixed to 1 mol. In addition, we fix
# pH of the aqueous solution with two given titrants $\mathrm{HCl}$ and $\mathrm{NaOH}$.

problem3 = EquilibriumInverseProblem(system)
problem3.add("H2O", 1, "kg")
problem3.add("NaCl", 0.1, "mol")
problem3.pH(8.0, "HCl", "NaOH")
problem3.fixSpeciesAmount("Calcite", 1, "mol")

state3 = equilibrate(problem3)

state3.output('state3.txt')

# According to the above instructions, the pH is fixed to 8 in `state3.txt`.
# The obtained ionic strength is less then in the previous case, i.e.,
# 0.103845 $\mathrm{[molal]}$.

# As it was mention in the introduction of the
# [EquilibriumInverseProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html) class,
# it is also possible to fix the volume of a phase at equilibrium in the inverse equilibrium problem. This is done in
# similar to earlier considered problem for all the phase valumes, i.e., gaseous, aqueous, and calcite, provided
# the name of titrans. For the gaseous phase, we use $\mathrm{CO_2}$ for titration, for the aqueous one, 1 kg of water
# and 0.1 mol of sodium-chloride, and $\mathrm{CaCO_3}$ for the mineral.

problem4 = EquilibriumInverseProblem(system)
problem4.add("H2O", 1, "kg")
problem4.add("NaCl", 0.1, "mol")
problem4.fixPhaseVolume("Gaseous", 0.2, "m3", "CO2")
problem4.fixPhaseVolume("Aqueous", 0.3, "m3", "1 kg H2O; 0.1 mol NaCl")
problem4.fixPhaseVolume("Calcite", 0.5, "m3", "CaCO3")

state4 = equilibrate(problem4)

state4.output('state4.txt')

# Again, in the column *Phases* the volumes of *Aqueous*, *Gaseous*, *Calcite* is fixed to 0.3, 0.2, and 0.5,
# respectively. As the result, we obtain considerably big amount of an aqueous phase 300.764 kg and calcite 1354.94 kg.
