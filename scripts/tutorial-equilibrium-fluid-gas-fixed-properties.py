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

# # Equilibrium system with gaseous and aqueous phases with fixed properties

# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation
# on an equilibrium system with gaseous and aqueous phases provided that some of the
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
editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
editor.addGaseousPhaseWithElements("H O C")

# Here, both [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html) and
# [GaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1GaseousPhase.html) are created by specifying the
# list of chemical element names. The database will be searched for all species that could be formed out of those
# elements.

# > **Note:** This automatic selection of chemical species for a phase can result in a large number of them. This
# potentially increases the computing cost of the chemical reaction calculations. If you are using Reaktoro in
# demanding applications, you might want to manually specify the chemical species of each phase in your chemical
# system. This can be achieved by providing an explicit list of species names, e.g., `editor.addAqueousPhase("H2O(l)
# H+ OH- CO2(aq)")`. Note, however, that care is required here to ensure relevant species are not missing. The just
# given example is a bad one in fact, with important species such as `HCO3-` and `CO3--` missing in the list.

# ## Chemical system definition
#
# Construction of the chemical system is done by calling

system = ChemicalSystem(editor)

# ## Inverse equilibrium problems
#
# Generally, the inverse equilibrium problem is handled by the class [EquilibriumInverseProblem](
# https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html).
# In an inverse equilibrium problem, not all elements have known molar amounts. Their amount constraints are
# replaced by other equilibrium constraints such as fixed species amount or activity, or the volume or total amount
# of a phase. Since the amounts of elements are not known a priori, an inverse equilibrium calculation tries to
# determine amounts of titrants that can control the specified equilibrium constraints. The amount of the titrants
# is unknown, and its addition or removal is done over the calculation so that the equilibrium state is driven
# towards a state where all given equilibrium constraints are satisfied.

# ### Problem with fixed ph, species' amount and activity

# First problem, which we consider, is the problem, which we initialize with 1 kg of $\mathrm{H_2O}$, 0.1 mol of
# sodium-chloride $\mathrm{NaCl}$, 2 mmol of calcium-chloride $\mathrm{CaCl_2}$, and
# 4 mmol of magnium-chloride $\mathrm{MgCl_2}$. From the properties of this inverse equilibrium problem, we fix
# pH to 3.0, providing the titran $\mathrm{HCl}$. Besides, $\mathrm{CO_2(g)}$ species' amount is fixed to 1.0 mol, as
# well as the species' activity of $\mathrm{O_2(g)}$ is prescribed to be 0.2.

problem1 = EquilibriumInverseProblem(system)
problem1.add("H2O", 1, "kg")
problem1.add("NaCl", 0.1, "mol")
problem1.add("CaCl2", 2, "mmol")
problem1.add("MgCl2", 4, "mmol")
problem1.pH(3.0, "HCl")
problem1.fixSpeciesAmount("CO2(g)", 1.0, "mol")
problem1.fixSpeciesActivity("O2(g)", 0.20)

# Using function [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#a8feea895affd774daaba5fbca8c4eeb3),
# we calculate the chemical equilibrium state of the system with given equilibrium conditions stored in the
# object `problem1`.

state1 = equilibrate(problem1)

# To evaluate the results of performed calculations, we save the chemical state to the txt-file `state1.txt`.

state1.output('state1.txt')

# According to the above instructions, the pH is fixed to 3.0 in `state1.txt`. We can also make sure that activity of
# $\mathrm{O_2(g)}$ (see sub-table with *Species*, column *Activity [-]*) is prescribed to 0.2, similarly to the
# species amount of $\mathrm{CO_2(g)}$ in the same table.
# The obtained ionic strength is 0.117995 molal, pE is 17.6023, the reduction potential is 1.04134 V, and
# the alkalinity is -0.0091649 eq/L.

# ### Problem with fixed ph

# The second equilibrium inverse problem has similar conditions to those used for the first problem. However,
# in this case, we fix pH to be equal to 4.0, providing $\mathrm{CO_2}$ as a titrant.

problem2 = EquilibriumInverseProblem(system)
problem2.add("H2O", 1, "kg")
problem2.add("NaCl", 0.1, "mol")
problem2.add("CaCl2", 2, "mmol")
problem2.add("MgCl2", 4, "mmol")
problem2.pH(4.0, "CO2")

state2 = equilibrate(problem2)

state2.output('state2.txt')

# Analogously to the `state1`, in `state2` the pH is fixed to 4.0 in `state1.txt`.
# The obtained ionic strength is slightly lower than in the previous case, i.e., 0.116759 molal,
# pE is also smaller, i.e., 13.6104, whereas the reduction potential is 0.805183 V.
# The final the alkalinity is -0.0078823 eq/L.