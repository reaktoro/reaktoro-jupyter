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

# # Inverse chemical equilibrium calculations
#
# This tutorial demonstrates how to use Reaktoro to perform *inverse chemical
# equilibrium calculations*.
#
# In an inverse chemical equilibrium problem, not all elements have known
# amounts. Instead of fully specifying their amounts, we impose one or more
# equilibrium constraints such as fixed species amount, fixed species activity,
# fixed volume of a phase, fixed pH of an aqueous solution; the possibilities
# are many.
#
# Since the amounts of elements are not known a priori, an inverse equilibrium
# calculation tries to determine amounts of titrants that can control the
# specified equilibrium constraints. The amounts of the titrants are unknown,
# and its addition or removal is done over the calculation so that the
# equilibrium state is driven towards a state where all given equilibrium
# constraints are satisfied.
#
# Below we create a chemical system containing an aqueous phase, a gaseous
# phase, and a mineral phase.

# +
from reaktoro import *

db = Database("supcrt98.xml")

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C");
editor.addGaseousPhaseWithElements("H O C");
editor.addMineralPhase("Calcite")

system = ChemicalSystem(editor)
# -

# > Check previous tutorials to learn the  steps above! For example,
# > [Basics of equilibrium calculation](eq.equilibrium-basics.ipynb),
# > [Equilibrium calculation of carbonate species](eq.equilibrium-carbonates.ipynb),
# > [Equilibrium calculations using equilibrium solver](eq.co2-brine-using-equilibrium-solver.ipynb), or
# > [Custom activity model for equilibrium calculations](eq.custom-activity-models.ipynb).

# #### Equilibrium problem with fixed pH, species amount, and species activity
#
# In this problem, we fix the pH of the aqueous solution to 3, and specify that
# HCl is the titrant that should be added or removed so that this constraint is
# satisfied. Therefore, the amounts of H and Cl are not known in advance.
#
# To obtain an aqueous solution saturated with CO<sub>2</sub>, we also add a
# constraint that enforces 1 mol for the amount of gaseous species
# CO<sub>2</sub>(g) at equilibrium.
#
# Finally, we add an activity constraint for gaseous species O<sub>2</sub>(g).
# This is equivalent to enforcing the fugacity of the gas.

# +
problem1 = EquilibriumInverseProblem(system)
problem1.add("H2O", 1, "kg")
problem1.add("NaCl", 0.1, "mol")
problem1.add("CaCl2", 2, "mmol")
problem1.add("MgCl2", 4, "mmol")
problem1.pH(3.0, "HCl")
problem1.fixSpeciesAmount("CO2(g)", 1.0, "mol")
problem1.fixSpeciesActivity("O2(g)", 0.20)

state1 = equilibrate(problem1)

state1.output("state1.txt")
# -

# > Check [EquilibriumInverseProblem](
# https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html) to learn more about all possible
# equilibrium constraints currently supported.
#
# > Check the output file `state1.txt` and check whether the constraints were
# > successfully satisfied. We plan to output these files in HTML format so that it can be more conveniently inspected.


# #### Equilibrium problem with fixed pH controlled by CO<sub>2</sub>
#
# The second equilibrium inverse problem has similar conditions to those used
# for the first problem. However, in this case, we fix pH to be equal to 4.0,
# with CO<sub>2</sub> being the titrant whose amount to be added into the system
# is unknown in advance.

# +
problem2 = EquilibriumInverseProblem(system)
problem2.add("H2O", 1, "kg")
problem2.add("NaCl", 0.1, "mol")
problem2.add("CaCl2", 2, "mmol")
problem2.add("MgCl2", 4, "mmol")
problem2.pH(4.0, "CO2")

state2 = equilibrate(problem2)

state2.output("state2.txt")
# -

# > Check the output file `state2.txt` and make sure whether the constraints were successfully satisfied.


# #### Equilibrium problem with fixed pH controlled by either HCl or NaOH
#
# In this equilibrium problem, we also fix the pH of the aqueous solution, but we
# specify that this constraint is to be satisfied by titrating either HCl or
# NaOH (not both!).

# +
problem3 = EquilibriumInverseProblem(system)
problem3.add("H2O", 1.0, "kg")
problem3.add("NaCl", 0.1, "mol")
problem3.add("Calcite", 1.0, "mol")
problem3.pH(8.0, "HCl", "NaOH")

state3 = equilibrate(problem3)

state3.output("state3.txt")
# -

# > Check the output file `state3.txt` and control whether the constraints were
# > successfully satisfied.
