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

# # Performing a chemical equilibrium calculation using equilibrium solver
#
# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation with the help of
# the [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html) class. First, we import
# everything from the `reaktoro` package by

from reaktoro import *

# ### Initializing chemical editor

# We start from creating an object of [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html) class to
# use in the initialization of the chemical system.

database = Database("supcrt98.xml")

# Next, we define which phases and species the chemical system should have. This is done using an instance of
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class.

editor = ChemicalEditor(database)
editor.addAqueousPhaseWithElementsOf("H2O NaCl CO2")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Halite")

# Here, [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html) is created by specifying the
# list of compound or substance names, i.e., H<sub>2</sub>O, NaCl, and CO<sub>2</sub>, that might not
# necessarily represent names of species in the database.
# Function [addAqueousPhaseWithElementsOf](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#a23e44f994b87c650a949226ddc195710) will brake the list
# of compounds into a list of element names, and the database will be similarly searched for all species that could
# be formed out of those elements.

# Then, [GaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1GaseousPhase.html) is composed of the names of
# the provided gaseous species H<sub>2</sub>O(g) and CO<sub>2</sub>(g), using function
# [addGaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#a8ec1e3a057d794df0dd6988c84bb5d3d).
# These names must conform to those used in the database that was specified  during the initialization of the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object, otherwise, an exception
# will be thrown.

# The [MineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1MineralPhase.html) object is created by specifying the
# names of the mineral species one by one. Analogously to the gaseous species, provided names must
# coincide with those used in the database (specified during the initialization of
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object) ,
# otherwise, an exception will be thrown. In this case, method
# [addMineralPhase](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#a05b263aa9d797a105feb9b83e05e1b86)
# is used to create single pure mineral phases with halite.

# ### Chemical system definition

# To initialize the chemical system, we use class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html), which requires the
# instance of [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) defined earlier.

system = ChemicalSystem(editor)

# ### Equilibrium problem definition

# The equilibrium problem is described by the class
# [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html). Here, different properties,
# such as temperature, pressure, and amounts of compounds, can be provided.

problem = EquilibriumProblem(system)
problem.setTemperature(60, "celsius")
problem.setPressure(300, "bar")
problem.add("H2O", 1, "kg")
problem.add("CO2", 100, "g")
problem.add("NaCl", 0.1, "mol")

# In particular, we set temperature to 60 &deg;C and pressure to 300 bar. To equilibrate the chemical problem,
# we also add 1kg of water, 100 g of carbon dioxide, and 0.1 mol of sodium chloride.

# ### Equilibrium solver definition

# The temperature, pressure, and mole amounts of the elements can be obtained from the definition of equilibrium
# problem.

T = problem.temperature()
P = problem.pressure()
b = problem.elementAmounts()

# Next, we create an object of [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)
# class that can be reused many times for equilibrium simulations.

solver = EquilibriumSolver(system)

# Similarly, an object of [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) class,
# which stores the equilibrium state of the system, can be created.

state = ChemicalState(system)

# ## Equilibration at 60 &deg;C

# Using method [solve](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html#ab01c678651bacb079f8f436c9a3a5148),
# the equilibrium state with given *(T, P, b)* inputs is generated and stored in the object `state`.

solver.solve(state, T, P, b)

# Here, the object `state` serves as the initial guess and the final state of the equilibrium calculation. If known,
# the temperature must be provided in Kelvin, and the pressure is expected in Pascal. Vector `b` provides the molar
# amounts of the elements in the equilibrium partition. Alternatively, one can call this method with given
# equilibrium problem.

# To save the calculated chemical equilibrium state into the text-file, we use method
# [output](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#ae5f2706f5be6e6856360a2f1073931e2).

state.output('state-T-60.txt')

# In the saved file, one can note that the amount of halite is of the order $10^{-21}$, which indicates its
# dissolution in sodium chloride brine.

# ## Equilibration at 70 &deg;C

# Calculate the new equilibrium state at temperature increased by 10 &deg;C. For that, we use previous
# equilibrium state as the initial guess for improved performance.

solver.solve(state, T + 10.0, P, b)

# Save newly calculated chemical equilibrium state with T = 70 &deg;C into the text-file.

state.output('state-T-70.txt')

# In comparison to chemical speciation obtained for equilibrium calculation at 60 &deg; C,
# the ionic strength, reduction potential, and alkalinity slightly decrease, whereas ph slightly increases.
