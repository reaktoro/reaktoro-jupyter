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

# # Using class EquilibriumSolver

# This tutorial demonstrates the use of class
# [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)
# for repeated equilibrium calculations.
#
# For this, we model the dissolution of mineral calcite (CaCO<sub>3</sub>) in a
# saline aqueous solution as we increase the amount of CO<sub>2</sub> into the system.

# > **Note:** Prefer class `EquilibriumSolver` instead of method `equilibrate`
# when many chemical equilibrium calculations are needed.

# We follow the usual initialization steps, such as importing the `reaktoro`
# package,

from reaktoro import *

# creating an object of class
# [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html),

db = Database("supcrt98.xml")

# and using this database object to initialize the object of class
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html)
# to specify the phases of interest in the chemical system:

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Na Cl C Ca")
editor.addGaseousPhase(["H2O(g)", "CO2(g)", "O2(g)"])
editor.addMineralPhase("Halite")
editor.addMineralPhase("Calcite")

# We now use the chemical system configuration stored in `editor` to create an
# instance of class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html),
# which will be constructed containing four phases: an aqueous phase, a gaseous
# phase, and two pure mineral phases.

system = ChemicalSystem(editor)

# We want to perform a sequence of equilibrium calculations. For this, we create
# first a numpy array containing evenly spaced numbers ranging from 0 to 4 (a
# total of 101 values). This array represents the amount of CO<sub>2</sub> in
# each equilibrium problem we shall solve later.

x = numpy.linspace(0.0, 4.0, num=101)

# We now create the list of equilibrium problems, one for each `amount` in array
# `x`. In each problem, we have 1 mol of calcite, 1 kg of H<sub>2</sub>O, 1 mol
# of NaCl, and `amount` of CO<sub>2</sub>:

problems = []

for amount in x:
    problem = EquilibriumProblem(system)
    problem.setTemperature(70.0, "celsius")
    problem.setPressure(100.0, "bar")
    problem.add("Calcite", 1.0, "mol")
    problem.add("H2O", 1.0, "kg")
    problem.add("NaCl", 1.0, "mol")
    problem.add("CO2", amount, "mol")

    problems.append(problem)  # append the new problem into problems

# Solving this sequence of equilibrium problems can be done with a single
# instance of class
# [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html):

solver = EquilibriumSolver(system)

# Each chemical equilibrium calculation results into an individual chemical
# state. Below we create a list of
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html)
# objects that will contain the solution of each problem.

states = [ChemicalState(system) for _ in range(len(problems))]

# We now solve the sequence of chemical equilibrium problems:

for i in range(len(problems)):
    solver.solve(states[i], problems[i])

# We want to plot the amount of calcite versus the amount of CO<sub>2</sub>
# added into the system. We collect below the amount of calcite in each computed
# equilibrium state.

y = [state.speciesAmount("Calcite") for state in states]

# We now use
# [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples)
# to do the plotting.

from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook
output_notebook()

# Create a new plot with given title and axis labels

fig1 = figure(title="Calcite dissolution",
    x_axis_label='Amount of CO2 [mol]',
    y_axis_label='Amount of Calcite [mol]')

# Plot the amount of calcite in `y` versus the amount of CO<sub>2</sub> in `x`:
fig1.line(x, y, line_width=4);

show(fig1)

# To understand why calcite does not dissolve further after
# about 0.76 mol of CO<sub>2</sub> is added into the system,
# we plot the amount of gaseous species CO<sub>2</sub>(g):

z = [state.speciesAmount("CO2(g)") for state in states]

fig2 = figure(title="CO2(g)",
    x_axis_label='Amount of CO2(g) [mol]',
    y_axis_label='Amount of CO2 added [mol]')

# Plot the amount of calcite in `y` versus the amount of CO<sub>2</sub> in `x`:
fig2.line(x, z, line_width=4);

show(fig2)

# Note above that after about 0.76 mol of CO<sub>2</sub> is added, the aqueous
# solution becomes CO<sub>2</sub>-saturated, and the amount of
# CO<sub>2</sub>(g), initially zero, starts to increase.
