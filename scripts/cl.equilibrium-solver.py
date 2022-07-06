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

# # Using class EquilibriumSolver
#
# This tutorial demonstrates the use of class
# [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)
# for repeated equilibrium calculations.
#
# For this, we model the dissolution of mineral calcite (CaCO<sub>3</sub>) in a
# saline aqueous solution as we increase the amount of CO<sub>2</sub> into the system.
#
# > **Note:** Prefer class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)
# instead of method [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e)
# when many chemical equilibrium calculations are needed.
#
# We follow the usual initialization steps, such as importing the `reaktoro`
# package,

from reaktoro import *
import numpy

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
# instance of the class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html),
# which will be constructed containing four phases: an aqueous phase, a gaseous
# phase, and two pure mineral phases.

system = ChemicalSystem(editor)

# We want to perform a sequence of equilibrium calculations. For this, we create
# first a numpy array containing evenly spaced numbers ranging from 0 to 4 (a
# total of 101 values). This array represents the amount of CO<sub>2</sub> in
# each equilibrium problem we shall solve later.

co2_amounts = numpy.linspace(0.0, 4.0, num=101)

# We now create the list of equilibrium problems stored in the list `problems`, one for each `amount` in array
# `x`. In each problem, we have 1 mol of calcite, 1 kg of H<sub>2</sub>O, 1 mol
# of NaCl, and `amount` of CO<sub>2</sub>:

# +
problems = []

for amount in co2_amounts:
    problem = EquilibriumProblem(system)
    problem.setTemperature(70.0, "celsius")
    problem.setPressure(100.0, "bar")
    problem.add("Calcite", 1.0, "mol")
    problem.add("H2O", 1.0, "kg")
    problem.add("NaCl", 1.0, "mol")
    problem.add("CO2", amount, "mol")

    problems.append(problem)  # append the new problem into problems
# -

# Solving this sequence of equilibrium problems can be done with a single
# instance of the class
# [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html):

solver = EquilibriumSolver(system)

# Each chemical equilibrium calculation results in an individual chemical
# state. Below, we create a list of
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html)
# objects that will contain the solution to each problem.

states = [ChemicalState(system) for _ in range(len(problems))]

# We now solve the sequence of chemical equilibrium problems:

for i in range(len(problems)):
    solver.solve(states[i], problems[i])

# We want to plot the amount of calcite versus the amount of CO<sub>2</sub>
# added into the system. We collect below the amount of calcite in each computed
# equilibrium state in a list.

calcite_amounts = [state.speciesAmount("Calcite") for state in states]

# We now use [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples)
# to do the plotting.

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

# Besides, we define a custom function that would generate figure of a size 600 x 300 with a label
# `Amount of CO2 [mol]` on the x-axis:

def custom_figure(title, y_axis_label): return figure(plot_width=600, plot_height=300,
                                                      title=title,
                                                      x_axis_label='Amount of CO2 added [mol]',
                                                      y_axis_label=y_axis_label)

# Create a new plot with the given title and y-axis label, where the amount of calcite in `calcite_amounts` versus
# the amount of CO<sub>2</sub> in `co2_amounts`:

fig1 = custom_figure(title="Calcite dissolution", y_axis_label='Amount of Calcite [mol]')
fig1.line(co2_amounts, calcite_amounts, line_width=4)
show(fig1)

# To understand why calcite does not dissolve further after
# about 0.76 mol of CO<sub>2</sub> is added into the system,
# we plot the amount of gaseous species CO<sub>2</sub>(g):

co2gas_amount = [state.speciesAmount("CO2(g)") for state in states]
fig2 = custom_figure(title="CO2(g)", y_axis_label='Amount of CO2(g) [mol]')
fig2.line(co2_amounts, co2gas_amount, line_width=4)
show(fig2)

# Note above that after about 0.76 mol of CO<sub>2</sub> is added, the aqueous
# solution becomes CO<sub>2</sub>-saturated, and the amount of
# CO<sub>2</sub>(g), initially zero, starts to increase.
