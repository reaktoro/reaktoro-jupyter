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
#       jupytext_version: 1.3.2
# ---

# # Carbon dioxide gas solubility in water
#
# This tutorial demonstrates how to simulate the solubility of CO2 gas in water or, in simpler words, the effect
# happening when one opens the bottle with soda.
#
# First, we defined considered chemical system:

from reaktoro import *
db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O C")
editor.addGaseousPhase(["CO2(g)"])
system = ChemicalSystem(editor)
print(system)

# Then we defined the pressure in the bottle before adn after opening (as two end states). To generates pressure values
# in between, we use function `linspace()` of the **numpy** library:

p_closed_bottle = 3.79 # in bars
p_open_bottle = 1.01325 # in bars
pressures = numpy.linspace(p_open_bottle, p_closed_bottle, num=100)
print(pressures)

# Next, we run through the generated pressure list and create a chemical problem corresponding to each pressure.
#
# > **Note**: A typical carbonated soft drink contains approximately 3–4 volumes (6–8 g/L) CO2.
# > To obtain amount of mol of CO2 we need to add to each problem, we do the following calculations:
# > 8 g/L = 8 / 44.01 =  0.18$, where 44.01 g/mol is the CO2 molar mass.

problems = []
for P in pressures:
    problem = EquilibriumProblem(system)
    problem.setTemperature(20.0, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 0.5, "kg")   # add ~ half a liter of water
    problem.add("C02", 0.18, "mol") # add calculated amount of gas
    problems.append(problem)        # append the new problem into the list of problems problems

# Equilibrate the list of generated chemical problems:

solver = EquilibriumSolver(system)
states = [ChemicalState(system) for _ in range(len(problems))]
for i in range(len(problems)):
    solver.solve(states[i], problems[i])

# To visualize the changes in the CO2 amount in the bottle, we export
# [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples) python plotting package.

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

def custom_figure(title, y_axis_label): return figure(plot_width=600, plot_height=300,
                                                      title=title,
                                                      x_axis_label='Pressures',
                                                      y_axis_label=y_axis_label)

co2gas_amount = [state.speciesAmount("CO2(g)") for state in states]
fig = custom_figure(title="CO2(g)", y_axis_label='Amount of CO2(g) [mol]')
fig.line(pressures, co2gas_amount, line_width=4, line_color='teal')
show(fig)

# From the plot, we see that by decreasing the pressure in the bottle, we also reduce the amount of CO2(g).
