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

# # Opening bottle with soda experiment

from reaktoro import *
db = Database("supcrt98.xml")
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O C")
editor.addGaseousPhase(["CO2(g)"])

system = ChemicalSystem(editor)
print(system)

p_closed_bottle = 3.79 # in bars
p_open_bottle = 1.01325 # in bars
pressures = numpy.linspace(p_open_bottle, p_closed_bottle, num=100)
print(pressures)

problems = []
for P in pressures:
    problem = EquilibriumProblem(system)
    problem.setTemperature(20.0, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 0.5, "kg")
    # A typical carbonated soft drink contains approximately 3–4 volumes (6–8 g/L) CO2.
    # 8 g/L = 8 / 44.01 =  0.18, where 44.01 g/mol is the CO2 molar mass
    problem.add("C02", 0.18, "mol")
    problems.append(problem)  # append the new problem into problems

solver = EquilibriumSolver(system)
states = [ChemicalState(system) for _ in range(len(problems))]
for i in range(len(problems)):
    solver.solve(states[i], problems[i])

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

def custom_figure(title, y_axis_label): return figure(plot_width=600, plot_height=300,
                                                      title=title,
                                                      x_axis_label='Pressures',
                                                      y_axis_label=y_axis_label)

co2gas_amount = [state.speciesAmount("CO2(g)") for state in states]
fig2 = custom_figure(title="CO2(g)", y_axis_label='Amount of CO2(g) [mol]')
fig2.line(pressures, co2gas_amount, line_width=4, line_color='teal')
show(fig2)
