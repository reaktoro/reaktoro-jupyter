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

# # Change of the carbon species in the reaction path vs. pH
#
# This tutorial demonstrates how to calculate a reaction path between two different chemical states in
# equilibrium, referred to as the *initial state* and *final state*.
# These states can have different temperatures, pressures, and/or molar amounts of elements. If we gradually adjust
# temperature, pressure, and elemental amounts in the system to bring the initial state to the final state, slowly
# enough so that **every intermediate state is in equilibrium**, the system would trace a co-called *reaction path*.
#
# Let the initial state have 0.5 mol of carbon-dioxide (CO<sub>2</sub>) and 1 mol of hydrogen chloride (HCl) mixed
# with 1 kg of water. We want to see how the addition of 2 mol of sodium hydroxide (NaOH) and removal of hydrogen
# chloride contributes to the chemical system. Thus, our initial and final states for a reaction path calculation can
# be described as follows:
#
# | Initial state    | Final state      |
# |------------------|------------------|
# | 1 kg of H2O      | 1 kg of H2O      |
# | 0.5 mol of CO2   | 0.5 mol of CO2   |
# | 1 mol of HCl     | 2 mol of NaOH    |
#
# As usual, we start by importing the `reaktoro` package:

from reaktoro import *

# Note that the object `editor` from class
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) was not initialized with a given
# [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html) object. Instead, it is initialized using the
# default built-in database file
# [supcrt98.xml](https://github.com/reaktoro/reaktoro/blob/master/databases/supcrt/supcrt98.xml).

editor = ChemicalEditor()

# For the aqueous phases, we list the chemical elements composing the phase instead of the species' exact names.
# Class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) searches for all species in
# the database those elements can form. Only species corresponding to the phase-type are selected
# (e.g., only aqueous species are searched in the current case).

# +
editor.addAqueousPhaseWithElements("H O C Cl Na")
editor.addGaseousPhase("CO2(g)")

system = ChemicalSystem(editor)
# -

# In the code below, two instances of the class
# [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html) are created:
# `initial_problem` describes the initial state, and `final_problem` corresponds to the final state.

# +
initial_problem = EquilibriumProblem(system)
initial_problem.setTemperature(30.0, "celsius")
initial_problem.setPressure(1.0, "bar")
initial_problem.add("H2O", 1, "kg")
initial_problem.add('CO2', 0.5, 'mol')
initial_problem.add('HCl', 1, 'mol')

final_problem = EquilibriumProblem(system)
final_problem.setTemperature(30.0, "celsius")
final_problem.setPressure(1.0, "bar")
final_problem.add("H2O", 1, "kg")
final_problem.add('CO2', 0.5, 'mol')
final_problem.add('NaOH', 2, 'mol')
# -

# Two instances of the class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) are created
# to store the initial and final equilibrium states calculated by the method
# [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e).

initial_state = equilibrate(initial_problem)
final_state = equilibrate(final_problem)

# Once the initial and final equilibrium states have been calculated, it is now time to trace the reaction path
# between them, with each intermediate state in chemical equilibrium. For this, we use the class
# [EquilibriumPath](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumPath.html). Note that its initialization
# requires a [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) instance:

path = EquilibriumPath(system)

# Before calling the method
# [EquilibriumPath::solve](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumPath.html#a008b74301618ed186caa95ec059eb204)
# , one can configure output-file to be generated during the calculation.
# To output quantities to a file or terminal during the calculation, use method
# [EquilibriumPath::output](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumPath.html#ac700ac6f939acbd4a8524e1346a1e588),
# which returns an instance of class [ChemicalOutput](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html):

output = path.output()
output.filename("result.txt")
output.add("speciesMolality(CO2(aq) units=mmolal)", "CO2(aq) [mmolal]")
output.add("speciesMolality(H+ units=mmolal)", "H+ [mmolal]")
output.add("speciesMolality(Cl- units=mmolal)", "Cl- [mmolal]")
output.add("speciesMolality(Na+ units=mmolal)", "Na+ [mmolal]")
output.add("speciesMolality(OH- units=mmolal)", "OH- [mmolal]")
output.add('speciesMolality(HCO3- units=mmolal)', 'HCO3- [mmolal]')
output.add('speciesMolality(CO3-- units=mmolal)', 'CO3-- [mmolal]')
output.add('chemicalPotential(CO2(aq) units=kJ/mol)', 'mu(CO2(aq)) [kJ/mol]')
output.add('chemicalPotential(HCO3- units=kJ/mol)', 'mu(HCO3-) [kJ/mol]')
output.add('chemicalPotential(CO3-- units=kJ/mol)', 'mu(CO3--) [kJ/mol]')
output.add("pH")

# The method [ChemicalOutput::add](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html#af3b5a7d6b0fbbc870664d6ad100b10dd)
# adds a quantity, which we want to be output to the file `result.txt`. The latter filename is specified in the call of
# the method [ChemicalOutput::filename](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html#ac5cc9d0f90cfe5c6e0972a55b7f7bf5d).
# Each call to [ChemicalOutput::add](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html#af3b5a7d6b0fbbc870664d6ad100b10dd)
# results in a new column of data in the output file.
#
# > **Note**: When two arguments are provided to the method
# [ChemicalOutput::add](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html#a54b0e4fd28823c4d1d1884c32eed1cf3),
# the first one is the name of the quantity to be output (e.g.,
# `time`, `elementAmount(Cl)`, `ionicStrength`). The second one is a label used as the heading of the column of data
# in the output file. When only one argument is provided, this single argument is both the label and the quantity name.
#
# Finally, after output files have been configured, the equilibrium path can be calculated using:

path = path.solve(initial_state, final_state)

# ### Plotting the results of equilibrium path calculation
#
# We now use [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples)
# to do the plotting.

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

# Besides, we define a custom function that would generate figure of a size 600 x 300:

def custom_figure(x_axis_label, y_axis_label):
    return figure(plot_width=600, plot_height=300,
                  x_axis_label=x_axis_label,
                  y_axis_label=y_axis_label)


# To load results from the outputfile, we use `loadtxt` function provided by the `numpy` package:

import numpy
filearray = numpy.loadtxt("result.txt", skiprows=1)
data = filearray.T
[co2aq_indx, h_indx, cl_indx, na_indx, oh_indx, hco3_indx, co3_indx,
 mu_co2aq_indx, mu_hco3_indx, mu_co3_indx, ph_indx] = numpy.arange(11)

# The first plot depicts the amount of species Cl<sup>-</sup> and H<sup>+</sup> in units of mmolal (on the *y*-axis)
# and the pH of the aqueous phase (on the *x*-axis). Below, we see how the concentration of both ions decrease,
# which is expected, as the initial state contains 1 mol of HCl and the final one none of it.

fig1 = custom_figure(x_axis_label="pH", y_axis_label="Species molality [mmolal]")
fig1.line(data[ph_indx], data[cl_indx], legend_label="Cl-", line_width=4, color="green")
fig1.line(data[ph_indx], data[h_indx], legend_label="H+", line_width=4, color="blue")
show(fig1)

# The second plot sets the *x*-axis to the pH level and
# the *y*-axis to the molalities of species Na<sup>+</sup> and OH<sup>-</sup>,
# i.e., the molar amount of Na<sup>+</sup> and OH<sup>-</sup> **in the aqueous
# phase**, divided by the mass of solvent water H<sub>2</sub>O(l). As expected, we see an increase of the ions amount
# as pH grows.

fig2 = custom_figure(x_axis_label="pH", y_axis_label="Species molality [mmolal]")
fig2.line(data[ph_indx], data[na_indx], legend_label="Na+", line_width=4, color="coral")
fig2.line(data[ph_indx], data[oh_indx], legend_label="OH-", line_width=4, color="gray")
show(fig2)

# The third plot sets the *x*-axis to pH, but the *y*-axis now contains three quantities: the molality of
# species CO<sub>2</sub>(aq), HCO<sub>3</sub><sup>-</sup>, and CO<sub>3</sub><sup>2-</sup>, all in units of mmolal.

fig3 = custom_figure(x_axis_label="pH", y_axis_label="Species molality [mmolal]")
fig3.line(data[ph_indx], data[co2aq_indx], legend_label="CO2(aq)", line_width=4, color="rosybrown")
fig3.line(data[ph_indx], data[hco3_indx], legend_label="HCO3-", line_width=4, color="cadetblue")
fig3.line(data[ph_indx], data[co3_indx], legend_label="CO3--", line_width=4, color="olivedrab")
show(fig3)

# We see here growing molalities of the HCO3<sup>-</sup> species until approximately pH = 7, after which it
# decreases again. Exactly around the point, where molality of HCO3<sup>-</sup> begins
# to decrease, molality of CO3<sup>2-</sup> starts increasing. Let us also plot the molality of CO<sub>2</sub>(aq)
# with respect to pH.

fig4 = custom_figure(x_axis_label="pH", y_axis_label="Species molality [mmolal]")
fig4.line(data[ph_indx], data[co2aq_indx], legend_label="CO2(aq)", line_width=4, color="rosybrown")
show(fig4)

# The plot above illustrated expected behaviour, where CO<sub>2</sub>(aq) molality in the brine remains over 25 mmolal
# until pH = 7 and rapidly drops after that value.
#
# Final fourth figure plots how the chemical potential of CO<sub>2</sub>(aq), HCO<sub>3</sub><sup>-</sup>,
# and CO<sub>3</sub><sup>2-</sup> depends on pH in the considered chemical path.

fig4 = custom_figure(x_axis_label="pH", y_axis_label="Chemical potential [kJ/mol]")
fig4.line(data[ph_indx], data[mu_co2aq_indx], legend_label="CO2(aq)", line_width=4, color="darkorange")
fig4.line(data[ph_indx], data[mu_hco3_indx], legend_label="HCO3-", line_width=4, color="palevioletred")
fig4.line(data[ph_indx], data[mu_co3_indx], legend_label="CO3--", line_width=4, color="indigo")
show(fig4)
