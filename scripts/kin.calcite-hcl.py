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
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Dissolution of calcite in an acidic HCl-solution
#
# This tutorial demonstrates how Reaktoro can be used for modeling the dissolution of calcite in an acidic
# HCl-solution at temperature 30 &deg;C and pressure 1 bar using chemical kinetics. A partial equilibrium
# assumption is considered here so that aqueous species react using a chemical equilibrium model, while calcite
# reacts with the aqueous solution using a chemical kinetics model.
#
# We start with again import the reaktoro Python package so that we can use its classes and methods for performing the
# chemical reaction calculations.

from reaktoro import *

# The class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) is used to conveniently
# create instances of classes [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html). In particular, we specify aqueous
# and mineral phases that should be considered in the chemical system. The aqueous phase is defined by the mixing of
# H2O, HCl, and CaCO<sub>3</sub> (effectively, collecting all aqueous species in the database that contains elements H,
# O, C, Cl,
# and Ca, which are the elements in this list of compounds). There is only one pure mineral phase: the calcite phase.

editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf("H2O HCl CaCO3")
editor.addMineralPhase("Calcite")

# We set the reaction equation using
# [setEquation](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#af3c6212c6edb42c0e6f110b493ece45c)
# method of the class [MineralReaction](https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html).
# Then, we add two mineral kinetic mechanisms for the reaction: neutral and acidic. This is done with the
# [addMechanism](https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html#a1bdeff51e51f42e4635208241cd54027)
# method, where we set, for example, `logk`, the kinetic rate constant of the reaction in log scale, and `Ea`, the
# Arrhenius activation energy. The values shown for `logk` and `Ea` were collected from:
#
# *Palandri, J.L., Kharaka, Y.K. (2004). A compilation of rate parameters of water-mineral interaction kinetics for
# application to geochemical modeling. U.S. Geological Survey Open File Report (Vol. 2004–1068). Menlo Park,
# California.*
#
# Finally, we provide the specific surface area of the mineral using method
# [setSpecificSurfaceArea](https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html#a9ea2feb68af0beddc856d6a60b863181)
# of class [MineralReaction](https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html), which can be specified
# in units of m<sup>2</sup>/g or m<sup>2</sup>/m<sup>3</sup>. Compatible units are allowed, such as cm<sup>2</sup>/mg or
# m<sup>2</sup>/dm<sup>3</sup>, and combinations.

editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

# Create instances of [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html).
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) is a class that represents a
# system and its attributes and properties, such as phases (in our case aqueous and mineral ones), species,
# elements (5 in total, i.e., H, O, Ca, C, Cl), formula matrix, as well as chemical and thermodynamical model. Class
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) serves as a system of the chemical
# reaction by a collection of [Reaction](https://reaktoro.org/cpp/classReaktoro_1_1Reaction.html) class instances.
# It provides convenient methods that calculate the equilibrium constants, reaction quotients, and rates of the
# reactions.

system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

# ### Specifying the equilibrium and kinetic species
#
# For the partition of a chemical system into equilibrium and kinetic species, we use the class
# [Partition](https://reaktoro.org/cpp/classReaktoro_1_1Partition.html). We only
# need to specify which species are kinetic species, and all others will be equilibrium species by default.
# We set species Calcite (the only species in the mineral phase also called Calcite!) to be the only kinetic species.
# This will allow us to model the dissolution of calcite using chemical kinetics, while all other species (the
# aqueous species) are modeled using chemical equilibrium (i.e., their amounts are updated over time using chemical
# equilibrium calculations).

partition = Partition(system)
partition.setKineticSpecies(["Calcite"])


# ### Defining the initial state of the equilibrium species
#
# After constructing the chemical system and specifying the partitioning of the species, we proceed to a chemical
# equilibrium calculation to set the initial state of the equilibrium species. For this, we use class
# [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html), as shown below:

problem = EquilibriumProblem(system)
problem.setPartition(partition)
problem.setTemperature(30, "celsius")
problem.setPressure(1, "bar")
problem.add("H2O", 1, "kg")
problem.add("HCl", 1, "mmol")

# We specify the equilibrium/kinetic partitioning of the chemical system using the method
# [setPartition](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html#a53a9496c9d4ffc72a85903146b390e44)
# of class [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html). We then prescribe
# what should be the initial state of the equilibrium species (the aqueous species in this case), before we start the
# chemical kinetics calculation that will simulate the dissolution of calcite in this aqueous fluid.
#
# By mixing 1 kg of H<sub>2</sub>O and 1 mmol of HCl at 30 &deg;C and 1 bar, we should produce a
# chemical equilibrium state that corresponds to an acidic aqueous fluid. The species in this fluid will be in
# disequilibrium with Calcite (our single kinetic species in this setup) since only equilibrium species
# (i.e., the aqueous species) are considered during the next chemical equilibrium calculation.
#
# ### Calculating the initial chemical equilibrium state of the fluid
#
# We now use the function equilibrate to calculate the chemical equilibrium state of the equilibrium partition,
# not the entire chemical system.

state0 = equilibrate(problem)

# For this calculation, Reaktoro uses an efficient Gibbs energy minimization algorithm to determine the amounts of
# the equilibrium species that correspond to a state of minimum Gibbs energy in the equilibrium partition only,
# at given conditions of temperature, pressure, and element amounts in the equilibrium partition. The result is
# stored in the object state0 of class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html),
# a computational representation of the state of a multiphase
# chemical system defined by its temperature (*T*), pressure (*P*), and vector of species amounts (*n*).
#
# To simulate the kinetic dissolution of calcite in the aqueous fluid we defined before, we need to specify its
# initial amount. Below, we set the initial mass of species Calcite to 100 g.

state0.setSpeciesMass("Calcite", 100, "g")

# ### Performing the kinetic path calculation
#
# To be able to simulate the chemical kinetic path, we use class
# [KineticPath](https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html). Note that here again, we need to
# specify the partitioning of the chemical system into equilibrium, kinetic, and inert species.

path = KineticPath(reactions)
path.setPartition(partition)

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result.txt`:

output = path.output()
output.filename("results.txt")
output.add("time(units=minute)")
output.add("elementMolality(Ca units=mmolal)", "Ca [mmolal]")
output.add("phaseMass(Calcite units=g)", "Calcite [units=g]")
output.add("speciesMolality(Ca++ units=mmolal)", "Ca++ [mmol]")
output.add("speciesMolality(HCO3- units=mmolal)", "HCO3- [mmol]")
output.add("pH")

# > **Note**: A list of all possible quantities that can be plotted is shown in the class
# > [ChemicalQuantity](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalQuantity.html),
# > which provides an interface for convenient ways of their retrieval.
#
# ### Solving the chemical kinetics problem
#
# Finally, we solve the kinetic path problem.
# This step executes the method solve of class
# [KineticPath](https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html), which requires the initial state of the
# system (100 g of calcite in disequilibrium with a 1 mmolal HCl aqueous solution at 30 &deg;C and 1 bar,
# represented with the object state0), the initial and final time of the kinetic path calculation (`t0` and `t1`,
# respectively), and the time unit of the specified time parameters (e.g., s, minute, day, year, etc.).

t0, t1 = 0.0, 5.0
path.solve(state0, t0, t1, "minute")

# ### Plotting the results of equilibrium path calculation
#
# To load results from the outputfile, we use [loadtxt](https://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html)
# function provided by the *numpy* package:

filearray = numpy.loadtxt("results.txt", skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, ca_elem_indx, calcite_indx, ca_species_indx, hco3_indx, ph_indx] = numpy.arange(0, 6) # assign indices of the corresponding data

# To visually analyze the obtained reaction path is with plots. For that, we export
# [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples) python plotting package.

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

# Below, we define a custom function that would generate figure of certain size (in this case, 600 by 300) with label
# `time` on the x-axis:

def custom_figure(title, y_axis_label):
    return figure(plot_width=600, plot_height=300,
                  title=title,
                  x_axis_label='time',
                  y_axis_label=y_axis_label)

# The plots below depict different chemical properties (x-axis) with respect to the time interval of the kinetic
# simulation (y-axis). We start from the behavior of the amount of element Ca with respect to time:

time = data[time_indx, :]  # fetch time from the data matrix
fig1 = custom_figure(title="Amount of Ca w.r.t. time", y_axis_label='Amount of Ca [mmolal]')
fig1.line(time, data[ca_elem_indx], line_width=4, color="coral")
show(fig1)

# The increase of the amount of Ca element is happening along the dissolution of calcite on the
# plot below:

fig2 = custom_figure(title="Mass of Calcite w.r.t. time", y_axis_label='Mass of Calcite [g]')
fig2.line(time, data[calcite_indx], line_width=4, color="blue")
show(fig2)

# As calcite dissolves, the molallities of species Ca<sup>2+</sup> and HCO3<sup>-</sup> are growing too:

fig3 = custom_figure(title="Species molality w.r.t. time", y_axis_label='Molality [mmolal]')
fig3.line(time, data[ca_species_indx], line_width=4, legend_label="Ca++", color="orange")
fig3.line(time, data[hco3_indx], line_width=4, legend_label="HCO3-", color="green")
show(fig3)

# Finally, the pH of the overall chemical system is increasing as well:

fig4 = custom_figure(title="pH w.r.t. time", y_axis_label='pH [-]')
fig4.line(time, data[ph_indx], line_width=4, color="darkviolet")
show(fig4)

