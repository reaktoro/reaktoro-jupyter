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

# # Precipitation of barite as a result of waterflooding
#
# This tutorial demonstrates how Reaktoro can be used for modeling barite precipitation as a result of the waterflooding
# technique often used in oil and gas industry. Waterflooding in this tutorial will be modeled by mixing chemical
# state representing the formation water (FW) and seawater (SW).
#
# We start with again import the reaktoro Python package (so that we can use its classes and methods for performing the
# chemical reaction calculations), initializing auxiliary variables, and thermodynamic condition.

# +
# Import reaktoro package
from reaktoro import *

# Define time related constants
second = 1
minute = 60 * second
hour = 60 * minute
day = 60 * hour

# Thermodynamic conditions
T = 60.0        # temperature (in units of Celsius)
P = 200.0       # pressure (in units of atm)
water_kg = 1.00 # water mass

# -

# Next, we construct the chemical system with its phases and species and fetch Debye-Huckel activity model parameters.

# +
# Define database
db = Database('supcrt07.xml')

# Fetch Debye-Huckel activity model parameters
dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()
# -

# For simulation below, we construct [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html)
# with its phases and species using [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html)
# class.

# +
# Define different phases of the chemical system
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H Cl S O Ba Ca Sr Na K Mg C Si"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase('Barite')    # BaSiO4

# Define barite mineral reaction and its parameters
eq_str_barite = "Barite = SO4-- + Ba++"
min_reaction_barite = editor.addMineralReaction("Barite") \
    .setEquation(eq_str_barite) \
    .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol") \
    .setSpecificSurfaceArea(0.006, "m2/g")
# -

# After initializing the system's phases, we create instances of
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html). Instance `system` represents the
# considered chemical system and its attributes and properties, such as phases (in our case aqueous and mineral ones),
# species, elements, formula matrix, as well as chemical and thermodynamical model. Class
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) serves as a set of the chemical
# reaction by a collection of [Reaction](https://reaktoro.org/cpp/classReaktoro_1_1Reaction.html) class instances.
# It provides convenient methods that calculate the equilibrium constants, reaction quotients, and rates of the
# reactions.

# Initialize chemical system
system = ChemicalSystem(editor)
# Initialize reaction system
reactions = ReactionSystem(editor)

# ### Specifying the equilibrium and kinetic species
#
# For the partition of a chemical system into equilibrium and kinetic species, we use the class
# [Partition](https://reaktoro.org/cpp/classReaktoro_1_1Partition.html). We only need
# to specify which species comprise kinetic species (barite, in this case), and all others will be equilibrium species
# by default. This will allow us to model the precipitation of barite using chemical kinetics, while all other species
# modeled by equilibrium.

partition = Partition(system)
partition.setKineticSpecies(["Barite"])

# ### Defining the formation water chemical state
#
# The chemical state corresponding to the formation water (FW) is taken from the manuscript of Bethke 2008, i.e.,
# Table 30.1 (Miller analysis). We recite it in the table below:
#
# | Aqueous species | Amount (mg / kg) |
# |-----------------|------------------|
# | Na<sup>+</sup>  | 27250            |
# | K<sup>+</sup>   | 1730             |
# | Mg<sup>2+</sup> | 110              |
# | Ca<sup>2+</sup> | 995              |
# | Sr<sup>2+</sup> | 105              |
# | Ba<sup>2+</sup> | 995              |
# | Cl<sup>-</sup>  | 45150            |
# | HCO<sub>3</sub><sup>-</sup> | 1980 |
# | SO<sub>4</sub><sup>2-</sup> | 10 · 10<sup>-3</sup>|
#
# The essential characteristics of the formation water are high concentration of the Ba<sup>2+</sup> and low
# concentrations of the SO<sub>4</sub><sup>2-</sup>. We also, the small initial seed of the barite into this state,
# to increase the precipitation scale.

# +
problem_fw = EquilibriumInverseProblem(system)
problem_fw.setTemperature(T, "celsius")
problem_fw.setPressure(P, "atm")
problem_fw.add('H2O', water_kg, 'kg')
problem_fw.add("SO4", 10 * water_kg, "ug")
problem_fw.add("Ca", 995 * water_kg, "mg")
problem_fw.add("Ba", 995 * water_kg, "mg")
problem_fw.add("Sr", 105 * water_kg, "mg")
problem_fw.add("Na", 27250 * water_kg, "mg")
problem_fw.add("K", 1730 * water_kg, "mg")
problem_fw.add("Mg", 110 * water_kg, "mg")
problem_fw.add("Cl", 45150 * water_kg, "mg")
problem_fw.add("HCO3", 1980 * water_kg, "mg")
problem_fw.pH(7.0, "HCl", "NaOH")

# Calculate the equilibrium states for the initial conditions
state_fw = equilibrate(problem_fw)
state_fw.setSpeciesAmount("Barite", 0.1, "mcmol")
state_fw.scaleVolume(1.0, "m3")
# -

# ### Defining the seawater chemical state

# Next, we define the seawater (SW) composition taken from Bethke 2008, i.e., Table 30.1 (Seawater). Further down, we
# list the quantities of the aqueous species:
#
# | Aqueous species | Amount (mg / kg) |
# |-----------------|------------------|
# | Na<sup>+</sup>  | 10760            |
# | K<sup>+</sup>   | 399              |
# | Mg<sup>2+</sup> | 1290             |
# | Ca<sup>2+</sup> | 411              |
# | Sr<sup>2+</sup> | 8                |
# | Ba<sup>2+</sup> | 0.01             |
# | Cl<sup>-</sup>  | 19350            |
# | HCO<sub>3</sub><sup>-</sup> | 142  |
# | SO<sub>4</sub><sup>2-</sup> | 2710 |
#
# Typically, seawater is rich in sulfate, poor in Ca<sup>2+</sup>, and nearly depleted in Sr<sup>2+</sup> and
# Ba<sup>2+</sup>. The pH of the seawater is fixed to 8.1.

# +
problem_sw = EquilibriumInverseProblem(system)
problem_sw.setTemperature(T, "celsius")
problem_sw.setPressure(P, "atm")
problem_sw.add('H2O', water_kg, 'kg')
problem_sw.add("SO4--", 2710 * water_kg, "mg")
problem_sw.add("Ca++", 411 * water_kg, "mg")
problem_sw.add("Ba++", 0.01 * water_kg, "mg")
problem_sw.add("Sr++", 8 * water_kg, "mg")
problem_sw.add("Na+", 10760 * water_kg, "mg")
problem_sw.add("K+", 399 * water_kg, "mg")
problem_sw.add("Mg++", 1290 * water_kg, "mg")
problem_sw.add("Cl-", 19350 * water_kg, "mg")
problem_sw.add("HCO3-", 142 * water_kg, "mg")
problem_sw.pH(8.1, "HCl", "NaOH")

# Calculate the equilibrium states for the initial conditions
state_sw = equilibrate(problem_sw)
state_sw.scaleVolume(1.0, "m3")
# -

# ### Solving kinetic barite precipitation problem
#
# Next, we define the time interval of the kinetic simulation and corresponding to it file-name (for the obtained
# numerical results to be stored):

t0, tfinal = 0.0, 1.0
result_file_name = "kineticpath-barite-precipitation-tfinal-" + str(tfinal*day) + ".txt"

# To be able to simulate the chemical kinetic path, we use the class
# [KineticPath](https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html). Note that here again, we need to
# specify the partitioning of the chemical system into equilibrium, kinetic, and inert species.

path = KineticPath(reactions)
path.setPartition(partition)

# To analyse the result of kinetic simulations, we save the evolution of various characteristics of the chemical
# state (with respect to time) into file `result_file_name`:

output = path.output()
output.filename(result_file_name)
output.add("time(units=s)")
output.add("speciesAmount(Barite units=mol)", "Barite")
output.filename(result_file_name)

# We solve the kinetic path problem with the method `solve()` of class
# [KineticPath](https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html). It requires the formation water state
# `state_fw` perturbed with the seawater `state_sw`, the initial and final time of the kinetic path calculation
# (`t0` and `tfinal`, respectively), and the time unit of the specified time parameters (days, in this case).

state_fw = state_fw + state_sw
path.solve(state_fw, t0, tfinal, "days")

# For plotting of the results of equilibrium path calculation, we load the results into the `data` array:

filearray = numpy.loadtxt(result_file_name, skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, barite_indx] = numpy.arange(0, 2)

# To visually analyze the obtained reaction path, we export
# [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples) python plotting package.

# +
from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

def custom_figure(title, y_axis_label, y_axis_type='auto'):
    return figure(plot_width=400, plot_height=200,
                  title=title,
                  x_axis_label='time',
                  y_axis_label=y_axis_label,
                  y_axis_type=y_axis_type,
                  background_fill_color="#fafafa")

time = data[time_indx, :]  # fetch time from the data matrix

fig = custom_figure(title="Minerals concentration w.r.t. time", y_axis_label='Concentration [mol/m3]')
fig.line(time, data[barite_indx], line_width=4, color="darkviolet", legend_label="Barite")
fig.legend.location = 'bottom_right'
show(fig)
# -

# From the plot, we see that the mineral precipitates right after addition of the seawater. It happens due to the
# contrasting compositions of the formation water (FW) and seawater (SW). In particular, SW
# is low on Ba<sup>2+</sup> and high on SO<sub>4</sub><sup>2-</sup> concentrations, whereas FW, on the opposite, is
# high on Ba<sup>2+</sup> and low on SO<sub>4</sub><sup>2-</sup>. During the mixing, both of these ions react creating
# barite according to the reaction Ba<sup>2+</sup> + SO<sub>4</sub><sup>2-</sup> &#8594; BaSO<sub>4</sub>(s).
# Such side effect of waterflooding (as part of the oil recovery techniques) reduces the near-wellbore permeability and
# hampers well productivity/injectivity.
# > **Note**: Other minerals that have the potential for scaling (when specific conditions are provided) are
# CaCO<sub>3</sub> (calcite), CaSO<sub>4</sub> (calcium sulfate), FeCO<sub>3</sub> (siderite).
