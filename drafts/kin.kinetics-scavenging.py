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

from reaktoro import *

# Define time related constants:

second = 1
minute = 60

# Define the time interval of the kinetic simulation and corresponding to it file-name (for the future results to be
# stored):

t0, t1 = 0.0, 10
result_file_name = "kineticpath-scavenging-tfinal-" + str(t1*minute) + ".txt"

# Define chemical system:

# +
# Construct the chemical system with its phases and species
db = Database('supcrt07.xml')

# Fetch Debye-Huckel activity model parameters
dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("C Ca Cl Fe H K Mg Na O S"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase('Siderite')      # FeCO3
editor.addMineralPhase('Pyrite')        # FeS2
editor.addMineralPhase('Hematite')      # Fe2O3

system = ChemicalSystem(editor)
# -

# We set the reaction equations, and particular parameters for calcite, dolomite, halite, k-feldspar, quartz, and
# kaolinite using chemical editor:

# +
# editor.addMineralReaction("Hematite") \
#     .setEquation("Hematite + 4*H+ = 2*H2O(l) + 2*Fe++ + 0.5*O2(aq)") \
#     .addMechanism("logk = -14.60 mol/(m2*s); Ea = 66.2 kJ/mol") \
#     .addMechanism("logk = -9.39 mol/(m2*s); Ea = 66.2 kJ/mol; a[H+] = 1.0") \
#     .setSpecificSurfaceArea(10, "cm2/g")
editor.addMineralReaction("Hematite") \
    .setEquation("Hematite + 6*H+ = 3*H2O(l) + 2*Fe+++") \
    .addMechanism("logk = -14.60 mol/(m2*s); Ea = 66.2 kJ/mol") \
    .addMechanism("logk = -9.39 mol/(m2*s); Ea = 66.2 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

# .addMechanism("logk = -4.55 mol/(m2*s); Ea = 56.9 kJ/mol; a[H+] = 0.500") \
editor.addMineralReaction("Pyrite") \
    .setEquation("Pyrite + H2O(l) = 0.25*H+ + 0.25*SO4-- + Fe++ + 1.75*HS-") \
    .addMechanism("logk = -7.52 mol/(m2*s); Ea = 56.9 kJ/mol; a[H+] = 0.500") \
    .setSpecificSurfaceArea(10, "cm2/g")

reactions = ReactionSystem(editor)
# -

# Specifying the partition including the kinetic species:

partition = Partition(system)
#partition.setKineticSpecies(["Hematite", "Pyrite"])
partition.setKineticSpecies(["Hematite"])


T = 25.0 + 273.15       # temperature (in units of celsius)
P = 1 * 1.01325 * 1e5   # pressure (in units of atm)
problem_ic = EquilibriumInverseProblem(system)
problem_ic.setPartition(partition)
problem_ic.setTemperature(T)
problem_ic.setPressure(P)
problem_ic.add("H2O", 58.0, "kg")
problem_ic.add("Cl-", 1122.3e-3, "kg")
problem_ic.add("Na+", 624.08e-3, "kg")
problem_ic.add("SO4--", 157.18e-3, "kg")
problem_ic.add("Mg++", 74.820e-3, "kg")
problem_ic.add("Ca++", 23.838e-3, "kg")
problem_ic.add("K+", 23.142e-3, "kg")
problem_ic.add("HCO3-", 8.236e-3, "kg")
problem_ic.add("O2(aq)", 58e-12, "kg")
problem_ic.add("Siderite", 0.5, "mol")
problem_ic.add("Pyrite", 0.0, "mol")
problem_ic.add("Hematite", 0.0, "mol")
problem_ic.add("HS-", 0.0196504, "mol")
problem_ic.add("H2S(aq)", 0.167794, "mol")
problem_ic.pH(5.726)
problem_ic.pE(8.220)
#problem_ic.pH(8.951)
#problem_ic.pE(8.676)

# Calculating the initial chemical equilibrium state of the fluid

# Calculate the equilibrium states for the initial conditions
state_ic = equilibrate(problem_ic)
state_ic.output('shell-kinetics-benchmark-initial.txt')

# Adding 0.5 mol of each mineral
state_ic.setSpeciesMass("Siderite", 115.86 * 0.5, "g")
state_ic.setSpeciesMass("Hematite", 55.845 * 0.5, "g")

# -

# Performing the kinetic path calculation:

path = KineticPath(reactions)
path.setPartition(partition)

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result_file_name`:

output = path.output()
output.filename(result_file_name)
output.add("time(units=minute)")
output.add("pH")
output.add("speciesMolality(H+)")
output.add("speciesMolality(HS-)")
output.add("speciesMolality(S2--)")
output.add("speciesMolality(CO3--)")
output.add("speciesMolality(HSO4-)")
output.add("speciesMolality(H2S(aq))")
output.add("speciesMolality(Siderite)")
output.add("speciesMolality(Pyrite)")
output.add("speciesMolality(Hematite)")

# Solving the chemical kinetics problem:

path.solve(state_ic, t0, t1, "minute")

# For plotting of the results of equilibrium path calculation, we load the results into the `data` array:

filearray = numpy.loadtxt(result_file_name, skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, ph_indx, hs_speices_indx, s2_species_indx, co3_species_indx, hso4_species_indx, h2s_species_indx,
 siderite_indx, pyrite_indx, hematite_indx] = numpy.arange(0, 10)

# To visually analyze the obtained reaction path is with plots. For that, we export
# To visually analyze the obtained reaction path is with plots. For that, we export
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

fig0 = custom_figure(title="pH w.r.t. time", y_axis_label='pH [-]')
fig0.line(time, data[ph_indx], line_width=4, color="darkviolet")
show(fig0)

fig1_1 = custom_figure(title="Minerals concentrations w.r.t. time", y_axis_label='Concentrations [mol/m3]', y_axis_type="log")
fig1_1.line(time, data[siderite_indx], line_width=4, color="yellow", legend_label="Siderite")
show(fig1_1)

fig1_2 = custom_figure(title="Minerals concentrations w.r.t. time", y_axis_label='Concentrations [mol/m3]', y_axis_type="log")
fig1_2.line(time, data[pyrite_indx], line_width=4, color="red", legend_label="Pyrite")
show(fig1_2)

fig1_3 = custom_figure(title="Minerals concentrations w.r.t. time", y_axis_label='Concentrations [mol/m3]', y_axis_type="log")
fig1_3.line(time, data[hematite_indx], line_width=4, color="green", legend_label="Hematite")
show(fig1_3)

fig2_1 = custom_figure(title="Aqueous species concentrations w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_1.line(time, data[hs_speices_indx], line_width=4, legend_label="HS-", color="pink")
show(fig2_1)

fig2_2 = custom_figure(title="Aqueous species concentrations w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_2.line(time, data[s2_species_indx], line_width=4, legend_label="S2-", color="brown")
show(fig2_2)

fig2_3 = custom_figure(title="Aqueous species concentrations w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_3.line(time, data[co3_species_indx], line_width=4, legend_label="CO3--", color="gold")
show(fig2_3)

fig2_4 = custom_figure(title="Aqueous species concentrations w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_4.line(time, data[hso4_species_indx], line_width=4, legend_label="HSO4-", color="olive")
show(fig2_4)

fig2_5 = custom_figure(title="Aqueous species concentrations w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_5.line(time, data[h2s_species_indx], line_width=4, legend_label="H2S(aq)", color="darkblue")
show(fig2_5)

