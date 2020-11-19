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

t0, t1 = 0.0, 1.0
result_file_name = "kineticpath-silica-metasomatism-benchmark-tfinal-" + str(t1*minute) + ".txt"

# Define chemical system:

# +
# Construct the chemical system with its phases and species
db = Database('supcrt07-organics.xml')

# Fetch Debye-Huckel activity model parameters
dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Mg Si"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase('Forsterite')    # Mg2SiO4
editor.addMineralPhase('Brucite')       # Mg(OH)2
editor.addMineralPhase('Talc')          # Mg3Si4O10(OH)2
editor.addMineralPhase('Chrysotile')    # Mg3Si2O5(OH)4, Lizardite

system = ChemicalSystem(editor)
#print(system)
#input()
# -

# We set the reaction equations, and particular parameters for calcite, dolomite, halite, k-feldspar, quartz, and
# kaolinite using chemical editor:

# +
# Forsterite
# 	Mg2SiO4 + 4 H+ =  SiO2 + 2 H2O + 2 Mg+2
# 	log_k		27.8626
# 	-delta_H	-205.614	kJ/mol
# #	deltafH		-520		kcal/mol
# 	-analytic	-7.6195e1 -1.4013e-2 1.4763e4 2.5090e1 -3.0379e5
# #	Range		0-350
# 	-Vm		43.79
# #	Extrapol	supcrt92
# #	Ref		HDN+78
# .setEquation("Forsterite + 4*H+= 2*Mg++ + SiO2(aq) + 2*H2O(l)") \
# .setEquation("Forsterite = Mg++ + SiO4--") \
editor.addMineralReaction("Forsterite") \
    .setEquation("Forsterite + 4*H+= 2*Mg++ + SiO2(aq) + 2*H2O(l)") \
    .addMechanism("logk = -10.64 mol/(m2*s); Ea = 79.0 kJ/mol") \
    .addMechanism("logk = -6.85 mol/(m2*s); Ea = 67.2 kJ/mol; a[H+] = 0.470") \
    .setSpecificSurfaceArea(1.0, "cm2/g")

# Brucite
# 	Mg(OH)2 + 2 H+ = Mg+2 + 2 H2O
# 	log_k		16.2980
# 	-delta_H	-111.34		kJ/mol
# #	deltafH		-221.39		kcal/mol
# 	-analytic	-1.0280e2 -1.9759e-2 9.0180e3 3.8282e1 1.4075e2
# #	Range		0-350
# 	-Vm		24.63
# #	Extrapol	supcrt92
# #	Ref		HDN+78
# .setEquation("Brucite = Mg++ + 2*OH-") \
# .setEquation("Brucite + 2*H+ = Mg++ + 2*H2O(l)") \
editor.addMineralReaction("Brucite") \
    .setEquation("Brucite + 2*H+ = Mg++ + 2*H2O(l)") \
    .addMechanism("logk = -8.24 mol/(m2*s); Ea = 42.0 kJ/mol") \
    .addMechanism("logk = -4.73 mol/(m2*s); Ea = 59.0 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(10, "cm2/g")

# Mg3Si4O10(OH)2
# Talc
# 	Mg3Si4O10(OH)2 + 4 H2O + 6 H+ = 3 Mg+2 + 4 H4SiO4
# 	-log_k	21.399
# 	-delta_h -46.352 kcal
# 	-Vm 68.34
# Talc
# 	Mg3Si4O10(OH)2 + 6 H+ = 3 Mg+2 + 4 H2O + 4 SiO2
# 	log_k		21.1383
# 	-delta_H	-148.737	kJ/mol
# #	deltafH		-1410.92	kcal/mol
# 	-analytic	1.1164e1 2.4724e-2 1.9810e4 -1.7568e1 -1.8241e6
# #	Range		0-350
# 	-Vm		136.25
# #	Extrapol	supcrt92
# #	Ref		HDN+78, Wilson+06 match
editor.addMineralReaction("Talc") \
    .setEquation("Talc + 6*H+ = 3*Mg++ + 4*SiO2(aq) + 4*H2O(l)") \
    .addMechanism("logk = -0.21 mol/(m2*s); Ea = 7.4 kJ/mol") \
    .setSpecificSurfaceArea(10, "cm2/g")

# Mg3Si2O5(OH)4
# Chrysotile
# 	Mg3Si2O5(OH)4 + 6 H+ = 2 SiO2 + 3 Mg+2 + 5 H2O
# 	log_k		31.1254
# 	-delta_H	-218.041	kJ/mol
# #	deltafH		-1043.12	kcal/mol
# 	-analytic	-9.2462e1 -1.1359e-2 1.8312e4 2.9289e1 -6.2342e5
# #	Range		0-350
# 	-Vm		108.5
# #	Extrapol	supcrt92
# #	Ref		HDN+78
editor.addMineralReaction("Chrysotile") \
    .setEquation("Chrysotile + 6*H+ = 2*SiO2(aq) + 3*Mg++ + 5*H2O(l)") \
    .addMechanism("logk = -12.00 mol/(m2*s); Ea = 73.5 kJ/mol") \
    .setSpecificSurfaceArea(10, "cm2/g")

reactions = ReactionSystem(editor)
# -

# Specifying the partition including the kinetic species:
partition = Partition(system)
partition.setKineticSpecies(["Forsterite", "Brucite", "Talc", "Chrysotile"])


T = 300.0 + 273.15       # temperature (in units of celsius)
P = 8.58 * 1e5   # pressure (in units of atm)
problem_ic = EquilibriumProblem(system)
problem_ic.setTemperature(T)
problem_ic.setPressure(P)
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('SiO(aq)', 7.12 * 1e-2, "mol") # 7.12 * 1e-8 mol/cm3 = 7.12 * 1e-2 mol / m3  # MM(SiO2) = 60.08 g/mol
#problem_ic.add('Forsterite', 10.0, 'mol')
#problem_ic.add('Mg2SiO4', 10.0, 'mol')

# Calculating the initial chemical equilibrium state of the fluid

# Calculate the equilibrium states for the initial conditions
state_ic = equilibrate(problem_ic)
print(state_ic)
state_ic.output('kineticpath-silica-metasomatism-initial-state.txt')
# Set the minerals initial values:

# +
# Clay minerals have the following proportions:
# Quartz        SiO2            85 %
# Calcite       CaCO3           6 %
# Dolomite      CaMg(CO3)2      4 %
# K-feldspar    K(AlSi3)O8      3 %
# Kaolinite     Al2Si2O5(OH)4   2 %

# Scale the volumes of the phases in the initial condition
#state_ic.scalePhaseVolume('Aqueous', 0.37, 'm3') # 37%
#state_ic.scalePhaseVolume('Forsterite', 0.63, 'm3') # for 63%

state_ic.setSpeciesAmount("SiO2(aq)", 7.41, "mol") # 7.41 * 1e-6 mol/cm3 = 7.41 mol / m3  # MM(SiO2) = 60.08 g/mol
state_ic.setSpeciesAmount("Forsterite", 10, "mol") # 7.41 * 1e-6 mol/cm3 = 7.41 mol / m3  # MM(SiO2) = 60.08 g/mol

# -

# Performing the kinetic path calculation:

path = KineticPath(reactions)
path.setPartition(partition)

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result_file_name`:

output = path.output()
output.filename(result_file_name)
output.add("time(units=second)")
output.add("pH")
output.add("speciesMolality(Mg++ units=mmolal)", "Mg++ [mmolal]")
output.add("speciesMolality(H+ units=mmolal)", "H+ [mmolal]")
output.add("speciesMolality(SiO2(aq) units=mmolal)", "SiO2(aq) [mmolal]")
output.add("speciesMolality(Forsterite units=mmolal)", "Forsterite [mmolal]")
output.add("speciesMolality(Brucite units=mmolal)", "Brucite [mmolal]")
output.add("speciesMolality(Talc units=mmolal)", "Talc [mmolal]")
output.add("speciesMolality(Chrysotile units=mmolal)", "Chrysotile [mmolal]")

# Solving the chemical kinetics problem:

path.solve(state_ic, t0, t1, "minute")
print(state_ic)
state_ic.output('kineticpath-silica-metasomatism-after-kinetics.txt')

# For plotting of the results of equilibrium path calculation, we load the results into the `data` array:

filearray = numpy.loadtxt(result_file_name, skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, ph_indx, mg_species_indx, h_species_indx, sio2aq_species_indx,
 forsterite_indx, brucite_indx, talc_indx, chrysotile_indx] = numpy.arange(0, 9)

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

fig1_1 = custom_figure(title="Minerals concentration w.r.t. time", y_axis_label='Concentration [mmol/m3]', y_axis_type="log")
fig1_1.line(time, data[forsterite_indx], line_width=4, color="yellow", legend_label="Forsterite")
show(fig1_1)

fig1_2 = custom_figure(title="Minerals concentration w.r.t. time", y_axis_label='Concentration [mmol/m3]', y_axis_type="log")
fig1_2.line(time, data[brucite_indx], line_width=4, color="red", legend_label="Brucite")
show(fig1_2)

fig1_3 = custom_figure(title="Minerals concentration w.r.t. time", y_axis_label='Concentration [mmol/m3]', y_axis_type="log")
fig1_3.line(time, data[talc_indx], line_width=4, color="green", legend_label="Talc")
show(fig1_3)

fig1_4 = custom_figure(title="Minerals concentration w.r.t. time", y_axis_label='Concentration [mmol/m3]', y_axis_type="log")
fig1_4.line(time, data[chrysotile_indx], line_width=4, color="green", legend_label="Chrysotile")
show(fig1_4)

fig2_1 = custom_figure(title="Aqueous species concentration w.r.t. time", y_axis_label='Concentration [mmolal]', y_axis_type="log")
fig2_1.line(time, data[mg_species_indx], line_width=4, legend_label="Mg++", color="pink")
show(fig2_1)

fig2_2 = custom_figure(title="Aqueous species concentration w.r.t. time", y_axis_label='Concentration [mmolal]', y_axis_type="log")
fig2_2.line(time, data[h_species_indx], line_width=4, legend_label="H+", color="brown")
show(fig2_2)

fig2_3 = custom_figure(title="Aqueous species concentration w.r.t. time", y_axis_label='Concentration [mmolal]', y_axis_type="log")
fig2_3.line(time, data[sio2aq_species_indx], line_width=4, legend_label="SiO2", color="gold")
show(fig2_3)


