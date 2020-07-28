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

second = 1
minute = 60

t0, t1 = 0.0, 10
result_file_name = "shell-kinetics-benchmark-tfinal-" + str(t1*minute) + ".txt"

from reaktoro import *

# Construct the chemical system with its phases and species
db = Database('supcrt07.xml')

dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H Cl S O Ca Sr Na K Mg C Si Al"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase('Halite')        # NaCl
editor.addMineralPhase('Calcite')       # CaCO3
editor.addMineralPhase('Dolomite')      # CaMg(CO3)2
editor.addMineralPhase('Quartz')        # SiO2
editor.addMineralPhase('K-Feldspar')    # K(AlSi3)O8
editor.addMineralPhase('Kaolinite')     # Al2Si2O5(OH)4

system = ChemicalSystem(editor)

# We set the reaction equation using

reaction = editor.addMineralReaction("Calcite")
reaction.setEquation("Calcite = Ca++ + CO3--")
reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
reaction.setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Dolomite") \
    .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
    .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
    .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Halite") \
    .setEquation("Halite = Na+ + Cl-") \
    .addMechanism("logk = -0.21 mol/(m2*s); Ea = 7.4 kJ/mol") \
    .setSpecificSurfaceArea(10, "cm2/g")

# K-Feldspar: K(AlSi3)O8 # potassium feldspar, orthoclase
# KAlSi3O8 + 4*H2O(l) + 4*H+ = Al+++ + 3*H4SiO4 + K+
# KAlSi3O8 + 4*H2O(l) + 4*H+ = Al+++ + 3*(H2O(l) + H+ + HSiO3-) + K+
# KAlSi3O8 + H2O(l) + H+ = Al+++ + 3*HSiO3- + K+
editor.addMineralReaction("K-Feldspar") \
    .setEquation("K-Feldspar + H2O(l) + H+ = Al+++ + 3*HSiO3- + K+") \
    .addMechanism("logk = -12.41 mol/(m2*s); Ea = 38.0 kJ/mol") \
    .addMechanism("logk = -10.06 mol/(m2*s); Ea = 51.7 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(10, "cm2/g")
# PHREEQC: K-Feldspar + 4*H+  =  Al+++ + K+ + 2*H2O(l) + 3*SiO2(aq)
# K-Feldspar
#         KAlSi3O8 + 4 H+  =  Al+++ + K+ + 2 H2O + 3 SiO2
#         log_k           -0.2753
#         -delta_H        -23.9408        kJ/mol  # Calculated enthalpy of reaction  K-Feldspar
# #        Enthalpy of formation:        -949.188 kcal/mol
# #        -analytic -1.0684e+000 1.3111e-002 1.1671e+004 -9.9129e+000 -1.5855e+006
#         -analytic 1.02685E+02	-1.40881E-02	-1.15294E+04	-3.05046E+01	1.38991E+06
# #       -Range:  0-300

# SiO2(quartz) + 2*H2O(l) = H4SiO4(aq)
# SiO2(quartz) + H2O(l) = H+ + HSiO3-
#.addMechanism("logk = -16.29 mol/(m2*s); Ea = 108366 kJ/mol; a[H+] = -0.5") \
editor.addMineralReaction("Quartz") \
    .setEquation("Quartz + H2O(l) = H+ + HSiO3-") \
    .addMechanism("logk = -13.99 mol/(m2*s); Ea = 87.7 kJ/mol") \
    .setSpecificSurfaceArea(10, "cm2/g")
# PHREEQC version: Quartz = SiO2(aq)
# Quartz
#         SiO2  =  SiO2
#         log_k           -3.9993
#         -delta_H        32.949        kJ/mol  # Calculated enthalpy of reaction  Quartz
# #        Enthalpy of formation:        -217.65 kcal/mol
# #        -analytic 7.7698e-002 1.0612e-002 3.4651e+003 -4.3551e+000 -7.2138e+005
#         -analytic 1.54450E+02	1.78164E-02	-1.09003E+04	-5.42517E+01	6.48510E+05
# #       -Range:  0-300


# # Al2Si2O5(OH)4
# Al2Si2O5(OH)4 + 6*H+ = H2O(l) + 2*H4SiO4 + 2*Al+++
# Al2Si2O5(OH)4 + 6*H+ = H2O(l) + 2*(H2O(l) + H+ + HSiO3-) + 2*Al+++
# Al2Si2O5(OH)4 + 4*H+ = 3*H2O(l) + 2*HSiO3- + 2*Al+++
# .addMechanism("logk = -17.05 mol/(m2*s); Ea = 17.9 kJ/mol; a[H+] = -0.472") \
editor.addMineralReaction("Kaolinite") \
    .setEquation("Kaolinite + 4*H+ = 3*H2O(l) + 2*HSiO3- + 2*Al+++") \
    .addMechanism("logk = -13.18 mol/(m2*s); Ea = 22.2 kJ/mol") \
    .addMechanism("logk = -11.31 mol/(m2*s); Ea = 65.9 kJ/mol; a[H+] = 0.777") \
    .setSpecificSurfaceArea(10, "cm2/g")
# PHREEQC version: Kaolinite + 6*H+  =  2*Al+++ + 2*SiO2(aq) + 5*H2O(l)
# Kaolinite
#         Al2Si2O5(OH)4 + 6 H+  =  2 Al+++ + 2 SiO2 + 5 H2O
#         log_k           6.8101
#         -delta_H        -151.779        kJ/mol  # Calculated enthalpy of reaction  Kaolinite
# #        Enthalpy of formation:        -982.221 kcal/mol
# #        -analytic 1.6835e+001 -7.8939e-003 7.7636e+003 -1.2190e+001 -3.2354e+005
#         -analytic -2.58805E+03	-3.42587E-01	1.60365E+05	9.15337E+02	-9.48341E+06
# #       -Range:  0-300


reactions = ReactionSystem(editor)

# ### Specifying the equilibrium and kinetic species

partition = Partition(system)
partition.setKineticSpecies(["Calcite", "Dolomite", "Halite", "Quartz", "K-Feldspar", "Kaolinite"])


# ### Defining the initial state of the equilibrium species
"""
temp 61
pressure 1
pH    7
units mol / kgw
Na 6.27
Cl 6.27 charge
Mg 0.0001
Ca 0.0001
C 0.0001
Si 0.0001
Al 0.0001
K 0.000000001
"""
T = 61.0 + 273.15       # temperature (in units of celsius)
P = 1 * 1.01325 * 1e5   # pressure (in units of atm)
problem_ic = EquilibriumInverseProblem(system)
problem_ic.setPartition(partition)
problem_ic.setTemperature(T)
problem_ic.setPressure(P)
problem_ic.add('H2O', 1.0, 'kg')
problem_ic.add('Na', 6.27 , 'mol') # 6.27 mol / kgw
problem_ic.add('Cl', 6.27, 'mol')  # 6.27 mol / kgw
problem_ic.add('Mg', 0.0001, 'mol')  # 0.0001 mol / kgw
problem_ic.add('Ca', 0.0001, 'mol')  # 0.0001 mol / kgw
problem_ic.add('C', 0.0001, 'mol')  # 0.0001 mol / kgw
problem_ic.add('Si', 0.0001, 'mol')  # 0.0001 mol / kgw
problem_ic.add('Al', 0.0001, 'mol')  # 6.27 mol / kgw
problem_ic.add('K', 1e-9, 'mol')
problem_ic.pH(7.0)

# ### Calculating the initial chemical equilibrium state of the fluid

# Calculate the equilibrium states for the initial conditions
state_ic = equilibrate(problem_ic)
state_ic.output('shell-kinetics-benchmark-initial.txt')

# Clay minerals have the following proportions:
# Quartz        SiO2            85 %
# Calcite       CaCO3           6 %
# Dolomite      CaMg(CO3)2      4 %
# K-feldspar    K(AlSi3)O8      3 %
# Kaolinite     Al2Si2O5(OH)4   2 %

# Scale the volumes of the phases in the initial condition
# state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
# state_ic.scalePhaseVolume('Quartz', 0.9 * 0.85, 'm3') # for 85 %
# state_ic.scalePhaseVolume('Calcite', 0.9 * 0.06, 'm3') # for 6 %
# state_ic.scalePhaseVolume('Dolomite', 0.9 * 0.04, 'm3')  # for 4 %
# state_ic.scalePhaseVolume('K-Feldspar', 0.9 * 0.03, 'm3')  # for 3 %
# state_ic.scalePhaseVolume('Kaolinite', 0.9 * 0.02, 'm3')  # for 2 %

state_ic.setSpeciesMass("Calcite", 100.0869 * 10, "g") #  molar mass of CaCO3 = 100.0869 g/mol
state_ic.setSpeciesMass("Dolomite", 184.4008 * 10, "g") # molar mass of CaMg(CO3)2 = 184.4008 g/mol
state_ic.setSpeciesMass("Halite", 58.44 * 10, "g") # molar mass of NaCl = 58.44 g/mol
state_ic.setSpeciesMass("Quartz", 60.08 * 10, "g") # molar mass of SiO2 = 60.08 g/mol
state_ic.setSpeciesMass("Kaolinite", 258.1604 * 10, "g") # molar mass of Al2Si2O5(OH)4 =  258.1604 g/mol
state_ic.setSpeciesMass("K-Feldspar", 278.3315 * 10, "g") # molar mass of K(AlSi3)O8 = 278.3315 g/mol

# ### Performing the kinetic path calculation

path = KineticPath(reactions)
path.setPartition(partition)

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result.txt`:

output = path.output()
output.filename(result_file_name)
output.add("time(units=minute)")
output.add("pH")
output.add("speciesMolality(Na+ units=mmolal)", "Na+ [mmol]")
output.add("speciesMolality(Cl- units=mmolal)", "Cl- [mmol]")
output.add("speciesMolality(Mg++ units=mmolal)", "Mg++ [mmol]")
output.add("speciesMolality(Ca++ units=mmolal)", "Ca++ [mmol]")
output.add("speciesMolality(CO3-- units=mmolal)", "CO3-- [mmol]")
output.add("speciesMolality(K+ units=mmolal)", "K+ [mmol]")
output.add("speciesMolality(Al+++ units=mmolal)", "Al+++ [mmol]")
output.add("speciesMolality(Calcite units=molal)", "Calcite [mol]")
output.add("speciesMolality(Dolomite units=molal)", "Dolomite [mol]")
output.add("speciesMolality(Halite units=molal)", "Halite [mol]")
output.add("speciesMolality(Quartz units=molal)", "Quartz [mol]")
output.add("speciesMolality(K-Feldspar units=molal)", "K-Feldspar [mol]")
output.add("speciesMolality(Kaolinite units=molal)", "Kaolinite [mol]")


# ### Solving the chemical kinetics problem
path.solve(state_ic, t0, t1, "minute")

# ### Plotting the results of equilibrium path calculation

filearray = numpy.loadtxt(result_file_name, skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, ph_indx, na_speices_indx, cl_species_indx, mg_species_indx, ca_species_indx, co3_species_indx, k_species_indx,
 al_species_indx, calcite_indx, dolomite_indx, halite_indx, quartz_indx, kfeldspar_indx, kaolinite_indx] = numpy.arange(0, 15)

# assign indices of the corresponding data

# To visually analyze the obtained reaction path is with plots. For that, we export
# [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples) python plotting package.

from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

def custom_figure(title, y_axis_label, y_axis_type='auto'):
    return figure(plot_width=600, plot_height=300,
                  title=title,
                  x_axis_label='time',
                  y_axis_label=y_axis_label,
                  y_axis_type=y_axis_type,
                  background_fill_color="#fafafa")


time = data[time_indx, :]  # fetch time from the data matrix

# The increase of the amount of Ca element is happening along the dissolution of calcite on the
# plot below:

fig0 = custom_figure(title="pH w.r.t. time", y_axis_label='pH [-]')
fig0.line(time, data[ph_indx], line_width=4, color="darkviolet")
show(fig0)

fig1_1 = custom_figure(title="Minerals molality w.r.t. time", y_axis_label='Molality [molal]', y_axis_type="log")
fig1_1.line(time, data[halite_indx], line_width=4, color="yellow", legend_label="Halite")
show(fig1_1)

fig1_2 = custom_figure(title="Minerals molality w.r.t. time", y_axis_label='Molality [molal]', y_axis_type="log")
fig1_2.line(time, data[calcite_indx], line_width=4, color="red", legend_label="Calcite")
fig1_2.line(time, data[dolomite_indx], line_width=4, color="orange", legend_label="Dolomite")
fig1_2.line(time, data[quartz_indx], line_width=4, color="indigo", legend_label="Quartz")
show(fig1_2)


fig1_3 = custom_figure(title="Minerals molality w.r.t. time", y_axis_label='Molality [molal]', y_axis_type="log")
fig1_3.line(time, data[kfeldspar_indx], line_width=4, color="green", legend_label="K-Feldspar")
fig1_3.line(time, data[kaolinite_indx], line_width=4, color="purple", legend_label="Kaolinite")
show(fig1_3)

# As calcite dissolves, the molallities of species Ca<sup>2+</sup> and HCO3<sup>-</sup> are growing too:

fig2_1 = custom_figure(title="Aqueous species molality w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_1.line(time, data[na_speices_indx], line_width=4, legend_label="Na+", color="pink")
show(fig2_1)

fig2_2 = custom_figure(title="Aqueous species molality w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_2.line(time, data[cl_species_indx], line_width=4, legend_label="Cl-", color="brown")
show(fig2_2)


fig2_3 = custom_figure(title="Aqueous species molality w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_3.line(time, data[co3_species_indx], line_width=4, legend_label="CO3--", color="gold")
fig2_3.line(time, data[al_species_indx], line_width=4, legend_label="Al+++", color="darkseagreen")
show(fig2_3)

fig2_4 = custom_figure(title="Aqueous species molality w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_4.line(time, data[mg_species_indx], line_width=4, legend_label="Mg++", color="olive")
fig2_4.line(time, data[ca_species_indx], line_width=4, legend_label="Ca++", color="gray")
show(fig2_4)


fig2_5 = custom_figure(title="Aqueous species molality w.r.t. time", y_axis_label='Molality [mmolal]', y_axis_type="log")
fig2_5.line(time, data[k_species_indx], line_width=4, legend_label="K+", color="darkblue")
show(fig2_5)
