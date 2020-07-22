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

# # Simple 1D RTM model of scaling
#

# +
print('============================================================')
print('Make sure you have the following Python packages installed: ')
print('     numpy, natsort, bokeh')
print('These can be installed with pip:')
print('     pip install numpy natsort bokeh')
print('============================================================')
from reaktoro import *
import numpy as np
from natsort import natsorted
from tqdm.notebook import tqdm
import os

# Import components of bokeh library
from bokeh.io import show, output_notebook
from bokeh.layouts import column
from bokeh.plotting import figure
from bokeh.models import Range1d, ColumnDataSource
from bokeh.layouts import gridplot
# -


second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# +
# ## Defining parameters for the reactive transport simulation
#
# +
# Discretization parameters
xl = 0.0                # x-coordinate of the left boundary
xr = 25.0              # x-coordinate of the right boundary
ncells = 243            # number of cells in the discretization
dx = (xr - xl) / ncells # length of the mesh cells (in units of m)
dt = 1 * hour           # time step

nsteps_cb = 45
nsteps_sw = 255
t_cb = 45 * hour
t_sw = 255 * hour
nsteps = nsteps_cb + nsteps_sw # number of steps in the reactive transport simulation

water_kg = 58

# Physical parameters
D = 0               # diffusion coefficient (in units of m2/s)
v = 1.05e-5         # fluid pore velocity (in units of m/s)
T = 60.0            # temperature (in units of celsius)
P = 200             # pressure (in units of atm)
phi = 0.1           # the porosity
# -

# Next, we generate the coordinates of the mesh nodes (array `xcells`) by equally dividing the interval *[xr, xl]* with
# the number of cells `ncells`. The length between each consecutive mesh node is computed and stored in `dx` (the
# length of the mesh cells).

xcells = np.linspace(xl, xr, ncells + 1)  # interval [xl, xr] split into ncells

# The boolean variable `dirichlet` is set to `True` or `False` depending on which boundary condition is considered in
# the numerical calculation. `False` corresponds to imposing the flux of the injected fluid, otherwise, `True` means
# imposing the composition of the fluid on the left boundary.

dirichlet = False  # parameter that determines whether Dirichlet BC must be used

# To make sure that the applied finite-volume scheme is stable, we need to keep track of Courant–Friedrichs–Lewy (CFL)
# number, which should be less than 1.0.

CFL = v * dt / dx
assert CFL <= 1.0, f"Make sure that CFL = {CFL} is less that 1.0"

# ## Specifying the quantities and properties to be outputted
#
# Before running the reactive transport simulations, we specify the list of parameters we are interested in outputting.
# In this case, it is `pH`, molality of `H+`, `HS-`, `S2--`, `SO4--`, `H2S(aq)`, as well as a phase amount/volume
# of pyrrhotite and siderite.

output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Cl-)
    speciesMolality(SO4--)
    speciesMolality(Ba++)
    speciesMolality(Ca++)
    speciesMolality(Sr++)
    speciesMolality(Na+)
    phaseAmount(Barite)
    activity(Ba++)
    activity(SO4--)
    activityCoefficient(Ba++)
    activityCoefficient(SO4--)
""".split()

# output_quantities = """
#     pH
#     speciesMolality(H+)
#     speciesMolality(Cl-)
#     speciesMolality(SO4--)
#     speciesMolality(Ba++)
#     speciesMolality(Ca++)
#     speciesMolality(Sr++)
#     speciesMolality(Na+)
#     phaseAmount(Barite)
#     phaseAmount(Anhydrite)
#     phaseAmount(Celestite)
#     activity(Ba++)
#     activity(SO4--)
#     activityCoefficient(Ba++)
#     activityCoefficient(SO4--)
# """.split()

# K_sp = a(Ba++) * a(SO4--) / a(Barite) = a(Ba++) * a(SO4--)
# IAP = a_actual(Ba++) * a_actual(SO4--)
# SI = log(a_actual(Ba++) * a_actual(SO4--) / a(Ba++) / a(SO4--))
# 𝑎𝑖 = 𝛾𝑖 * 𝑚𝑖 - actualy activity? 

# Then, we define the list of names for the DataFrame columns. Note, that they must correspond
# to the order of the properties defined in the `output_quantities` list:

column_quantities = """
    pH
    Hcation
    Clanion
    SO4anion
    Bacation
    Cacation
    Srcation
    Nacation
    Barite_phase_amount
    Bacation_activity
    SO4anion_activity
    Bacation_activity_cofficient
    SO4anion_activity_cofficient
""".split()

# column_quantities = """
#     pH
#     Hcation
#     Clanion
#     SO4anion
#     Bacation
#     Cacation
#     Srcation
#     Nacation
#     Barite_phase_amount
#     Anhydrite_phase_amount
#     Celestite_phase_amount
#     Bacation_activity
#     SO4anion_activity
#     Bacation_activity_cofficient
#     SO4anion_activity_cofficient
# """.split()

# Create the list of columns stored in dataframes
columns = ['step', 'x'] + column_quantities
import pandas as pd

# Initialize dataframes with above defined columns
df = pd.DataFrame(columns=columns)

# ## Organization of the program
#
folder_results = 'results-rt-scaling'
def make_results_folders():
    os.system('mkdir -p ' + folder_results)

# ## Performing the reactive transport simulation

# Subsections below correspond to the methods responsible for each of the functional parts of `simulate()` method.
#
# ### Construction of the chemical system with its phases and species
#

def define_chemical_system():

    # Construct the chemical system with its phases and species
    db = Database('supcrt07.xml')
    
    dhModel = DebyeHuckelParams()
    dhModel.setPHREEQC()

    editor = ChemicalEditor(db)
    editor.addAqueousPhaseWithElements("H Cl S O Ba Ca Sr Na K Mg C Si").\
        setChemicalModelDebyeHuckel(dhModel)

    editor.addMineralPhase('Barite')
    #editor.addMineralPhase('Witherite')
    #editor.addMineralPhase('Anhydrite')
    #editor.addMineralPhase('Celestite')

    system = ChemicalSystem(editor)
    #print(system)

    return system

# ### Initial condition (IC) of the reactive transport problem

def define_initial_condition_fw(system):

    # Formation water at equilbrium:
    # contain bivalent cations in relative abundance 
    # little sulfate
    # the Miller analysis:
    # Na+ = 27250 mg/kg
    # K+ = 1730 mg/kg
    # Mg++ = 110 mg/kg
    # Ca++ = 995 mg/kg
    # Sr++ = 105 mg/kg
    # Ba++ = 995 mg/kg
    # Cl- = 45150 mg/kg
    # HCO3- = 1980 mg/kg
    # SO4-- = 10 ug/kg
    #
    #  FW → high Ba2+ and low SO42- concentration

    problem_ic = EquilibriumProblem(system)
    problem_ic.setTemperature(T, "celsius")
    problem_ic.setPressure(P, "atm")
    problem_ic.add("H2O", water_kg, "kg")
    problem_ic.add("SO4--", 1e-6 * water_kg, "kg") # SO4-- = 10 ug/kg
    problem_ic.add("Ca++", 0.995 * water_kg, "kg") # Ca++ = 995 mg/kg
    problem_ic.add("Ba++", 0.995 * water_kg, "kg") # Ba++ = 995 mg/kg
    problem_ic.add("Sr++", 0.105 * water_kg, "kg") # Sr++ = 105 mg/kg
    problem_ic.add("Na+", 27.250 * water_kg, "kg") # Na+ = 27250 mg/kg
    problem_ic.add("K+", 1.730 * water_kg, "kg") # K+ = 1730 mg/kg
    problem_ic.add("Mg++", 0.110 * water_kg, "kg") # Mg++ = 110 mg/kg
    problem_ic.add("Cl-", 45.150 * water_kg, "kg") # Cl- = 45150 mg/kg
    problem_ic.add("HCO3-", 1.980 * water_kg, "kg") # HCO3- = 1980 mg/kg

    # Calculate the equilibrium states for the initial conditions
    state_ic = equilibrate(problem_ic)

    # Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3') # 10% of porosity
    state_ic.scaleVolume(1.0, 'm3')

    props = state_ic.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)

    #print("state_ic = \n", state_ic)
    #print("state_ic = \n", state_ic)
    print("ph(FW) = ", pH.val)

    return state_ic

def define_initial_condition_fw_phreeqc(system):

    # Formation water at equilbrium:
    # contain bivalent cations in relative abundance
    # little sulfate
    #
    # Not correct yet!

    problem_ic = EquilibriumProblem(system)
    problem_ic.setTemperature(T, "celsius")
    problem_ic.setPressure(P, "atm")
    problem_ic.add("H2O", water_kg, "kg")
    problem_ic.add("SO4--", 1e-6 * water_kg, "kg") # SO4-- = 10 ug/kg
    problem_ic.add("Ca++", 0.995 * water_kg, "kg") # Ca++ = 995 mg/kg
    problem_ic.add("Ba++", 0.995 * water_kg, "kg") # Ba++ = 995 mg/kg
    problem_ic.add("Sr++", 0.105 * water_kg, "kg") # Sr++ = 105 mg/kg
    problem_ic.add("Na+", 27.250 * water_kg, "kg") # Na+ = 27250 mg/kg
    problem_ic.add("K+", 1.730 * water_kg, "kg") # K+ = 1730 mg/kg
    problem_ic.add("Mg++", 0.110 * water_kg, "kg") # Mg++ = 110 mg/kg
    problem_ic.add("Cl-", 45.150 * water_kg, "kg") # Cl- = 45150 mg/kg
    problem_ic.add("HCO3-", 1.980 * water_kg, "kg") # HCO3- = 1980 mg/kg

    # Calculate the equilibrium states for the initial conditions
    state_ic = equilibrate(problem_ic)

    # Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3') # 10% of porosity
    state_ic.scaleVolume(1.0, 'm3')

    #print("state_ic = \n", state_ic)

    props = state_ic.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(FW, PHREEQC) = ", pH.val)

    return state_ic

# ### Boundary condition (BC) of the reactive transport problem

def define_boundary_condition_cb(system):

    # Define the boundary condition of the reactive transport modeling problem
    problem_bc = EquilibriumProblem(system)
    problem_bc.setTemperature(T, "celsius")
    problem_bc.setPressure(P, "atm")
    problem_bc.add("H2O", water_kg, "kg")
    problem_bc.add("NaCl", 7, "mol")
    
    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    #print("state_bc_cb = \n", state_bc)

    props = state_bc.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(CB) = ", pH.val)

    return state_bc

def define_boundary_condition_sw(system):

    # Seawater: 
    # rich in sulfate > 2500 mg / kg
    # poor in Ca++ and
    # nearly depleted in Sr++ and Ba++
    # SW → low Ba2+ and high SO42- concentration

    # Seewater from Bethke, Table 30.1:
    # Cl - = 19350 mg / kg
    # Ca ++ = 411 mg / kg
    # Mg ++ = 1290 mg / kg
    # Na + = 10760 mg / kg
    # K + = 399 mg / kg
    # SO4 -- = 2710 mg / kg
    # HCO3 - = 142 mg / kg
    # SiO2(aq) = 6 mg / kg
    # Ca ++ = 411 mg / kg
    # Sr ++ = 8 mg / kg
    # Ba ++ = 0.01 mg / kg

    problem_bc = EquilibriumInverseProblem(system)
    problem_bc.setTemperature(T, "celsius")
    problem_bc.setPressure(P, "atm")
    problem_bc.add("H2O", water_kg, "kg")
    problem_bc.add("SO4--", 2.710 * water_kg, "kg") # 2710 mg / kg = 2.710 kg / kg * 58 kg = 157.18 kg
    problem_bc.add("Ca++", 0.411 * water_kg, "kg")  # 411 mg / kg = 0.411 kg / kg * 58 kg = 23.838 kg
    problem_bc.add("Ba++", 0.00001 * water_kg, "kg")    # 0.01 mg / kg = 0.00001 kg / kg * 58 kg = 0.00058 kg
    problem_bc.add("Sr++", 0.008 * water_kg, "kg")   # 8 mg / kg = 0.008 kg / kg * 58 kg = 0.464 kg
    problem_bc.add("Na+", 10.760 * water_kg, "kg") # 10760  mg / kg = 10.760 kg / kg * 58 kg = 624.08 kg
    problem_bc.add("K+", 0.399 * water_kg, "kg") # 399 mg / kg = 0.399 kg / kg * 58 kg = 23.142 kg
    problem_bc.add("Mg++", 1.29 * water_kg, "kg") # 1290 mg / kg = 1.29 kg / kg * 58 kg = 74.82 kg
    problem_bc.add("Cl-", 19.350 * water_kg, "kg") # 19350 mg / kg = 19.350 kg / kg * 58 kg = 1122.3 kg
    problem_bc.add("HCO3-", 0.142 * water_kg, "kg") # 142 mg / kg = 0.142 kg / kg * 58 kg = 8.236 kg
    #problem_bc.pH(8.1)
    problem_bc.pH(8.1, "NaCl(aq)")
    #problem_bc.pH(8.1, "HCl")
    #problem_bc.pH(8.1, "CO2")
    #problem_bc.pH(8.1, "HCl", "NaOH")

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    #print("state_bc_sw = \n", state_bc)

    props = state_bc.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(SW) = ", pH.val)

    return state_bc

def define_boundary_condition_sw_phreeqc(system):

    # Seawater:
    # rich in sulfate > 2500 mg / kg
    # poor in Ca++ and
    # nearly depleted in Sr++ and Ba++
    # SW → low Ba2+ and high SO42- concentration

    # Seawater:
    # Concentration is in ppm (parts per million)
    # Assuming that 1 ppm = 1 mg/kg
    #         pH      8.22
    #         pe      8.451
    #         density 1.023
    #         temp    25.0
    #         Ca              412.3
    #         Mg              1291.8
    #         Na              10768.0
    #         K               399.1
    #         Si              4.28
    #         Cl              19353.0
    #         Alkalinity      141.682 as HCO3
    #         SO4--           2712.0
    problem_bc = EquilibriumInverseProblem(system)
    problem_bc.setTemperature(T, "celsius")
    problem_bc.setPressure(P, "atm")
    problem_bc.add("H2O", water_kg, "kg")
    problem_bc.add("Ca++", 0.4123 * water_kg, "kg")  # 412.3 mg / kg = 0.4123 kg / kg => 0.4123 * 58 = 23.9134
    problem_bc.add("Mg++", 1.2918 * water_kg, "kg")  # 1291.8 mg / kg = 1.2918 kg / kg => 1.2918 * 58 = 74.9244
    problem_bc.add("Na+", 10.768 * water_kg, "kg")  # 10768.0  mg / kg = 10.768 kg / kg => 10.768 * 58 = 624.544
    problem_bc.add("K+", 0.3991 * water_kg, "kg")  # 399.1 mg / kg = 0.3991 kg / kg => 0.3991 * 58 = 23.1478
    problem_bc.add("Si", 0.00428 * water_kg, "kg")  # 4.28 mg / kg = 0.00428 kg / kg => 0.00428 * 58 = 0.24824
    problem_bc.add("Cl-", 19.353 * water_kg, "kg")  # 19353.0 mg / kg = 19.353 kg / kg => 19.353 * 58 = 1122.474
    problem_bc.add("HCO3-", 0.142682 * water_kg, "kg")  # 141.682 mg / kg = 0.142682 kg / kg => 0.142682 * 58 = 8.275556
    problem_bc.add("SO4--", 2.712 * water_kg, "kg") # 2712.0 mg / kg = 2.712 kg / kg => 2.712 * 58 = 157.296
    #problem_bc.pH(8.22)
    #problem_bc.pE(8.451)
    problem_bc.pH(8.22, "NaCl(aq)")
    #problem_bc.pH(8.1, "CO2")
    #problem_bc.pH(8.1, "HCl", "NaOH")

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    #print("state_bc_sw = \n", state_bc)

    props = state_bc.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(SW, PHREEQC) = ", pH.val)

    return state_bc

# ### Indices of partitioning fluid and solid species

def partition_indices(system):
    nelems = system.numElements()

    ifluid_species = system.indicesFluidSpecies()
    isolid_species = system.indicesSolidSpecies()

    return nelems, ifluid_species, isolid_species


# ### Partitioning fluid and solid species

def partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc_cb, state_bc_sw):
    # The concentrations of each element in each mesh cell (in the current time step)
    b = np.zeros((ncells, nelems))
    # Initialize the concentrations (mol/m3) of the elements in each mesh cell
    b[:] = state_ic.elementAmounts()

    # The concentrations (mol/m3) of each element in the fluid partition, in each mesh cell
    bfluid = np.zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the solid partition, in each mesh cell
    bsolid = np.zeros((ncells, nelems))

    # Initialize the concentrations (mol/m3) of each element on the boundary, while injecting completion brine
    b_bc_cb = state_bc_cb.elementAmounts()

    # Initialize the concentrations (mol/m3) of each element on the boundary, while injecting completion brine
    b_bc_sw = state_bc_sw.elementAmounts()

    return b, bfluid, bsolid, b_bc_cb, b_bc_sw


# ### Reactive transport cycle
#
# #### Transport

def transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species):
    # Collect the amounts of elements from fluid and solid partitions
    for icell in range(ncells):
        bfluid[icell] = states[icell].elementAmountsInSpecies(ifluid_species)
        bsolid[icell] = states[icell].elementAmountsInSpecies(isolid_species)

    # Get the porosity of the boundary cell
    bc_cell = 0
    phi_bc = states[bc_cell].properties().fluidVolume().val / states[bc_cell].properties().volume().val

    # Transport each element in the fluid phase
    for j in range(nelems):
        transport_fullimplicit(bfluid[:, j], dt, dx, v, D, phi_bc * b_bc[j])

    # Update the amounts of elements in both fluid and solid partitions
    b[:] = bsolid + bfluid

    return bfluid, bsolid, b


# ##### Transport calculation with finite-volume scheme

def transport_fullimplicit(u, dt, dx, v, D, ul):
    # Number of DOFs
    n = len(u)
    alpha = D * dt / dx ** 2
    beta = v * dt / dx

    # Upwind finite volume scheme
    a = np.full(n, -beta - alpha)
    b = np.full(n, 1 + beta + 2 * alpha)
    c = np.full(n, -alpha)

    # Set the boundary condition on the left cell
    if dirichlet:
        # Use Dirichlet BC boundary conditions
        b[0] = 1.0
        c[0] = 0.0
        u[0] = ul

    else:
        # Flux boundary conditions (implicit scheme for the advection)
        # Left boundary
        b[0] = 1 + alpha + beta
        c[0] = -alpha  # stays the same as it is defined -alpha
        u[0] += beta * ul  # = dt/dx * v * g, flux that we prescribe is equal v * ul

    # Right boundary is free
    a[-1] = - beta
    b[-1] = 1 + beta

    # Solve a tridiagonal matrix equation
    thomas(a, b, c, u)


# ##### Solving the system of equations obtained from finite volume discretization

def thomas(a, b, c, d):
    n = len(d)
    c[0] /= b[0]
    for i in range(1, n - 1):
        c[i] /= b[i] - a[i] * c[i - 1]
    d[0] /= b[0]
    for i in range(1, n):
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1])
    x = d
    for i in reversed(range(0, n - 1)):
        x[i] -= c[i] * x[i + 1]
    return x


# #### Reactive chemistry

def reactive_chemistry(solver, states, b):
    # Equilibrating all cells with the updated element amounts
    for icell in range(ncells):
        solver.solve(states[icell], T, P, b[icell])
    return states


# ### Results saving and analyzing
#
# Function `outputstate_df` is the auxiliary function to add data to the DataFrame at each time step.

def outputstate_df(step, system, states):
    # Define the instance of ChemicalQuantity class
    quantity = ChemicalQuantity(system)

    # Create the list with empty values to populate with chemical properties
    values = [None] * len(columns)
    for state, x in zip(states, xcells):

        # Populate values with number of reactive transport step and spacial coordinates
        values[0] = step
        values[1] = x

        # Update the
        quantity.update(state)
        for quantity_name, i in zip(output_quantities, range(2, len(states))):
            values[i] = quantity.value(quantity_name) * (100 / (1 - phi) if "phaseVolume" in quantity_name else 1)
        df.loc[len(df)] = values

# ### Plotting of the obtained results

def titlestr(t):
    t = t / minute  # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: %2dh %2dm' % (h, m)

# Routines `plot_figures_ph()`, `plot_figures_pyrrhotite_siderite_volume()`, `plot_figures_pyrrhotite_siderite_amount()`,
# and 'plot_figures_aqueous_species()' are dedicated to drawing the plots with chemical properties on the selected steps
# that are specified by the user below.

def plot_figures_ph(steps):
    # Plot ph on the selected steps
    plots = []
    for i in steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='pH', color='teal', line_width=2, legend_label='pH', source=source)
        p.x_range = Range1d(xl, xr-1)
        #p.y_range = Range1d(6.0, 26.0)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'pH'
        p.legend.location = 'bottom_right'
        p.title.text = titlestr(t)

        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figures_barite_phase_amount(steps):
    plots = []
    for i in steps:
        print("On barite figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='Barite_phase_amount', color='blue', line_width=2, legend_label='Barite',
               muted_color='blue', muted_alpha=0.2, source=source)
        p.x_range = Range1d(xl, xr-1)
        #p.y_range = Range1d(-0.001, 60.0)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Mineral Phase Amount [mol]'
        p.legend.location = 'center_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figures_aqueous_species(steps):
    plots = []
    for i in steps:
        print("On aqueous-species figure at time step: {}".format(i))
        source = ColumnDataSource(df[df['step'] == i])
        t = dt * i

        p = figure(plot_width=600, plot_height=300, y_axis_type = 'log',)
        #p.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
        p.line(x='x', y='Clanion', color='darkcyan', line_width=2, legend_label='Cl-', source=source)
        p.line(x='x', y='SO4anion', color='darkorange', line_width=2, legend_label='SO4--', source=source)
        p.line(x='x', y='Bacation', color='seagreen', line_width=2, legend_label='Ba++', source=source)
        #p.line(x='x', y='Cacation', color='indianred', line_width=2, legend_label='Ca++', source=source)
        #p.line(x='x', y='Srcation', color='darkblue', line_width=2, legend_label='Sr++', source=source)
        #p.line(x='x', y='Nacation', color='blue', line_width=2, legend_label='Na+', source=source)
        p.x_range = Range1d(xl, xr-1)
        #p.y_range = Range1d(1e-12, 1e-1)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Concentration [molal]'
        p.legend.location = 'top_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figures_aqueous_species(steps):
    plots = []
    for i in steps:
        print("On aqueous-species figure at time step: {}".format(i))
        source = ColumnDataSource(df[df['step'] == i])
        t = dt * i

        p = figure(plot_width=600, plot_height=300, y_axis_type = 'log',)
        #p.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
        p.line(x='x', y='Clanion', color='darkcyan', line_width=2, legend_label='Cl-', source=source)
        p.line(x='x', y='SO4anion', color='darkorange', line_width=2, legend_label='SO4--', source=source)
        p.line(x='x', y='Bacation', color='seagreen', line_width=2, legend_label='Ba++', source=source)
        #p.line(x='x', y='Cacation', color='indianred', line_width=2, legend_label='Ca++', source=source)
        #p.line(x='x', y='Srcation', color='darkblue', line_width=2, legend_label='Sr++', source=source)
        #p.line(x='x', y='Nacation', color='blue', line_width=2, legend_label='Na+', source=source)
        p.x_range = Range1d(xl, xr)
        #p.y_range = Range1d(1e-12, 1e-1)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Concentration [molal]'
        p.legend.location = 'top_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

# # Main parts of the tutorial
#
# First, we create folders for the results:

make_results_folders()

# # Run the reactive transport simulations:

# Construct the chemical system with its phases and species
system = define_chemical_system()

# Define the initial condition of the reactive transport modeling problem
state_ic = define_initial_condition_fw(system)

# # Define the boundary condition of the reactive transport modeling problem composed of two different stages

# Define the completion brine (CB)
state_bc_cb = define_boundary_condition_cb(system)

# Define the seawater (SW)
state_bc_sw = define_boundary_condition_sw(system)

# Define the seawater (SW)
state_bc_sw = define_boundary_condition_sw_phreeqc(system)

# Generate indices of partitioning fluid and solid species
nelems, ifluid_species, isolid_species = partition_indices(system)

# Partitioning fluid and solid species
b, bfluid, bsolid, b_bc_cb, b_bc_sw \
    = partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc_cb, state_bc_sw)

# Create a list of chemical states for the mesh cells (one for each cell, initialized to state_ic)
states = [state_ic.clone() for _ in range(ncells + 1)]
print(len(states))

# Create the equilibrium solver object for the repeated equilibrium calculation
solver = EquilibriumSolver(system)

# Running the reactive transport simulation loop
step = 0  # the current step number
t = 0.0  # the current time (in seconds)

# Output the initial state of the reactive transport calculation
outputstate_df(step, system, states)

df.shape

with tqdm(total=nsteps_cb, desc="45 hours of completion brine (CB) injection") as pbar:
    while step < nsteps_cb:
        # Perform transport calculations
        bfluid, bsolid, b = transport(states, bfluid, bsolid, b, b_bc_cb, nelems, ifluid_species, isolid_species)

        # Perform reactive chemical calculations
        states = reactive_chemistry(solver, states, b)

        # Increment time step and number of time steps
        t += dt
        step += 1

        # Output the current state of the reactive transport calculation
        outputstate_df(step, system, states)

        # Update a progress bar
        pbar.update(1)

# (rows, columns) = df.shape

print(f"time: {t / hour} hours")


with tqdm(total=nsteps_sw, desc="855 hours of seawater (SW) injection") as pbar:
    while step < nsteps_sw + nsteps_cb:
        # Perform transport calculations
        bfluid, bsolid, b = transport(states, bfluid, bsolid, b, b_bc_sw, nelems, ifluid_species, isolid_species)

        # Perform reactive chemical calculations
        states = reactive_chemistry(solver, states, b)

        # Increment time step and number of time steps
        t += dt
        step += 1

        # Output the current state of the reactive transport calculation
        outputstate_df(step, system, states)

        # Update a progress bar
        pbar.update(1)

# To inspect the collected data, one can run:

df.shape

df

# To save the results in csv-format, please execute:

df.to_csv(folder_results + '/rt.scaling.csv', index=False)

# Select the steps, on which results must plotted:

selected_steps_to_plot = [20, 45, 46, 60, 120, 260, 300]
assert all(step <= nsteps for step in selected_steps_to_plot), f"Make sure that selceted steps are less than " \
                                                               f"total amount of steps {nsteps}"

# Outputting the plots to the notebook requires the call of `output_notebook()` that specifies outputting the plot
# inline in the Jupyter notebook:

output_notebook()

# Plot ph on the selected steps:
plot_figures_ph(selected_steps_to_plot)

# Plot calcite and dolomite on the selected steps:
plot_figures_barite_phase_amount(selected_steps_to_plot)

# Plot aqueous species on the selected steps:
plot_figures_aqueous_species(selected_steps_to_plot)


step = 10

# The data streaming is looped, i.e., we will return to the initial time step when reaching the end of the reactive
# transport simulations.

def modify_doc(doc):
    # Initialize the data by the initial chemical state
    source = ColumnDataSource(df[df['step'] == 0])

    # Auxiliary function that returns a string for the title of a figure in the format Time: #h##m
    def titlestr(t):
        t = t / minute  # Convert from seconds to minutes
        h = int(t) / 60  # The number of hours
        m = int(t) % 60  # The number of remaining minutes
        return 'Time: %2dh %2dm' % (h, m)

    # Plot for ph
    p1 = figure(plot_width=600, plot_height=250)
    p1.line(x='x', y='pH', color='teal', line_width=2, legend_label='pH', source=source)
    p1.x_range = Range1d(xl, xr-1)
    #p1.y_range = Range1d(6.0, 26.0)
    p1.xaxis.axis_label = 'Distance [m]'
    p1.yaxis.axis_label = 'pH'
    p1.legend.location = 'bottom_right'
    p1.title.text = titlestr(0 * dt)

    # Plot for calcite and dolomite
    p2 = figure(plot_width=600, plot_height=250)
    p2.line(x='x', y='Barite_phase_amount', color='blue', line_width=2,
            legend_label='Barite', muted_color='blue', muted_alpha=0.2,
            source=source)
    p2.x_range = Range1d(xl, xr-1)
    p2.y_range = Range1d(-0.001, 60.0)
    p2.xaxis.axis_label = 'Distance [m]'
    p2.yaxis.axis_label = 'Phase Amount [mol]'
    p2.legend.location = 'center_right'
    p2.title.text = titlestr(0 * dt)
    p2.legend.click_policy = 'mute'

    p3 = figure(plot_width=600, plot_height=300, y_axis_type='log')
    #p3.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
    p3.line(x='x', y='Clanion', color='darkcyan', line_width=2, legend_label='Cl-', source=source)
    p3.line(x='x', y='SO4anion', color='darkorange', line_width=2, legend_label='SO4--', source=source)
    p3.line(x='x', y='Bacation', color='seagreen', line_width=2, legend_label='Ba++', source=source)
    #p3.line(x='x', y='Cacation', color='indianred', line_width=2, legend_label='Ca++', source=source)
    #p3.line(x='x', y='Srcation', color='darkblue', line_width=2, legend_label='Sr++', source=source)
    #p3.line(x='x', y='Nacation', color='blue', line_width=2, legend_label='Na+', source=source)
    p3.x_range = Range1d(xl, xr-1)
    p3.xaxis.axis_label = 'Distance [m]'
    p3.yaxis.axis_label = 'Concentration [molal]'
    p3.legend.location = 'top_right'
    p3.title.text = titlestr(0 * dt)
    p3.legend.click_policy = 'mute'

    layout = column(p1, p2, p3)

    # Function that return the data dictionary with provided index of the file
    def update():

        if source.data['step'][0] + 1 <= nsteps:
            step_number = source.data['step'][0] + step
        else:
            step_number = 0

        new_source = ColumnDataSource(df[df['step'] == step_number])
        new_data = dict(index=np.linspace(0, ncells, ncells + 1, dtype=int),
                        step=new_source.data['step'],
                        x=new_source.data['x'],
                        pH=new_source.data['pH'],
                        Barite_phase_amount=new_source.data['Barite_phase_amount'],
                        Hcation=new_source.data['Hcation'],
                        SO4anion=new_source.data['SO4anion'],
                        Clanion=new_source.data['Clanion'],
                        Bacation=new_source.data['Bacation'],
                        Cacation=new_source.data['Cacation'],
                        Srcation=new_source.data['Srcation'],
                        Nacation=new_source.data['Nacation'],
                        Bacation_activity=new_source.data['Bacation_activity'],
                        SO4anion_activity=new_source.data['SO4anion_activity'],
                        Bacation_activity_cofficient=new_source.data['Bacation_activity_cofficient'],
                        SO4anion_activity_cofficient=new_source.data['SO4anion_activity_cofficient'])

        p1.title.text = titlestr(step_number * dt)
        p2.title.text = titlestr(step_number * dt)
        p3.title.text = titlestr(step_number * dt)

        source.stream(new_data, rollover=ncells + 1)

    doc.add_periodic_callback(update, 500)
    doc.add_root(layout)

# Outputting the plots to the notebook requires the call of `output_notebook()` that specifies outputting the plot
# inline in the Jupyter notebook. Finally, the function `modify_doc()` must be passed to `show`, so that the app defined
# by it is displayed inline.
#
# > **Important:** If you run this tutorial in the *localhost*, make sure that number provided to the variable
# `notebook_url` below coincides with the number of the localhost you have in your browser.
#
# In the app below, we refresh the reactive time step in a loop, which automatically updates the data source for the
# plots for ph, volume phases of calcite and dolomite, and mollalities of aqueous species (in logarithmic scale).

output_notebook()
show(modify_doc, notebook_url="http://localhost:8888")
