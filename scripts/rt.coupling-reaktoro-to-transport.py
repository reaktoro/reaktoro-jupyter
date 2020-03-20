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

# # Coupling Reaktoro into other reactive transport codes
#
# In this tutorial, we show how Reaktoro can be used in other codes for reactive transport modeling. Here,
# you will find that we have split the mass transport and chemical reaction calculations so that you can see
# how a dedicated and advanced transport solver could be combined with Reaktoro's solvers for chemical reaction
# calculations. A basic transport solver is used here for the sake of the software coupling demonstration.
#
# The reactive transport problem we consider here is also relatively simple. We proceed with a step-by-step
# explanation.
#
# ## Importing python packages
#
# First, we need to import a few Python packages to enable us to perform the numerical calculations and plotting.

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
from bokeh.models.widgets import Slider

# -

# We import the **reaktoro** Python package so that we can use its classes and methods for performing chemical
# reaction calculations, **numpy** for working with arrays, **tqdm** for the progress bar functionality and **os**,
# to provide a portable way of using operating system dependent functionality. For plotting capabilities of obtained
# results, we use **bokeh** library.
#
# > **Note**: To make sure that all the widgets are working correctly, make sure to run:
# > `$ jupyter nbextension enable --py widgetsnbextension` and
# > `$jupyter labextension install @jupyter-widgets/jupyterlab-manager`.
#
# ## Initializing auxiliary time-related constants
# In this step, we initialize auxiliary time-related constants from seconds up to years used in the rest of the code.

second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# ## Defining parameters for the reactive transport simulation
#
# Next, we define reactive transport and numerical discretization parameters. In particular, we specify the considered
# rock domain by setting coordinates of its left and right boundaries to 0.0 m and 100.0 m, respectively. The
# discretization parameters, i.e., the number of cells and steps in time, are both set to 100. The reactive
# transport modeling procedure assumes a constant fluid velocity of 1 m/day (1.16 · 10<sup>-5</sup> m/s) and the same
# diffusion coefficient of 10<sup>-9</sup> m<sup>2</sup>/s for all fluid species (without dispersivity). The
# size of the time-step is set to 30 minutes. Temperature and pressure are set to 60 &deg;C and 100 bar, respectively,
# throughout the whole tutorial. 

# +
# Discretization parameters
xl = 0.0  # x-coordinate of the left boundary
xr = 1.0  # x-coordinate of the right boundary
ncells = 100  # number of cells in the discretization
nsteps = 100  # number of steps in the reactive transport simulation
dx = (xr - xl) / ncells  # length of the mesh cells (in units of m)
dt = 10 * minute  # time step

# Physical parameters
D = 1.0e-9  # diffusion coefficient (in units of m2/s)
v = 1.0 / day  # fluid pore velocity (in units of m/s)
T = 60.0 + 273.15  # temperature (in units of K)
P = 100 * 1e5  # pressure (in units of Pa)
phi = 0.1  # the porosity
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
# Before running the reactive transport simulations, we specify the list of parameters we are interested in
# outputting. In this case, it is pH, molality of `H+`, `Ca++`, `Mg++`, `HCO3-`, `CO2(aq)`, as well as a phase volume
# of calcite and dolomite.

output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Ca++)
    speciesMolality(Mg++)
    speciesMolality(HCO3-)
    speciesMolality(CO2(aq))
    phaseVolume(Calcite)
    phaseVolume(Dolomite)
""".split()

# Then, we define the list of name for the DataFrame columns. Note, that they must correspond
# to the order of the properties defined in the `output_quantities` list:

column_quantities = """
    pH
    Hcation
    Cacation
    Mgcation
    HCO3anion
    CO2aq
    calcite
    dolomite
""".split()

# We create a DataFrame storage `df` using python library **pandas**. The columns in it match include the number of
# reactive transport step, the space coordinates of the discretization cells, and the name of the properties we
# assigned to be output (see `output_quantities` list above).

# Create the list of columns stored in dataframes
columns = ['step', 'x'] + column_quantities
import pandas as pd

# Initialize dataframes with above defined columns
df = pd.DataFrame(columns=columns)

# ## Organization of the program
#
# The main part of the program (at the bottom of this tutorial) consists of three parts, each represented by a Python
# function and documented in the following sections:
# * creation of folders for the results (function `make_results_folders()`),
# * simulation of reactive transport problem (method `simulate()`), and
# * plotting of the obtained results.
#
# ## Creating folders for the outputted results
# Using **os** package, we create required folders for outputting the obtained results and for the plot and video
# files later.

folder_results = 'results-rt-coupling'
def make_results_folders():
    os.system('mkdir -p ' + folder_results)


# ## Performing the reactive transport simulation
# The reactive transport simulation is performed in the function `simulate()`, which consists of several building
# blocks:
# * initialization of the reactive transport problem and
# * performing the reactive transport simulation along defined time interval.
#
# The *preparatory initialization* step consists of the following sub-steps:
# * definition of chemical system with its phases and species using `define_chemical_system()`,
# * definition of the initial condition of the reactive transport problem in `define_initial_condition()`,
# * definition of the boundary condition of the reactive transport problem in `define_initial_condition()`,
# * generation of auxiliary indices to partition elements using `partition_indices()` and elements' partitioning 
# corresponding to fluid and solid species with function `partition_elements_in_mesh_cell()`, and finally
# * definition of instance of [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html).
#
# The *simulation of the reactive transport problem* is represented by the loop over discretized time interval until
# the final time is reached. On each step of this loop, the following functionality of performed:
# * transport calculations using method `transport()`,
# * reactive chemical calculations with `reactive_chemistry()` function, and
# * saving the obtained results by means of `outputstate()`.
#
# Performing of the transport and reactive chemistry sequentially is possible due to the *operator splitting
# procedure*, in which we first update the amounts of elements `b`. These updated amounts of
# elements in each cell are used to evaluate its new chemical equilibrium state, thus producing new amounts of the
# species in both the fluid and solid phases (available in the list `states`, the instances of
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) classes). This chemical reaction
# equilibrium calculation step, at each mesh cell, permits, for example, aqueous species and minerals to react,
# and thus causes mineral dissolution or precipitation, depending on how much the amount of mineral species changes.
# This can then be used, for example, to compute new porosity value for the cell.
#
# The structure of the function `simulate`, performing the simulation of reactive transport problem, can be inspected
# below:

def simulate():
    # Construct the chemical system with its phases and species
    system = define_chemical_system()

    # Define the initial condition of the reactive transport modeling problem
    state_ic = define_initial_condition(system)

    # Define the boundary condition of the reactive transport modeling problem
    state_bc = define_boundary_condition(system)

    # Generate indices of partitioning fluid and solid species
    nelems, ifluid_species, isolid_species = partition_indices(system)

    # Partitioning fluid and solid species
    b, bfluid, bsolid, b_bc = partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc)

    # Create a list of chemical states for the mesh cells (one for each cell, initialized to state_ic)
    states = [state_ic.clone() for _ in range(ncells + 1)]

    # Create the equilibrium solver object for the repeated equilibrium calculation
    solver = EquilibriumSolver(system)

    # Running the reactive transport simulation loop
    step = 0  # the current step number
    t = 0.0  # the current time (in seconds)
    # Save the initial state of the reactive transport calculation to the dataframe
    outputstate_df(step, system, states)

    with tqdm(total=nsteps, desc="Reactive transport simulations") as pbar:
        while step <= nsteps:
            # Perform transport calculations
            bfluid, bsolid, b = transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species)

            # Perform reactive chemical calculations
            states = reactive_chemistry(solver, states, b)

            # Increment time step and number of time steps
            t += dt
            step += 1

            # Save the current state of the reactive transport calculation to the dataframe
            outputstate_df(step, system, states)

            # Update a progress bar
            pbar.update(1)

    print("Finished!")


# Subsections below correspond to the methods responsible for each of the functional parts of `simulate()` method.
#
# ### Construction of the chemical system with its phases and species
# Method `define_chemical_system()` begins by defining the chemical system using [ChemicalEditor](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class. In particular, we specify aqueous and mineral
# phases that should be considered in the chemical system. For performance reasons, the aqueous phase is defined by
# manually specifying the chemical species, instead of automatic collection of species from the database. There are
# three pure mineral phase considered: quartz (SiO<sub>2</sub>), calcite (CaCO<sub>3</sub>), and dolomite
# (CaMg(CO<sub>3</sub>)<sub>2</sub>). Finally, using class
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html),
# we define the chemical system.

def define_chemical_system():
    # Define thermodynamical database to be used for chemical simulations
    db = Database('supcrt98.xml')

    # Define chemical editor to work with the chemical system
    editor = ChemicalEditor(db)

    # Define phases and corresponding to these phases species
    editor.addAqueousPhaseWithElements('H O Na Cl Ca Mg C Si Ca')
    editor.addMineralPhase('Quartz')
    editor.addMineralPhase('Calcite')
    editor.addMineralPhase('Dolomite')

    # Create the chemical system using the configured editor
    system = ChemicalSystem(editor)

    return system


# ### Initial condition of the reactive transport problem
# After constructing the chemical system of interest, we can proceed to the definition of a chemical equilibrium
# problem to set the *initial condition* for the fluid and rock composition that is consistent with our intention of
# modeling reactive transport in a porous rock column composed of quartz and calcite, at 60 &deg;C and 100 bar,
# with a 0.7 NaCl molal brine in equilibrium with the rock minerals. To ensure equilibrium with both quartz and
# calcite, the equilibrium problem setup below considers a relatively large amount for these minerals (10 moles each).
#
# Next, we use method [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e)
# to calculate the chemical equilibrium state of the system with the given initial conditions stored in the objects
# `problem_ic`. The numerical solution of each problem results in the objects `state_ic` of class
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html), which stores the temperature,
# pressure, and the amounts of every species in the system.
#
# The function ends with scaling the volumes of the aqueous and mineral phases so that they are consistent with a 10 %
# porosity and the required volume percentages of the rock minerals (98 %<sub>vol</sub> of quartz and 2 %<sub>vol</sub>
# of calcite).

def define_initial_condition(system):
    problem_ic = EquilibriumProblem(system)
    problem_ic.setTemperature(T)
    problem_ic.setPressure(P)
    problem_ic.add('H2O', 1.0, 'kg')
    problem_ic.add('NaCl', 0.7, 'mol')
    problem_ic.add('CaCO3', 10, 'mol')
    problem_ic.add('SiO2', 10, 'mol')

    # Calculate the equilibrium states for the initial conditions
    state_ic = equilibrate(problem_ic)

    # Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')
    state_ic.scalePhaseVolume('Quartz', 0.882, 'm3')
    state_ic.scalePhaseVolume('Calcite', 0.018, 'm3')

    return state_ic


# ### Boundary condition of the reactive transport problem
# For the boundary condition, we need to specify the composition of the fluid that is injected into the rock. This is
# done below, by defining an equilibrium problem that will later be solved to produce an aqueous fluid with 0.9 molal
# NaCl, 0.05 molal MgCl<sub>2</sub>, 0.01 CaCl<sub>2</sub>, and 0.75 molal CO<sub>2</sub>, in a state
# very close to CO<sub>2</sub> saturation.
#
# After equilibration, the obtained chemical state representing the boundary condition for the injected fluid
# composition, we scale its volume to 1 m<sup>3</sup>. This is done so that the amounts of the species in the fluid are
# consistent with a mol/m<sup>3</sup> scale.

def define_boundary_condition(system):
    # Define the boundary condition of the reactive transport modeling problem
    problem_bc = EquilibriumProblem(system)
    problem_bc.setTemperature(T)
    problem_bc.setPressure(P)
    problem_bc.add('H2O', 1.0, 'kg')
    problem_bc.add('NaCl', 0.90, 'mol')
    problem_bc.add('MgCl2', 0.05, 'mol')
    problem_bc.add('CaCl2', 0.01, 'mol')
    problem_bc.add('CO2', 0.75, 'mol')

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    return state_bc


# ### Indices of partitioning fluid and solid species
# Only species in fluid phases are mobile and transported by advection and diffusion mechanisms. The solid phases are
# immobile. The code below identifies the indices of the fluid and solid species. We use methods
# [indicesFluidSpecies](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html#ac2a8b713f46f7a66b2731ba63faa95ad)
# and [indicesSolidSpecies](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html#a8b0c237fff1d827f7bf2dbc911fa5bbf)
# of class [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) to get the indices of the
# fluid and solid species, which are stored in the lists `ifluid_species` and `isolid_species`, respectively.

def partition_indices(system):
    nelems = system.numElements()

    ifluid_species = system.indicesFluidSpecies()
    isolid_species = system.indicesSolidSpecies()

    return nelems, ifluid_species, isolid_species


# ### Partitioning fluid and solid species
# In this function, we create arrays to keep track of the amounts of elements in the fluid and solid partition
# (i.e., the amounts of elements among all fluid phases, here only an aqueous phase, and the amounts of elements among
# all solid phases, here the mineral phases). For that, we define the arrays `b`, `bfluid`, `bsolid`, that
# will store, respectively, the concentrations (mol/m<sup>3</sup>) of each element in the system, in the fluid
# partition, and in the solid partition at every time step.
#
# The array `b` is initialized with the concentrations of the elements at the initial chemical state, `state_ic`,
# using method [elementAmounts](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#a827457e68a90f89920c13f0cc06fda78)
# of class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html). The array `b_bc` stores
# the concentrations of each element on the boundary in mol/m<sup>3</sup><sub>fluid</sub>.

def partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc):
    # The concentrations of each element in each mesh cell (in the current time step)
    b = np.zeros((ncells, nelems))
    # Initialize the concentrations (mol/m3) of the elements in each mesh cell
    b[:] = state_ic.elementAmounts()

    # The concentrations (mol/m3) of each element in the fluid partition, in each mesh cell
    bfluid = np.zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the solid partition, in each mesh cell
    bsolid = np.zeros((ncells, nelems))

    # Initialize the concentrations (mol/m3) of each element on the boundary
    b_bc = state_bc.elementAmounts()

    return b, bfluid, bsolid, b_bc


# ### Reactive transport cycle
#
# #### Transport
# This step updates in the fluid partition `bfluid` using the transport equations (without reactions).
# The `transport_fullimplicit()` function below is responsible for solving an advection-diffusion equation, that is
# later applied to transport the concentrations mol/m<sup>3</sup> of elements in the fluid partition (*a
# simplification that is possible because of common diffusion coefficients and velocities of the fluid species,
# otherwise the transport of individual fluid species would be needed*).
#
# To match the units of concentrations of the elements in the fluid measured in mol/m<sup>3</sup><sub>bulk</sub> and the
# imposed concentration `b_bc[j]` mol/m<sup>3</sup><sub>fluid</sub>, we need to multiply it by the porosity `phi_bc`
# m<sup>3</sup><sub>fluid</sub>/m<sup>3</sup><sub>bulk</sub> on the boundary cell . We use function
# [properties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#ad3fa8fd9e1b948da7a698eb020513f3d)
# of the class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) to retrieve fluid volume
# m<sup>3</sup><sub>fluid</sub> and total volume m<sup>3</sup><sub>bulk</sub> in the inflow boundary cell.
#
# The updated amounts of elements in the fluid partition are then summed with the amounts of elements in the solid partition
# `bsolid`, which remained constant during the transport step), and thus updating the amounts of elements in the
# chemical system `b`. Reactive transport calculations involve the solution of a system of
# advection-diffusion-reaction equations.

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
# The function `transport()` expects a conservative property (argument `u`) (e.g., the concentration mol/m<sup>3</sup>
# of *j*th element in the fluid given by `bfluid[j]`), the time step (`dt`), the mesh cell length (`dx`),
# the fluid velocity (`v`), the diffusion coefficient (`D`), and the boundary condition of the conservative property
# (`g`) (e.g., the concentration of the *j*th element in the fluid on the left boundary).
#
# The transport equations are solved with a finite volume method, where diffusion and convection are treated implicitly.
# Its discretization in space and time (implicit) results in the constants `alpha` and `beta`. These correspond to
# the diffusion and advection terms in the equation: `D*dt/dx**2` and `v*dt/dx`, respectively.
#
# Arrays `a`, `b`, `c` are the diagonals in the tridiagonal matrix that results by writing all discretized equations
# in a matrix equation. This system of linear equations is solved by the tridiagonal matrix algorithm, also known
# as the Thomas algorithm.

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
# The tridiagonal matrix equation is solved using the Thomas algorithm (or the TriDiagonal Matrix Algorithm (TDMA)).
# It is a simplified form of Gaussian elimination that can be used to solve tridiagonal systems of equations.

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
# The chemical equilibrium calculations performed in each mesh cell, using *Gibbs energy minimisation* algorithm (
# provided by the class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)).

def reactive_chemistry(solver, states, b):
    # Equilibrating all cells with the updated element amounts
    for icell in range(ncells):
        solver.solve(states[icell], T, P, b[icell])
    return states

# ### Saving and analyzing of the results
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


# # Main parts of the tutorial
#
# First, we create folders for the results:

make_results_folders()

# Run the reactive transport simulations:

simulate()

# To inspect the collected data, one can run:

df

# ## Plotting of the obtained results
#
# The last block of the main routine is dedicated to the plotting of the results in a Jupyter app generated by the
# library **bokeh**. It is an interactive visualization library that provides elegant, concise construction of
# versatile graphics, and affords high-performance interactivity over large or streaming datasets.
#
# First, we have to collect files with results corresponding to each step:

print("Collecting files...")
files = [file for file in natsorted(os.listdir(folder_results))]


# To make a Bokeh application that helps us study the results of the reactive transport simulations, a function
# `modify_doc_df(doc)` that creates Bokeh content and adds it to the app.

def modify_doc_df(doc):
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
    p1.line(x='x', y='pH',
            color='teal', line_width=2, legend_label='pH', source=source)
    p1.x_range = Range1d(-0.001, 1.001)
    p1.y_range = Range1d(4.0, 9.0)
    p1.xaxis.axis_label = 'Distance [m]'
    p1.yaxis.axis_label = 'pH'
    p1.legend.location = 'bottom_right'
    p1.title.text = titlestr(0 * dt)

    # Plot for calcite and dolomite
    p2 = figure(plot_width=600, plot_height=250)
    p2.line(x='x', y='calcite', color='blue', line_width=2, legend_label='Calcite',
            muted_color='blue', muted_alpha=0.2,
            source=source)
    p2.line(x='x', y='dolomite', color='orange', line_width=2, legend_label='Dolomite',
            muted_color='orange', muted_alpha=0.2,
            source=source)
    p2.x_range = Range1d(-0.001, 1.001)
    p2.y_range = Range1d(-0.1, 2.1)
    p2.xaxis.axis_label = 'Distance [m]'
    p2.yaxis.axis_label = 'Mineral Volume [%vol]'
    p2.legend.location = 'center_right'
    p2.title.text = titlestr(0 * dt)
    p2.legend.click_policy = 'mute'

    # Plot for aqueous species
    p3 = figure(plot_width=600, plot_height=300, y_axis_type='log')
    p3.line(x='x', y='Cacation', color='blue', line_width=2, legend_label='Ca++', source=source)
    p3.line(x='x', y='Mgcation', color='orange', line_width=2, legend_label='Mg++', source=source)
    p3.line(x='x', y='HCO3anion', color='green', line_width=2, legend_label='HCO3-', source=source)
    p3.line(x='x', y='CO2aq', color='red', line_width=2, legend_label='CO2(aq)', source=source)
    p3.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
    p3.x_range = Range1d(-0.001, 1.001)
    p3.y_range = Range1d(0.5e-6, 1e0)
    p3.xaxis.axis_label = 'Distance [m]'
    p3.yaxis.axis_label = 'Concentration [molal]'
    p3.legend.location = 'top_right'
    p3.title.text = titlestr(0 * dt)
    p3.legend.click_policy = 'mute'

    step_slider = Slider(start=0,
                         end=nsteps,
                         value=0,
                         step=1,
                         title='Reactive Transport Step')

    # Callback function that changes the data and the title of the plots
    def update_step(attrname, old, new):
        step = step_slider.value
        new_source = ColumnDataSource(df[df['step'] == step])
        source.data = dict(x=new_source.data['x'],
                           pH=new_source.data['pH'],
                           calcite=new_source.data['calcite'],
                           dolomite=new_source.data['dolomite'],
                           Hcation=new_source.data['Hcation'],
                           Cacation=new_source.data['Cacation'],
                           Mgcation=new_source.data['Mgcation'],
                           HCO3anion=new_source.data['HCO3anion'],
                           CO2aq=new_source.data['CO2aq'])
        p1.title.text = titlestr(step * dt)
        p2.title.text = titlestr(step * dt)
        p3.title.text = titlestr(step * dt)

    # Put the update function on the slider
    step_slider.on_change('value', update_step)

    layout = column(step_slider,
                    column(p1, p2, p3))

    doc.add_root(layout)


# Outputting the plots to the notebook requires the call of `output_notebook()` that specifies outputting the plot
# inline in the Jupyter notebook. Finally, the function `modify_doc()` must be passed to `show`, so that the app defined
# by it is displayed inline.
#
# > **Important:** If you run this tutorial in the *localhost*, make sure that number provided to the variable
# `notebook_url` below  coincides with the number of the localhost you have in your browser.
#
# In the app below, we provide an interactive *slider* witget that defines the number of the reactive transport time
# step. This automatically updates the plots with ph, volume phases of calcite and dolomite, and mollalities of aqueous
# species (in logarithmic scale).

output_notebook()
show(modify_doc_df, notebook_url="http://localhost:8890")
