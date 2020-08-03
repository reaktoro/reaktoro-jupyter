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

# # One-dimensional reactive transport modeling of scaling (without oil)
#
# This tutorial demonstrates sequential reactive transport calculations of the barite
# scaling resulting from the waterflooding of the oil reservoirs.
#
# We start by importing Python packages (including **reaktoro**) to enable all the necessary numerical calculations
# as well as the analysis and visualization of the obtained results.

# +
from reaktoro import *
import numpy as np
from tqdm.notebook import tqdm
import os
import pandas as pd

# Import components of bokeh library
from bokeh.io import show, output_notebook
from bokeh.plotting import figure
from bokeh.models import Range1d, ColumnDataSource
from bokeh.layouts import gridplot
# -

# ## Initializing auxiliary time-related constants, discretization, and physical parameters
#
# We start with the initialization of the auxiliary time-related constants from seconds up to hours used in the rest of
# the code.

second = 1
minute = 60
hour = 60 * minute

# Next, we define reactive transport and numerical discretization parameters. In particular, we specify the rock size
# by setting coordinates of its left and right boundaries to 0.0 m and 25.0 m, respectively. The discretization
# parameters, i.e., the number of cells and steps in time, are set to 243 and 900, respectively. Considering that
# the time-step `dt` is fixed to one hour, the simulation lasts 900 hours. Out of these 900 hours, the first
# 45 hours the completion brine (CB) is injected, whereas the rest of the time seawater (SW) is pumped.
# The reactive transport modeling procedure assumes a constant fluid velocity of 0.8e-5 · 10<sup>-5</sup> m/s and
# the zero diffusion coefficient for all fluid species. Temperature and pressure are set to 60 &deg;C and 1 atm,
# respectively, throughout the tutorial. The porosity of the rock is assumed to be 10%. Finally, we use 1kg of water
# in all the mixtures.

# +
# Discretization parameters
xl = 0.0                # x-coordinate of the left boundary
xr = 25.0               # x-coordinate of the right boundary
ncells = 243            # number of cells in the discretization
dx = (xr - xl) / ncells # length of the mesh cells (in units of m)
dt = 1 * hour           # time step

t_cb = 45 * hour                # duration of the completion brine (CB) injection
t_sw = 855 * hour               # duration of the seawater (SW) injection
nsteps_cb = 45                  # number of time steps of the completion brine (CB) injection
nsteps_sw = 855                 # number of time steps of the seawater (SW) injection
nsteps = nsteps_cb + nsteps_sw  # the total number of steps in the reactive transport simulation

# Physical parameters
D = 0               # diffusion coefficient (in units of m2/s)
v = 0.8e-5          # fluid pore velocity (in units of m/s)
T = 60.0            # temperature (in units of celsius)
P = 200             # pressure (in units of atm)
phi = 0.1           # the porosity
water_kg = 1        # amount of water in the initial and injected chemical states
# -

# Next, we generate the mesh nodes (array `xcells`) by equally dividing the interval *[xr, xl]* with the number of cells
# `ncells`. The length between each consecutive mesh nodes (mesh-size) is computed and stored in `dx`.

xcells = np.linspace(xl, xr, ncells + 1)  # interval [xl, xr] split into ncells

# The boolean variable `dirichlet` is set to `True` or `False`, depending on which boundary condition, is considered in
# the numerical calculation. `False` corresponds to imposing the flux of the injected fluid; otherwise, `True` means
# imposing fixed fluid's composition on the left boundary.

dirichlet = False  # parameter that determines whether Dirichlet BC must be used

# To ensure that the applied finite-volume scheme is stable, we need to keep track of the Courant–Friedrichs–Lewy (CFL)
# number, which should be less than 1.0.

CFL = v * dt / dx
print("CFL = ", CFL)
assert CFL <= 1.0, f"Make sure that CFL = {CFL} is less that 1.0"

# ## Specifying the quantities and properties to be outputted
#
# Before running the reactive transport simulations, we provide a list of parameters we are interested in outputting.
# In this case, it is `pH`, molality of `H+`, `Cl-`, `SO4--`, `Ba++`, `Ca++`, `Sr++`, `Na+`, as well as the concentration,
# the phase amount, and the volume of barite (BaSO<sub>4</sub>) mineral.

output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Cl-)
    speciesMolality(SO4--)
    speciesMolality(Ba++)
    speciesMolality(Ca++)
    speciesMolality(Sr++)
    speciesMolality(Na+)
    speciesMolality(Barite)
    phaseAmount(Barite)
    phaseVolume(Barite)
""".split()

# Then, we define the list of names for the `DataFrame` columns. Note, that they must correspond to the order of the
# properties defined in the `output_quantities` list:

column_quantities = """
    pH
    Hcation
    Clanion
    SO4anion
    Bacation
    Cacation
    Srcation
    Nacation
    Barite
    Barite_phase_amount
    Barite_phase_volume
    Barite_SI
""".split()

# > **Note**: All the properties mentioned above (except barite's saturation/equilibrium index) are available using
# the class [ChemicalOutput](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html). To retrieve properties
# that are not included in `ChemicalOutput` functionality, we add several code lines to the `outputstate_df()` function.

# Create the list of columns stored in the dataframe structure and initialize the `DataFrame` instance with it.

columns = ['step', 'x'] + column_quantities
df = pd.DataFrame(columns=columns)

# Finally, we create required folders for outputting the obtained results:

folder_results = 'results-rt-scaling'
def make_results_folders():
    os.system('mkdir -p ' + folder_results)

# ## Performing the reactive transport simulation

# Subsections below correspond to the methods responsible for each of the functional parts of the reactive transport
# simulation performed in the main part of the tutorial.
#
# ### Construction of the chemical system with its phases and species
#
# To define the chemical system, we need to initialize the class
# [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html)
# that provides operations to retrieve physical and thermodynamic data of chemical species. Here,
# [supcrt07.xml](https://github.com/reaktoro/reaktoro/blob/master/databases/supcrt/supcrt07.xml) database file is used.
# In addition to that, we initialize parameters in the Debye-Huckel activity model used for aqueous mixtures. Method
# `setPHREEQC` allows setting parameters *&#229;* and *b* of the ionic species according to those used in PHREEQC v3.
#
# To specify how the chemical system should be defined, i.e., defining all phases and the chemical species they
# contain, one must use [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class (see the
# method `define_chemical_system()` below). Here, we specify two phases, an *aqueous* and a *mineral*. The aqueous phase
# is defined by providing the list of elements (so that all possible combinations of these elements are considered when
# creating a set of chemical species). Function `setChemicalModelDebyeHuckel()` helps to set the chemical model of
# the phase with the Debye-Huckel equation of state, providing specific parameters `dhModel` defined earlier.
# The mineral phase corresponds to mineral barite (BaSO<sub>4</sub>).

def define_chemical_system():

    # Construct the chemical system with its phases and species
    db = Database('supcrt07.xml')

    # Initialize parameters in the Debye-Huckel activity model and set the corresponding parameters
    dhModel = DebyeHuckelParams()
    dhModel.setPHREEQC()

    # Define the phases of the chemical system
    editor = ChemicalEditor(db)
    editor.addAqueousPhaseWithElements("H Cl S O Ba Ca Sr Na K Mg C Si").\
        setChemicalModelDebyeHuckel(dhModel)
    editor.addMineralPhase('Barite')

    # Define the chemical system via the chemical editor
    system = ChemicalSystem(editor)

    return system

# ### Initial condition (IC) of the reactive transport problem
#
# The chemical state corresponding to the **initial condition** of the reactive transport simulation is defined in
# the function `define_initial_condition_fw()`. The composition of the formation water (FW) is taken from the manuscript
# of Bethke 2008, i.e., Table 30.1 (Miller analysis). We recite it in the table below:
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
# The main characteristics of the formation water are high concentration of the Ba<sup>2+</sup> and low concentrations
# of the SO<sub>4</sub><sup>2-</sup>.

def define_initial_condition_fw(system):

    # Define composition of the initial chemical problem
    problem_ic = EquilibriumInverseProblem(system)
    problem_ic.setTemperature(T, "celsius")
    problem_ic.setPressure(P, "atm")
    problem_ic.add("H2O", water_kg, "kg")
    problem_ic.add("SO4", 10 * water_kg, "ug")
    problem_ic.add("Ca", 995 * water_kg, "mg")
    problem_ic.add("Ba", 995 * water_kg, "mg")
    problem_ic.add("Sr", 105 * water_kg, "mg")
    problem_ic.add("Na", 27250 * water_kg, "mg")
    problem_ic.add("K", 1730 * water_kg, "mg")
    problem_ic.add("Mg", 110 * water_kg, "mg")
    problem_ic.add("Cl", 45150 * water_kg, "mg")
    problem_ic.add("HCO3", 1980 * water_kg, "mg")
    problem_ic.pH(7.0, "HCl", "NaOH")

    # Calculate the equilibrium states for the initial conditions
    state_ic = equilibrate(problem_ic)

    # Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3') # 10% of porosity
    state_ic.scaleVolume(1.0, 'm3')

    # Fetch teh value of the ph in the initial chemical state
    props = state_ic.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(FW) = ", pH.val)

    return state_ic

# ### Boundary condition (BC) of the reactive transport problem

# Next, we define the **boundary condition** of the constructed chemical system applied the first 45 hours of simulations.
# The 7-mol sodium chloride brine below defines completion brine (CB).

def define_boundary_condition_cb(system):

    # Define the boundary condition of the reactive transport modeling problem corresponding to completion brine
    problem_bc = EquilibriumProblem(system)
    problem_bc.setTemperature(T, "celsius")
    problem_bc.setPressure(P, "atm")
    problem_bc.add("H2O", water_kg, "kg")
    problem_bc.add("NaCl", 7, "mol")

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    # Fetch ph of the evaluated chemical state
    props = state_bc.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(CB) = ", pH.val)

    return state_bc

# For the rest of the simulation time, we inject the seawater (SW). Similarly, the seawater (SW) composition
# is taken from Bethke 2008, i.e., Table 30.1 (Seawater). Below, we list the quantities of the aqueous species:
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
# Normally, seawater is rich in sulfate, poor in Ca<sup>2+</sup>, and nearly depleted in Sr<sup>2+</sup> and
# Ba<sup>2+</sup>. The pH of the seawater is fixed to 8.1.

def define_boundary_condition_sw(system):

    # Define the boundary condition of the reactive transport modeling problem corresponding to seawater
    problem_bc = EquilibriumInverseProblem(system)
    problem_bc.setTemperature(T, "celsius")
    problem_bc.setPressure(P, "atm")
    problem_bc.add("H2O", water_kg, "kg")
    problem_bc.add("SO4--", 2710 * water_kg, "mg")
    problem_bc.add("Ca++", 411 * water_kg, "mg")
    problem_bc.add("Ba++", 0.01 * water_kg, "mg")
    problem_bc.add("Sr++", 8 * water_kg, "mg")
    problem_bc.add("Na+", 10760 * water_kg, "mg")
    problem_bc.add("K+", 399 * water_kg, "mg")
    problem_bc.add("Mg++", 1290 * water_kg, "mg")
    problem_bc.add("Cl-", 19350 * water_kg, "mg")
    problem_bc.add("HCO3-", 142 * water_kg, "mg")
    problem_bc.pH(8.1, "HCl", "NaOH")

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    props = state_bc.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(SW) = ", pH.val)

    return state_bc

# Alternative (but very similar to the previous) definition of the SW chemical state can be found from the official
# PHREEQC tutorials, i.e., [Example 1](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-70.html).

def define_boundary_condition_sw_phreeqc(system):

    # Define the boundary condition of the reactive transport modeling problem corresponding to seawater
    problem_bc = EquilibriumInverseProblem(system)
    problem_bc.setTemperature(T, "celsius")
    problem_bc.setPressure(P, "atm")
    problem_bc.add("H2O", water_kg, "kg")
    problem_bc.add("Ca++", 412.3 * water_kg, "mg")  # 412.3 mg / kg = 0.4123 kg / kg => 0.4123 * 58 = 23.9134
    problem_bc.add("Mg++", 1291.8 * water_kg, "mg")  # 1291.8 mg / kg = 1.2918 kg / kg => 1.2918 * 58 = 74.9244
    problem_bc.add("Na+", 10768.0 * water_kg, "mg")  # 10768.0  mg / kg = 10.768 kg / kg => 10.768 * 58 = 624.544
    problem_bc.add("K+", 399.1 * water_kg, "mg")  # 399.1 mg / kg = 0.3991 kg / kg => 0.3991 * 58 = 23.1478
    problem_bc.add("Si", 4.28 * water_kg, "mg")  # 4.28 mg / kg = 0.00428 kg / kg => 0.00428 * 58 = 0.24824
    problem_bc.add("Cl-", 19353 * water_kg, "mg")  # 19353.0 mg / kg = 19.353 kg / kg => 19.353 * 58 = 1122.474
    problem_bc.add("HCO3-", 141.682 * water_kg, "mg")  # 141.682 mg / kg = 0.142682 kg / kg => 0.142682 * 58 = 8.275556
    problem_bc.add("SO4--", 2712 * water_kg, "mg") # 2712.0 mg / kg = 2.712 kg / kg => 2.712 * 58 = 157.296
    problem_bc.pH(8.22, "HCl")
    problem_bc.pE(8.451)

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    props = state_bc.properties()
    evaluate_pH = ChemicalProperty.pH(system)
    pH = evaluate_pH(props)
    print("ph(SW, PHREEQC) = ", pH.val)

    return state_bc

# ### Indices of partitioning fluid and solid species
#
# We use methods
# [indicesFluidSpecies](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html#ac2a8b713f46f7a66b2731ba63faa95ad)
# and [indicesSolidSpecies](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html#a8b0c237fff1d827f7bf2dbc911fa5bbf)
# of class [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) to get the indices of the
# fluid and solid species, to separate mobile from immobile ones.

def partition_indices(system):

    # Get number of elements
    nelems = system.numElements()

    # Get indices of species in fluids and solid partitions
    ifluid_species = system.indicesFluidSpecies()
    isolid_species = system.indicesSolidSpecies()

    return nelems, ifluid_species, isolid_species


# ### Partitioning fluid and solid species
#
# Next, we create arrays to track the amounts of elements in the fluid and solid partition.
# We define the arrays `b`, `bfluid`, `bsolid`, storing the concentrations (mol/m<sup>3</sup>) of each element in the
# system, and those in the fluid partition and in the solid partition at every time step.
#
# The array `b` is initialized with the concentrations of the elements at the initial chemical state, `state_ic`, using
# method [elementAmounts](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#a827457e68a90f89920c13f0cc06fda78)
# of class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html). The array `b_bc` stores the
# concentrations of each element on the boundary in mol/m<sup>3</sup><sub>fluid</sub> and is obtained similarly from
# the `state_bc`.

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
#
# This step updates in the fluid partition `bfluid` using the transport equations (without reactions).
# The `transport_fullimplicit()` function below is responsible for solving an advection-diffusion equation, that is
# later applied to transport the concentrations of elements in the fluid partition. This is a
# simplification that is possible because of common diffusion coefficients and velocities of the fluid species,
# otherwise, the transport of individual fluid species would be needed.
#
# To match the units of concentrations of the elements in the fluid measure in mol/m<sup>3</sup><sub>bulk</sub> and the
# imposed concentration `b_bc[j]` mol/m<sup>3</sup><sub>fluid</sub>, we need to scale it by the porosity `phi_bc`
# on the boundary cell m<sup>3</sup><sub>fluid</sub>/m<sup>3</sup><sub>bulk</sub>. We use function
# [properties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#ad3fa8fd9e1b948da7a698eb020513f3d)
# of the class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) to retrieve fluid volume
# m<sup>3</sup><sub>fluid</sub> and total volume m<sup>3</sup><sub>bulk</sub> in the inflow boundary cell.
#
# Finally, the updated amounts of elements in the fluid partition are summed with the amounts of elements in the solid
# partition (remained constant during the transport step), and thus updating the amounts of elements in the chemical
# system `b`. Reactive transport calculations involve the solution of a system of advection-diffusion-reaction equations.

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


# #### Transport calculation with the finite-volume scheme
#
# The function `transport()` expects a conservative property (argument `u`) (e.g., the concentration mol/m<sup>3</sup>
# of the *j*th element in the fluid given by `bfluid[j]`), the time step (`dt`), the mesh cell length (`dx`),
# the fluid velocity (`v`), the diffusion coefficient (`D`), and the boundary condition of the conservative property
# (`g`) (e.g., the concentration of the *j*th element in the fluid on the left boundary).
#
# The transport equations are solved with a finite volume method, where diffusion and convection are treated implicitly.
# Its discretization in space and time (implicit) results in the constants `alpha` and `beta`. These correspond to
# the diffusion and advection terms in the equation: `D*dt/dx**2` and `v*dt/dx`.
#
# Arrays `a`, `b`, `c` are the diagonals in the tridiagonal matrix that results by writing all discretized equations
# in a matrix equation. This linear equation system is solved by the tridiagonal matrix algorithm, also known
# as the [Thomas algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).

def transport_fullimplicit(u, dt, dx, v, D, ul):

    # Number of DOFs
    n = len(u)
    # Fetch the coefficients bearing the diffusion and advection terms
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


# #### Solving the system of equations obtained from finite volume discretization
#
# The tridiagonal matrix equation is solved using the
# [Thomas algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) (or the TriDiagonal Matrix Algorithm (TDMA)).

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
#
# The chemical equilibrium calculations performed in each mesh cell, using the *Gibbs energy minimization* algorithm
# (provided by the class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)).
# **Note:** Before providing temperature and pressure to the `solve()` method, we need to convert Celsius and atmosphere
# to kelvin and bars, respectively.

def reactive_chemistry(solver, states, b):
    # Equilibrating all cells with the updated element amounts
    for icell in range(ncells):
        T_kelvin = T + 273.15
        P_bar = P * 1.01325 * 1e5
        solver.solve(states[icell], T_kelvin, P_bar, b[icell])
    return states


# ### Saving and analyzing of the results
#
# Function `outputstate_df` is the auxiliary function to add data to the `DataFrame` instance at each time step.

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

        # Fetch Barite's saturation index
        phase_SI = state.phaseStabilityIndices()
        barite_phase_index = system.indexPhase("Barite")
        values[-1] = phase_SI[barite_phase_index]

        # Add values into the dataframe
        df.loc[len(df)] = values

# ### Plotting of the obtained results

# The library **bokeh** enables the plotting of the results in a Jupyter app. Below, we list auxiliary functions that
# we use in plotting. Function `titlestr` returns a string for the title of a figure in the  format Time: #h##m

def titlestr(t):
    t = t / minute  # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: %2dh %2dm' % (h, m)

# Routines `plot_figures_ph()`, `plot_figures_barite_phase_amount()`, `plot_figures_barite_concetration()`,
# `plot_figues_barite_saturation_index()`, and 'plot_figures_aqueous_species()' drawing the plots with
# chemical states or properties on the selected steps.

def plot_figures_ph(steps):

    plots = []
    for i in steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='pH', color='teal', line_width=2, legend_label='pH', source=source)
        p.x_range = Range1d(xl-0.1, xr+0.1)
        p.y_range = Range1d(6.8, 10.0)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'pH'
        p.legend.location = 'top_right'
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
        p.line(x='x', y='Barite_phase_amount', color='steelblue', line_width=2, legend_label='Barite',
               muted_color='steelblue', muted_alpha=0.2, source=source)
        p.x_range = Range1d(xl-0.1, xr+0.1)
        p.y_range = Range1d(-0.0001, 0.7)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Mineral Phase Amount [mol]'
        p.legend.location = 'center_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figures_barite_concetration(steps):
    plots = []
    for i in steps:
        print("On barite figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='Barite', color='darkorchid', line_width=2, legend_label='Barite',
               muted_color='darkorchid', muted_alpha=0.2, source=source)
        p.x_range = Range1d(xl-0.1, xr+0.1)
        p.y_range = Range1d(-1e-5, 8e-4)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Concentration [mol/m3]'
        p.legend.location = 'center_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figues_barite_saturation_index(steps):
    plots = []
    for i in steps:
        print("On barite's SI figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='Barite_SI', color='indianred', line_width=2, legend_label='SI (Barite)',
               muted_color='teal', muted_alpha=0.2, source=source)
        p.x_range = Range1d(xl-0.1, xr+0.1)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'SI [-]'
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

        p = figure(plot_width=600, plot_height=250, y_axis_type = 'log',)
        p.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
        p.line(x='x', y='Clanion', color='darkcyan', line_width=2, legend_label='Cl-', source=source)
        p.line(x='x', y='SO4anion', color='darkorange', line_width=2, legend_label='SO4--', source=source)
        p.line(x='x', y='Bacation', color='seagreen', line_width=2, legend_label='Ba++', source=source)
        p.line(x='x', y='Cacation', color='indianred', line_width=2, legend_label='Ca++', source=source)
        p.line(x='x', y='Srcation', color='darkblue', line_width=2, legend_label='Sr++', source=source)
        p.x_range = Range1d(xl-0.1, xr+0.1)
        p.y_range = Range1d(1e-10, 1e2)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Concentration [molal]'
        p.legend.location = 'top_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

# # Main part of the tutorial
#
# First, we create folders for the result files.

make_results_folders()

# Construct the chemical system with its phases and species.

system = define_chemical_system()

# Define the initial condition of the reactive transport modeling problem.

state_ic = define_initial_condition_fw(system)

# Define the boundary condition of the reactive transport modeling problem composed of two different stages.

# +
# Define the completion brine (CB)
state_bc_cb = define_boundary_condition_cb(system)

# Define the seawater (SW)
state_bc_sw = define_boundary_condition_sw(system)
# -

# Generate indices of partitioning fluid and solid species, as well as the vectors of elements, including those in
# fluid and solid partition, and vector of element amounts corresponding to the injected completion brine and seawater.

nelems, ifluid_species, isolid_species = partition_indices(system)
b, bfluid, bsolid, b_bc_cb, b_bc_sw \
    = partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc_cb, state_bc_sw)

# Create a list of chemical states for each mesh cell initialized `to state_ic`.

states = [state_ic.clone() for _ in range(ncells + 1)]

# Create the equilibrium solver object for the repeated equilibrium calculation.

solver = EquilibriumSolver(system)

# Running the reactive transport simulation loop. We start with the completion brine injection.

# +
step = 0  # the current step number
t = 0.0  # the current time (in seconds)

# Output the initial state of the reactive transport calculation
outputstate_df(step, system, states)

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
print(f"time: {t / hour} hours")
# -

# Simulation continues with the seawater injection.

with tqdm(total=nsteps_sw, desc="855 hours of seawater (SW) injection") as pbar:
    while step < nsteps_cb + nsteps_sw:
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

# To save the results in csv-format, the line below can be executed.

df.to_csv(folder_results + '/rt.scaling.csv', index=False)


# Outputting the plots to the notebook requires the call of `output_notebook()` that specifies outputting the plot
# inline in the Jupyter notebook:

output_notebook()

# Select the steps on which results with pH must be output.

selected_steps_to_plot = [45, 60, 300]
assert all(step <= nsteps for step in selected_steps_to_plot), f"Make sure that selceted steps are less than " \
                                                               f"total amount of steps {nsteps}"
plot_figures_ph(selected_steps_to_plot)

# Select the steps on which the rest of the species must be demonstrated.

selected_steps_to_plot = [120, 300, 600]
assert all(step <= nsteps for step in selected_steps_to_plot), f"Make sure that selceted steps are less than " \
                                                               f"total amount of steps {nsteps}"

# After the execution of functions `plot_figures_barite_concetration()` and `plot_figues_barite_saturation_index()`,
# we see that the mineral starts to precipitate right after we initiate the injection of seawater. It happens due to
# the mixing of the formation water (FW) and seawater (SW), which have contrasting compositions. In particular, SW
# is low on Ba<sup>2+</sup> and high on SO<sub>4</sub><sup>2-</sup> concentrations, whereas FW, on the opposite, is
# high on Ba<sup>2+</sup> and low on SO<sub>4</sub><sup>2-</sup>. During the mixing, both of these ions react creating
# barite according to the reaction Ba<sup>2+</sup> + SO<sub>4</sub><sup>2-</sup> &#8594; BaSO<sub>4</sub>(s).
# Such side effect of waterflooding (as part of the oil recovery techniques) reduces the near-wellbore permeability and
# hampers well productivity/injectivity. Other mineral that have potential to scaling are CaCO<sub>3</sub> (calcite),
# CaSO<sub>4</sub> (calcium sulfate), FeCO<sub>3</sub> (siderite).

# Plot barite's concetration on the selected steps:
plot_figures_barite_concetration(selected_steps_to_plot)

# Plot barite's saturation index on the selected steps:
plot_figues_barite_saturation_index(selected_steps_to_plot)

# In the plots of the aqueous species below, we see that we have a decrease of Ba<sup>2+</sup> concentrations during
# mixing, consumed by the barite. The concentration of SO<sub>4</sub><sup>2-</sup> in SW is considerably
# higher than in FW, which explains the monotonically decreasing curve the left to the right
# boundary. Finally, we observe a slight increase in the Cl<sup>-</sup> concentration due to the initial NaCl-brine
# injection for 45 hours.

# Plot aqueous species on the selected steps:
plot_figures_aqueous_species(selected_steps_to_plot)

# To study the time-dependent behavior of the chemical properties, we create a Bokeh application using the function
# `modify_doc(doc)`. It creates Bokeh content and adds it to the app. The speed of streaming of reactive transport
# data can be  controlled by the parameter `step` defined below (bigger the step, faster we will run through available
# data set):

step = 10
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
    p1 = figure(plot_width=500, plot_height=250)
    p1.line(x='x', y='pH', color='teal', line_width=2, legend_label='pH', source=source)
    p1.x_range = Range1d(xl-0.1, xr+0.1)
    p1.y_range = Range1d(6.8, 10.0)
    p1.xaxis.axis_label = 'Distance [m]'
    p1.yaxis.axis_label = 'pH'
    p1.legend.location = 'bottom_right'
    p1.title.text = titlestr(0 * dt)

    p2 = figure(plot_width=500, plot_height=250)
    p2.line(x='x', y='Barite', color='darkorchid', line_width=2, legend_label='Barite',
           muted_color='darkorchid', muted_alpha=0.2, source=source)
    p2.x_range = Range1d(xl, xr - 1)
    p2.y_range = Range1d(0.0, 8e-4)
    p2.xaxis.axis_label = 'Distance [m]'
    p2.yaxis.axis_label = 'Concentration [mol/m3]'
    p2.legend.location = 'center_right'
    p2.title.text = titlestr(0 * dt)
    p2.legend.click_policy = 'mute'

    p3 = figure(plot_width=500, plot_height=250, y_axis_type='log')
    p3.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
    p3.line(x='x', y='Clanion', color='darkcyan', line_width=2, legend_label='Cl-', source=source)
    p3.line(x='x', y='SO4anion', color='darkorange', line_width=2, legend_label='SO4--', source=source)
    p3.line(x='x', y='Bacation', color='seagreen', line_width=2, legend_label='Ba++', source=source)
    p3.line(x='x', y='Cacation', color='indianred', line_width=2, legend_label='Ca++', source=source)
    p3.line(x='x', y='Srcation', color='darkblue', line_width=2, legend_label='Sr++', source=source)
    p3.x_range = Range1d(xl-0.1, xr+0.1)
    p3.y_range = Range1d(1e-8, 1e1)
    p3.xaxis.axis_label = 'Distance [m]'
    p3.yaxis.axis_label = 'Concentration [molal]'
    p3.legend.location = 'top_right'
    p3.title.text = titlestr(0 * dt)
    p3.legend.click_policy = 'mute'

    p4 = figure(plot_width=500, plot_height=250)
    p4.line(x='x', y='Barite_SI', color='indianred', line_width=2, legend_label='SI (Barite)',
           muted_color='teal', muted_alpha=0.2, source=source)
    p4.x_range = Range1d(xl, xr - 1)
    p4.xaxis.axis_label = 'Distance [m]'
    p4.yaxis.axis_label = 'SI [-]'
    p4.legend.location = 'center_right'
    p4.title.text = titlestr(0 * dt)
    p4.legend.click_policy = 'mute'

    layout =  gridplot([[p1, p2], [p3, p4]])

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
                        Barite=new_source.data['Barite'],
                        Barite_phase_volume=new_source.data['Barite_phase_volume'],
                        Barite_SI=new_source.data['Barite_SI'])

        p1.title.text = titlestr(step_number * dt)
        p2.title.text = titlestr(step_number * dt)
        p3.title.text = titlestr(step_number * dt)
        p4.title.text = titlestr(step_number * dt)

        source.stream(new_data, rollover=ncells + 1)

    doc.add_periodic_callback(update, 500)
    doc.add_root(layout)

# Outputting the plots to the notebook requires the call of `output_notebook()` that specifies outputting the plot
# inline in the Jupyter notebook. Finally, the function `modify_doc()` must be passed to `show`, so that the app defined
# by it is displayed inline.
#
# > **Important:** If this tutorial is executed in the *localhost*, make sure that the number provided to the variable
# `notebook_url` below coincides with the number of the localhost of the browser.
#
# In the app below, we refresh the reactive time step in a loop, automatically updating the data source for the
# plots for ph, volume phases of calcite and dolomite, and mollalities of aqueous species (in logarithmic scale).

output_notebook()
show(modify_doc, notebook_url="http://localhost:8888")
