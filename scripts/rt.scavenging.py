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

# # Reactive transport modeling of the H<sub>2</sub>S scavenging process along a rock core
#
# In this tutorial, we show how Reaktoro can be used for sequential reactive transport calculations of the injected
# H<sub>2</sub>S-rich brine into siderite bearing reservoir with subsequent FeS precipitation.
#
# Using Reaktoro in Python requires first an import of the python package **reaktoro**. From this point on,
# we can use the library components of Reaktoro (classes, methods, constants), which are needed to define our
# chemical system and chemical reaction modeling problems.
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
from bokeh.layouts import gridplot
# -

# We import the **reaktoro** Python package so that we can use its classes and methods for performing chemical
# reaction calculations, **numpy** for working with arrays, **tqdm** for the progress bar functionality and **os**,
# to provide a portable way of using operating system dependent functionality. For plotting capabilities of obtained
# results, we use **bokeh** library.
#
# > **Note**: To simplify the tutorials, we use `from reaktoro import *`, which imports all components of the
# > **reaktoro** package into the default Python namespace. We note that this can potentially create name conflicts
# > when used in bigger projects. For your applications, consider using `import reaktoro as rkt` instead,
# > and call classes and methods as `rkt.Database`, `rkt.ChemicalSystem`, `rkt.equilibrate`, etc.
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
# transport modeling procedure assumes a constant fluid velocity of 1.05 · 10<sup>-5</sup> m/s and
# the zero diffusion coefficient for all fluid species. The size of the time-step is set to 0.05 days (1.2 hours).
# Temperature and pressure are set to 25 &deg;C and 1.01325 bar, respectively, throughout the whole tutorial.
# The porosity fo the rock is set to 10%.

# +
# Discretization parameters
xl = 0.0                # x-coordinate of the left boundary
xr = 100.0              # x-coordinate of the right boundary
ncells = 100            # number of cells in the discretization
nsteps = 1000            # number of steps in the reactive transport simulation
dx = (xr - xl) / ncells # length of the mesh cells (in units of m)
dt = 0.05*day           # time step

# Physical parameters
D = 0               # diffusion coefficient (in units of m2/s)
v = 1.05e-5         # fluid pore velocity (in units of m/s)
T = 25.0 + 273.15   # temperature (in units of K)
P = 1.01325 * 1e5   # pressure (in units of Pa)
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
    speciesMolality(HS-)
    speciesMolality(S2--)
    speciesMolality(CO3--)
    speciesMolality(HSO4-)
    speciesMolality(H2S(aq))
    phaseAmount(Pyrrhotite)
    phaseAmount(Siderite)
    phaseVolume(Pyrrhotite)
    phaseVolume(Siderite)
    speciesMolality(Fe++)
""".split()

# Then, we define the list of names for the DataFrame columns. Note, that they must correspond
# to the order of the properties defined in the `output_quantities` list:

column_quantities = """
    pH
    Hcation
    HSanion
    S2anion
    CO3anion
    HSO4anion
    H2Saq
    pyrrhotite_phase_amount
    siderite_phase_amount
    pyrrhotite_phase_volume
    siderite_phase_volume
    Fe2cation
""".split()

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
#
# Using **os** package, we create required folders for outputting the obtained results and for the plot and video
# files later.

folder_results = 'results-rt-scaveging'
def make_results_folders():
    os.system('mkdir -p ' + folder_results)

# ## Performing the reactive transport simulation
#
# The reactive transport simulation is performed in the function `simulate`, which consists of several building
# blocks (functions):
# * initialization of the reactive transport problem and
# * performing the reactive transport simulation along with a defined time interval.
#
# The preparatory initialization step consists of the following sub-steps:
# * definition of the chemical system with its phases and species using `define_chemical_system()`,
# * definition of the initial condition of the reactive transport problem in `define_initial_condition()`,
# * definition of the boundary condition of the reactive transport problem in `define_initial_condition()`,
# * generation of auxiliary indices to partition elements using `partition_indices()` and elements' partitioning
# corresponding to fluid and solid species with function `partition_elements_in_mesh_cell()`, and finally
# * definition of instance of [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html).
#
# The simulation of the reactive transport problem is represented by the loop over discretized time interval until
# the final time is reached. On each step of this loop, the following functionality of performed:
# * transport calculations using method `transport()`,
# * reactive chemical calculations with `reactive_chemistry()` function, and
# * saving the obtained results by means of `outputstate()`.
#
# Performing of the transport and reactive chemistry sequentially is possible due to the *operator splitting
# procedure*, in which we first update the amounts of elements `b`. These updated amounts of
# elements in the cell are used to evaluate its new chemical equilibrium state, thus producing new amounts of the
# species in both the fluid and solid phases (available in the list `states` of
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) objects). This chemical reaction
# equilibrium calculation step, at each mesh cell, permits, for example, aqueous species and minerals to react,
# and thus causes mineral dissolution or precipitation, depending on how much the amount of mineral species changes.
# This can then be used, for example, to compute new porosity value for the cell.

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

    # Output the initial state of the reactive transport calculation
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

            # Output the current state of the reactive transport calculation
            outputstate_df(step, system, states)

            # Update a progress bar
            pbar.update(1)


# Subsections below correspond to the methods responsible for each of the functional parts of `simulate()` method.
#
# ### Construction of the chemical system with its phases and species
#
# To define the chemical system, we need to initialize the class
# [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html)
# that provides operations to retrieve physical and thermodynamic data of chemical species. To achieve that, we use
# [supcrt07.xml](https://github.com/reaktoro/reaktoro/blob/master/databases/supcrt/supcrt07.xml) database file.
#
# **Note:** If filename does not point to a valid database file or the database file is not found, then a default
# built-in database with the same name will be tried. If no default built-in database exists with a given name,
# an exception will be thrown.
#
# In addition to the database, we also need to initialize parameters in the Debye-Huckel activity model used for aqueous
# mixtures. Method `setPHREEQC` allows to set parameters *&#229;* and *b* of the ionic species according to those used
# in PHREEQC v3.
#
# Reaktoro is a general-purpose chemical solver that avoids as much as possible presuming specific assumptions about
# your problems. Thus, you need to specify how your chemical system should be defined. This encompasses the
# specification of all phases in the system as well as the chemical species that compose each phase. By using the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class, you can conveniently achieve
# this as shown below in method `define_chemical_system()`.
#
# In this step, we create an object of class [ChemicalEditor](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) and specify two phases, an *aqueous* and a
# *mineral*, should be considered in the chemical system. The aqueous phase is defined by specifying the list of
# chemical species. Function `setChemicalModelDebyeHuckel()` helps to set the chemical model of the phase with
# the Debye-Huckel equation of state, providing specific parameters `dhModel` defined earlier. The mineral phases
# are defined as two mineral species: pyrrhotite (FeS) and siderite (FeCO<sub>3</sub>).
#
# Finally, we create an object of class [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html)
# using the chemical system definition details stored in the object `editor`.
#
# > **Note**: [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) is perhaps the main class
# > in Reaktoro. An object of this class stores the phases, species, and elements in our defined chemical system,
# > as well as provides the means to compute many types of thermodynamic properties, such as *standard thermodynamic
# > properties* (e.g., standard Gibbs energies, standard enthalpies, and standard volumes of species),
# > and *thermo-chemical properties* (e.g., activity and activity coefficients of species; density, enthalpy and
# > internal energy of phases). As you learn more about other Reaktoro's classes, you will note that an object of class
# > [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) is almost always needed for their
# > initialization.
#
# The *activity coefficients* of the aqueous species in this tutorial are calculated using the
# *Debye-Huckel model* for solvent water and ionic species.
# The *standard chemical potentials* of the species are calculated using the equations of state of Helgeson and
# Kirkham (1974), Helgeson et al. (1974), Tanger and Helgeson (1988), Shock and Helgeson (1988), and Shock et al. (
# 1992). The database file [slop07.dat](https://github.com/reaktoro/reaktoro/blob/master/databases/supcrt/slop07.dat)
# from the software SUPCRT07 is used to obtain the parameters for the equations of state.
# The equation of state of Wagner and Pruss (2002) is used to calculate the *density of water* and its temperature and
# pressure derivatives. Kinetics of *dissolution* and *precipitation* of both pyrrhotite and siderite is neglected, i.e.,
# the local equilibrium assumption is employed.

def define_chemical_system():

    # Construct the chemical system with its phases and species
    db = Database('supcrt07.xml')
    #db = Database('phreeqc.dat')

    dhModel = DebyeHuckelParams()
    dhModel.setPHREEQC()

    editor = ChemicalEditor(db)
    editor.addAqueousPhase(['Ca(HCO3)+', 'CO3--', 'CaCO3(aq)', 'Ca++', 'CaSO4(aq)', 'CaOH+', 'Cl-',
                            'FeCl++', 'FeCl2(aq)', 'FeCl+', 'Fe++', 'FeOH+', 'FeOH++', 'Fe+++',
                            'H2S(aq)', 'H2(aq)', 'HS-', 'H2O(l)', 'H+', 'OH-', 'HCO3-', 'HSO4-',
                            'KSO4-',  'K+',
                            'Mg++', 'Mg(HCO3)+', 'MgCO3(aq)', 'MgSO4(aq)', 'MgOH+',
                            'Na+', 'NaSO4-',
                            'O2(aq)',
                            'S5--', 'S4--', 'S3--', 'S2--', 'SO4--']).\
        setChemicalModelDebyeHuckel(dhModel)

    editor.addMineralPhase('Pyrrhotite')
    editor.addMineralPhase('Siderite')

    system = ChemicalSystem(editor)

    return system

# ### Initial condition (IC) of the reactive transport problem
#
# We have now defined and constructed the chemical system of interest, enabling us to move on to the next step in
# Reaktoro's modeling workflow: *defining our chemical reaction problems*. Below, we define its **initial condition**
# with already prescribed equilibrium conditions for *temperature*, *pressure*, and *amounts of elements* that are
# consistent to model reactive transport of injected hydrogen sulfide brine into the rock-fluid composition of siderite
# at 25 &deg;C and 1.01325 bar. The resident fluid in the rock is obtained by the mixture of the aqueous species
# summarized in the following table:
#
# | Aqueous species | Amount (kg)                          |
# |-----------------|--------------------------------------|
# | H<sub>2</sub>O  | 58.0                                 |
# | Cl<sup>-</sup>O | 1122.3 · 10<sup>-3</sup>             |
# | Na<sup>+</sup>O | 624.08 · 10<sup>-3</sup>             |
# | SO<sub>4</sub><sup>2-</sup>O | 157.18 · 10<sup>-3</sup>|
# | Mg<sup>2+</sup>O | 74.820 · 10<sup>-3</sup>            |
# | Ca<sup>2+</sup>O | 23.838 · 10<sup>-3</sup>            |
# | K<sup>+</sup>O | 23.142 · 10<sup>-3</sup>              |
# | HCO<sub>3</sub><sup>-</sup>O | 8.236 · 10<sup>-3</sup> |
# | O<sub>2</sub>(aq) | 58 · 10<sup>-12</sup>              |
#
# For that purpose, the class
# [EquilibriumInverseProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html) is used, where
# specific fixed pH and pE are be prescribed to 8.951 and 8.676, respectively.

def define_initial_condition(system):

    problem_ic = EquilibriumInverseProblem(system)
    problem_ic.setTemperature(T)
    problem_ic.setPressure(P)
    problem_ic.add("H2O", 58.0, "kg");
    problem_ic.add("Cl-", 1122.3e-3, "kg");
    problem_ic.add("Na+", 624.08e-3, "kg");
    problem_ic.add("SO4--", 157.18e-3, "kg");
    problem_ic.add("Mg++", 74.820e-3, "kg");
    problem_ic.add("Ca++", 23.838e-3, "kg");
    problem_ic.add("K+", 23.142e-3, "kg");
    problem_ic.add("HCO3-", 8.236e-3, "kg");
    problem_ic.add("O2(aq)", 58e-12, "kg");
    problem_ic.add("Pyrrhotite", 0.0, "mol");
    problem_ic.add("Siderite", 0.5, "mol");
    problem_ic.pH(8.951);
    problem_ic.pE(8.676);

    # Calculate the equilibrium states for the initial conditions
    state_ic = equilibrate(problem_ic)

    # Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3') # 10% of porosity
    state_ic.scaleVolume(1.0, 'm3')

    return state_ic

# > **Note**: After providing the amounts of substances H<sub>2</sub>O, aqueous species, pyrrhotite, and siderite in
# > the above code, Reaktoro parses these chemical formulas (using the thermodynamic database) and determines the
# > elements and their coefficients. Once this is done, the amount of each element stored inside the object
# > `problem_ic` is incremented according to the given amount of substance and its coefficient in the formula. The
# > amounts of elements you provide are then used as constraints for the Gibbs energy minimization calculation when
# > computing the state of chemical equilibrium (i.e., when we try to find the amounts of all species in the system
# > that corresponds to a state of minimum Gibbs energy and at the same time satisfying the *element amounts
# > constraints*).
#
# > **Note**: Please note that we are not condemning the input form shown above in terms of element amounts,
# > but only telling you to be attentive with the values you input. If you are using Reaktoro as a chemical reaction
# > solver in a reactive transport simulator, for example, you'll most likely need to work directly with given amounts
# > of elements, which shows that this input form is required in certain cases. For such time-dependent modeling
# > problems, you often only need to ensure that the initial conditions for elements amounts result in feasible initial
# > species amounts.
#
# To calculate the chemical equilibrium state of the system with the given initial conditions, we use the method
# [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e), the numerical
# solution of which is written in the objects `problem_ic`. It is an instance of the class
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) that stores the temperature,
# pressure, and the amounts of every species in the system.
# For this calculation, Reaktoro uses an efficient **Gibbs energy minimization** computation to determine the species
# amounts that correspond to a state of minimum Gibbs energy in the system, while satisfying the prescribed amount
# conditions for temperature, pressure, and element amounts. In an inverse equilibrium problem, however, not all
# elements have known molar amounts. Their amount's constraints are replaced by other equilibrium constraints such as
# fixed pH and pE.
#
# > *See tutorials [**EquilibriumInverseProblem**](eq.inverse-equilibrium.ipynb) for more detailed explanation of
# > capabilities of this class*.
#
# The function ends with scaling the volume to 1 m<sup>3</sup>. Moreover, we specify the 10% porosity of the rock
# by calling `state_ic.scalePhaseVolume('Aqueous', 0.1, 'm3')`.
#
#
# ### Boundary condition (BC) of the reactive transport problem
#
# Next, we define the **boundary condition** of the constructed chemical system with its *temperature*, *pressure*,
# and *amounts of elements*. We prescribe the amount of injected hydrogen sulfide brine, in particular,
# 0.0196504 mol of hydrosulfide ion (HS<sup>-</sup>) and 0.167794 mol of aqueous hydrogen sulfide (H2S(aq)).
# Here, the ph is lowered in comparison to the initial state to 5.726.
#
# After equilibration, the obtained chemical state representing the boundary condition for the injected fluid
# composition, we scale its volume to 1 m<sup>3</sup>. This is done so that the amounts of the species in the fluid are
# consistent with a mol/<sup>3</sup> scale.

def define_boundary_condition(system):

    # Define the boundary condition of the reactive transport modeling problem
    problem_bc = EquilibriumInverseProblem(system)
    problem_bc.setTemperature(T)
    problem_bc.setPressure(P)
    problem_bc.add("H2O", 58.0, "kg")
    problem_bc.add("Cl-", 1122.3e-3, "kg")
    problem_bc.add("Na+", 624.08e-3, "kg")
    problem_bc.add("SO4--", 157.18e-3, "kg")
    problem_bc.add("Mg++", 74.820e-3, "kg")
    problem_bc.add("Ca++", 23.838e-3, "kg")
    problem_bc.add("K+", 23.142e-3, "kg")
    problem_bc.add("HCO3-", 8.236e-3, "kg")
    problem_bc.add("O2(aq)", 58e-12, "kg")
    problem_bc.add("Pyrrhotite", 0.0, "mol")
    problem_bc.add("Siderite", 0.0, "mol")
    problem_bc.add("HS-", 0.0196504, "mol")
    problem_bc.add("H2S(aq)", 0.167794, "mol")
    problem_bc.pH(5.726)
    problem_bc.pE(8.220)

    # Calculate the equilibrium states for the boundary conditions
    state_bc = equilibrate(problem_bc)
    # Scale the boundary condition state to 1 m3
    state_bc.scaleVolume(1.0, 'm3')

    return state_bc


# ### Indices of partitioning fluid and solid species
#
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
#
# In this function, we create arrays to keep track of the amounts of elements in the fluid and solid partition
# (i.e., the amounts of elements among all fluid phases, here only an aqueous phase, and the amounts of elements among
# all solid phases, here the mineral phases). For that, we define the arrays `b`, `bfluid`, `bsolid`, that
# will store, respectively, the concentrations (mol/m<sup>3</sup>) of each element in the system, in the fluid
# partition, and in the solid partition at every time step.
#
# The array `b` is initialized with the concentrations of the elements at the initial chemical state, `state_ic`,
# using method
# [elementAmounts](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#a827457e68a90f89920c13f0cc06fda78)
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
#
# This step updates in the fluid partition `bfluid` using the transport equations (without reactions).
# The `transport_fullimplicit()` function below is responsible for solving an advection-diffusion equation, that is
# later applied to transport the concentrations (mol/m<sup>3</sup>) of elements in the fluid partition (*a
# simplification that is possible because of common diffusion coefficients and velocities of the fluid species,
# otherwise the transport of individual fluid species would be needed*).
#
# To match the units of concentrations of the elements in the fluid measure in mol/m<sup>3</sup><sub>bulk</sub> and the
# imposed concentration `b_bc[j]` mol/m<sup>3</sup><sub>fluid</sub>, e need to multiply it by the porosity `phi_bc`
# on the boundary cell m<sup>3</sup><sub>fluid</sub>/m<sup>3</sup><sub>bulk</sub>. We use function
# [properties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#ad3fa8fd9e1b948da7a698eb020513f3d)
# of the class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) to retrieve fluid volume
# m<sup>3</sup><sub>fluid</sub> and total volume m<sup>3</sup><sub>bulk</sub> in the inflow boundary cell.
#
# The updated amounts of elements in the fluid partition are then summed with the amounts of elements in the solid
# partition `bsolid`, which remained constant during the transport step), and thus updating the amounts of elements
# in the chemical system `b`. Reactive transport calculations involve the solution of a system of
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
#
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
#
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
#
# The chemical equilibrium calculations performed in each mesh cell, using *Gibbs energy minimisation* algorithm (
# provided by the class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)).

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
#
# The last block of the main routine is dedicated to the plotting of the results in a Jupyter app generated by the
# library **bokeh**. It is an interactive visualization library that provides elegant, concise construction of
# versatile graphics, and affords high-performance interactivity over large or streaming datasets.
#
# Below, we list auxiliary functions that we use in plotting. Function `titlestr` returns a string for the title
# of a figure in the  format Time: #h##m

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
        p.x_range = Range1d(-1, 101)
        p.y_range = Range1d(4, 9.0)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'pH'
        p.legend.location = 'bottom_right'
        p.title.text = titlestr(t)

        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figures_pyrrhotite_siderite_volume(steps):
    plots = []
    for i in steps:
        print("On pyrrhotite-siderite figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='pyrrhotite_phase_volume', color='blue', line_width=2, legend_label='Pyrrhotite',
               muted_color='blue', muted_alpha=0.2, source=source)
        p.line(x='x', y='siderite_phase_volume', color='orange', line_width=2, legend_label='Siderite',
               muted_color='orange', muted_alpha=0.2, source=source)
        p.x_range = Range1d(-1, 101)
        p.y_range = Range1d(-0.001, 0.018)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Mineral Volume [%vol]'
        p.legend.location = 'center_right'
        p.title.text = titlestr(t)
        p.legend.click_policy = 'mute'
        plots.append([p])

    grid = gridplot(plots)
    show(grid)

def plot_figures_pyrrhotite_siderite_amount(steps):
    plots = []
    for i in steps:
        print("On pyrrhotite-siderite figure at time step: {}".format(i))
        t = i * dt
        source = ColumnDataSource(df[df['step'] == i])

        p = figure(plot_width=600, plot_height=250)
        p.line(x='x', y='pyrrhotite_phase_amount', color='blue', line_width=2, legend_label='Pyrrhotite',
               muted_color='blue', muted_alpha=0.2, source=source)
        p.line(x='x', y='siderite_phase_amount', color='orange', line_width=2, legend_label='Siderite',
               muted_color='orange', muted_alpha=0.2, source=source)
        p.x_range = Range1d(-1, 101)
        p.y_range = Range1d(-0.5, 5.5)
        p.xaxis.axis_label = 'Distance [m]'
        p.yaxis.axis_label = 'Mineral Amount [mol]'
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
        p.line(x='x', y='HSanion', color='darkcyan', line_width=2, legend_label='HS-', source=source)
        p.line(x='x', y='S2anion', color='darkorange', line_width=2, legend_label='S2--', source=source)
        p.line(x='x', y='CO3anion', color='seagreen', line_width=2, legend_label='CO3--', source=source)
        p.line(x='x', y='HSO4anion', color='indianred', line_width=2, legend_label='HSO4-', source=source)
        p.line(x='x', y='H2Saq', color='gray', line_width=2, legend_label='H2S(aq)', source=source)
        p.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
        p.line(x='x', y='Fe2cation', color='darkblue', line_width=2, legend_label='Fe++', source=source)
        p.x_range = Range1d(-1, 101)
        p.y_range = Range1d(1e-12, 1e-1)
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

# Run the reactive transport simulations:

simulate()

# To inspect the collected data, one can run:

df

# To save the results in csv-format, please execute:

df.to_csv(folder_results + '/rt.scaveging.csv', index=False)

# Select the steps, on which results must plotted:

selected_steps_to_plot = [60, 120, 960]
assert all(step <= nsteps for step in selected_steps_to_plot), f"Make sure that selceted steps are less than " \
                                                               f"total amount of steps {nsteps}"

# Outputting the plots to the notebook requires the call of `output_notebook()` that specifies outputting the plot
# inline in the Jupyter notebook:

output_notebook()

# Plot ph on the selected steps:

plot_figures_ph(selected_steps_to_plot)

# Plot calcite and dolomite on the selected steps:

plot_figures_pyrrhotite_siderite_amount(selected_steps_to_plot)

# One can also call `plot_figures_pyrrhotite_siderite_volume(selected_steps_to_plot)` instead of the function
# `plot_figures_pyrrhotite_siderite_amount()`.

# Plot aqueous species on the selected steps:

plot_figures_aqueous_species(selected_steps_to_plot)

# We see on the plots above that the main chemical reactions can be divided into two parts.
# On the first step, the iron ions Fe<sup>2+</sup> are being released by the siderite (FeCO<sub>3</sub>) reacting with
# the H<sub>2</sub>S-brine. On the second one, free Fe<sup>2+</sup> ions are reacting with hydrogen sulfide
# precipitating iron sulfide (pyrrhotite): Fe<sup>2+</sup> + H<sub>2</sub>S &#8594; FeS<sub>2</sub>.
#
# The minerals' dissolution and precipitation are accompanied by the formation and dilution of aqueous species.
# For instance, we see the sharp increase of Fe<sup>2+</sup> at the point of the phase transformation from the siderite
# phase to pyrrhotite and followed by the gradual decrease as it gets used during the formation of FeS<sub>2</sub>.
# Both curves representing HS<sup>-</sup> and H<sub>2</sub>S(aq) have two points of the sharp decrease.
# The first one is where both species are involved in the dissolution of FeCO<sub>3</sub>, and the second one is where
# they are being exhausted by the reaction with iron ions to form iron sulfide. The CO<sub>3</sub><sup>2-</sup> anion
# line has an interesting shape: it first locally drops but then increases due to continuous dissolution of siderite.
#
# To study the time-dependent behavior of the chemical properties, we create a Bokeh application using the function
# `modify_doc(doc)`. It creates Bokeh content and adds it to the app. The speed of streaming of the reactive transport
# data can be  controlled by the parameter `step` defined below (bigger the step, faster we will run through available
# data set):

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
    p1.x_range = Range1d(-1, 101)
    p1.y_range = Range1d(4.0, 9.0)
    p1.xaxis.axis_label = 'Distance [m]'
    p1.yaxis.axis_label = 'pH'
    p1.legend.location = 'bottom_right'
    p1.title.text = titlestr(0 * dt)

    # Plot for calcite and dolomite
    p2 = figure(plot_width=600, plot_height=250)
    p2.line(x='x', y='pyrrhotite_phase_volume', color='blue', line_width=2,
            legend_label='Pyrrhotite', muted_color='blue', muted_alpha=0.2,
            source=source)
    p2.line(x='x', y='siderite_phase_volume', color='orange', line_width=2,
            legend_label='Siderite', muted_color='orange', muted_alpha=0.2,
            source=source)
    p2.x_range = Range1d(-1, 101)
    p2.y_range = Range1d(-0.001, 0.018)
    p2.xaxis.axis_label = 'Distance [m]'
    p2.yaxis.axis_label = 'Phase Volume [%vol]'
    p2.legend.location = 'center_right'
    p2.title.text = titlestr(0 * dt)
    p2.legend.click_policy = 'mute'
    p3 = figure(plot_width=600, plot_height=300, y_axis_type='log')
    p3.line(x='x', y='HSanion', color='darkcyan', line_width=2, legend_label='HS-', source=source)
    p3.line(x='x', y='S2anion', color='darkorange', line_width=2, legend_label='S2--', source=source)
    p3.line(x='x', y='CO3anion', color='seagreen', line_width=2, legend_label='CO3--', source=source)
    p3.line(x='x', y='HSO4anion', color='indianred', line_width=2, legend_label='HSO4-', source=source)
    p3.line(x='x', y='H2Saq', color='gray', line_width=2, legend_label='H2S(aq)', source=source)
    p3.line(x='x', y='Hcation', color='darkviolet', line_width=2, legend_label='H+', source=source)
    p3.line(x='x', y='Fe2cation', color='darkblue', line_width=2, legend_label='Fe++', source=source)
    p3.x_range = Range1d(-1, 101)
    p3.y_range = Range1d(1e-12, 1e-1)
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
                        pyrrhotite_phase_volume=new_source.data['pyrrhotite_phase_volume'],
                        siderite_phase_volume=new_source.data['siderite_phase_volume'],
                        pyrrhotite_phase_amount=new_source.data['pyrrhotite_phase_amount'],
                        siderite_phase_amount=new_source.data['siderite_phase_amount'],
                        HSanion=new_source.data['HSanion'],
                        S2anion=new_source.data['S2anion'],
                        CO3anion=new_source.data['CO3anion'],
                        HSO4anion=new_source.data['HSO4anion'],
                        H2Saq=new_source.data['H2Saq'],
                        Hcation=new_source.data['Hcation'],
                        Fe2cation=new_source.data['Fe2cation'])

        p1.title.text = titlestr(step_number * dt)
        p2.title.text = titlestr(step_number * dt)
        p3.title.text = titlestr(step_number * dt)

        source.stream(new_data, rollover=ncells+1)

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
