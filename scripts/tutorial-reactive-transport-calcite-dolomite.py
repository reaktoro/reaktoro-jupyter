# -*- coding: utf-8 -*-
# # Reactive transport modelling along a rock core after injection of the fluid-rock composition

# In this tutorial, we show how Reaktoro can be used for sequential calculations of the reactive transport along a
# rock column after injecting the fluid and rock composition at temperature 60 &deg; C and pressure 100 bar.

# Using Reaktoro in Python requires first an import of the python package **reaktoro**. From this point on,
# we can use the library components of Reaktoro (classes, methods, constants), which are needed to define our
# chemical system and chemical reaction modeling problems.

# ## Importing python packages

# First, we need to import a few Python packages to enable us to perform the numerical calculations and plotting.

print('============================================================')
print('Make sure you have the following Python packages installed: ')
print('     numpy, matplotlib, joblib, ffmpeg')
print('These can be installed with pip:')
print('     pip install numpy matplotlib joblib ffmpeg')
print('============================================================')
from reaktoro import *
from numpy import *
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import os

# We import the **reaktoro** Python package so that we can use its classes and methods for performing chemical
# reaction calculations, **numpy** for working with arrays, **matplotlib** for plotting capabilities, **joblib** for
# simple parallel computing, **os**, to provide a portable way of using operating system dependent
# functionality. Finally, *ffmpeg* must be installed for handling video, audio, and other multimedia files and streams.

# > **Note**: To simplify the tutorials, we use `from reaktoro import *`, which imports all components of the
# > **reaktoro** package into the default Python namespace. We note that this can potentially create name conflicts
# > when used in bigger projects. For your applications, consider using `import reaktoro as rkt` instead,
# > and call classes and methods as `rkt.Database`, `rkt.ChemicalSystem`, `rkt.equilibrate`, and etc.

# ## Initializing auxiliary time-related constants
# In this step, we initialize auxiliary time-related constants from seconds up to years used in the rest of the code.

second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

# ## Defining parameters for the reactive transport simulation
# Next, we define reactive transport and numerical discretization parameters. In particular, we specify the considered
# rock domain by setting coordinates of its left and right boundaries to 0.0 m and 100.0 m, respectively. The
# discretization parameters, i.e., the number of cells and steps in time, are both set to 100. The reactive
# transport modeling procedure assumes a constant fluid velocity of 1 m/week (1.65 · $10^{-6}$ m/s) and
# the same diffusion coefficient of $10^{-9}$ m2/s for all fluid species (without dispersivity). The size of the
# time-step is set to 30 minutes. Temperature and pressure are set to 60 &deg; C and 100 bar, respectively,
# throughout the whole tutorial. 

xl = 0.0  # x-coordinate of the left boundary
xr = 1.0  # x-coordinate of the right boundary
ncells = 100  # number of cells in the discretization
nsteps = 10  # number of steps in the reactive transport simulation
D = 1.0e-9  # diffusion coefficient (in units of m2/s)
v = 1.0 / week  # fluid pore velocity (in units of m/s)
dx = (xr - xl) / ncells  # length of the mesh cells (in units of m)
dt = 30 * minute  # time step
T = 60.0 + 273.15  # temperature (in units of K)
P = 100 * 1e5  # pressure (in units of Pa)

# Next, we generate the coordinates of the mesh nodes (array `x`) by equally
# dividing the interval *[xr, xl]* into `ncells`. The length between each
# consecutive mesh nodes is computed and stored in `dx` (the length of the mesh
# cells).

x = linspace(xl, xr, ncells + 1)  # interval [xl, xr] split into ncells

# The boolean variable `dirichlet` is set to `True` or `False` depending on which boundary condition is considered in
# the numerical calculation. Set to `False` for imposing the flux of the injected fluid, otherwise, set to
# `True` for imposing the composition of the fluid on the left boundary.

dirichlet = False  # parameter that determines whether Dirichlet BC must be used

# Another auxiliary parameter is the number of digits in the number of steps (e.g., 100 has 3 digits). It is needed for
# generation of the names of the files, where chemical states are saved, as well as the creating of the videos

ndigits = len(str(nsteps))

# To make sure that the applied finite-volume scheme is stable, we need to keep track of Courant–Friedrichs–Lewy (CFL)
# number, which should be less than 1.0.

CFL = v * dt / dx
print(f"Make sure that CFL = {v * dt / dx} is less that 1.0")

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


# ## Organization of the main function
# The main function (at the bottom of this tutorial) consists of three parts, each represented by a Python function and
# documented in the following sections:
# * creation of folders for the results (function `make_results_folders()`),
# * simulation of reactive transport problem (method `simulate()`), and
# * plotting of the obtained results (function `plotting()`).

# ## Creating folders for the outputted results
#
# Using **os** package, we create required folders for outputting the obtained results and for the plot and video
# files later.

def make_results_folders():
    os.system('mkdir -p results')
    os.system('mkdir -p figures/ph')
    os.system('mkdir -p figures/aqueous-species')
    os.system('mkdir -p figures/calcite-dolomite')
    os.system('mkdir -p videos')


# ## Performing the reactive transport simulation
#
# The reactive transport simulation is performed in the function `simulate`, which consists from the several building
# blocks (functions):
# * initialization of the reactive transort problem and 
# * performing the reactive transport simulation along defined time interval.
#
# The preparitory initialization step consists of the following substeps:
# * definition of chemical system with its phases and species using `define_chemical_system()`,
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
    states = [state_ic.clone() for _ in range(ncells)]

    # Create the equilibrium solver object for the repeated equilibrium calculation
    solver = EquilibriumSolver(system)

    # Running the reactive transport simulation loop
    step = 0  # the current step number
    t = 0.0  # the current time (in seconds)

    while step <= nsteps:
        # Print the progress of the simulation
        print("Progress: {}/{} steps, {} min".format(step, nsteps, t / minute))

        # Perform transport calculations
        bfluid, bsolid, b = transport(states, bfluid, bsolid, b, b_bc, nelems, ifluid_species, isolid_species)

        # Perform reactive chemical calculations
        states = reactive_chemistry(solver, states, b)

        # Output the current state of the reactive transport calculation
        outputstate(step, system, states)

        # Increment time step and number of time steps
        t += dt
        step += 1

    print("Finished!")


# Subsections below correspond to the methods responsible for each of the functional parts of `simulate()` method.

# ### Construction of the chemical system with its phases and species
#
# Reaktoro is a general-purpose chemical solver that avoids as much as possible presuming specific assumptions about
# your problems. Thus, you need to specify how your chemical system should be defined. This encompasses the
# specification of all phases in the system as well as the chemical species that compose each phase. By using the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class, you can conveniently achieve
# this as shown below in method `define_chemical_system()`.

# In this step, we create an object of class [ChemicalEditor](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) and specify two phases, an *aqueous* and a
# *mineral*, should be considered in the chemical system. The aqueous phase is defined by using a list of
# compounds, which will be broken into a list of element names, and the database will then be searched for all
# species that could be formed out of those elements. The mineral phase is defined as three mineral species:
# quartz ($\mathrm{SiO_2}$), calcite ($\mathrm{CaCO3}$), and dolomite ($\mathrm{CaMg(CO_3)_2}$).
#
# > **Note**: An automatic search for chemical species can result in a large number of species in the phase,
# > potentially causing the chemical reaction calculations to be more computationally expensive. If you are using
# > Reaktoro in demanding applications (e.g., as a chemical solver in a reactive transport simulator), you might want
# > to manually specify the chemical species of each phase in your chemical system. This can be achieved by providing a
# > list of species names as `editor.addAqueousPhaseWithElements([ 'H2O(l)', 'H+', 'OH-', 'Na+', 'Cl-', 'Ca++', 'Mg++',
# > 'HCO3-', 'CO2(aq)', 'CO3--' ]).

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
# > **Note**: The *activity coefficients* of the aqueous species are calculated using the *HKF extended Debye-Hückel
# > model* for solvent water and ionic species, except for the aqueous species $\mathrm{CO_2(aq)}$, for which the
# >  *Drummond model* is used.
# > The *standard chemical potentials* of the species are calculated using the equations of state of Helgeson and
# > Kirkham (1974), Helgeson et al. (1974), Tanger and Helgeson (1988), Shock and Helgeson (1988), and Shock et al. (
# > 1992). The database file [slop98.dat](https://github.com/reaktoro/reaktoro/blob/master/databases/supcrt/slop98.dat)
# > from the software SUPCRT92 is used to obtain the parameters for the equations of state.
# > The equation of state of Wagner and Pruss (2002) is used to calculate the *density of water* and its temperature and
# > pressure derivatives. Kinetics of *dissolution* and *precipitation* of both calcite and dolomite is neglected, i.e.,
# > the local equilibrium assumption is employed.

def define_chemical_system():
    # Step 7.1: Construct the chemical system with its phases and species
    db = Database('supcrt98.xml')

    editor = ChemicalEditor(db)

    editor.addAqueousPhaseWithElements('H O Na Cl Ca Mg C Si Ca')
    editor.addMineralPhase('Quartz')
    editor.addMineralPhase('Calcite')
    editor.addMineralPhase('Dolomite')

    system = ChemicalSystem(editor)

    return system


# ### Initial condition (IC) of the reactive transport problem
#
# We have now defined and constructed our chemical system of interest, enabling us to move on to the next step in
# Reaktoro's modeling workflow: *defining our chemical reaction problems*. Below we define its **initial condition**
# with already prescribed equilibrium conditions for *temperature*, *pressure*, and *amounts of elements* that are
# consistent with an intention to model reactive transport of injected $\mathrm{NaCl-MgCl_2-CaCl_2}$ brine into the
# rock-fluid composition of quartz and calcite at 60 &deg;C and 100 bar. In particular, we consider resident fluid is a
# 0.7 molal $\mathrm{NaCl}$ brine in equilibrium with the rock minerals with a calculated pH of 10.0.
#

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

# > **Note**: Our **0.7 molal NaCl brine** is here represented by the mixture of 1 kg of $\mathrm{H_20}$ and 0.7 mol of
# > $\mathrm{NaCl}$.
#
# > **Note**: After providing the amounts of substances $\mathrm{H_20}$, $\mathrm{NaCl}$, quartz $\mathrm{SiO_2}$,
# > and calcite $\mathrm{CaCO_3}$ in the above code, Reaktoro parses these chemical formulas and determines the
# > elements and their coefficients. Once this is done, the amount of each element stored inside the object
# > `problem_ic` is incremented according to the given amount of substance and its coefficient in the formula. The
# > amounts of elements you provide are then used as constraints for the Gibbs energy minimization calculation when
# > computing the state of chemical equilibrium (i.e., when we try to find the amounts of all species in the system
# > that corresponds to a state of minimum Gibbs energy and at the same time satisfying the *element amounts
# > constraints*).
#
# > **Note**: Please note that we are not condemning the input form shown above in terms of element amounts,
# but only telling you to be attentive with the values you input. If you are using Reaktoro as a chemical reaction
# solver in a reactive transport simulator, for example, you'll most likely need to work directly with given amounts
# of elements, which shows that this input form is required in certain cases. For such time-dependent modeling
# problems, you often only need to ensure that the initial conditions for elements amounts result in feasible initial
# species amounts.
#
# Next, we use method [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e)
# to calculate the chemical equilibrium state of the system with the given initial conditions stored in the objects
# `problem_ic`. The numerical solution of each problem results in the objects `state_ic` of class
# [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html), which stores the temperature,
# pressure, and the amounts of every species in the system.
# For this calculation, Reaktoro uses an efficient **Gibbs energy minimization** computation to determine the species
# amounts that correspond to a state of minimum Gibbs energy in the system, while satisfying the prescribed amount
# conditions for temperature, pressure, and element amounts.

# The function ends with scaling the volumes of the aqueous and mineral phases so that they are consistent with a 10 %
# porosity and the required volume percentages of the rock minerals (98 $\%_{\rm vol}$ of quartz and 2 $\%_{\rm vol}$
# of calcite).

# ### Boundary condition (BC) of the reactive transport problem
#
# Next, we define the **boundary condition** of the constructed chemical system with its *temperature*, *pressure*,
# and *amounts of elements*. In particular, we prescribe the amount of injected aqueous fluid resulting from
# the mixture of 1 kg of water with 0.90 moles of $\mathrm{NaCl}$, 0.05 molal $\mathrm{MgCl_2}$, 0.01 $\mathrm{
# CaCl_2}$, and 0.75 molal $\mathrm{CO_2}$, in a state very close to $\mathrm{CO_2}$ saturation.
# The temperature and the pressure stays the same, i.e., 60 &deg; C and 100 bar, respectively.
#
# After equilibration, the obtained chemical state representing the boundary condition for the injected fluid
# composition, we scale its volume to 1 $\mathrm{m^3}$. This is done so that the amounts of the species in the fluid are
# consistent with a \mathrm{mol/m^3} scale.

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
# will store, respectively, the concentrations ($\mathrm{mol/m^3}$) of each element in the system, in the fluid
# partition, and in the solid partition at every time step.
#
# The array `b` is initialized with the concentrations of the elements at the initial chemical state, `state_ic`,
# using method [elementAmounts](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#a827457e68a90f89920c13f0cc06fda78)
# of class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html). The array `b_bc` stores
# the concentrations of each element on the boundary in $\mathrm{[mol/m^3_{fluid}]}$.

def partition_elements_in_mesh_cell(ncells, nelems, state_ic, state_bc):
    # The concentrations of each element in each mesh cell (in the current time step)
    b = zeros((ncells, nelems))
    # Initialize the concentrations (mol/m3) of the elements in each mesh cell
    b[:] = state_ic.elementAmounts()

    # The concentrations (mol/m3) of each element in the fluid partition, in each mesh cell
    bfluid = zeros((ncells, nelems))

    # The concentrations (mol/m3) of each element in the solid partition, in each mesh cell
    bsolid = zeros((ncells, nelems))

    # Initialize the concentrations (mol/m3) of each element on the boundary
    b_bc = state_bc.elementAmounts()

    return b, bfluid, bsolid, b_bc


# ### Reactive transport cycle
#
# #### Transport
#
# This step updates in the fluid partition `bfluid` using the transport equations (without reactions).
# The `transport_fullimplicit()` function below is responsible for solving an advection-diffusion equation, that is
# later applied to transport the concentrations ($\mathrm{mol/m^3}$) of elements in the fluid partition (*a
# simplification that is possible because of common diffusion coefficients and velocities of the fluid species,
# otherwise the transport of individual fluid species would be needed*).
#
# To match the units of concentrations of the elements in the fluid measure in $\mathrm{[mol/m^3_{bulk}]}$ and the
# imposed concentration `b_bc[j]` $\mathrm{[mol/m^3_{fluid}]}$, e need to multiply it by the porosity `phi_bc` on the
# boundary cell $\mathrm{[m^3_{fluid} / m^3_{bulk}]}$. We use function
# [properties](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html#ad3fa8fd9e1b948da7a698eb020513f3d)
# of the class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) to retrieve fluid volume
# $\mathrm{[m^3_{fluid}]}$ and total volume $\mathrm{[m^3_{bulk}]}$ in the inflow boundary cell.
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
    phi_bc = states[bc_cell].properties().fluidVolume().val / states[bc_cell].properties().volume().val;

    # Transport each element in the fluid phase
    for j in range(nelems):
        transport_fullimplicit(bfluid[:, j], dt, dx, v, D, phi_bc * b_bc[j])

    # Update the amounts of elements in both fluid and solid partitions
    b[:] = bsolid + bfluid

    return bfluid, bsolid, b


# ##### Transport calculation with finite-volume scheme

# The function `transport()` expects a conservative property (argument `u`) (e.g., the concentration $\mathrm{[mol/m^3]}$
# of *j*th element in the fluid given by `bfluid[j]`), the time step (`dt`), the mesh cell length (`dx`),
# the fluid velocity (`v`), the diffusion coefficient (`D`), and the boundary condition of the conservative property
# (`g`) (e.g., the concentration of the *j*th element in the fluid on the left boundary).
#
# The transport equations are solved with a finite volume method, where diffusion and convection are treated implicitly.
# Its discretization in space and time (implicit) results in the constants `alpha` and `beta`. These correspond to
# the diffusion and advection terms in the equation: `D*dt/dx**2` and `v*dt/dx`, respectively.

# Arrays `a`, `b`, `c` are the diagonals in the tridiagonal matrix that results by writing all discretized equations
# in a matrix equation. This system of linear equations is solved by the tridiagonal matrix algorithm, also known
# as the Thomas algorithm.

def transport_fullimplicit(u, dt, dx, v, D, ul):
    # Number of DOFs
    n = len(u)
    alpha = D * dt / dx ** 2
    beta = v * dt / dx

    # Upwind finite volume scheme
    a = full(n, -beta - alpha)
    b = full(n, 1 + beta + 2 * alpha)
    c = full(n, -alpha)

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
#
# The chemical equilibrium calculations performed in each mesh cell, using *Gibbs energy minimisation* algorithm (
# provided by the class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html)).

def reactive_chemistry(solver, states, b):
    # Equilibrating all cells with the updated element amounts
    for icell in range(ncells):
        solver.solve(states[icell], T, P, b[icell])
    return states


# ### Results saving and analyzing

# Function `outputstate` is the auxiliary function to create an output file each time step.

def outputstate(step, system, states):
    # Create the instance of ChemicalOutput class
    output = ChemicalOutput(system)

    # Provide the output file name, which will correspond
    output.filename('results/{}.txt'.format(str(step).zfill(ndigits)))

    # We define the columns' tags filled with the name of the quantities
    # The first column has a tag 'x' (which corresponds to the center coordinates of the cells )
    output.add('tag', 'x')  # The value of the center coordinates of the cells

    # The rest of the columns correspond to the requested properties
    for quantity in output_quantities:
        output.add(quantity)

    # We update the file with states that correspond to the cells' coordinates stored in x
    output.open()
    for state, tag in zip(states, x):
        output.update(state, tag)
    output.close()


# ## Plotting of the obtained results
#
# The last block of the main routine is dedicated to plotting of the results and generating a video from the plots to
# illustrate the time-dependent behavior of chemical properties.

def plot():
    # Plot all result files
    files = sorted(os.listdir('results'))
    Parallel(n_jobs=16)(delayed(plotfile)(file) for file in files)

    # Create videos for the figures
    ffmpegstr = 'ffmpeg -y -r 30 -i figures/{0}/%0' + str(
        ndigits) + 'd.png -codec:v mpeg4 -flags:v +qscale -global_quality:v 0 videos/{0}.mp4'
    os.system(ffmpegstr.format('calcite-dolomite'))
    os.system(ffmpegstr.format('aqueous-species'))
    os.system(ffmpegstr.format('ph'))


# Generate figures for a result file
def plotfile(file):
    # Fetch the step number from name of the file
    step = int(file.split('.')[0])
    # Calculate corresponding to this step time moment
    t = step * dt

    # Auxiliary function that return a string for the title of a figure in the format Time: #h##m
    def titlestr(t):
        t = t / minute  # Convert from seconds to minutes
        h = int(t) / 60  # The number of hours
        m = int(t) % 60  # The number of remaining minutes
        return 'Time: %2dh %2dm' % (h, m)

    print('Plotting figure', step, '...')

    # Loading data from he analyzed file
    filearray = loadtxt('results/' + file, skiprows=1)
    data = filearray.T

    # Plotting ph w.r.t. rock length
    plt.figure()
    plt.xlim(left=xl - 0.02, right=xr + 0.02)
    plt.ylim(bottom=2 - 0.1, top=9.0)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('pH')
    plt.plot(data[0], data[1])
    plt.tight_layout()
    plt.savefig('figures/ph/{}.png'.format(str(step).zfill(ndigits)))

    # Plotting phase volumes of calcite and dolomite w.r.t. rock length
    plt.figure()
    plt.xlim(left=xl - 0.02, right=xr + 0.02)
    plt.ylim(bottom=-0.1, top=2.1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    plt.plot(data[0], data[7] * 100, label='Calcite')
    plt.plot(data[0], data[8] * 100, label='Dolomite')
    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig('figures/calcite-dolomite/{}.png'.format(str(step).zfill(ndigits)))

    # Plotting concentrations of aqueous species w.r.t. rock length
    plt.figure()
    plt.yscale('log')
    plt.xlim(left=xl - 0.02, right=xr + 0.02)
    plt.ylim(bottom=0.5e-5, top=2)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Concentration [molal]')
    plt.plot(data[0], data[3], label='Ca++')
    plt.plot(data[0], data[4], label='Mg++')
    plt.plot(data[0], data[5], label='HCO3-')
    plt.plot(data[0], data[6], label='CO2(aq)')
    plt.plot(data[0], data[2], label='H+')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('figures/aqueous-species/{}.png'.format(str(step).zfill(ndigits)))

    plt.close('all')


# Step 5: Define the main function
if __name__ == '__main__':
    # Create folders for the results
    make_results_folders()

    # Run the reactive transport simulations
    simulate()

    # Plotting the result
    plot()
