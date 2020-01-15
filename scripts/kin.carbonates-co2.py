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

# # Dissolution of carbonate minerals in a CO<sub>2</sub>-saturated brine
#
# In this tutorial, we demonstrate how Reaktoro can be used to kinetically model the dissolution of carbonate
# minerals (calcite, magnesite, and dolomite) in a CO<sub>2</sub>-saturated brine.

# ### Importing the reaktoro Python package
#
# Here, we import the **reaktoro** Python package to enable us to use all its library components (classes, methods,
# constants).

from reaktoro import *

# ### Defining the chemical system

# In this simulation, we consider an *aqueous phase* (to model the brine solution), a *gaseous phase* (to model the
# CO<sub>2</sub>-rich phase with water vapor), and four pure *mineral phases*: halite, calcite, magnesite,
# and dolomite. These are the phases that will either exist initially during the simulation, or that could
# potentially appear later as the calculations proceed.

# All potential phases that could appear in a reactive process should ideally be considered when defining the
# chemical system. If one or more of these phases are ignored, then the equilibrium and kinetics calculations cannot
# identify if they should be present or not with positive amounts. Unrealistic results may occur, such as,
# for example, an aqueous phase containing more CO<sub>2</sub> dissolved than it could, because a gaseous phase,
# which should contain the excess of CO<sub>2</sub>, was not considered in the chemical system.

# The code below uses class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) to define
# our chemical system with the phases of interest and their species:

editor = ChemicalEditor()
editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3 MgCO3")
editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Magnesite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Halite")

# The aqueous phase is defined by considering all aqueous species in the database that could form once the substances
# H<sub>2</sub>O, NaCl, CaCO<sub>3</sub>, and MgCO<sub>3</sub> are mixed. The gaseous phase is defined
# so that only the gaseous species H<sub>2</sub>O(g) and CO<sub>2</sub>(g) are considered. There are four pure
# mineral phases: calcite (CaCO<sub>3</sub>), magnesite (MgCO<sub>3</sub>), dolomite (CaMg(CO<sub>3</sub>)<sub>2</sub>),
# and halite (NaCl).

# ### Defining the kinetically-controlled reactions

# A *partial equilibrium assumption* is considered in this simulation. This simplifies the problem by assuming that
# those species that react at much faster rates do so *continually in chemical equilibrium*. They are referred to as
# *equilibrium species*. The remaining species (those reacting at relatively slower rates) are referred to as
# *kinetic species*.

# Because aqueous and gaseous species, as well as halite, react at relatively fast rates, they are reasonable
# candidates in this problem for being equilibrium species. The kinetic species are thus the carbonate minerals:
# calcite, magnesite, and dolomite.

# With this partial equilibrium assumption, there is no need to specify kinetic rate laws for the fast-reacting
# equilibrium species. For these species, chemical equilibrium calculations are performed to update their amounts as the
# amounts of each kinetic species change over time (i.e., the equilibrium species react instantaneously to a new
# state of equilibrium as the kinetic species react and perturb the current partial equilibrium state).

# Thus, we only need to define kinetic rates for the relatively slow-reacting carbonate minerals (our kinetic species
# in this simulation), which is shown below

editor.addMineralReaction("Calcite") \
    .setEquation("Calcite = Ca++ + CO3--") \
    .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Magnesite") \
    .setEquation("Magnesite = Mg++ + CO3--") \
    .addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
    .addMechanism("logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
    .setSpecificSurfaceArea(10, "cm2/g")

editor.addMineralReaction("Dolomite") \
    .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
    .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
    .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
    .setSpecificSurfaceArea(10, "cm2/g")

# We set the equation of each mineral reaction using method `setEquation` of class
# [MineralReaction](https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html).

# We then prescribe neutral and acidic mechanisms for each mineral reaction using the method ``addMechanism`` of class
# [MineralMechanism](https://reaktoro.org/cpp/structReaktoro_1_1MineralMechanism.html). We set values for ``logk``,
# the kinetic rate constant of the reaction in log scale, and ``Ea``, the Arrhenius activation energy. The units of
# both parameters must be provided as shown in the example, and other compatible units are allowed.

# Finally, we define the specific surface area of the mineral using method `setSpecificSurfaceArea` of class
# [MineralReaction](https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html). Any units compatible to `m2/kg`
# or `m2/m3` are allowed (e.g., ``cm2/g``, ``mm2/mm3``).

# ### Creating the chemical and reaction systems

# Create instances of
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html).

system = ChemicalSystem(editor)
reactions = ReactionSystem(editor)

# Class [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) is used to represent a
# chemical system, containing one or more phases  (in this problem, aqueous, gaseous, and mineral phases). Each
# phase can contain one or more species (the aqueous phase includes many aqueous species, the gaseous phase two,
# and each mineral phase containing a single mineral species with the same name as the name of the phase). This class
# is also used to calculate the thermodynamic properties of the phases and species (standard thermodynamic properties,
# activities, chemical potentials, phase molar volume, etc.).

# Class [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) is a collection of
# [Reaction](https://reaktoro.org/cpp/classReaktoro_1_1Reaction.html) objects used to represent a system of chemical
# reactions that are controlled by chemical kinetics. These classes provide convenient methods for the calculation of
# equilibrium constants, reaction quotients, and rates of the reactions.

# ### Specifying equilibrium and kinetic species

# In Reaktoro, the species in a chemical system can be partitioned into groups of:
# *equilibrium species*, *kinetic species*, and *inert species*.

# For the *equilibrium species*, their amounts are governed by chemical equilibrium (calculated via a Gibbs energy
# minimization). The amount of *kinetic species* are controlled by chemical kinetics (by solving a system of ordinary
# differential equations that models the kinetics of a system of reactions). The *inert species* maintain their
# amounts constant in all chemical equilibrium or kinetics calculations.

# This classification of species can be done using class
# [Partition](https://reaktoro.org/cpp/classReaktoro_1_1Partition.html). By default, all species are considered to
# be equilibrium species in this class. Thus, we only need to specify which are kinetic ones:

partition = Partition(system)
partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])

# In this case, the mineral species calcite, magnesite, and dolomite are specified to be *kinetic species* using
# method `setKineticSpecies` of class [Partition](https://reaktoro.org/cpp/classReaktoro_1_1Partition.html).

# Method `setKineticPhases` of class [Partition](https://reaktoro.org/cpp/classReaktoro_1_1Partition.html) could also
# be used here. This method sets all species in the given phases to be kinetic species, and it is more convenient if
# a phase has many species. However, since each of the mineral phases considered here only contains a single mineral
# species, the method `setKineticSpecies` is a convenient alternative.

# ### Defining the initial state of the equilibrium species

# In a chemical kinetics calculation, an *initial condition* is needed for the amounts of both equilibrium and
# kinetic species.

# The equilibrium problem formulated below, using class
# [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html), is done so that the
# initial condition for the amounts of each equilibrium species result from the solution of a chemical equilibrium
# problem in which 1 kg of water is mixed with 1 mol of NaCl and 1 mol of CO<sub>2</sub> at 60 &deg;C
# and 100 bar. This amount of CO<sub>2</sub> is sufficient to saturate the brine solution. The excess will exist in
# the CO<sub>2</sub>-rich gaseous phase.

problem = EquilibriumProblem(system)
problem.setPartition(partition)
problem.setTemperature(60, "celsius")
problem.setPressure(100, "bar")
problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CO2", 1, "mol")

# *Note:* To ensure that the equilibrium calculation performed in the next step ignores the kinetic species in the
# system so that we maintain a disequilibrium state between equilibrium and kinetic species, it is important not to
# forget to call `setPartition(partition)` method. Ignoring this step will produce an initial condition for the
# amounts of equilibrium and kinetic species that correspond to a complete equilibrium state in the system so that no
# kinetics calculation makes sense afterwards.

# ### Calculating the initial amounts of the equilibrium species

# We now use the convenient `equilibrate` function to calculate the amounts of the equilibrium species by minimizing
# the Gibbs energy of the equilibrium partition only, and not of the entire system. The result is stored in the
# object `state0` of class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html),
# a computational representation of the state of a multiphase chemical system defined by its temperature (*T*),
# pressure (*P*), and vector of species amounts (*n*). We then output this chemical state to a file.

# ### Setting the initial mass of minerals

# We have now to prescribe the initial amounts of the kinetic species (i.e., the carbonate minerals). This is done
# below by setting the initial mass of calcite to 100 g and of dolomite to 50 g. The initial amount of magnesite is
# zero.

state0 = equilibrate(problem)
state0.output('demo-kineticpath-carbonates-co2-before-kinetics.txt')

state0.setSpeciesMass("Calcite", 100, "g")
state0.setSpeciesMass("Dolomite", 50, "g")

# ### Setting the kinetic path calculations

# To be able to simulate of the kinetic process of mineral dissolution/precipitation, we introduce the
# instance of the class [KineticPath](https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html) that enables this
# functionality. The instance of kinetic path solver is provided with the partition to the equilibrium, kinetic,
# and inert species defined above.

path = KineticPath(reactions)
path.setPartition(partition)

# **Note**: For repeated chemical kinetics calculations (e.g., in a reactive transport simulation, where kinetics is
# performed for each mesh cell/node), consider using the class
# [KineticSolver](https://reaktoro.org/cpp/classReaktoro_1_1KineticSolver.html) instead for avoiding some overhead of
# class [KineticPath](https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html). For equilibrium calculations,
# consider also the class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html),
# instead of method `equilibrate`, for similar reasons.

# To analyse the result of kinetic simulations, we save the evolution of different properties of the chemical system
# into file `result.txt`:

output = path.output()
output.filename("results.txt")
output.add("time(units=minute)")
output.add("pH")
output.add("elementMolality(Ca units=mmolal)", "Ca [mmolal]")
output.add("elementMolality(Mg units=mmolal)", "Mg [mmolal]")
output.add("phaseMass(Calcite units=g)", "Calcite [units=g]")
output.add("phaseMass(Dolomite units=g)", "Dolomite [units=g]")

# **Note**: A list of all possible quantities that can be plotted is shown in the class
# [ChemicalQuantity](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalQuantity.html), which provides an interface
# for convenient ways of their retrieval.

# To perform the calculation of kinetic path, we need to provide time interval, on which this path must be recovered. In
# this case, we choose 25 hours of simulations. Initial chemical state provided by the instance `state0`, as well as
# the time interval with units, in which time is measured, are given with a call of function `solve`.

t0, t1 = 0.0, 25.0
path.solve(state0, t0, t1, "hours")
state0.output('demo-kineticpath-carbonates-co2-after-kinetics.txt')

# ### Plotting the results of equilibrium path calculation
#
# The best way to visually analyze the obtained reaction path is with plots. For that, we export python plotting
# package `matplotlib` and `numpy`, the fundamental package for scientific computing with Python.

import matplotlib.pyplot as plt
import numpy as np

# To load results from the outputfile, we use `loadtxt` function provided by the `numpy` package:

filearray = np.loadtxt("results.txt", skiprows=1)
data = filearray.T
time = data[0, :]
[ph_indx, ca_elem_indx, mg_elem_indx, calcite_indx, dolomite_indx] = np.arange(1, 6)

# The plots below depict different chemical properties (x-axis) with respect to the time interval of the kinetic
# simulation (y-axis):

plt.figure()
plt.plot(time, data[ph_indx], label="pH")
plt.xlabel("time")
plt.ylabel("pH")
plt.legend(loc='center right')
plt.savefig("ph-vs-time.png")

plt.figure()
plt.plot(time, data[ca_elem_indx], label="Ca")
plt.plot(time, data[mg_elem_indx], label="Mg")
plt.xlabel("time")
plt.ylabel("Molality [mmolal]")
plt.legend(loc='center right')
plt.savefig("molality-ca-mg-vs-time.png")

plt.figure()
plt.plot(time, data[calcite_indx], label="Calcite")
plt.xlabel("time")
plt.ylabel("Mass of Calcite [g]")
plt.legend(loc='center right')
plt.savefig("mass-calcite-vs-time.png")

plt.figure()
plt.plot(time, data[dolomite_indx], label="Dolomite")
plt.xlabel("time")
plt.ylabel("Mass of Dolomite [g]")
plt.legend(loc='center right')
plt.savefig("mass-dolomite-vs-time.png")
