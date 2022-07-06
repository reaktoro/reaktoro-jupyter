# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Performing a chemical equilibrium calculation with customized activity models
#
# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation with customized activity
# models. We use an example of equilibration of  H<sub>2</sub>O-NaCl-CO<sub>2</sub> system when 1 kg of H<sub>2</sub>O,
# 100 g of CO<sub>2</sub>, and 1 mol of NaCl are mixed at temperature 60 &deg;C and pressure 300 bar.
# First, we import everything from the `reaktoro` package by

from reaktoro import *

# To indicate phases and corresponding to them species, as well as models used to evaluate activities of the aqueous
# species, we use the class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html). Here,
# the default database SUPCRT98 is used.

editor = ChemicalEditor()

# ## Specifying different phases and corresponding models
#
# Specifying the phases and their species is not enough to fully describe a chemical system in the *computational
# sense*. Every phase in Reaktoro has two associated models: a *thermodynamic model* and a *chemical model*. These
# denominations are not standard in the literature, but they are useful in the differentiation of two needed types
# of models for a phase.
#
# * A *thermodynamic model* is a model for the calculation of *standard thermodynamic properties* of the species.
# Examples include standard Gibbs energies, or standard chemical potentials, standard molar volumes, standard heat
# capacities, standard enthalpies, and so forth. These models are functions of *temperature* and *pressure* only.
#
# * A *chemical model* is a model that describes the *non-ideal behavior* of phases. These models not only depend on
# temperature and pressure, like the thermodynamic models, but also on the amounts of the species in the phase. To be
# more precise, on the concentrations of these species, which can be calculated from the amounts of the species.
#
# One can define a chemical system with many phases, each phase containing one or more species. This does not mean
# that all phases and their species exist at positive amounts! To be precise, this means that the chemical
# calculations, equilibrium or kinetics, are capable of deciding if a phase in the chemical system should exist at
# positive amounts for some given conditions (e.g., temperature, pressure, overall composition).
#
# By selecting as many phases as possible, with the possibilities constrained by the *thermodynamic database* being
# used, one can increase the confidence level of the estimated chemical states. Note, however, that accurate and
# realistic estimates depend on many more factors than just the selection of potential phases, such as the *choice of
# thermodynamic models for non-ideal phases*. Furthermore, note that adding too many phases and species to the
# definition of the chemical system can result in *more computationally expensive* chemical calculations. In critical
# performance applications, for instance, when combining chemical reactions and fluid flow and species transport
# modeling, restricting the number of phases and species might be necessary for achieving feasible simulation times.
# The modeler is responsible to decide to which extent the number of phases and species can be compromised for
# efficiency reasons at the expense of chemical realism.

# ### Specifying aqueous phase
#
# To initialize [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html#details),
# we can add the list of exact names of aqueous species we wish to be simulated in computations.
#
# > **Note**: Alternative methods to specify species in aqueous phase are described in documentation of the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) class or tutorial
# [**ChemicalEditor**](cl.chemical-editor.ipynb).
#
# The default model to calculate the activities of solvent water and ionic species is the HKF model.
# Note that activity models are also needed for *neutral species* then an ideal model is used, in which their activity
# coefficients are one. For some neutral aqueous species, such as CO<sub>2</sub>(aq), we specify the Drummond model

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setActivityModelDrummondCO2()

# Alternatively, one can select ideal activity model by
# [AqueousPhase::setActivityModelIdeal](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html#aef6fbb39539c771554b628de53edeab7)
# method for neutral aqueous species by providing the corresponding name, e.g.,

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setActivityModelIdeal("CO2(aq)")

# Finally, we can set the activity model of neutral aqueous species to be the Setschenow one, where not only name of
# the species but also the Setschenow constant must be provided in method
# [AqueousPhase::setActivityModelSetschenow](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html#ab1d6cc44d10a01c7bc1380d6cedfff79).

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--", "NaCl(aq)"]) \
    .setChemicalModelDebyeHuckel() \
    .setActivityModelSetschenow("NaCl(aq)", 0.1)

# To choose, for instance, the Pitzer equation of state for aqueous species, we can use

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setChemicalModelPitzerHMW() \
    .setActivityModelDrummondCO2()

# The instance of the [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html#details) class can
# also be constructed explicitly and have its default *chemical model* changed to the Debye-Huckel model, for example:

aqueous_phase = editor.addAqueousPhaseWithElementsOf("H2O NaCl CO2")
aqueous_phase.setChemicalModelDebyeHuckel()

# Method
# [ChemicalEditor::addAqueousPhaseWithElementsOf](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html#a23e44f994b87c650a949226ddc195710)
# used above is practical in those occasions, where it would be preferable to
# just specify a few compound or substance names, *not necessarily named as in the database*, and then let
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktorobl_1_1ChemicalEditor.html) to automatically select all chemical
# species that could be formed out of the combination of those compounds.

# ### Specifying gaseous phase

# Similarly, [GaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1GaseousPhase.html) can be defined either from
# list of provided elements or from substance names that are parsed by the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktorobl_1_1ChemicalEditor.html) to generate all possible
# species that can be combined from elements used in those lists. However, for this particular example,
# the water vapor, H<sub>2</sub>O(g), and gaseous/supercritical carbon dioxide, CO<sub>2</sub>(g), suffices to
# represent the gaseous phase:

editor.addGaseousPhase(["H2O(g)", "CO2(g)"]) \
    .setChemicalModelSpycherPruessEnnis()

# Here, method [GaseousPhase::setChemicalModelSpycherPruessEnnis](https://reaktoro.org/cpp/classReaktoro_1_1FluidPhase.html#a7fc8a1bc87f24b31939ce3295cd92e8f)
# sets Spycher et al. (2003) equation of state. This model only
# supports the gaseous species H<sub>2</sub>O(g) and CO<sub>2</sub>(g). Any other species will result in a runtime
# error. Alternately, the Spycher and Reed (1988) equation of state can be set (only for H<sub>2</sub>O(g),
# CO<sub>2</sub>(g), and CH<sub>4</sub>(aq)).
# If no model is explicitly specified, the Peng-Robinson equation of state is chosen by default to calculate the
# thermodynamic and chemical properties of this GaseousPhase object.

# ### Specifying mineral phase
#
# In most geochemical modeling applications, one or more mineral phases are needed. In most cases, these mineral
# phases are *pure mineral phases*, i.e., they contain only one mineral species. If more than one minerals are
# present, they are often called *solid solutions*. Defining a pure mineral phase or a solid solution phase is
# similar to defining any other phase type. The code below demonstrates addition of pure mineral phases halite
# NaCl:

editor.addMineralPhase("Halite")

# Definition of solid solution can performed as follows:

editor.addMineralPhase(["Calcite", "Magnesite"])

# ### Creating chemical system
#
# Next, we create chemical system using the class [ChemicalSystem](
# https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html), which helps to specify system's attributes and
# properties:

system = ChemicalSystem(editor)

# Once system is set, the equilibrium problem using [EquilibriumProblem](
# https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html#details) class must be initialized by

problem = EquilibriumProblem(system)

# We use class [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html) to specify
# the conditions at which our system should be in equilibrium.
# For the equilibrium calculation, we have set temperature and pressure with optional units.
#
# > **Note**: The default values are 25 &deg;C for the temperature and 1 bar for pressure.

problem.setTemperature(60, "celsius")
problem.setPressure(300, "bar")

# Additionally, we have to add amount of the solutions used in the equilibrium calculations. For NaCl-brine,
# we mix 1 kg of water with 1 mol of salt. Plus, we take 100 g of CO<sub>2</sub>:

problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CO2", 100, "g")

# The units above can be changed, or even suppressed. If not provided, default units are used, such as K for
# temperatures, Pa for pressures, and mol for amounts. The
# [EquilibriumProblem::add](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html#a18d8ff0e8b8fc66f1eec0127057c7d54)
# method supports both amount and mass units, such as `mmol`,  `umol`, `g`, `mg`, etc.

# To provide computational representation of the state of a multiphase chemical system resulting from  equilibration
# process, class [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) must be used.
# Function [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e)
# equilibrates a chemical state instance with an equilibrium problem.

# The code below uses the definition of the equilibrium problem stored in the object `problem` to perform the
# equilibrium calculation with utility method [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e). The result of the calculation is the object `state`,
# an instance of [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) class, which is used
# to store the chemical state (i.e., the temperature, pressure, and molar amounts of all species) of the system at
# prescribed equilibrium conditions. The [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html)
# class also provides methods for querying thermodynamic properties of the system.

state = equilibrate(problem)

# **Note:** Method [equilibrate](https://reaktoro.org/cpp/namespaceReaktoro.html#af2d3b39d3e0b8f9cb5a4d9bbb06b697e)
# is not the optimal method for performing a sequence of equilibrium calculations
# (e.g., when coupling Reaktoro with other codes for simulating fluid flow, species transport, etc.). In situations,
# where many equilibrium calculations need to be performed and sufficient initial guesses are available each time
# (e.g., the equilibrium state at a previous time step serving as an initial guess for the new equilibrium
# calculation), use class [EquilibriumSolver](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html).

# The chemical state can be printed to console by `print(state)` command or saved to a file as follows:

state.output("state.txt")

# The output information contains details about the equilibrium state of the defined chemical system for the
# given equilibrium conditions, including, for example, the *amounts*, *masses*, *mole fractions*, *activities*,
# *activity coefficients*, and *chemical potentials* of the chemical species in each phase. It will also list
# properties of each phase, such as *density*, *molar volume*, *volume fraction*, as well specific properties of some
# phases (e.g., *ionic strength*, *pH*, *pe* for the aqueous phase).

# ### Options for the equilibrium calculation

# To customize options for the equilibrium calculations, class
# [EquilibriumOptions](https://reaktoro.org/cpp/structReaktoro_1_1EquilibriumOptions.html) can be used.
# For instance, [NonlinearOptions](https://reaktoro.org/cpp/structReaktoro_1_1NonlinearOptions.html) contains
# information about the nonlinear solver used.

options = EquilibriumOptions()
options.optimum.output.active = True
options.epsilon = 1e-50

# Here, we set the field of [OutputterOptions](https://reaktoro.org/cpp/structReaktoro_1_1OutputterOptions.html)
# class to be `True` to determine whether the intermediate values of equilibrium simulations must be output
# to the console.
# Then, we set the parameter *epsilon* (used for the numerical representation of a zero molar amount) to be equal to
# 1e-50. The molar amount of the *i*th species is considered zero if *n[i] < epsilon min b*, where *b* is
# the vector of element molar amounts.

# To encounter provided above options, we use function [equilibrate](
# https://reaktoro.org/cpp/namespaceReaktoro.html#a908245bfa7d236d8c556241dc87f489e)
# that accepts not only the instance of equilibrium problem but also the specified `options`.

state = equilibrate(problem, options)
