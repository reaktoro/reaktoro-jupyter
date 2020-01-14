# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: py:light,../notebooks//ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Performing a chemical equilibrium calculation with customized activity models
#
# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation with customized activity
# models. We use an example of equilibration of  $\mathrm{H_2O–NaCl–CO_2}$ system when 1 kg of $\mathrm{H_2O}$,
# 100 g of $\mathrm{CO_2}$, and 1 mol of $\mathrm{NaCl}$ are mixed at temperature 60 $^\circ$C and pressure 300 bar.
# First, we import everything from the `reaktoro` package by

from reaktoro import *

# To indicate phases and corresponding to them species, as well as models used to evaluate activities of the aqueous
# species, we use `ChemicalEditor`. Here, the default database SUPCRT98 is used.

editor = ChemicalEditor()

# ## Specifying different phases and corresponding models
#
# Specifying the phases and their species is not enough to fully describe a chemical system in the *computational
# sense*. Every phase in Reaktoro has two associated models: a *thermodynamic model* and a *chemical model*. These
# denominations are not standard in the literature, but they are useful in the differentiation of two needed types
# of models for a phase.

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

# To initialize [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html#details),
# we add the list of exact names of aqueous species we wish to be simulated in computations.
# The default model to calculate the activities of solvent water and ionic species is the HKF model.
# Note that activity models are also needed for *neutral species* then an ideal model is used, in which their activity
# coefficients are one. For some neutral aqueous species, such as $\mathrm{CO_2(aq)}$, we specify the Drummond model

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setActivityModelDrummondCO2()

# Alternatively, one can select ideal activity model by `setActivityModelIdeal()` method for neutral aqueous species by
# providing the corresponding name, e.g.,

editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"]) \
    .setActivityModelIdeal("CO2(aq)")

# Finally, we can set the activity model of neutral aqueous species to be the Setschenow one, where not only name of
# the species but also the Setschenow constant must be provided in method `setActivityModelSetschenow()`.

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

# Method `addAqueousPhaseWithElementsOf()` used above is practical in those occasions, where it would be preferable to
# just specify a few compound or substance names, *not necessarily named as in the database*, and then let
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktorobl_1_1ChemicalEditor.html) to automatically select all chemical
# species that could be formed out of the combination of those compounds.

# ### Specifying gaseous phase

# Similarly, [GaseousPhase](https://reaktoro.org/cpp/classReaktoro_1_1GaseousPhase.html) can be defined either from
# list of provided elements or from substance names that are parsed by the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktorobl_1_1ChemicalEditor.html) to generate all possible
# species that can be combined from elements used in those lists. However, for this particular example,
# the water vapor, $\mathrm{H_2O(g)}$, and gaseous/supercritical carbon dioxide, $\mathrm{CO_2(g)}$, suffices to
# represent the gaseous phase:

editor.addGaseousPhase(["H2O(g)", "CO2(g)"]) \
    .setChemicalModelSpycherPruessEnnis()

# Here, method `setChemicalModelSpycherPruessEnnis()` sets Spycher et al. (2003) equation of state. This model only
# supports the gaseous species $\mathrm{H_2O(g)}$ and $\mathrm{CO_2(g)}$. Any other species will result in a runtime
# error. Alternately, the Spycher and Reed (1988) equation of state can be set (only for $\mathrm{H_2O(g)}$,
# $\mathrm{CO_2(g)}$, and $\mathrm{CH_4(g)}$).
# If no model is explicitly specified, the Peng-Robinson equation of state is chosen by default to calculate the
# thermodynamic and chemical properties of this GaseousPhase object.

# ### Specifying mineral phase
#
# In most geochemical modeling applications, one or more mineral phases are needed. In most cases, these mineral
# phases are *pure mineral phases*, i.e., they contain only one mineral species. If more than one minerals are
# present, they are often called *solid solutions*. Defining a pure mineral phase or a solid solution phase is
# similar to defining any other phase type. The code below demonstrates addition of pure mineral phases halite
# $\mathrm{NaCl}$:

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
# https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html#details) class must be initialized:

problem = EquilibriumProblem(system)

# We use class [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html) to specify
# the conditions at which our system should be in equilibrium.
# For the equilibrium calculation, we have set temperature and pressure with optional units.
#
# > **Note**: The default values are 25 $^\circ$C for the temperature and 1 bar for pressure.

problem.setTemperature(60, "celsius")
problem.setPressure(300, "bar")

# Additionally, we have to add amount of the solutions used in the equilibrium calculations. For $\mathrm{
# NaCl}$-brine, we mix 1 kg of water with 1 mol of salt. Plus, we take 100 g of $\mathrm{CO_2}$:

problem.add("H2O", 1, "kg")
problem.add("NaCl", 1, "mol")
problem.add("CO2", 100, "g")

# The units above can be changed, or even suppressed. If not provided, default units are used, such as K for
# temperatures, Pa for pressures, and mol for amounts. The `add` method in
# [EquilibriumProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html) supports both amount and
# mass units, such as `mmol`,  `umol`, `g`, `mg`, etc.

# To provide computational representation of the state of a multiphase chemical system resulting from  equilibration
# process, class [ChemicalState]() must be used. Function `equilibrate()` equilibrates a chemical state instance with
# an equilibrium problem.

# The code below uses the definition of the equilibrium problem stored in the object `problem` to perform the
# equilibrium calculation with utility method `equilibrate()`. The result of the calculation is the object `state`,
# an instance of [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) class, which is used
# to store the chemical state (i.e., the temperature, pressure, and molar amounts of all species) of the system at
# prescribed equilibrium conditions. The [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html)
# class also provides methods for querying thermodynamic properties of the system.

state = equilibrate(problem)

# **Note:** Method `equilibrate` is not the optimal method for performing a sequence of equilibrium calculations
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

