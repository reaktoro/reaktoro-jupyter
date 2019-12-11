# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: notebooks//ipynb,scripts//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Performing a chemical equilibrium calculation
#
# This tutorial demonstrates how to use Reaktoro to perform a chemical equilibrium calculation.

# %% [markdown]
# We start by importing everything from the `reaktoro` package:

# %% pycharm={"is_executing": false}
from reaktoro import *

# %% [markdown]
# Next, we need a thermodynamic database that enable us to compute thermodynamic properties of species and reactions. 
# To initialize it, we execute

# %% pycharm={"is_executing": false}
db = Database('supcrt98.xml')

# %% [markdown]
# To indicate phases as well as species, which may potentially exist in them at equilibrium, we use 
# `ChemicalEditor`. 

# %% pycharm={"is_executing": false}
editor = ChemicalEditor(db)

# %% [markdown]
# We consider an aqueous phase with species that can be formed with given chemical element symbols:
#  

# %% pycharm={"is_executing": false}
editor.addAqueousPhaseWithElements('H O Na Cl C Ca Mg Si');

# %% [markdown]
# An automatic search for chemical species can result in a large number of species in the phase, potentially causing 
# the chemical reaction calculations to be more computationally expensive. If you are using Reaktoro in demanding  
# applications , you might want to manually specify the chemical species of each phase in your chemical system. This 
# can be achieved by providing explicit list of species  names as it is done for gaseous species below:

# %% pycharm={"is_executing": false}
editor.addGaseousPhase('H2O(g) CO2(g)');  

# %% [markdown] pycharm={"name": "#%% md\n"}
# Finally, we consider all possible pure minerals that could exist in our equilibrium calculations:

# %% pycharm={"is_executing": false}
editor.addMineralPhase('Halite');
editor.addMineralPhase('Calcite');
editor.addMineralPhase('Magnesite');
editor.addMineralPhase('Dolomite');
editor.addMineralPhase('Quartz');

# %% [markdown] pycharm={"name": "#%% md\n"}
# Next follows an important step, where we create the chemical system with the information so far collected in the `ChemicalEditor` object `editor`:

# %% pycharm={"is_executing": false}
system = ChemicalSystem(editor) 

# %% [markdown] pycharm={"name": "#%% md\n"}
# We use `EquilibriumProblem` to specify the conditions, at which our system should be in equilibrium. 

# %% pycharm={"is_executing": false}
problem = EquilibriumProblem(system) 

# %% [markdown]
# In particular, we can specify temperature and pressure as well as the initial condition for substance amounts

# %% pycharm={"is_executing": false}
problem.setTemperature(70, 'celsius');
problem.setPressure(100, 'bar');

problem.add('H2O', 1.0, 'kg');
problem.add('CO2', 2.0, 'mol');
problem.add('NaCl', 1.0, 'mol');
problem.add('CaCO3', 10.0, 'g');
problem.add('MgCO3', 5.0, 'g');
problem.add('Quartz', 1.0, 'mol');

# %% [markdown] pycharm={"name": "#%% md\n"}
# To perform a fast Gibbs energy minimization calculation and calculate the chemical equilibrium `state`, we have to 
# equilibrate the object `problem` with similarly named method:

# %% pycharm={"is_executing": false}
state = equilibrate(problem)  

# %% [markdown]
# We can access the amounts of species at equilibrium with method `speciesAmounts()`, i.e., 

# %% pycharm={"is_executing": false}
state.speciesAmounts()  

# %% [markdown]
# To print the name and the amount (in mol) of each species in the obtained chemical system, we execute the following 
# loop:

# %% pycharm={"name": "#%%\n", "is_executing": false}
for species in system.species():
    print(f'{species.name():>15} = {state.speciesAmount(species.name())}')

# %% [markdown]
# You can also output the chemical state to a file (check the results folder on the left pane of jupyter lab!)

# %% pycharm={"name": "#%%\n", "is_executing": false}
state.output('results/equilibrium-basics-state.txt')

# %% [markdown]
# Use `ChemicalProperties` to obtain all chemical and thermodynamic properties of a chemical system at some 
# equilibrium state by

# %% pycharm={"is_executing": false}
properties = ChemicalProperties(system)

# %% [markdown]
# and compute the chemical properties of the system at given state

# %% pycharm={"is_executing": false}
properties.update(state.temperature(), state.pressure(), state.speciesAmounts())  

# %% [markdown] pycharm={"name": "#%% md\n"}
# Alternatively, one can also use the shorter notation below, at the expense of creating a new `ChemicalProperties` object
#  at each call, instead of updating an existing one

# %% pycharm={"is_executing": false}
properties = state.properties()

# %% [markdown] pycharm={"name": "#%% md\n"}
# To evaluate some of the obtained properties, let's print the activities of each species. First, we get the values of 
# $\ln$ activities of the species by a call `.val`. **Note!** If you use `.ddT`, `.ddP`, or `.ddn` instead of `.val`, 
# you'll get, respectively, the T, P, n derivatives. 

# %% pycharm={"is_executing": false}
lna = properties.lnActivities().val

# %% [markdown] pycharm={"name": "#%% md\n"}
# To compute the actual activities, not the natural logo of them, we use

# %% pycharm={"name": "#%%\n", "is_executing": false}
a = numpy.exp(lna)  

for i, species in enumerate(system.species()):
    print(f'{species.name():>15} = {a[i]}')

# %% [markdown] pycharm={"name": "#%% md\n"}
# Let's create a pH function that determines the pH of the aqueous solution given the chemical properties of the system 
# (this will be soon simplified!)

# %% pycharm={"is_executing": false}
evaluate_pH = ChemicalProperty.pH(system)
pH = evaluate_pH(properties)
print(f'The pH of the aqueous phase is {pH.val}.')
print(f'Its sensitivity with respect to speciation, ∂(pH)/∂n, is:')
for i, species in enumerate(system.species()):
    print(f'{species.name():>15} = {pH.ddn[i]}')
