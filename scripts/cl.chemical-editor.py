# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
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

# # Functionality of ChemicalEditor class
#
# In this tutorial, we provide clarification of functionality of class
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) that is used to
# conveniently create chemical and reaction systems.
#
# ### Importing the reaktoro Python package
#
# Using **Reaktoro** in Python requires an import of the python package `reaktoro`:

from reaktoro import *

# The default thermodynamic databases embedded into Reaktoro is SUPCRT92, so you do not have to initialize the
# database `db = Database('supcrt98.xml')`, unless an alternative database must be used.

# ### Initializing chemical editor
#
# Class [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html)
# provides convenient operations to initialize
# [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) and
# [ReactionSystem](https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html) instances.

# Define the editor of the chemical system from the default database SUPCRT92
editor = ChemicalEditor()

# Alternatively, the editor can be initialized by the database instance:

# Initialize a thermodynamic database with supcrt98.xml
db = Database("supcrt98.xml")
# Define the editor of the chemical system
editor = ChemicalEditor(db)

# ### Preparation of chemical system definition
#
# Before definition of chemical system, aqueous, gaseous, and mineral phases must be added. It can be done in various
# ways. Let us consider, first, definition of aqueous species:
#
# * With method `addAqueousPhase()`, the [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html)
# can be created by specifying the names of the species one by
# one. These species names must conform to those used in the database that was specified during the initialization of
# the [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object, otherwise, an exception
# will be thrown.

editor.addAqueousPhase(
    ["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO3--", "CO2(aq)"]
)

# * Alternatively, instead of listing the names of the species one by one, which might require prior knowledge of the
# species names in the database, we can use method `addAqueousPhaseWithElements()`. It permits the `AqueousPhase`
# object to be constructed by using a list of chemical element names. The database will be searched for all
# species that could be formed out of those elements. These species will then be used to construct the `AqueousPhase`
# object.

editor.addAqueousPhaseWithElements("H O C Ca Cl Mg")

# * Finally, [AqueousPhase](https://reaktoro.org/cpp/classReaktoro_1_1AqueousPhase.html)
# object can be also constructed by using a list of compound or substance names that might not
# necessarily represent names of species in the database. The list of compounds will be broken into a list of element
# names, and the database will be similarly searched for all species that could be formed out of those elements.

editor.addAqueousPhaseWithElementsOf("H2O NaCl CO2")

# > **Note:** The call of `addAqueousPhaseWithElements()` will result in larger chemical system than call of
# > `addAqueousPhase()`, where the specific aqueous species are provided.  In more demanding applications (e.g., as a
# > chemical solver in a reactive transport simulator), you might want to manually specify the chemical species of each
# > phase in your chemical system.

# To add gaseous phase, similar methods can be used:

editor.addGaseousPhase(["H2O(g)", "CO2(g)", "H2(g)", "O2(g)", "CH4(g)"])
editor.addGaseousPhaseWithElements(["H", "O", "C"])
editor.addGaseousPhaseWithElementsOf(["H2O", "CO2"])

# The `MineralPhase` object is created by specifying the names of the species one by one. These species names must
# conform to those used in the database that was specified during the initialization of the
# [ChemicalEditor](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html) object,
# otherwise, an exception will be thrown. The example below describes the usage of this method for the creation of two
# pure mineral phases and one solid solution with two mineral species.

editor.addMineralPhase("Calcite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase(["Dolomite", "Calcite"])

# There two more alternatives to add [MinerialPhase](https://reaktoro.org/cpp/classReaktoro_1_1MineralPhase.html), i.e.,

editor.addMineralPhaseWithElements(["Ca", "C", "O"])
editor.addMineralPhaseWithElementsOf(["CaCO3", "MgCO3"])

# Besides phases, one could also set temperatures and pressures of considered chemical system for constructing
# interpolation tables of thermodynamic properties., e.g.,

editor.setTemperatures([60, 80, 100, 120, 140, 160], "celsius")
editor.setPressures([1, 10, 100], "bar")

# ### Chemical system definition
#
# Finally, definition of chemical system is done by calling

# Construct the chemical system
system = ChemicalSystem(editor)
