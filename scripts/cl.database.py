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

# # Functionality of Database class
#
# In this tutorial, we provide an explanation on the functionality of the class
# [Database](https://reaktoro.org/cpp/classReaktoro_1_1Database.html)
# that provides operations to retrieve physical and thermodynamic data of chemical species.
#
# ### Importing the reaktoro Python package
#
# Using **Reaktoro** in Python requires first an import of the python package `reaktoro`:

from reaktoro import *

# From this point on, we are able to use the library components of Reaktoro (classes, methods, constants), which are
# needed to define our chemical system and chemical reaction modeling problems.
#
# > **Note:** To simplify the tutorials, we
# use `from reaktoro import *`, which imports all components of the `reaktoro` package into the default Python
# namespace. In a more complex project, this can potentially create name conflicts. Thus, for bigger applications,
# consider using `import reaktoro as rkt`, and then refer to Reaktoro’s classes and methods as `rkt.Database`,
# `rkt.ChemicalSystem`, `rkt.equilibrate`, etc.
#
# Reaktoro currently supports the following thermodynamic databases:
#
# * SUPCRT92;
# * PHREEQC; and
# * GEMS.
#
# More information can be found on the web-page with the
# [Reaktoro's Documentation](https://reaktoro.org/thermodynamic-databases.html).
#
# ### Initializing a thermodynamic database
#
# To initialize a thermodynamic database, we must provide xml-file, i.e.,

db = Database("supcrt98.xml")

# Here, we use [supcrt98.xml](https://github.com/reaktoro/reaktoro/blob/master/databases/supcrt/supcrt98.xml)
# database file generated from the original **SUPCRT92** database file slop98.dat.
#
# > **Note:** If filename does not point
# to a valid database file or the database file is not found, then a default built-in database with the same name
# will be tried. If no default built-in database exists with a given name, an exception will be thrown.
#
# ### Accessing the content of thermodynamic database
#
# It is possible to print all the aqueous species contained in the database SUPCRT92, i.e.,

print("List of all aqueous species in database SUPCRT92:")
for aqueous_species in db.aqueousSpecies():
    print(aqueous_species.name())

# Similar output can be written for the gaseous species (with exception of the function that return the list of such
# species), i.e.,

print("List of all gaseous species in database SUPCRT92:")
for gaseous_species in db.gaseousSpecies():
    print(gaseous_species.name())

# All the minerals included in the database can be accessed as well:

print("List of all minerals in database SUPCRT92:")
for minerals in db.mineralSpecies():
    print(minerals.name())

# To check if certain species is present in the database, one can use

print("Is Zn++ present in the database? ", db.containsAqueousSpecies("Zn++"))
print("Is Calcite present in the database? ", db.containsMineralSpecies("Calcite"))
print("Is CO2(g) present in the database? ", db.containsGaseousSpecies("CO2(g)"))
print("Is Zn+ present in the database? ", db.containsAqueousSpecies("Zn+"))

# Besides, we can output all, say, aqueous species containing a particular element (e.g., hydrogen):

print("List of all aqueous species in database with hydrogen:")
for species in db.aqueousSpeciesWithElements(["H"]):
    print(species.name())

# A thermodynamic database contains model parameters for the evaluation of standard thermodynamic properties of
# species and/or reactions (e.g., standard Gibbs energies, equilibrium constants). These properties can be accessed by
# fetching a particular element from the database:

species = db.aqueousSpecies("CaCl2(aq)")

# After obtaining a particular species, its properties can be obtained by calling corresponding function, e.g., we can
# output charge and dissociation of CaCl<sub>2</sub>:

print("Charge of CaCl2(aq): ", species.charge())
print("Dissociation CaCl2(aq): ", species.dissociation())

# Alternatively, one can also use more recent SUPCRT databased from 2017:

db = Database("supcrt07.xml")

# To list all the mineral species, we use:

print("List of all minerals in database SUPCRT07:")
for minerals in db.mineralSpecies():
    print(minerals.name())

# Assume that we want to study the properties of pyrrhotite (FeS). To fetch it from the database,
# access its formula, molar mass, charge, as well as elements, which Pyrrhotite, is constructed from,
# the following code it used:

# +
species = db.mineralSpecies("Pyrrhotite")

print("Formula of Pyrrhotite:", species.formula())
print(f"Molar mass of Pyrrhotite: {species.molarMass():4.4e} kg/mol")
print("Charge of Pyrrhotite:", species.charge())
print("Name : Molar mass ")
for elements in species.elements():
    print(f"{elements.name():>4} : {elements.molarMass()}")
# -

# Similar output can be done, for instance, for Siderite (FeCO3):

# +
species = db.mineralSpecies("Siderite")

print("Formula of Siderite:", species.formula())
print(f"Molar mass of Siderite: {species.molarMass():4.4e} kg/mol")
print("Charge of Siderite:", species.charge())
print("Name : Molar mass ")
for elements in species.elements():
    print(f"{elements.name():>4} : {elements.molarMass()}")
