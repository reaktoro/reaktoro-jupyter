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
# ---

# #
#
from reaktoro import *

# # Analysis of the Evian water (using PHREEQC script and backend)
#
# The tutorial below aims at contracting the problem dedicated to checking the quality of the Evian water.
#
# The string below defines a PHREEQC script problem.
ex1 = r'''(
title exemple 1 Evian

SOLUTION 1
    temp      25
    pH        7.2
    pe        4
    redox     pe
    units     mg/kgw
    density   1
    C(4)      350
    Ca        80
    Cl        6.8
    K         1
    Mg        26
    N(5)      3.7
    Na        6.5
    S(6)      12.6
    Si        15
    -water    1 # kg
end
)'''

# Initialize a Phreeqc instance with the official phreeqc.dat database file
phreeqc = Phreeqc('../databases/phreeqc/phreeqc.dat')

# Execute a PHREEQC script defining a geochemical problem.
# Here this script is actually embedded into a string named `ex1`.
# However, `ex1` could also be a string containing the path to a script file.
# Method execute will automatically identify when the contents are embedded in
# the string and when the string is actually a path to a script file.
phreeqc.execute(ex1)

# Initialize a ChemicalSystem instance using the current state of the Phreeqc
# instance. This will allow the use of both PHREEQC thermodynamic data and
# PHREEQC activity models in the subsequent equilibrium calculations using
# Reaktoro's algorithms.
system = ChemicalSystem(phreeqc)

# Initialize an ChemicalState instance using the current state of the
# Phreeqc instance.
state = phreeqc.state(system)

# Output the equilibrium state calculated by PHREEQC to a file.
state.output('evian-water.txt')
