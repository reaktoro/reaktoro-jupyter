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

# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *

# This string defines a PHREEQC script problem.
# This problem was taken from the official PHREEQC example named ex1.
ex1 = r'''(
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
phreeqc = Phreeqc('databases/phreeqc/phreeqc.dat')

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
state.output('state-water-analysis-with-phreeqc.txt')
