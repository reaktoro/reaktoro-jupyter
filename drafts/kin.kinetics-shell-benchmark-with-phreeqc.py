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

# This string defines a PHREEQC input file. The contents of the input string
# below was taken from the official PHREEQC example named ex1.
ex1 = r'''(
DATABASE /shell-kinetics-benchmark/DatabaseLLNLbased.dat
title shell-kinetics-benchmark

SOLUTION 1-25
      temp  61
      pressure 1
      pH    7 
      units mol/kgw
      Na	  6.27 
      Cl	  6.27 charge
      Mg 0.0001
      Ca 0.0001
      C  0.0001
      Si 0.0001
      Al 0.0001
	  K  0.000000001
END

INCREMENTAL_REACTIONS true

#EQUILIBRIUM_PHASES 1-25
#Halite         0 10
#Calcite        0 10
#Dolomite-dis   0 10
#K-Feldspar     0 10
#Quartz         0 10
#Kaolinite      0 10

KINETICS 1-25
Halite        
  -m0 10 
  -parms 1.0   0 0 1 1
Calcite       
  -m0 10 
  -parms 1.0   0 0 1 1
Dolomite-dis  
  -m0 10 
  -parms 1.0   0 0 1 1
K-Feldspar    
  -m0 10 
  -parms 1.0   0 0 1 1
Quartz         
  -m0 10 
  -parms 1.0   0 0 1 1
Kaolinite 
  -m0 10 
  -parms 1.0   0 0 1 1
-rk 6

RUN_CELLS
-cells 1-25
-start_time 0 sec
-time_step 0.1 sec
)'''

# Initialize a Phreeqc instance with the official phreeqc.dat database file
# Make sure the file phreeqc.dat is located in the given path below relative to
# the path where you execute this script!
phreeqc = Phreeqc('databases/phreeqc/phreeqc.dat')

# Execute a PHREEQC script defining a geochemical problem. Here this script is
# actually embedded into a string named `ex1`. However, `ex1` could also be a
# string containing the path to a PHREEQC input file. The method `execute`
# below will automatically identify if you are providing an input string (as
# shown above) or a path string (e.g., 'path/to/input-file').
phreeqc.execute(ex1)

# Initialize a ChemicalSystem instance using the initialized Phreeqc instance.
# This will allow the use of thermodynamic data and activity models of PHREEQC
# in subsequent equilibrium calculations, which are performed by Reaktoro's
# numerical algorithms.
system = ChemicalSystem(phreeqc)

# Initialize a ChemicalState instance using the current chemical state stored
# in the Phreeqc instance `phreeqc`.
state = phreeqc.state(system)

# Output this chemical state to a file.
state.output('phreeqc-ex1-result.txt')
