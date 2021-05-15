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

# # Analysis of the Evian water
#
# The tutorial considers the problem dedicated to checking the quality of the Evian water using
# Reaktoro with custom-defined chemical composition and the PHREEQC backend.
#
# ![title](../images/evian-chemical-water-composition.png)
#
# ## Using Reaktoro backend
#
# First, we use purely Reaktoro functionality to define the chemical state . We start from setting a chemical system.

# +
from reaktoro import *
from math import *

db = Database("supcrt98.xml")

editor = ChemicalEditor(db)

editor.addAqueousPhaseWithElements("C Ca Cl H K Mg N Na O S Si")
editor.addMineralPhase("Calcite")
editor.addMineralPhase("Dolomite")
editor.addMineralPhase("Quartz")

system = ChemicalSystem(editor)
# -

# To evaluate the saturation indices of the mineral phases later, we define the corresponding reactions:

# +
# Define reaction equations for each phase
reaction_equation_calcite = ReactionEquation("Calcite + H+ = Ca++ + HCO3-")
reaction_equation_dolomite = ReactionEquation("Dolomite + 2*H+ = Ca++ + Mg++ + 2*HCO3-")
reaction_equation_quartz = ReactionEquation("Quartz + H2O(l) = H+ + HSiO3-")

# Define chemical reactions based on the reaction equation and chemical system
reaction_calcite = Reaction(reaction_equation_calcite, system)
reaction_dolomite = Reaction(reaction_equation_dolomite, system)
reaction_quartz = Reaction(reaction_equation_quartz, system)
# -

# To define the equilibrium problem with fixed pH, we use
# [EquilibriumInverseProblem](https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumInverseProblem.html) class.

# +
problem = EquilibriumInverseProblem(system)

problem.setTemperature(25, "celsius")
problem.setPressure(1.01325, "bar") # 1 atm
problem.add("H2O", 1.0, "kg")
# -

# Since the concentration of the species in Evian water (see picture above) are defined in milligrams per liter (mg/L),
# we need to convert these values to moles using species molar masses. The latter can be obtained using function
# `system.species("H+").molarMass()`. Also, to convert `mg` into `kg`, we multiply the values by 1e-6.

# Calcium, Ca++
problem.add("Ca++", 80 * 1e-6 / system.species("Ca++").molarMass(), "mol")  # 80 * 1e-6 kg / (MW(Ca++) kg / mol)
# Chloride, Cl-
problem.add("Cl-", 6.8 * 1e-6 / system.species("Cl-").molarMass(), "mol")  # 6.8 * 1e-6 kg / (MW(Cl-) kg / mol)
# Bicarbonate, HCO3-
problem.add("HCO3-", 350 * 1e-6 / system.species("HCO3-").molarMass(), "mol")  # 350 * 1e-6 kg / (MW(HCO3-) kg / mol)
# Magnesium, Mg++
problem.add("Mg++", 26 * 1e-6 / system.species("Mg++").molarMass(), "mol")   # 26 * 1e-6 kg / (MW(Mg++) kg / mol)
# Nitrate, NO3-
problem.add("NO3-", 3.7 * 1e-6 / system.species("NO3-").molarMass(), "mol")  # 3.7 * 1e-6 kg / (MW(NO3-) kg / mol)
# Potassium, K
problem.add("K+", 1 * 1e-6 / system.species("K+").molarMass(), "mol")  # 1 * 1e-6 kg / (MW(K+) kg / mol)
# Sodium, Na+
problem.add("Na+", 6.5 * 1e-6 / system.species("Na+").molarMass(), "mol")  # 6.5 * 1e-6 kg / (MW(Na+) kg / mol)
# Sulfates, SO4
problem.add("SO4--", 12.6 * 1e-6 / system.species("S2O4--").molarMass(), "mol")  # 12.6 * 1e-6 kg / (MW(SO4--) kg / mol)

# We also add silica amount manually, using molar weight 28 g/mol:

# Silica, Si
problem.add("Si", 15 * 1e-3 / 28.0855, "mol")  # 15 * 1e-3 g / (28.0855 g) * mol, where Molar Weight of Si = 28.0855

# Initial calcite, dolomite, and quartz phases are set to zero and pH is set to 7.2 according to the picture above.

problem.add("Calcite", 0.0, "mol")
problem.add("Dolomite", 0.0, "mol")
problem.add("Quartz", 0.0, "mol")
problem.pH(7.2, "HCl", "NaOH")

# Finally, we equilibrate the above-defined problem.

state = equilibrate(problem)
print(state)

# In order to obtain saturation indices of the carbonates and quartz, we need to access chemical properties of the
# calculated chemical state. The saturation index is defined as a log10 of the ratio of equilibrium constant and
# reaction quotient. It is 0 for minerals that are precipitated (i.e., in equilibrium with the solution), SI > 0 for
# supersaturated minerals, and SI < 0 for undersaturated minerals.

# +
# Calculate calcite's saturation index

# Fetch chemical properties
props = state.properties()
# Calculate equilibrium constant
lnK_calcite = reaction_calcite.lnEquilibriumConstant(props)
# Calculate reaction quotient
lnQ_calcite = reaction_calcite.lnReactionQuotient(props)
# Calculate saturation ratio
lnSR_calcite = lnQ_calcite.val - lnK_calcite.val
SI_calcite = lnSR_calcite / log(10)

# Calculate dolomite's saturation index
lnK_dolomite = reaction_dolomite.lnEquilibriumConstant(props)
lnQ_dolomite = reaction_dolomite.lnReactionQuotient(props)
lnSR_dolomite = lnQ_dolomite.val - lnK_dolomite.val
SI_dolomite = lnSR_dolomite / log(10)

# Calculate quartz's saturation index
lnK_quartz = reaction_quartz.lnEquilibriumConstant(props)
lnQ_quartz = reaction_quartz.lnReactionQuotient(props)
lnSR_quartz = lnQ_quartz.val - lnK_quartz.val
SI_quartz = lnSR_quartz / log(10)

print("Saturation Index (Calcite) = ", SI_calcite)
print("Saturation Index (Dolomite) = ", SI_dolomite)
print("Saturation Index (Quartz) = ", SI_quartz)
# -

# Based on the obtained results, we see that the water is saturated with dolomite and quartz. Whereas calcite is
# undersaturated. Similarly, we can analyze the stability indices of the calcite,
# dolomite, and quartz phases. It is 0 if the phase is stable, bigger than 0 if the phase is supersaturated, and less
# than 0 if the phase is undersaturated. Stability indices confirm the results we have obtained for saturation indices.

# +
# Fetch the list of stability indices
phase_stability_indices = state.phaseStabilityIndices()

# Get the index of calcite, dolomite, and quartz
calcite_phase_index = system.indexPhase("Calcite")
dolomite_phase_index = system.indexPhase("Dolomite")
quartz_phase_index = system.indexPhase("Quartz")

# Output
print("Stability Index (Calcite)  = ", phase_stability_indices[calcite_phase_index])
print("Stability Index (Dolomite) = ", phase_stability_indices[dolomite_phase_index])
print("Stability Index (Quartz)   = ", phase_stability_indices[quartz_phase_index])
# -

# ## Using PHREEQC backend
#
# Similarly, PHREEQC simulations can be executed using the PHREEQC backend. For this, we first define the problem using
# the following script:

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
EQUILIBRIUM_PHASES 1
    Dolomite 0 0
    Calcite 0 0
    Quartz 0 0
end
)'''

# Next, we initialize a Phreeqc instance with the official phreeqc.dat database file:

phreeqc = Phreeqc('../databases/phreeqc/phreeqc.dat')

# To define a geochemical problem, we execute a PHREEQC script `ex1`. Here, `ex1` could also be a string containing
# the path to a script file. The method `execute()` automatically identifies if the content is embedded in the string or
# if the string is a path to a script file.

phreeqc.execute(ex1)

# Next, we initialize a [ChemicalSystem](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html) instance using
# the current state of the Phreeqc instance `phreeqc`. It will allow using both PHREEQC thermodynamic data and
# activity models in the subsequent equilibrium calculations using Reaktoro's algorithms.

system = ChemicalSystem(phreeqc)

# Finally, to initialize a [ChemicalState](https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html) instance,
# using the current state of `phreeqc`.

state = phreeqc.state(system)

# The final equilibrium state calculated by PHREEQC is output into the file `state-water-analysis-with-phreeqc.txt` and
# printed out to the console:

print(state)
