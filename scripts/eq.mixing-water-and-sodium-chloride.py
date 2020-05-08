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

# # Mixing of the table salt with water
#
# Importing the 'reaktoro' package:

from reaktoro import *

# Initialize a thermodynamic database:

db = Database("supcrt98.xml")

# Initializing chemical editor:

editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H O Na Cl")
editor.addMineralPhase("Halite") # sodium chloride (solid)

# Initializing chemical system:

system = ChemicalSystem(editor)

# To evaluate the content of the system, we type:

print("system = ", system)

# We see that all the discussed during the lecture species are there, i.e.,
# H<sub>2</sub>O, H<sup>+</sup>. OH<sup>-</sup>, O<sub>2</sub>(aq),
# H<sub>2</sub>(aq), Na<sup>+</sup>, Cl<sup>-</sup>, NaCl(aq), HCl(aq), NaOH(aq) among others.

# Initializing chemical problem:

problem = EquilibriumProblem(system)
problem.setTemperature(25, "celsius")
problem.setPressure(1, "bar")
problem.add("H2O", 1.0, "kg") # water
problem.add("NaCl", 1.0, "mg") # sodium chloride / table salt

# After mixing water and salt, the following reactions start:
#
# \begin{alignat}{4}
# {\rm NaCl(s)} & \rightleftharpoons- {\rm Na}^+ + {\rm Cl}^- &\qquad& (1) \\
# {\rm Na}^+ + {\rm Cl}^- &  \rightleftharpoons {\rm NaCl(aq)} &\qquad& (2) \\
# {\rm H_2O} & \rightleftharpoons {\rm H}^+ + {\rm 0H}^- &\qquad& (3) \\
# \end{alignat}
#
# Reaction (1) corresponds to dissolution of solid sodium chloride provides ions of sodium and chloride.
# In (2), ions of sodium and chloride combine aqueous sodium chloride, NaCl(aq).
# In reaction (3), water provides solvent species for other solutions.
# Another species that can be forming are: HCl(aq), NaOH(aq) by the following reactions:
# $${\rm H}^{+} + {\rm Cl}^- \rightleftharpoons {\rm HCl(aq)}$$
#
# Let us equilibrate the problem:

state = equilibrate(problem)

# To evaluate the result of equilibration, we use `print` function to output the `state`
# to the console

print("state = ", state)

# We see that the Halite mineral has practically zero amount and `under stable` state
# in the resulting mixture, which indicates that it was completely dissolved in water.
# Let us add some more salt

problem.add("NaCl", 1.0, "kg") # additional 1kg of salt
state = equilibrate(problem)
print("state = ", state)

# To evaluate the amount of a single (or a particular species'), one can use:

print(f"Amount of NaCl(s) : {state.speciesAmount('Halite')}")
