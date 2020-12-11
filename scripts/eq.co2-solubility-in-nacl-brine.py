from reaktoro import *

import numpy as np
import matplotlib.pyplot as plt


def solubility_co2(system, T, P, mNaCl):

    n0CO2g = 10.0

    problem = EquilibriumProblem(system)
    problem.setTemperature(T, "celsius")
    problem.setPressure(P, "bar")
    problem.add("H2O", 1.0, "kg")
    problem.add("CO2", n0CO2g, "mol")
    problem.add("NaCl", mNaCl, "mol")

    state = equilibrate(problem)

    return state.speciesAmount("CO2(g)")

database = Database("supcrt98.xml")

editor = ChemicalEditor(database)
editor.addAqueousPhaseWithElements("H O C Na Cl")
editor.addGaseousPhase(["CO2(g)"])
editor.addMineralPhase("Halite")

system = ChemicalSystem(editor)

T = np.arange(25.0, 90.0, 5.0)
P = 1.0
#P = 100.0

mCO2_1 = [solubility_co2(system, x, P, mNaCl=1.0) for x in T]  # [0] is needed to get the value of autodiff.real
mCO2_2 = [solubility_co2(system, x, P, mNaCl=2.0) for x in T]  # [0] is needed to get the value of autodiff.real
mCO2_3 = [solubility_co2(system, x, P, mNaCl=4.0) for x in T]  # [0] is needed to get the value of autodiff.real

fig, ax = plt.subplots()

ax.plot(T, mCO2_1, label=f"1 NaCl molal")
ax.plot(T, mCO2_2, label=f"2 NaCl molal")
ax.plot(T, mCO2_3, label=f"4 NaCl molal")

ax.legend(loc="upper right")
ax.grid(True)

ax.set(xlabel='Temperature [°C]')
ax.set(ylabel='Solubility [mol/kgw]')
ax.set(title='Solubility of CO2 in NaCl brine, P = ' + str(P) + ' bar')

fig.savefig('co2-solubility-nacl-h2o-'+ str(P) + 'bar.png')
