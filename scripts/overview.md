---
jupyter:
  jupytext:
    formats: notebooks//ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.3.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<div>
<a href="https://reaktoro.org"><img src="https://reaktoro.org/_images/reaktoro-header.svg"></a>
</div>


## Tutorial Overview


The tutorial is broken into several sections, which are each presented in their own notebook.
Tutorials, explaining the functionality of the most important classes in Reaktoro:

1.  [Database class](cl.database.ipynb)
2.  [ChemicalEditor class](cl.chemical-editor.ipynb)
3.  [ChemicalSystem class](cl.chemical-system.ipynb)
4.  [EquilibriumSolver class](cl.equilibrium-solver.ipynb)
5.  [Thermo class](cl.thermo.ipynb)

Tutorials on functionality related to modeling equilibrium problems:

1.  [Basics of equilibrium calculation](eq.1.equilibrium-basics.ipynb)
2.  [Equilibrium calculation of carbonate species](eq.2.equilibrium-carbonates.ipynb)
3.  [Equilibrium calculation till full evaporation of the water](eq.3.co2-brine-full-water-evaporation.ipynb)
4.  [Equilibrium calculations using equilibrium solver](eq.4.co2-brine-using-equilibrium-solver.ipynb)
5.  [Custom activity model for equilibrium calculations](eq.5.custom-activity-models.ipynb)
6.  [Inverse chemical equilibrium calculation](eq.6.inverse-chemical-equilibrium-calculations.ipynb)
7.  [Equilibrium path calculation](eq.7.equilibriumpath.ipynb)

Tutorials on modelling chemical paths of kinetically controlled reactions:

1.  [Dissolution of calcite in an acidic HCl-solution](kin.1.calcite-hcl.ipynb)
2.  [Dissolution of carbonate minerals in a carbon dioxide saturated brine](kin.2.carbonates-co2.ipynb)

Finally, the last group is tutorials on using Reaktoro in reactive-transport simulations:

1. [Reactive transport of CO<sub>2</sub>-saturated brine along a porous rock column](rt.1.calcite-brine.ipynb)
2. [Reactive transport modeling along a rock core after injection of the fluid-rock composition](rt.2.calcite-dolomite.ipynb)
3. [Coupling Reaktoro into other reactive transport codes](rt.3.coupling-reaktoro-to-transport.ipynb)
