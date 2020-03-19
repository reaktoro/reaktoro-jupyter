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


# Tutorial

Find below a list of tutorials explaining how to use Reaktoro for solving
chemical equilibrium, chemical kinetics, and reactive transport problems.

> If you need assistance, please don't hesitate to chat with us on
[Gitter](https://gitter.im/reaktoro/community)!

### Tutorials on performing chemical equilibrium calculations:

* [Basics of equilibrium calculation](eq.equilibrium-basics.ipynb)
* [Equilibrium calculation of carbonate species](eq.equilibrium-carbonates.ipynb)
* [Equilibrium calculations using equilibrium solver](eq.co2-brine-using-equilibrium-solver.ipynb)
* [Custom activity model for equilibrium calculations](eq.custom-activity-models.ipynb)
* [Inverse chemical equilibrium calculation](eq.inverse-equilibrium.ipynb)
* [Equilibrium path calculation](eq.equilibriumpath.ipynb)

### Tutorials on performing chemical kinetics calculations:

* [Dissolution of calcite in an acidic HCl-solution](kin.calcite-hcl.ipynb)
* [Dissolution of carbonate minerals in a carbon dioxide saturated brine](kin.carbonates-co2.ipynb)

### Tutorials on performing reactive transport calculations:

* [Reactive transport of CO<sub>2</sub>-saturated brine along a porous rock column (using transport solver embedded
 into Reaktoro)](rt.calcite-brine.ipynb)
* [Reactive transport modeling along a rock core after injection of the fluid-rock composition](rt.calcite-dolomite.ipynb)
* [Coupling Reaktoro into other reactive transport codes](rt.coupling-reaktoro-to-transport.ipynb)

### Tutorials on explaining the most important classes in Reaktoro:

For a comprehensive documentation of all Reaktoro’s classes and methods, please check [Reaktoro’s C++ API programming reference](https://reaktoro.org/cpp/index.html).
Below you find a more interactive documentation of a few selected classes.

* [**Database**](cl.database.ipynb)
* [**ChemicalEditor**](cl.chemical-editor.ipynb)
* [**ChemicalSystem**](cl.chemical-system.ipynb)
* [**EquilibriumSolver**](cl.equilibrium-solver.ipynb)
* [**Thermo**](cl.thermo.ipynb)

