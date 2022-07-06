---
jupyter:
  jupytext:
    formats: notebooks//ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
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
* [Solubility of sodium chloride salt in water](eq.eq.sodium-chloride-solubility-in-water.ipynb)
* [Solubility of carbon dioxide in water](eq.co2-solubility-in-water.ipynb)
* [Analysis of the Evian water (with Reaktoro and PHREEQC backend)](eq.evian-water-analysis.ipynb)
* [Solubility of sodium chloride salt in water (using Setschenow activity model)](eq.equilibrium-setschenow-activity-model.ipynb)
* [Carbon dioxide gas solubility in the NaCl-brine](eq.co2-solubility-in-nacl-brine.ipynb)
* [Dependence of the pH on the CO2(g) amount in seawater](eq.ph-dependence-on-co2-in-seawater.ipynb)
* [Dependence of the pH on added contaminant in water](eq.ph-dependence-on-contaminants-in-water.ipynb)
* [Change of the carbon species in the reaction path vs pH](eq.equilibriumpath-carbon-species-vs-ph.ipynb)
* [Calcite solubility in water and CO2-saturated rainwater](eq.calcite-solubility.py)
* [Gypsum/anhydrite solubility in water](eq.anhydrite-solubility.ipynb)

### Tutorials on performing chemical kinetics calculations:

* [Dissolution of calcite in an acidic HCl-solution](kin.calcite-hcl.ipynb)
* [Dissolution of carbonate minerals in a carbon dioxide saturated brine](kin.carbonates-co2.ipynb)
* [Precipitation of barite as a result or water-flooding](kin.barite-precipitation.ipynb)

### Tutorials on performing reactive transport calculations:

* [Reactive transport of CO<sub>2</sub>-saturated brine along a porous rock column (using transport solver embedded
 into Reaktoro)](rt.calcite-brine.ipynb)
* [Reactive transport modeling along a rock core after injection of the fluid-rock composition](rt.calcite-dolomite.ipynb)
* [Coupling Reaktoro into other reactive transport codes](rt.coupling-reaktoro-to-transport.ipynb)
* [Reactive transport modeling of the H<sub>2</sub>S scavenging process along a rock core](rt.scavenging.ipynb)
* [Reactive transport modeling of the H<sub>2</sub>S scavenging process along a rock core (with hematite and pyrite)](rt.scavenging-with-hematite-and-pyrite.ipynb)
* [Scaling of barite as a result of the water-flooding](rt.scaling.ipynb)
* [Reactive transport in granite simulation](rt.granite.ipynb)

### Tutorials on basic properties related to chemical equilibrium calculations:

* [Calculation of ionic strength and activity coefficients of aqueous species](eq.ionic-strength-and-activity-coefficients.ipynb)
* [Mass balance and mass action equations and related chemical properties](eq.mass-balance-mass-action.ipynb)

### Tutorials on explaining the most important classes in Reaktoro:

For a comprehensive documentation of all Reaktoro’s classes and methods, please check [Reaktoro’s C++ API programming reference](https://reaktoro.org/cpp/index.html).
Below you find a more interactive documentation of a few selected classes.

* [**Database**](cl.database.ipynb)
* [**ChemicalEditor**](cl.chemical-editor.ipynb)
* [**ChemicalSystem**](cl.chemical-system.ipynb)
* [**EquilibriumSolver**](cl.equilibrium-solver.ipynb)
* [**Thermo**](cl.thermo.ipynb)
