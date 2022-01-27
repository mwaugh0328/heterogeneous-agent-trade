# HAT: Heterogenous Agent Trade Models

This repository contain codes code to teach and compute HAT-models.

### What is HAT?

HAT models are comprised of a heterogenous agent block within each country and a neoclassical trade block where demand and supply of goods across countries must equilibrate.

- The heterogenous agent block is characterized by variants on the "standard incomplete markets model" where agents face to idiosyncratic productivity shocks and have access to a risk free asset and the ability vary labor supply to smooth these fluctuations. Aggregate demand and supply in this economy are endogenous (i.e. not pinned down by the production function) and are determined by the micro-level responses of these agents.

- The supply block is characterized by variants on neoclassical trade models like Armington or Eaton and Korum (2002). In general, these models take demand and supply of goods in each country and determine the pattern of trade. In standard trade models, however, both domestic demand and domestic supply are technologically determined unlike in the HA block.

- HAT combines the two with the HA block and the trade block interacting.

---

### How?

The goal here is to provide code and informative notebooks to illustrate how things work. Code will be in [julia](https://github.com/JuliaLang) and (maybe python) with the goal of implementing things **fast** using transparent and well developed methods.

---

### Want to know more?

- Star and keep watching this repository.

- [Lyon and Waugh (2019)](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/lw_quant_losses.pdf) and [Lyon and Waugh (2018) JIE](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/lw_tax.pdf) are precursors to this work.

- Waugh (2022) (an evolution of [Waugh (2019)](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/waugh_consumption.pdf)) is an example as well.
