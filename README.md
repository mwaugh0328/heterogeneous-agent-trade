# HAT: Heterogenous Agent Trade

<p align="center">
<img src="./notes/figures/micro-elasticity.png">
</p>

This repository contain code associated with the paper Heterogenous Agent Trade.

**This is a work in progress. Please let me know of any problems. Functionality may change rapidly.**

---

### What is HAT?

The HAT model is a model where aggregete trade arises from the explicit aggregation accross households. 

Households are characterized by the "standard incomplete markets model" where agents face to idiosyncratic productivity and taste shocks and have access to a risk free asset. It the model's current formulation, trade in goods follows the Armington tradition (a Ricardian version is in progress) with producers in each country producing a national variety. The twist is that households have random utility over these varieties and they make a discrete choice over the varieties to consume in addition to their savings decisions. The explicit aggregation of household-level decisions then determines aggregate trade flows, trade elasticities, and the gains from trade.

HAT models provide a foundation to systematically think about the distributional affects of trade reforms, their dynamics, and complementary policies to mitigate the downsides of globalization.  

---

### How?

The goal here is to provide code and informative notebooks to illustrate how things work. Code is in [julia](https://github.com/JuliaLang) and python (needs to be updated) with the goal of implementing things **fast** using transparent and well developed methods. Most the notation in the code tries to closely follow the paper. 

---

### Want to know more?

Much of this is a continuation of my thinking across several papers:

- [Lyon and Waugh (2019)](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/lw_quant_losses.pdf) and [Lyon and Waugh (2018) JIE](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/lw_tax.pdf) are precursors to this work.

- Waugh (2023) (an evolution of [Waugh (2019)](https://www.waugheconomics.com/uploads/2/2/5/6/22563786/waugh_consumption.pdf)) is an example as well.
