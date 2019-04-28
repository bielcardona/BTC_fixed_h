# Supplementary files for "On binary tree-child phylogenetic networks with a fixed number of hybrid nodes"
## By G. Cardona

If you have `pipenv` installed, the requirements can be installed by executing `pipenv install`. Notice however that `graphviz` must be installed previously (only needed to run the Jupyter notebook).
Once installed, the jupyter 
notebooks can be run with `pipenv run jupyter-notebook`, and the scripts with `pipenv run python script.py`.

The jupyter notebook can also be run without installing anything: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bielcardona/BTC_fixed_h/master?filepath=BTC_nh-Generator.ipynb)

The `sage` scripts need to be executed using SageMath (http://www.sagemath.org/)

The included files are:

* `BTC_fixed_h.py`: Python module that implements the generation of BTC networks with a fixed number of hybrid nodes.
* `BTC_nh-Generator.ipynb`: Jupyter notebook that demonstrates the sequential generation of BTC networks with fixed number of hybrid nodes using the module `BTC_fixed_h.py`.
* `formula_maker.sage`: Sage script that implements the obtention of formulas for the number of BTC networks with a fixed number of hybrid nodes, and uses these formulas to compute these numbers.
* `counts_BTC_nh.txt`: Results showing the number of BTC networks with n <= 8 leaves and h < n hybrid nodes.
* `formulas.txt`: Formulas for the number of BTC networks with h <= 7 hybrid nodes.
