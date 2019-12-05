.. image:: https://travis-ci.org/eng-tools/liquepy.svg?branch=master
   :target: https://travis-ci.org/eng-tools/liquepy
   :alt: Testing Status

.. image:: https://circleci.com/gh/eng-tools/liquepy.svg?style=svg
    :target: https://circleci.com/gh/eng-tools/liquepy

.. image:: https://img.shields.io/pypi/v/liquepy.svg
   :target: https://pypi.python.org/pypi/liquepy
   :alt: PyPi version

.. image:: https://coveralls.io/repos/github/eng-tools/liquepy/badge.svg
   :target: https://coveralls.io/github/eng-tools/liquepy

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://github.com/eng-tools/liquepy/blob/master/LICENSE
    :alt: License

.. image:: https://zenodo.org/badge/133526042.svg
    :target: https://zenodo.org/badge/latestdoi/133526042
    :alt: DOI

*******
Liquepy
*******

Python tools for solving problems related to soil liquefaction

Features
========

* Triggering:
    * Calculate liquefaction factor of safety according to Boulanger and Idriss (2014)
    * Calculate expected lateral strain and volumetric strain according to Zhang et al (2002)
* Settlement:
    * Calculate foundation settlement on liquefiable deposit according to Karamitros et al. (2013) or Bray et al. (2018)
* Element tests:
    * Calculate the dissipated energy required to liquefy
* Site response analysis
    * Perform equivalent linear site response analysis using the `pysra` package

Installation
============

.. code:: bash

    pip install liquepy

Optional modules (`sra`, `fig`, `spatial`) require large third-party dependencies and therefore do
not import unless dependencies are satisfied.

to install all dependencies for these optional modules run (example for installing `sra` dependencies)

.. code:: bash

    pip install liquepy[sra]

Contributors
============

How do I get set up?
--------------------

1. Run ``pip install -r requirements.txt``


Code suggestions
----------------

* Implementations of published liquefaction methods should be written as `calc_<property>_<first_author>[_et_al]_<date>()` for two authors include both.

* Plotting and visualisation should be not included in computation objects. Instead plotting functions or objects should receive computation objects as inputs.

Formatting
----------

* Follow `pep8 formatting standard <https://www.python.org/dev/peps/pep-0008/>`_ (except for line character limit not strictly observed)

* Docstrings written in `numpy format <https://numpydoc.readthedocs.io/en/latest/format.html>`_

* Tabulated data stored as comma separated or semi-colon separated files (not xlsx or xls)

Testing
-------

Tests are run with pytest

* Locally run: ``pytest`` on the command line.

* Tests are run on every push using travis, see the ``.travis.yml`` file

To test the docs:

1. Install the check docs package: Run ``pip install collective.checkdocs``

2. Run the check docs package and fix the errors: Run ``python setup.py checkdocs``


Deployment
----------

To deploy the package to pypi.com you need to:

 1. Push to the *pypi* branch. This executes the tests on circleci.com

 2. Create a git tag and push to github, run: ``trigger_deploy.py`` or manually:

 .. code:: bash

    git tag 0.5.2 -m "version 0.5.2"
    git push --tags origin pypi
