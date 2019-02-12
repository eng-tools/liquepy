.. image:: https://travis-ci.org/eng-tools/liquepy.svg?branch=master
   :target: https://travis-ci.org/eng-tools/liquepy
   :alt: Testing Status

.. image:: https://img.shields.io/pypi/v/liquepy.svg
   :target: https://pypi.python.org/pypi/liquepy
   :alt: PyPi version

.. image:: https://coveralls.io/repos/github/eng-tools/liquepy/badge.svg
   :target: https://coveralls.io/github/eng-tools/liquepy

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://github.com/eng-tools/liquepy/blob/master/LICENSE
    :alt: License

*******
Liquepy
*******

Python tools for solving problems related to soil liquefaction

Installation
============

.. code:: bash

    pip install liquepy


Contributors
============

How do I get set up?
--------------------

1. Run ``pip install -r requirements.txt``

Testing
-------

Tests are run with pytest

* Locally run: ``pytest`` on the command line.

* Tests are run on every push using travis, see the ``.travis.yml`` file


Deployment
----------

To deploy the package to pypi.com you need to:

 1. Push to the *pypi* branch. This executes the tests on circleci.com

 2. Create a git tag and push to github, run: ``trigger_deploy.py`` or manually:

 .. code:: bash

    git tag 0.5.2 -m "version 0.5.2"
    git push --tags origin pypi
