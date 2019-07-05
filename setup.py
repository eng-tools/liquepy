from setuptools import setup, find_packages

about = {}
with open("liquepy/__about__.py") as fp:
    exec(fp.read(), about)

with open('README.rst') as readme_file:
    readme = readme_file.read()

# with open('HISTORY.rst') as history_file:
#     history = history_file.read()

setup(name=about['__project__'],
      version=about['__version__'],
      description='Tools for soil liquefaction analysis',
      long_description=readme + '\n\n',# + history,
      url='https://github.com/eng-tools/liquepy',
      author=about['__author__'],
      author_email='mmi46@uclive.ac.nz',
      keywords='geotechnical engineering foundations settlement',
      license='MIT',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      install_requires=['sfsimodels>=0.9.1',
                        'geofound',
                        'numpy',
                        'pyexcel',
                        'pysra>=0.3.2',
                        'eqsig>=0.5.35',
                        'geographiclib',
                        'matplotlib'],
      # List additional groups of dependencies here (e.g. development
      # dependencies). You can install these using the following syntax,
      # for example:
      # $ pip install -e .[dev,test]
      extras_require={
          'test': ['pytest'],
      },
      python_requires='>=3',
      package_data={},
      zip_safe=False)


# From python packaging guides
# versioning is a 3-part MAJOR.MINOR.MAINTENANCE numbering scheme,
# where the project author increments:

# MAJOR version when they make incompatible API changes,
# MINOR version when they add functionality in a backwards-compatible manner, and
# MAINTENANCE version when they make backwards-compatible bug fixes.
