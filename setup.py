from setuptools import setup, find_packages

about = {}
with open("liquepy/__about__.py") as fp:
    exec(fp.read(), about)

setup(name='liquepy',
      version=about['__version__'],
      description='Tools for soil liquefaction analysis',
      url='https://github.com/eng-tools/liquepy',
      author='Julieth Quintero and Maxim Millen',
      author_email='mmi46@uclive.ac.nz',
      keywords='geotechnical engineering foundations settlement',
      license='MIT',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      install_requires=['sfsimodels>=0.5.6', 'geofound', 'numpy', 'eqsig'],
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
