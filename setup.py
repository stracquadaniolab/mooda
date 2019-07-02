import os
from setuptools import find_packages, setup

# determining the directory containing setup.py
setup_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(setup_path, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

setup(
    # package information
    name = 'mooda',
    packages = find_packages(),
    version = '0.0.1-dev',
    description = 'Multi-Objective Optimisation algorithm for sequence Design and Assembly',
    long_description = readme,
    license = 'MIT',
    url='git@github.com:stracquadaniolab/mooda.git',
    keywords='',

    # author information
    author = 'Angelo Gaeta Giovanni Stracquadanio',
    author_email = 'a.gaeta@sms.ed.ac.uk',

    # installation info and requirements
    install_requires=[
        "biopython == 1.73",
        "numpy == 1.16.2",
        "pandas == 0.24.2",
        "PyYAML == 5.1",
    ],
    setup_requires=[],

    # test info and requirements
    test_suite='tests',
    tests_require=[

    ],

    # package deployment info
    include_package_data=True,
    zip_safe=False,

    # all tools have cli interface
    entry_points={
        'console_scripts': [
            'mooda=mooda.cli:main',
        ],
    },
)
