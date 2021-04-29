import os
from setuptools import find_packages, setup, Extension

# determining the directory containing setup.py
setup_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(setup_path, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

setup(
    # package information
    name = 'mooda-dna',
    packages = find_packages(),
    version = '0.9.5',
    description = 'A Multi-Objective algorithm for DNA Design and Assembly',
    long_description = readme,
    long_description_content_type="text/markdown",
    license = 'MIT',
    url='https://github.com/stracquadaniolab/mooda',
    keywords='',

    # author information
    author = 'Angelo Gaeta, Giovanni Stracquadanio',
    author_email = 'a.gaeta@sms.ed.ac.uk',

    # installation info and requirements
    install_requires=[
        "biopython==1.73",
        "numpy==1.16.2",
        "pandas==0.24.2",
        "PyYAML==5.4",
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

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
