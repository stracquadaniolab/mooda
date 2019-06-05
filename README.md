# mooda

Current version: 0.0.1-dev

Multi-Objective Optimisation algorithm for sequence Design and Assembly

## Authors

* Angelo Gaeta Giovanni Stracquadanio, a.gaeta@sms.ed.ac.uk

## Overview

### Step 1: Cloning the package

    git clone git@github.com:stracquadaniolab/mooda.git

### Step 2: Coding
The package follows standard python packaging convention.
**PEP8** is the coding standard.

#### Command line interface
The entry point for the command line interface is the file `cli.py`.
Commands provided by the CLI must be written in `commands.py` and imported into `cli.py` as needed.

#### Documentation
Code must be documented with `autodoc` in-line in the code and docs updated with `sphinx-autodoc`.
Before release, the *installation* and *example usage* sections in the `tutorial.rst` file must be completed.

### Step 3: Versioning
Packages should use the following semantic versioning scheme:
```
{MAJOR}.{MINOR}.{PATCH}-{RELEASE}
```
where:
- MAJOR: the version introduces incompatible changes with the previous one or new code/functioning structure.
- MINOR: the version introduces compatible changes with the previous one (e.g. new features).
- PATCH: bug fixes.
- RELEASE: package is in `dev` stage up to the point when it is released to the public, where it changes into `stable`.

Version should be tracked using `bumpversion`. You can bump all parts of the version scheme as follows:

- MAJOR: `bumpversion major`
- MINOR: `bumpversion minor`
- PATCH: `bumpversion patch`
- RELEASE: `bumpversion release`

### Step 4: Sharing changes
Changes should be commited to GIT by running:
```
git push
git push --tags
```

