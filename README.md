# CTGoMartini-tests

This repository contains test cases for the CTGoMartini package, which appears to be a molecular dynamics simulation framework that implements the Martini coarse-grained force field with various extensions.

## Project Structure

The test suite is organized into the following categories:

### API Tests (`tests/api/`)

These tests validate the core API functionality of CTGoMartini:

- **Classic Martini Tests**: Verify the implementation of the standard Martini force field by comparing OpenMM and GROMACS simulation results for energy and forces.
- **Multiple Basin Go-Martini Tests**: Test the multiple basin potential implementation that allows simulation of conformational changes between different states.
- **Contacts Tests**: Validate contact map implementations.
- **Energy Item Comparison**: Compare energy terms between different implementations.

### Functional Tests (`tests/func/`)

These tests focus on specific functional components:

- **WriteItp Tests**: Verify the generation of GROMACS topology (.itp) files.
- **ConvertLongShortElasticBonds Tests**: Test the conversion between different elastic network representations.

### Run Tests (`tests/run/`)

These tests validate the execution of simulation workflows:

- **GenMBItp Tests**: Test the generation of multiple basin topology files.
- **RunMBGoMartini Tests**: Verify the execution of Multiple Basin Go-Martini simulations.

## Data Directory

The `tests/data/` directory contains reference data for various test cases, organized by test category.

## Running Tests

To run all tests:

```bash
python -m tests.tests
```

This will execute the pytest suite on all test files in the project.

To run specific test categories:

```bash
python -m pytest tests/api/  # Run only API tests
python -m pytest tests/func/  # Run only functional tests
python -m pytest tests/run/   # Run only run tests
```

## Requirements

The test suite requires:

- Python 3.12
- CTGoMartini
- pytest
