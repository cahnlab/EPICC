# Unit Tests for Direct Methylation (dmC) Feature

This directory contains unit tests for the direct methylation (dmC) analysis feature added to the EpigeneticButton pipeline. The dmC feature enables processing of native base modification data from long-read sequencing platforms (Oxford Nanopore, PacBio) without bisulfite conversion.

## Overview

The tests cover:

1. **validate_dmc_input.py** - Validation of dmC input files (modBAM with MM/ML tags and bedMethyl format)
2. **mC.smk helpers** - Helper functions for sample type detection, input format determination, and parameter selection

## Test Files

- `test_validate_dmc_input.py` - Tests for input validation script
- `test_mC_helpers.py` - Tests for Snakemake helper functions

## Prerequisites

Install test dependencies:

```bash
# Activate your conda environment
conda activate epicc

# Install pytest and dependencies
pip install pytest pytest-cov pytest-mock
```

## Running Tests

### Run All Unit Tests

```bash
# From the repository root
pytest tests/unit/ -v

# Or from the tests directory
cd tests
pytest unit/ -v
```

### Run Specific Test Files

```bash
# Test only the validation script
pytest tests/unit/test_validate_dmc_input.py -v

# Test only the helper functions
pytest tests/unit/test_mC_helpers.py -v
```

### Run Specific Test Classes or Functions

```bash
# Run all tests in a class
pytest tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl -v

# Run a single test function
pytest tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl::test_valid_bedmethyl_11_columns -v
```

### Run Tests with Coverage

```bash
# Generate coverage report
pytest tests/unit/ --cov=workflow/scripts --cov-report=html --cov-report=term

# View HTML coverage report
firefox htmlcov/index.html  # or your preferred browser
```

### Run Tests in Parallel

```bash
# Install pytest-xdist
pip install pytest-xdist

# Run tests in parallel (faster)
pytest tests/unit/ -n auto
```

## Test Structure

### test_validate_dmc_input.py

Tests are organized into classes by function:

- **TestValidateBedMethyl** - Tests for bedMethyl file validation
  - Valid format tests (11 columns, 10 columns, with headers)
  - Invalid format tests (wrong columns, bad coordinates, invalid values)
  - Edge cases (empty files, gzipped files, chromosome validation)

- **TestValidateModBAM** - Tests for modBAM file validation
  - Valid modBAM tests (with MM/ML tags)
  - Invalid modBAM tests (missing tags, empty files)
  - Reference alignment validation tests

- **TestValidationMain** - Tests for CLI behavior
  - Command-line argument parsing
  - Exit codes

### test_mC_helpers.py

Tests are organized by function:

- **TestParseSampleName** - Tests for sample name parsing
- **TestIsDmcSample** - Tests for dmC sample detection
- **TestGetDmcInputType** - Tests for input type determination
- **TestParametersForMc** - Tests for parameter set selection
- **TestEdgeCases** - Edge cases and error conditions
- **TestIntegrationScenarios** - Realistic workflow scenarios

## Mock Data

Test fixtures are defined in:
- `tests/conftest.py` - Shared fixtures for all tests
- `tests/data/` - Mock data files for integration testing

### Available Fixtures

- `temp_dir` - Temporary directory for test files
- `sample_chrom_sizes` - Mock chromosome sizes file
- `valid_bedmethyl_content` - Valid 11-column bedMethyl data
- `valid_bedmethyl_10col_content` - Valid 10-column bedMethyl data
- `bedmethyl_with_header` - bedMethyl with track header
- `mock_bam_header` - Mock BAM file header
- `mock_bam_read_with_mm_ml` - Mock BAM read with methylation tags
- `mock_bam_read_without_mm` - Mock BAM read without methylation tags

## Writing New Tests

### Test Guidelines

1. **Arrange-Act-Assert Pattern**
   ```python
   def test_example(self, fixture):
       # Arrange - set up test data
       input_data = "test_value"

       # Act - call the function
       result = function_under_test(input_data)

       # Assert - verify the outcome
       assert result == expected_value
   ```

2. **Use Descriptive Names**
   - Test names should describe the scenario and expected result
   - Format: `test_<function>_<scenario>_<expected>`
   - Example: `test_validate_bedmethyl_negative_coverage`

3. **Test Both Success and Failure**
   - Always test both valid and invalid inputs
   - Test edge cases (empty, null, boundary values)
   - Test error handling

4. **Use Fixtures for Reusable Data**
   - Add new fixtures to `conftest.py`
   - Keep fixtures focused and single-purpose
   - Document fixture purpose in docstring

5. **Mock External Dependencies**
   - Use `@patch` for subprocess calls, file I/O
   - Mock only what you need to control
   - Verify mocks are called correctly

### Example Test

```python
def test_validate_bedmethyl_valid_file(self, temp_dir, valid_bedmethyl_content):
    """Test validation passes for valid bedMethyl file."""
    # Arrange
    bedmethyl_path = temp_dir / "test.bedmethyl"
    bedmethyl_path.write_text(valid_bedmethyl_content)

    # Act
    is_valid, message = validate_bedmethyl(str(bedmethyl_path))

    # Assert
    assert is_valid is True
    assert "Valid bedMethyl format" in message
```

## Continuous Integration

These tests can be integrated into CI/CD pipelines:

```yaml
# Example GitHub Actions workflow
name: Unit Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      - name: Install dependencies
        run: |
          pip install pytest pytest-cov
      - name: Run tests
        run: |
          pytest tests/unit/ --cov=workflow/scripts --cov-report=xml
      - name: Upload coverage
        uses: codecov/codecov-action@v2
```

## Troubleshooting

### Import Errors

If you get import errors when running tests:

```bash
# Make sure you're in the repository root
cd /path/to/EPICC

# Run tests with PYTHONPATH set
PYTHONPATH=. pytest tests/unit/
```

### Samtools Not Found (for modBAM tests)

The modBAM validation tests use mocked subprocess calls, so samtools doesn't need to be installed for unit tests. However, if you want to run integration tests:

```bash
# Make sure samtools is in your PATH
conda activate epicc
which samtools

# Or install it
conda install -c bioconda samtools
```

### Test Failures

If tests fail:

1. Check the test output for specific error messages
2. Run with `-v` flag for verbose output
3. Run with `--pdb` flag to drop into debugger on failure
4. Check that fixtures are being used correctly

```bash
# Verbose output
pytest tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl -v

# Drop into debugger on failure
pytest tests/unit/test_validate_dmc_input.py --pdb
```

## Test Coverage Goals

Target coverage metrics:
- **validate_dmc_input.py**: >90% line coverage
- **mC.smk helper functions**: >95% line coverage

Current coverage can be checked with:
```bash
pytest tests/unit/ --cov=workflow/scripts/validate_dmc_input --cov-report=term-missing
```

## Contributing

When adding new features to the dmC (direct methylation) workflow:

1. Write tests FIRST (TDD approach)
2. Ensure all existing tests pass
3. Add tests for new functionality
4. Maintain or improve coverage percentage
5. Update this README if adding new test categories

## References

- [pytest documentation](https://docs.pytest.org/)
- [pytest fixtures](https://docs.pytest.org/en/stable/fixture.html)
- [unittest.mock documentation](https://docs.python.org/3/library/unittest.mock.html)
