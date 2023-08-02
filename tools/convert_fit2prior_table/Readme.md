# Fit2Prior

This provides a tool to convert posterior from an analysis to a 1D or a 2D tabulated prior, that can be used with TAMCMC v1.85.0 and higher.

The main function is `make_prior_table()` within `make_priors_tables.py`. But you can also use the more convenient function dopriors() within do_many_processes.py.
This function allows you to scan a directory for processes and automatically identifies the parameters for which a tabulation is requested by their names. It has also
some (limited) conversion capabilities between some parameters that are equivalent in various version of TAMCMC, but that could be named differently or parametrised differently.

You must have compiled the main TAMCMC code in order to get working `bin2txt` binary. And you need to specify the path of that binary (default is cpp_prg="../bin/", assumes that you run the program from this tools subdirectory).

Unit tests are provided in `unittest_priors_tables.py`.

## Usage

```python
binning = 30
work_dir = os.getcwd()
process_name = "10280410_Gaussfit"
dir_tamcmc_outputs = work_dir + "/test_data/" #+ process_name + "/outputs/"
output_dir = work_dir + "/test_data/Outputs/"

print("Generate a 1D tabulated prior")
index1 = 7  # numax parameter of the test file
make_prior_table(
    index1,
    binning,
    output_dir,
    dir_tamcmc_outputs,
    process_name,
    phase='A',
    chain=0,
    first_index=0,
    last_index=-1,
    period=5,
    index2=None,
    cpp_prg="../../bin/"
)

print("Generate a 2D tabulated prior")
index1 = 8  # numax parameter of the test file in the x-axis
index2 = 7  # Amax parameter of the test file in the y-axis
make_prior_table(
    index1,
    binning,
    output_dir,
    dir_tamcmc_outputs,
    process_name,
    phase='A',
    chain=0,
    first_index=0,
    last_index=-1,
    period=5,
    index2=index2,
    cpp_prg="../../bin/"
)
