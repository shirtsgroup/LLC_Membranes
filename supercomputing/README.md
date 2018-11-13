# Setting up super computing resources

Reference this document if you are setting up jobs on a super computing resource for the first time

## Requirements for Bridges
### Python

Load the appropriate python module: `module load python/3.6.4_gcc5_np1.14.5`

pip install the following modules with the command:  `pip install --user module_name`: \
`tqdm` \
`mdtraj`

Use bridges_slurm_gpu.sh as a guide for creating job submissions

## Gromacs
Bridges has pre-compiled versions of gromacs: \
`gromacs/2018_cpu` -- compiled for use with MPI \
`gromasc/2018_gpu` -- compiled for use with GPU and MPI

You can compile your own version of gromacs using these [instructions](https://github.com/shirtsgroup/shirts-group-wiki/wiki/Installing-MPI-and-GPU-enabled-GROMACS#on-bridges)

*IMPORTANT*: GROMACS will not parallelize properly unless you load the correct MPI module. Make sure to include the following line in your batch submission scripts:\
`module load mpi/intel_mpi` \
Note that this may cause issues with the some python modules so load the module _after_ any python scripts.
