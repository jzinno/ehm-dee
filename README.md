# ehmdee: Molecular Dynamics Simulation with OpenMM

This program is a Python script that runs a Molecular Dynamics (MD) simulation using the OpenMM library. It is designed to be run from the command line and accepts several arguments that allow you to customize the simulation.

## Getting Started

To run the simulation, navigate to the directory containing the script and run:

python ehmdee.py --input input.pdb --output output.pdb --use_gpu --solvent --minimize --timesteps 10000 --report_interval 100

## Command-line Arguments

- `--input`: Path to the input file in PDB format.
- `--output_path`: Path to save the output file.
- `--use_gpu`: Use GPU for the simulation.
- `--solvent`: Add solvent to the simulation.
- `--minimize`: Minimize the energy of the system before running the simulation.
- `--timesteps`: Number of time steps to run in the simulation.
- `--report_interval`: Interval at which to report simulation data.

## Simulation Demonstration

![Simulation GIF](assets/GLH-1_MIP-1_seperated-MD_bounce.gif)
