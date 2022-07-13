"""
A small and simple Molecular Dynamics tool using OpenMM
"""
import os
import sys
import argparse
from posixpath import basename
from openmm.app import *
from openmm import *
from openmm.unit import *

# I/O
parser = argparse.ArgumentParser(description="Runs a Molecular Dynamics Simulation")
parser.add_argument("-i", "--input", dest="input", required=True, help="Input pdb file")
parser.add_argument(
    "-o",
    "--output_path",
    dest="output_path",
    type=str,
    default="./output",
    help="Path to save output files",
)
parser.add_argument(
    "-gpu",
    "--use_gpu",
    dest="gpu",
    type=bool,
    default=True,
    help="Enables running simulation on GPU",
)
parser.add_argument(
    "-s",
    "--solvent",
    dest="solvent",
    type=bool,
    default=True,
    help="Adds solvent to the simulation",
)
parser.add_argument(
    "-m",
    "--minimize",
    dest="minimize",
    type=bool,
    default=True,
    help="Performs Energy Minimization of the structure before simulation",
)
parser.add_argument(
    "-t",
    "--time_step",
    dest="time_step",
    type=int,
    default=5000,
    help="Number of time-steps to perform",
)
parser.add_argument(
    "-r",
    "--report_interval",
    dest="report_interval",
    type=int,
    default=1000,
    help="Frequency of simulation reports in time-steps",
)
args = parser.parse_args()

# Arguments
inputfile = args.input
output_path = args.output_path
gpu = args.gpu
solvent = args.solvent
minimize = args.minimize
time_step = args.time_step
report_interval = args.report_interval

# Setting up ID name
sim_name = basename(inputfile).rsplit(".", 1)[0]
print("Preparing simulation of: " + sim_name)

# Setting up output directory
output_dir = os.path.join(output_path, sim_name)
os.makedirs(output_dir)

# Setting up platform
if gpu is True:
    platform = Platform.getPlatformByName("CUDA")
else:
    platform = Platform.getPlatformByName("CPU")
print("Simulation platform: " + platform.getName())

# Loading input
pdb = PDBFile(inputfile)
print("Input file loaded")

modeller = Modeller(pdb.topology, pdb.positions)

forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

if solvent is True:
    print("Filling the box with water")
    modeller.addSolvent(forcefield, padding=1.0 * nanometers)

    print("Saving solvated structure")
    with open(
        os.path.join(output_dir, sim_name) + "_post-solvent.pdb", "w"
    ) as postsolvent:
        PDBFile.writeFile(modeller.topology, modeller.positions, postsolvent)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=NoCutoff,
    nonbondedCutoff=1 * nanometer,
    constraints=HBonds,
)

integrator = LangevinMiddleIntegrator(
    293.15 * kelvin, 1 / picosecond, 0.004 * picoseconds
)

simulation = Simulation(modeller.topology, system, integrator, platform)

simulation.context.setPositions(modeller.positions)

if minimize is True:
    print("Running Energy Minimization")
    simulation.minimizeEnergy()

    state = simulation.context.getState(getEnergy=True)

    print("Minimization complete E: " + str(state.getPotentialEnergy()))

    print("Saving minimized structure")
    with open(os.path.join(output_dir, sim_name) + "_minimized.pdb", "w") as minimized:
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, minimized)


simulation.reporters.append(
    PDBReporter(os.path.join(output_dir, sim_name) + "_out.pdb", report_interval)
)

simulation.reporters.append(
    StateDataReporter(
        sys.stdout, report_interval, step=True, potentialEnergy=True, temperature=True
    )
)

simulation.reporters.append(
    StateDataReporter(
        os.path.join(output_dir, sim_name) + "_log.csv",
        report_interval,
        step=True,
        potentialEnergy=True,
        temperature=True,
    )
)

print("Beginning Time Steps")
simulation.step(time_step)
