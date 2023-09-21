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


class MolecularDynamicsSimulation:
    def __init__(self, args):
        self.inputfile = args.input
        self.output_path = args.output_path
        self.gpu = args.gpu
        self.solvent = args.solvent
        self.minimize = args.minimize
        self.time_step = args.time_step
        self.report_interval = args.report_interval
        self.sim_name = basename(self.inputfile).rsplit(".", 1)[0]
        self.output_dir = os.path.join(self.output_path, self.sim_name)
        os.makedirs(self.output_dir, exist_ok=True)
        self.platform = (
            Platform.getPlatformByName("CUDA")
            if self.gpu
            else Platform.getPlatformByName("CPU")
        )

    def load_input(self):
        self.pdb = PDBFile(self.inputfile)
        self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
        self.forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        self.modeller.addHydrogens(self.forcefield)

    def add_solvent(self):
        if self.solvent:
            self.modeller.addSolvent(self.forcefield, padding=1.0 * nanometers)
            with open(
                os.path.join(self.output_dir, self.sim_name) + "_post-solvent.pdb", "w"
            ) as postsolvent:
                PDBFile.writeFile(
                    self.modeller.topology, self.modeller.positions, postsolvent
                )

    def create_system(self):
        self.system = self.forcefield.createSystem(
            self.modeller.topology,
            nonbondedMethod=NoCutoff,
            nonbondedCutoff=1 * nanometer,
            constraints=HBonds,
        )
        self.integrator = LangevinMiddleIntegrator(
            293.15 * kelvin, 1 / picosecond, 0.004 * picoseconds
        )
        self.simulation = Simulation(
            self.modeller.topology, self.system, self.integrator, self.platform
        )
        self.simulation.context.setPositions(self.modeller.positions)

    def energy_minimization(self):
        if self.minimize:
            self.simulation.minimizeEnergy()
            state = self.simulation.context.getState(getEnergy=True)
            with open(
                os.path.join(self.output_dir, self.sim_name) + "_minimized.pdb", "w"
            ) as minimized:
                positions = self.simulation.context.getState(
                    getPositions=True
                ).getPositions()
                PDBFile.writeFile(self.simulation.topology, positions, minimized)

    def add_reporters(self):
        self.simulation.reporters.append(
            PDBReporter(
                os.path.join(self.output_dir, self.sim_name) + "_out.pdb",
                self.report_interval,
            )
        )
        self.simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                self.report_interval,
                step=True,
                potentialEnergy=True,
                temperature=True,
                speed=True,
                elapsedTime=True,
                progress=True,
                remainingTime=True,
                totalSteps=self.time_step,
            )
        )
        self.simulation.reporters.append(
            StateDataReporter(
                os.path.join(self.output_dir, self.sim_name) + "_log.csv",
                self.report_interval,
                step=True,
                potentialEnergy=True,
                temperature=True,
                speed=True,
                elapsedTime=True,
                progress=True,
                remainingTime=True,
                totalSteps=self.time_step,
            )
        )

    def run_simulation(self):
        self.simulation.step(self.time_step)


def main():
    parser = argparse.ArgumentParser(description="Runs a Molecular Dynamics Simulation")
    parser.add_argument(
        "-i", "--input", dest="input", required=True, help="Input pdb file"
    )
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

    simulation = MolecularDynamicsSimulation(args)
    simulation.load_input()
    simulation.add_solvent()
    simulation.create_system()
    simulation.energy_minimization()
    simulation.add_reporters()
    simulation.run_simulation()


if __name__ == "__main__":
    main()
