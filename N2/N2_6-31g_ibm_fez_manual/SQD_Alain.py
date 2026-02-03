# Sample-based Quantum Diagonalization (SQD) by Alain Chancé

## MIT License

# MIT_License Copyright (c) 2025 Alain Chancé
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the \"Software\"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.'
# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.'

#------------------------------------------------------------------------------------------------------------------------
## Credit
# This Python code is derived from the tutorial "Improving energy estimation of a chemistry Hamiltonian with SQD"
# https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb
# which is distributed under the "Apache License 2.0", https://github.com/Qiskit/qiskit-addon-sqd/blob/main/LICENSE.txt.
#
## Additions by Alain Chancé
# The SQD class integrates the eco2AI tracking feature, a python library which accumulates statistics about power
# consumption and CO2 emission during running code. The Eco2AI is licensed under a Apache licence 2.0.
# https://github.com/sb-ai-lab/Eco2AI
#------------------------------------------------------------------------------------------------------------------------

# Import common packages first
import psutil
import sys
import os
import time
import numpy as np
import pandas as pd
from math import comb
import warnings
import pyscf
import matplotlib.pyplot as plt
import pickle
from functools import partial

#------------------------------------------------------------------------------------------------------------------------
# Import qiskit.aer
# Additional circuit methods. On import, Aer adds several simulation-specific methods to QuantumCircuit for convenience. 
# These methods are not available until Aer is imported (import qiskit_aer). 
# https://qiskit.github.io/qiskit-aer/apidocs/circuit.html
#------------------------------------------------------------------------------------------------------------------------
import qiskit_aer

#---------------------
# Import AerSimulator
#---------------------
from qiskit_aer import AerSimulator

#---------------------------
# Import StatevectorSampler
#---------------------------
from qiskit.primitives import StatevectorSampler

#-----------------
# Import BitArray
#-----------------
from qiskit.primitives import BitArray

#-------------------------------------
# Import Fake provider for backend V2
#-------------------------------------
# from qiskit_ibm_runtime import fake_provider
# https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/fake-provider-fake-provider-for-backend-v2
# https://github.com/Qiskit/qiskit-ibm-runtime/blob/stable/0.40/qiskit_ibm_runtime/fake_provider/fake_provider.py
from qiskit_ibm_runtime.fake_provider import FakeProviderForBackendV2
from qiskit_ibm_runtime.fake_provider import FakeTorino

#-----------------------------------------------------------------------
# Import from Qiskit Aer noise module (kept for optional fake backends)
#-----------------------------------------------------------------------
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    ReadoutError,
    depolarizing_error,
    pauli_error,
    thermal_relaxation_error,
)

#-----------------------
# Import qiskit classes
#-----------------------
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit.visualization import plot_gate_map

# https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
from qiskit_addon_sqd.fermion import SCIResult, diagonalize_fermionic_hamiltonian, solve_sci_batch
from qiskit_addon_sqd.counts import generate_bit_array_uniform, counts_to_arrays

# Import qiskit ecosystems
import ffsim
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
from qiskit_ibm_runtime import SamplerOptions

# Import the function get_zigzag_physical_layout
from zigzag import get_zigzag_physical_layout

warnings.filterwarnings("ignore")

#-----------------------------------------------------------------------------------------
# If the code is running in an IPython terminal, then from IPython.display import display
# else if it is running in a plain Python shell or script: 
# we assign display = print and array_to_latex = identity
#-----------------------------------------------------------------------------------------
try:
    shell = get_ipython().__class__.__name__
except NameError:
    shell = None

if shell == 'TerminalInteractiveShell':
    # The code is running in an IPython terminal
    from IPython.display import display
    
elif shell == None:
    # The code is running in a plain Python shell or script: 
    # we assign display = print and array_to_latex = identity
    display = print
    array_to_latex = lambda x: x

#--------------------------------------------------------
# Print version of Qiskit, Qiskit Aer and Qiskit runtime
#--------------------------------------------------------
import qiskit
print(f"Qiskit version: {qiskit.__version__}")

import qiskit_aer
print(f"Qiskit Aer version: {qiskit_aer.__version__}")

import qiskit_ibm_runtime
print(f"Qiskit runtime version: {qiskit_ibm_runtime.__version__}")

#------------------------------------------------------------------------------------------------
# If the import "from SQD_Alain import SQD" fails and file SQD_Alain.py is in the same directory 
# as your python or Jupyter notebook, try adding the following lines:
# import sys
# import os
# cwd = os.getcwd()
# _= (sys.path.append(cwd))
#------------------------------------------------------------------------------------------------
class SQD:
    def __init__(self,
                 #-------------
                 # Run options
                 #-------------
                 config_filename = None,                           # Configuration file name (.json)
                 backend_name = None,                              # IBM cloud backend name
                 job_id = None,                                    # job_id of a previously run job
                 do_plot_gate_map = True,                          # Whether to plot the gate map 
                 load_bit_array_file = None,                       # If provided, function step_3 will load samples from this file
                 save_bit_array_file = None,                       # If provided, function step_3 will save samples into this file
                 n_ancillary_qubits = 0,                           # Number of ancillary qubits
                 run_on_QPU = False,                               # Whether to run the quantum circuit on the target hardware
                 nshots = 1000,                                    # Number of shots
                 opt_level = 1,                                    # Optimization level
                 poll_interval = 5,                                # Poll interval in seconds for job monitor
                 timeout = 600,                                    # Time out for job monitor
                 nsamples = 10000,                                 # Number of samples to be drawn from the uniform distribution
                 #-------------------------------------
                 # eco2AI Tracker options
                 # https://github.com/sb-ai-lab/Eco2AI
                 #-------------------------------------
                 do_eco2ai = False,                                # Whether to use the eco2AI Tracker
                 project_name = "SQD_Alain",                       # Project name
                 experiment_description = "SQD",                   # Experiment description
                 eco2ai_file_name = "SQD.csv",                     # eco2AI file name
                 #---------------------------------------------------------------------------------
                 # Ballpark figure (kW) for the power consumption of the IBM cloud backend
                 # "The power consumption of a quantum computer is about 15-25kW"
                 # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
                 #---------------------------------------------------------------------------------
                 power_QPU = 25,                                   # Ballpark figure (kW) for the power consumption of the IBM cloud backend
                 #-------------------------------------------------------------------------------
                 # PySCF options
                 # https://pyscf.org/user/gto.html#initializing-a-molecule
                 # https://pyscf.org/pyscf_api_docs/pyscf.gto.html#pyscf.gto.mole.MoleBase.build
                 #-------------------------------------------------------------------------------
                 basis = "6-31g",                                  # Basis set
                 atom = [["N", (0, 0, 0)], ["N", (1.0, 0, 0)]],    # Atom configuration, PySCF Initializing a molecule
                                                                   # https://pyscf.org/user/gto.html#initializing-a-molecule
                 spin = 0,                                         # The number of unpaired electrons of a molecule
                 charge = None,                                    # Charge of a molecule
                 symmetry = "Dooh",                                # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry
                 n_frozen = 2,                                     # Number of frozen orbitals
                 #---------------------
                 # Molecular integrals
                 #---------------------
                 compute_exact_energy = True,                      # Whether to compute the exact energy with PySCF, cas.run().e_tot
                 exact_energy = None,
                 #----------------------------------------------------------------------------------
                 # SQD options
                 # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
                 #----------------------------------------------------------------------------------
                 chem_accuracy = 1e-3,                             # Chemical accuracy (+/- 1 milli-Hartree)
                 energy_tol = 1e-4,                                # Numerical tolerance for convergence of the energy
                 occupancies_tol = 1e-3,                           # Numerical tolerance for convergence of the average orbital occupancies
                 max_iterations = 10,                              # Limit on the number of configuration recovery iterations
                 #-----------------------------------------------------------------------------------
                 # Eigenstate solver options
                 # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
                 #-----------------------------------------------------------------------------------
                 num_batches = 5,                                  # The number of batches to subsample in each configuration recovery iteration
                 samples_per_batch = 300,                          # The number of bitstrings to include in each subsampled batch of bitstrings
                 symmetrize_spin = True,                           # Whether to always merge spin-alpha and spin-beta CI strings into a single list
                                                                   # Spin symmetrization is only possible if the numbers of alpha and beta electrons are equal
                 carryover_threshold = 1e-4,                       # Threshold for carrying over bitstrings with large CI weight from one iteration to the next
                 max_cycle = 200,                                  # The maximum number of Davidson cycles run by the eigenstate solver
                 seed = 24,                                        # A seed for the pseudorandom number generator
                 spin_sq = 0.0                                     # Spin square
              ):

        # Initialize self.backend to None
        self.backend = None
        
        #-------------------
        # Print run options
        #-------------------
        print("\nRun options")
        print("Configuration file name (.json):", config_filename)
        print("Backend name:", backend_name)
        
        if job_id != None:
            print("job_id:", job_id)
        
        print("do plot gate map: ", do_plot_gate_map)
        
        if load_bit_array_file != None:
            print("load_bit_array_file:", load_bit_array_file)
            
        if save_bit_array_file != None:
            print("save_bit_array_file:", save_bit_array_file)
            
        print("Run on QPU:", run_on_QPU)

        if run_on_QPU:
            nshots = min(max(1000, nshots), 200000)
            print("Number of shots:", nshots)
        else:
            nsamples = min(max(10000, nsamples), 200000)
            print("Number of samples drawn from the uniform distribution:", nsamples)

        #-------------------------------------
        # Print eco2AI Tracker options
        # https://github.com/sb-ai-lab/Eco2AI
        #-------------------------------------
        if do_eco2ai:
            print("\neco2AI Tracker options")
            print("project_name:", project_name)
            print("experiment_description:", experiment_description)
            print("eco2ai_file_name:", eco2ai_file_name)
            
        #---------------------------------------------------------------------------------
        # Print Ballpark figure (kW) for the power consumption of the IBM cloud backend
        # "The power consumption of a quantum computer is about 15-25kW"
        # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
        #---------------------------------------------------------------------------------
        power_QPU = min(max(15.0, power_QPU), 25.0)
        print("\nBallpark figure (kW) for QPU power consumption:", power_QPU)
        
        #---------------------
        # Print PySCF options
        #---------------------
        print("\nPySCF options")
        print("Basis:", basis)
        print("Atom configuration: ", atom)
        print("spin: ", spin)
        if charge != None:
            print("charge: ", charge)
        print("Symmetry: ", symmetry)
        print("Number of frozen orbitals: ", n_frozen)

        #-----------------------------------
        # Print molecular integrals options
        #-----------------------------------
        print("Compute exact energy: ", compute_exact_energy)

        if not compute_exact_energy and exact_energy is not None:
            print("exact energy: ", exact_energy)

        #-------------------
        # Print SQD options
        #-------------------
        print("\nSQD options")
        print("Chemical accuracy: ", chem_accuracy)
        print("Numerical tolerance for convergence of the energy: ", energy_tol)
        print("Numerical tolerance for convergence of the average orbital occupancies: ", occupancies_tol)
        print("Limit on the number of configuration recovery iterations: ", max_iterations)

        #---------------------------
        # Eigenstate solver options
        #---------------------------
        print("\nEigenstate solver options")
        print("Number of batches to subsample in each configuration recovery iteration: ", num_batches)
        print("Number of bitstrings to include in each subsampled batch of bitstrings: ", samples_per_batch)
        print("Symmetrize spin: ", symmetrize_spin)
        print("Threshold for carrying over bitstrings with large CI weight from one iteration to the next: ", carryover_threshold)
        print("Maximum number of Davidson cycles run by the eigenstate solver: ", max_cycle)
        print("Seed for the pseudorandom number generator: ", seed)
        print("Spin square: ", spin_sq) 
        
        #-------------------------
        # Set up param dictionary
        #-------------------------
        self.param = {
            "job": None,                                     # job = service.job(job_id) if job_id is provided
            "circuit": None,                                 # Quantum circuit returned by function step_1
            "circuit_depth": 0,                              # Depth of the quantum circuit returned by function step_1
            "n_qubits": 2,                                   # Number of qubits set by function step_1
            "isa_circuit": None,                             # ISA circuit returned by function step_2
            "isa_circuit_depth": 0,                          # Depth of the ISA quantum circuit returned by function step_2
            "bit_array": None,                               # Bit array returned by function step_3
            "result": None,                                  # result returned by post_process function
            "result_history": None,                          # result_history returned by post_process function
            "SQD_energy": None,                              # SQD energy returned by post_process function
            "Absolute_error": None,                          # Absolute error returned by post_process function
            "QPU_usage": None,                               # Qiskit Runtime usage in seconds returned by function get_QPU_usage
            "QPU_power_consumption": None,                   # QPU power consumption returned by function get_QPU_usage
            "eco2ai_tracker": None,                          # eco2AI tracker
            "duration": None,                                # duration of classical processing returned by get_classical_power_usage
            "classical_power_usage": None,                   # classical power usage returned by get_classical_power_usage
            "plt_energy": None,                              # Plot energy and occupancy returned by function plot_energy_and_occupancy
            #-------------------------------
            # Layout set by setup_zig_zag()
            #-------------------------------
            "initial_layout": [],                            # Initial_layout
            #------------------------------------------
            # Files containing token (API key) and CRN
            #------------------------------------------
            "token_file": "Token.txt",                       # Token file
            "CRN_file": "CRN.txt",                           # CRN file
            #---------------------------
            # Bit array file (optional)
            #---------------------------
            "load_bit_array_file": load_bit_array_file,      # If provided, function step_3 will load samples from this file
            "save_bit_array_file": save_bit_array_file,      # If provided, function step_3 will save samples into this file
            #-------------
            # Run options
            #-------------
            "config_filename": config_filename,              # Configuration file name (.json)
            "backend_name": backend_name,                    # IBM cloud backend name
            "job_id": job_id,                                # job_id
            "do_plot_gate_map": do_plot_gate_map,            # Whether to plot the gate map 
            "n_ancillary_qubits": n_ancillary_qubits,        # Number of ancillary qubits
            "run_on_QPU": run_on_QPU,                        # Whether to run the quantum circuit on the target hardware
            "nshots": nshots,                                # Number of shots
            'opt_level':opt_level,                           # Optimization level
            "poll_interval": poll_interval,                  # Poll interval in seconds for job monitor
            "timeout": timeout,                              # Time out in seconds for gob monitor
            "nsamples": nsamples,                            # Number of samples to be drawn from the uniform distribution
            #-------------------------------------
            # eco2AI Tracker options
            # https://github.com/sb-ai-lab/Eco2AI
            #-------------------------------------
            "do_eco2ai": do_eco2ai,                           # Whether to use the eco2AI Tracker
            "project_name": project_name,                     # Project name
            "experiment_description": experiment_description, # Experiment description
            "eco2ai_file_name": eco2ai_file_name,             # File name
            #---------------------------------------------------------------------------------
            # Ballpark figure (kW) for the power consumption of the IBM cloud backend
            # "The power consumption of a quantum computer is about 15-25kW"
            # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
            #---------------------------------------------------------------------------------
            "power_QPU": power_QPU,                           # Ballpark figure (kW) for the power consumption of the IBM cloud backend
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            "mol": pyscf.gto.Mole(),                          # Molecule
            "basis": basis,                                   # Basis set
            "spin": spin,                                     # Spin
            "charge": charge,                                 # Charge
            "atom": atom,                                     # Atom configuration, PySCF Initializing a molecule
            "symmetry": symmetry,                             # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry 
            "n_frozen": n_frozen,                             # Number of frozen orbitals
            #---------------------
            # Molecular integrals
            #---------------------
            "num_orbitals": 2,                                # The number of spatial orbitals
            "hcore": 0.0,                                     # The one-body tensor of the Hamiltonian
            "nuclear_repulsion_energy": 0.0,                  # Nuclear repulsion energy
            "eri": 0.0,                                       # The two-body tensor of the Hamiltonian
            "nelec": 1,                                       # The numbers of alpha and beta electrons
            "compute_exact_energy": compute_exact_energy,     # Whether to compute the exact energy with PySCF, cas.run().e_tot
            "exact_energy": exact_energy,                     # Exact energy
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            "chem_accuracy": chem_accuracy,                   # Chemical accuracy (+/- 1 milli-Hartree)
            "energy_tol": energy_tol,                         # Numerical tolerance for convergence of the energy
            "occupancies_tol": occupancies_tol,               # Numerical tolerance for convergence of the average orbital occupancies
            "max_iterations": max_iterations,                 # Limit on the number of configuration recovery iterations
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            "num_batches": num_batches,                    # The number of batches to subsample in each configuration recovery iteration
            "samples_per_batch": samples_per_batch,        # The number of bitstrings to include in each subsampled batch of bitstrings
            "symmetrize_spin": symmetrize_spin,            # Whether to always merge spin-alpha and spin-beta CI strings into a single list
                                                           # Spin symmetrization is only possible if the numbers of alpha and beta electrons are equal
            "carryover_threshold": carryover_threshold,    # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            "max_cycle": max_cycle,                        # The maximum number of Davidson cycles run by the eigenstate solver
            "seed": seed,                                  # A seed for the pseudorandom number generator
            "spin_sq": spin_sq                             # spin square
            }

        #--------------------------------------------------------
        # If do_eco2ai is True, try importing the eco2AI library
        #--------------------------------------------------------
        if do_eco2ai:
            try:
                import eco2ai
            except Exception as e:
                print(f"Error importing the eco2AI library: {e}")
                print(f"Install eco2AI with the command 'pip install eco2ai'")
                do_eco2ai = False
                self.param['do_eco2ai'] = False

        #-----------------------------------------
        # Start an instance of the eco2AI Tracker
        # https://github.com/sb-ai-lab/Eco2AI
        #-----------------------------------------
        if do_eco2ai:

            # Delete eco2ai_file_name to reset duration counter
            if os.path.isfile(eco2ai_file_name): 
                try: 
                    os.remove(eco2ai_file_name) 
                    print(f"\nDeleted: {eco2ai_file_name}") 
                except Exception as e: 
                    print(f"\nError deleting file: {e}")
            
            tracker = eco2ai.Tracker(
                project_name = project_name, 
                experiment_description = experiment_description,
                file_name = eco2ai_file_name
            )
            self.param['eco2ai_tracker'] = tracker

            tracker.start()
            print(f"eco2AI tracker started")

        #-------------------------------------------------------------------------------------------
        # Spin symmetrization is only possible if the numbers of alpha and beta electrons are equal
        #-------------------------------------------------------------------------------------------
        if spin != 0:
            self.param['symmetrize_spin'] = False
            print("\nSpin symmetrization will not be performed because the number of alpha and beta electrons is unequal") 

        #---------------------------------------------------------
        # Define and initialize a molecule with PySCF
        # PySCF Initializing a molecule
        # https://pyscf.org/user/gto.html#initializing-a-molecule
        #---------------------------------------------------------
        print("\nDefining and initializing a molecule with PySCF")
        self.build_molecule()
        
        #---------------------------------------------------
        # Get Qiskit Runtime account credentials from files
        #---------------------------------------------------
        token_file = self.param['token_file']
        CRN_file = self.param['CRN_file']

        if os.path.isfile(token_file):
            f = open(token_file, "r") 
            token = f.read() # Read token from token_file
            print("\nToken read from file: ", token_file)
            f.close()

            if os.path.isfile(CRN_file):
                f = open(CRN_file, "r") 
                crn = f.read() # Read CRN code from token_file
                print("CRN code read from file: ", CRN_file, "\n")
                f.close()

                #---------------------------------------------------------------------------------------------------------------------
                # Save the Qiskit Runtime account credentials
                # Save your account on disk
                # The credentials are saved in the $HOME/.qiskit/qiskit-ibm.json file, where $HOME is your home directory.
                # https://github.com/Qiskit/qiskit-ibm-runtime?tab=readme-ov-file#save-your-account-on-disk
                #---------------------------------------------------------------------------------------------------------------------
                QiskitRuntimeService.save_account(channel="ibm_cloud", token=token, instance=crn, set_as_default=True, overwrite=True)

        return

    #--------------------------------------------------------------------------------------
    # Define function build_molecule() which defines and initializes a molecule with pySCF
    # PySCF Initializing a molecule
    # https://pyscf.org/user/gto.html#initializing-a-molecule
    # https://pyscf.org/pyscf_api_docs/pyscf.gto.html#pyscf.gto.mole.MoleBase.build
    #--------------------------------------------------------------------------------------
    def build_molecule(self):
        
        if self.param is None:
            print("build_molecule: missing parameter param")
            return None 
        param = self.param
        
        mol = pyscf.gto.Mole()
        mol.build(
            atom = param['atom'],
            basis = param['basis'],
            spin = param['spin'],
            charge = param['charge'],
            symmetry = param['symmetry']
        )
        param['mol'] = mol
        return

    #-------------------------------------
    # Define the function setup_backend() 
    #-------------------------------------
    def setup_backend(self):

        if self.param is None:
            print("setup_backend: missing parameter param")
            return None 
        param = self.param
        
        #-------------------------------------------------------------------------------------------
        # Instantiate the service
        # Once the account is saved on disk, you can instantiate the service without any arguments:
        # https://docs.quantum.ibm.com/api/migration-guides/qiskit-runtime
        #-------------------------------------------------------------------------------------------
        try:
            self.service = QiskitRuntimeService()
        except:
            print(f"Error creating an instance of QiskitRuntimeService(): {e}")
            self.service = None

        service = self.service

        #---------------------------------------------------------------
        # If a valid job_id is provided, then get the corresponding job
        #---------------------------------------------------------------
        job_id = param['job_id']
        if job_id is not None:
            try:
                param['job'] = service.job(job_id)
                print(f"Successfully retrieved job with job_id: {job_id}")
            except Exception as e:
                print(f"Error retrieving job with job_id {job_id}: {e}")
        
        backend_name = param['backend_name']
        n_qubits = param['n_qubits']
        opt_level = param['opt_level']
        
        print("backend_name:", backend_name)
        
        if service is None or backend_name == "AerSimulator noiseless":
            # Use AerSimulator(method='statevector')
            # https://docs.quantum.ibm.com/migration-guides/local-simulators#aersimulator
            self.backend = AerSimulator(method='statevector')
            print("\nUsing AerSimulator with method statevector and noiseless")

            self.sampler = StatevectorSampler()

        else:
            self.backend = None
            if backend_name is None or backend_name == "None":
                # Assign least busy device to backend
                # https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/qiskit-runtime-service#least_busy
                try:
                    self.backend = service.least_busy(min_num_qubits=n_qubits, simulator=False, operational=True)
                    backend_name = self.backend.name
                    param['backend_name'] = self.backend.name

                    # Print the least busy device
                    print(f"The least busy device: {self.backend}")
                    
                except Exception as e:
                    print(f"No suitable backend found with minimum: {n_qubits} qubits - Default to 'fake_torino'")
                    backend_name = "fake_torino"

            if not backend_name[:4] in ['None', 'ibm_', 'fake']:
                print(f"Unknown backend name: {backend_name} - Default to 'fake_torino'")
                backend_name = "fake_torino"

            if backend_name[:4] == "ibm_":

                if backend_name == "ibm_torino":
                    self.backend = service.backend(backend_name)
                    opts = SamplerOptions()
                    opts.dynamical_decoupling.enable = True
                    opts.twirling.enable_measure = True
                    self.sampler = Sampler(mode=self.backend, options=opts)

                else:
                    # Get the operational real backends that have a minimum of n_qubits
                    # https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/qiskit-runtime-service
                    try:
                        backends = service.backends(backend_name, min_num_qubits=n_qubits, simulator=False, operational=True)
                        self.backend = backends[0]
                        backend_name = self.backend.name
                        param['backend_name'] = self.backend.name
                        self.sampler = Sampler(mode=self.backend)
                    except Exception as e:
                        print(f"No suitable backend found with name {backend_name} and minimum: {n_qubits} qubits - Default to 'fake_torino'")
                        backend_name = "fake_torino"

            if backend_name[:4] == "fake":
                # https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/fake-provider-fake-provider-for-backend-v2
                # https://github.com/Qiskit/qiskit-ibm-runtime/blob/stable/0.40/qiskit_ibm_runtime/fake_provider/fake_provider.py
                fake_provider = FakeProviderForBackendV2()
                try:
                    backend = fake_provider.backend(backend_name)
                except Exception as e:
                    print(f"Unknown fake backend name: {backend_name} - Default to 'fake_torino'")
                    backend_name = "fake_torino"
                    backend = FakeTorino()
                
                noise_model = NoiseModel.from_backend(backend)
                self.backend = AerSimulator(method='statevector', noise_model=noise_model)
                param['backend_name'] = self.backend.name
                self.sampler = StatevectorSampler()
                print("\nUsing AerSimulator with method statevector and noise model from", backend_name)

        #----------------------------------------------------------------------------------------------
        # In function step2() - Optimize problem for quantum execution generates a staged pass manager
        #----------------------------------------------------------------------------------------------
        if isinstance(self.backend, AerSimulator):
            # Check that there is enough memory to perform a simulation with AerSimulator
            self.check_size()
                    
        else:
            # Setup the zig-zag pattern for qubit interactions
            self.setup_zig_zag()
            
            print(f"Backend name: {self.backend.name}\n"
                  f"Version: {self.backend.version}\n"
                  f"Number of qubits: {self.backend.num_qubits}"
                 )
        return

    #---------------------------------------------------------------------------------------
    # Setup the zig-zag pattern for qubit interactions
    #
    # Sample-based quantum diagonalization of a chemistry Hamiltonia
    # https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization
    #
    # Step 2: Optimize problem for quantum hardware execution
    # def create_lucj_zigzag_layout()
    # def get_zigzag_physical_layout()
    # initial_layout, _ = get_zigzag_physical_layout(num_orbitals, backend=backend)
    #---------------------------------------------------------------------------------------
    def setup_zig_zag(self):
        if isinstance(self.backend, AerSimulator):
            return

        if self.param is None:
            print("setup_zig_zag: missing parameter param")
            return None 
        param = self.param

        backend = self.backend
        num_orbitals = param['num_orbitals']
        
        initial_layout, num_alpha_beta_qubits = get_zigzag_physical_layout(num_orbitals, backend=backend)
        param['initial_layout'] = initial_layout
        print("initial_layout:", initial_layout)

        self.my_plot_gate_map()
        
        return

    #--------------------------------------------------------------------------------------------------
    # Plot gate map for the backend with the spin-a layout shown in red, and the spin-b layout in blue
    # The 'Graphviz' library is required to use 'plot_coupling_map'.  
    # To install, follow the instructions at https://graphviz.org/download/
    # On Debian, Ubuntu, run the following command: 
    # sudo apt install graphviz
    #--------------------------------------------------------------------------------------------------
    def my_plot_gate_map(self):

        if self.param is None:
            print("my_plot_gate_map: missing parameter param")
            return None 
        param = self.param

        if not param['do_plot_gate_map']:
            return

        initial_layout = param['initial_layout']

        if initial_layout != []:
            try:
                qubit_color = []
                for i in range(self.backend.num_qubits):
                    if i in initial_layout:
                        qubit_color.append("#FFA07A") # light red
                    else:
                        qubit_color.append("#888888") # gray
                
                print("\nGate map for :", self.backend.name, ", the initial layout is shown in red\n")

                # Call the plot_gate_map function which is part of Qiskit’s visualization tools
                fig = plot_gate_map(self.backend, qubit_color=qubit_color, plot_directed=False)

                # Save the figure to a file
                fig.savefig("plot_gate_map.png", bbox_inches='tight')

                # Display the gate map
                display(fig)
            
            except:
                print("The 'Graphviz' library is required to use 'plot_coupling_map'")
                print(sys.exc_info())
            
        return

    #---------------------------------------------------------------------------------------------------------------------
    # Define the function check_size() which checks that there is enough memory to perform a simulation with AerSimulator
    #---------------------------------------------------------------------------------------------------------------------
    def check_size(self):

        if self.param is None:
            print("check_size: missing parameter param")
            return 
        param = self.param
        
        #-----------------------------------------------------------------------------
        # Return if run_on_QPU is False or backend is not an instance of AerSimulator
        #-----------------------------------------------------------------------------
        if not param['run_on_QPU'] or not isinstance(self.backend, AerSimulator):
            return

        n_qubits = param['n_qubits']
        
        # Statevector simulator requires 2**instruction.num_qubits data of type complex:
        # Let's compute the amount of memory required to store 2**n_qubits numbers of data type complex128
        # Amount of memory required to store one quantum circuit that simulates a permutation operation
        # https://numpy.org/doc/stable/reference/arrays.dtypes.html
        mem_circuit = np.dtype(np.complex128).itemsize*8*2**n_qubits/10**9
    
        # Get available memory for processes
        # https://www.geeksforgeeks.org/how-to-get-current-cpu-and-ram-usage-in-python/
        mem_avail = psutil.virtual_memory()[1]/10**9

        # Check that there is enough available memory for the statevector simulation
        if mem_circuit > mem_avail:
        
            print(f"Available memory for processes (GB): {mem_avail} < Amount of memory required by AerSimulator: {mem_circuit}")
            print("Simulating the run by reading bit array file or generating 10k samples drawn from the uniform distribution")
        
            # Set run_on_QPU to False
            self.param['run_on_QPU'] = False

        return
    
    #------------------------------------------------------------------------------------------------------------------
    # Define function step_1() - Map classical inputs to a quantum problem
    #
    # We perform a CCSD calculation. The t1 and t2 amplitudes from this calculation will be used to initialize
    # the parameters of the `LUCJ` ansatz circuit.
    # https://en.wikipedia.org/wiki/Coupled_cluster#Cluster_operator
    #
    # We use ffsim to create the ansatz circuit. ffsim is a software library for simulating fermionic quantum circuits
    # that conserve particle number and the Z component of spin. This category includes many quantum circuits used for
    # quantum chemistry simulations. By exploiting the symmetries and using specialized algorithms, ffsim can simulate
    # these circuits much faster than a generic quantum circuit simulator.
    # https://github.com/qiskit-community/ffsim
    #
    # Depending on the number of unpaired electrons, we use either:
    #
    #   - the spin-balanced variant of the unitary cluster Jastrow (UCJ) ansatz, ffsim class UCJOpSpinBalanced
    #     https://qiskit-community.github.io/ffsim/api/ffsim.html#ffsim.UCJOpSpinBalanced
    #
    #   - the spin-unbalanced (local) unitary cluster Jastrow ansatz
    #     https://qiskit-community.github.io/ffsim/api/ffsim.html#ffsim.UCJOpSpinUnbalanced
    #------------------------------------------------------------------------------------------------------------------
    def step_1(self):

        if self.param is None:
            print("step_1: missing parameter param")
            return None 
        param = self.param
    
        #-------------------------------------------
        # Retrieve parameters from param dictionary
        #-------------------------------------------
        n_frozen = param['n_frozen']
        mol = param['mol']
        n_ancillary_qubits = param['n_ancillary_qubits']

        #---------------------
        # Define active space
        #---------------------
        active_space = range(n_frozen, mol.nao_nr())

        #-------------------------
        # Get molecular integrals
        #-------------------------
        scf = pyscf.scf.RHF(mol).run()

        num_orbitals = len(active_space)
        self.param['num_orbitals'] = num_orbitals
        print("Number of spatial orbitals: ", num_orbitals)

        print("Number of spin orbitals: ", num_orbitals * 2)

        n_electrons = int(sum(scf.mo_occ[active_space]))
        print("Number of electrons: ", n_electrons)

        num_elec_a = (n_electrons + mol.spin) // 2
        print("Number of α-spin electrons: ", num_elec_a)

        num_elec_b = (n_electrons - mol.spin) // 2
        print("Number of β-spin electrons: ", num_elec_b)

        cas = pyscf.mcscf.CASCI(scf, num_orbitals, (num_elec_a, num_elec_b))

        mo = cas.sort_mo(active_space, base=0)

        hcore, nuclear_repulsion_energy = cas.get_h1cas(mo)
        param['hcore'] = hcore
        param['nuclear_repulsion_energy'] = nuclear_repulsion_energy

        eri = pyscf.ao2mo.restore(1, cas.get_h2cas(mo), num_orbitals)
        param['eri'] = eri

        #----------------------------------------------------------------
        # Get CCSD calculation t2 amplitudes for initializing the ansatz
        # https://en.wikipedia.org/wiki/Coupled_cluster#Cluster_operator
        #----------------------------------------------------------------
        ccsd = pyscf.cc.CCSD(scf, frozen=[i for i in range(mol.nao_nr()) if i not in active_space]).run()
        t1 = ccsd.t1
        t2 = ccsd.t2
        print("\n")

        print(f"Total CCSD energy: {ccsd.e_tot}")
        print(f"CCSD correlation energy: {ccsd.e_corr}")
        print(f"SCF energy = Total CCSD energy - CCSD correlation energy: {ccsd.e_tot - ccsd.e_corr}")

        #-----------------------------------------------------------
        # Added by Alain Chancé
        # Use the total CCSD energy as default for the exact energy
        #-----------------------------------------------------------
        if param['exact_energy'] is None:
            param['exact_energy'] = ccsd.e_tot     # The total CCSD energy
            
            if not param['compute_exact_energy']:
                print(f"The total CCSD energy {ccsd.e_tot} is used as default exact energy")

        #--------------------------
        # Compute the exact energy
        #--------------------------
        if param['compute_exact_energy']:
            try:
                param['exact_energy'] = cas.run().e_tot
                print(f"exact energy: {param['exact_energy']}") 
            except:
                print(f"\nFailed to compute exact energy - using the total CCSD energy {ccsd.e_tot}")
                param['exact_energy'] = ccsd.e_tot 

        n_reps = 1

        if mol.spin == 0:
            #---------------------------------------------------------------------------------
            # Use the spin-balanced variant of the unitary cluster Jastrow (UCJ) ansatz
            # https://qiskit-community.github.io/ffsim/api/ffsim.html#ffsim.UCJOpSpinBalanced
            #---------------------------------------------------------------------------------
            #---------------------------------------------------------------------
            # Create a pair of lists, for alpha-alpha and alpha-beta interactions
            #---------------------------------------------------------------------
            alpha_alpha_indices = [(p, p + 1) for p in range(num_orbitals - 1)]
            alpha_beta_indices = [(p, p) for p in range(0, num_orbitals, 4)]
            
            ucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
                t2=t2, # The t2 amplitudes.
                t1=t1, # The t1 amplitudes
                n_reps=n_reps, # The number of ansatz repetitions
                #---------------------------------------------------------------------------------------------------------
                # Interaction_pairs should be a pair of lists, for alpha-alpha and alpha-beta interactions, in that order
                #---------------------------------------------------------------------------------------------------------
                interaction_pairs=(alpha_alpha_indices, alpha_beta_indices),
            )
        else:
            #-----------------------------------------------------------------------------------
            # Use the spin-unbalanced (local) unitary cluster Jastrow ansatz
            # https://qiskit-community.github.io/ffsim/api/ffsim.html#ffsim.UCJOpSpinUnbalanced
            #-----------------------------------------------------------------------------------
            ucj_op = ffsim.UCJOpSpinUnbalanced.from_t_amplitudes(
                t2=t2, # The t2 amplitudes.
                t1=t1, # The t1 amplitudes
                n_reps=n_reps, # The number of ansatz repetitions
                #-------------------------------------------------------------------------------------------------------------------------
                # interaction_pairs: Optional restrictions on allowed orbital interactions for the diagonal Coulomb operators.
                # If specified, `interaction_pairs` should be a tuple of 3 lists, for alpha-alpha, alpha-beta, and beta-beta interactions, 
                # in that order. Any list can be substituted with ``None`` to indicate no restrictions on interactions.
                # Each list should contain pairs of integers representing the orbitals that are allowed to interact. 
                # These pairs can also be interpreted as indices of diagonal Coulomb matrix entries that are allowed to be nonzero.
                # For the alpha-alpha and beta-beta interactions, each integer pair must be upper triangular, that is, of the form 
                # :math:`(i, j)` where :math:`i \leq j`.
                #-------------------------------------------------------------------------------------------------------------------------
                interaction_pairs=None,
            )

        #--------------------------------------------------------------------------
        # Create the ansatz circuit with Hartree-Fock state as the reference state
        #--------------------------------------------------------------------------
        nelec = (num_elec_a, num_elec_b)   # The numbers of alpha and beta electrons
        self.param['nelec'] = nelec

        #---------------------------------------------------------
        # Compute the number of qubits including ancillary qubits
        #---------------------------------------------------------
        n_qubits = 2*num_orbitals + n_ancillary_qubits
        param['n_qubits'] = n_qubits
        print("\nstep_1 - Number of qubits:", n_qubits)

        #---------------------------------
        # Create an empty quantum circuit
        #---------------------------------
        qubits = QuantumRegister(n_qubits, name="q")
        circuit = QuantumCircuit(qubits)

        #----------------------------------------------------------------------------------------
        # Prepare Hartree-Fock state as the reference state and append it to the quantum circuit
        #----------------------------------------------------------------------------------------
        circuit.append(ffsim.qiskit.PrepareHartreeFockJW(num_orbitals, nelec), qubits)

        #-----------------------------------------------
        # Apply the UCJ operator to the reference state
        #-----------------------------------------------
        if mol.spin == 0:
            circuit.append(ffsim.qiskit.UCJOpSpinBalancedJW(ucj_op), qubits)
        else:
            circuit.append(ffsim.qiskit.UCJOpSpinUnbalancedJW(ucj_op), qubits)
        circuit.measure_all()

        # Decompose the quantum circuit
        circuit = circuit.decompose().decompose()
        self.param['circuit'] = circuit

        # Get the depth of the quantum circuit
        circuit_depth = circuit.depth()
        self.param['circuit_depth'] = circuit_depth
        print("Depth of the quantum_circuit:", circuit_depth)

        # Draw the circuit using matplotlib
        fig = circuit.draw("mpl", fold=-1)

        # Save the figure to a file
        fig.savefig("quantum_circuit.png", bbox_inches='tight')

        # Display the figure
        display(fig)
        
        # Return the figure
        return fig

    #--------------------------------------------------------------------------------------------------------------
    # Define function step_2() - Optimize problem for quantum execution
    #
    # Next, we optimize the circuit for a target hardware. 
    # We generate a staged pass manager using the function generate\_preset\_pass\_manager,
    # https://docs.quantum.ibm.com/api/qiskit/transpiler_preset#generate_preset_pass_manager) 
    # from qiskit with the specified `backend` and `initial_layout`.
    #
    # We set the `pre_init` stage of the staged pass manager to `ffsim.qiskit.PRE_INIT`.
    # It includes qiskit transpiler passes that decompose gates 
    # into orbital rotations and then merges the orbital rotations, resulting in fewer gates in the final circuit.
    #--------------------------------------------------------------------------------------------------------------
    def step_2(self):

        if self.param is None:
            print("step_2: missing parameter param")
            return None
        param = self.param

        #-------------------------------
        # Return if run_on_QPU is False
        #-------------------------------
        if not param['run_on_QPU']:
            return None

        #---------------
        # Setup backend
        #---------------
        self.setup_backend()
        
        if not param['run_on_QPU']:
            return None
        
        #--------------------------------
        # Generate a staged pass manager
        #--------------------------------
        pass_manager = generate_preset_pass_manager(optimization_level=3, backend=self.backend, initial_layout=param['initial_layout'])

        #-----------------------------------------------------------------------
        # Use the circuit generated by this pass manager for hardware execution
        #-----------------------------------------------------------------------
        pass_manager.pre_init = ffsim.qiskit.PRE_INIT

        # Run the pass manager
        try:
            isa_circuit = pass_manager.run(self.param['circuit'])
            print(f"Gate counts (w/ pre-init passes): {isa_circuit.count_ops()}")
        except:
            print("Simulating the run by reading bit array file or Generating 10k samples drawn from the uniform distribution")
            param['run_on_QPU'] = False
            isa_circuit = None

        self.param['isa_circuit'] = isa_circuit

        # Get the depth of the ISA quantum circuit
        if isa_circuit is not None:
            isa_circuit_depth = isa_circuit.depth()
            self.param['isa_circuit_depth'] = isa_circuit_depth
            print("Depth of the ISA quantum circuit:", isa_circuit_depth)

        return isa_circuit

    #----------------------------------------------------------------------------------------------------
    # Define the function get_bit_array() which loads a bit_array from a file or generates random samples 
    # drawn from the uniform distribution
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb
    #-----------------------------------------------------------------------------------------------------
    def get_bit_array(self):

        if self.param is None:
            print("get_bit_array: missing parameter param")
            return None
        param = self.param
        
        load_bit_array_file = param['load_bit_array_file']
        nsamples = param['nsamples']
        
        if load_bit_array_file != None and os.path.isfile(load_bit_array_file):
            print("Reading in samples from bit array file: ", load_bit_array_file)
            # https://numpy.org/doc/2.2/reference/generated/numpy.load.html
            bit_array = np.load(load_bit_array_file, allow_pickle=True).item()
            
        else:
            print("Generating samples drawn from the uniform distribution")
            rng = np.random.default_rng(param['seed'])
            bit_array = generate_bit_array_uniform(nsamples, param['num_orbitals'] * 2, rand_seed=rng)
            
        return bit_array

    #---------------------------------
    # Define function get_QPU_usage()
    #---------------------------------
    def get_QPU_usage(self):

        if self.param is None:
            print("get_QPU_usage: missing parameter param")
            return None, None
        param = self.param

        if param['job'] is None:
            #---------------------
            # Retrieve the job_id
            #---------------------
            job_id = param['job_id']
            if job_id is None:
                return None, None

            #--------------------------
            # Instantiate the service
            #-------------------------
            if self.service is None:
                try:
                    self.service = QiskitRuntimeService()
                except:
                    print(f"Error creating an instance of QiskitRuntimeService(): {e}")
                    return None, None
                    
            service = self.service

            #-----------------------------------------
            # Get the job corresponding to the job_id
            #-----------------------------------------
            try:
                param['job'] = service.job(job_id)
                print(f"Successfully retrieved job with job_id: {job_id}")
            except Exception as e:
                print(f"Error retrieving job with job_id {job_id}: {e}")
                return None, None
        
        job = param['job']
        power_QPU = param['power_QPU']

        try:
            metrics = job.metrics()           # Fetch metrics from server
            usage = metrics.get("usage", {})
            QPU_usage = usage.get("quantum_seconds")
            param['QPU_usage'] = QPU_usage
        except Exception as e:
            print(f"Error retrieving job metrics: {e}")
            return None, None

        if QPU_usage is None:
            return None, None
        print(f"\nQiskit Runtime usage (s): {QPU_usage:.2f}")

        try:
            QPU_power_consumption = power_QPU*QPU_usage/3600.0
            param['QPU_power_consumption'] = QPU_power_consumption
            
            print("\nA rough estimate for QPU power consumption is computed as power QPU (kW) * QPU usage (h)")
            print(f"power QPU (kW): {power_QPU}")
            print(f"QPU usage (s): {QPU_usage:.2f}, (h): {QPU_usage/3600.0:.4f}")
            print(f"Rough estimate for QPU power consumption (kWh): {QPU_power_consumption:.4f}")
            
        except Exception as e:
            print(f"Error computing a rough estimate for QPU power consumption: {e}")
            return QPU_usage, None

        return QPU_usage, QPU_power_consumption

    #---------------------------------------------
    # Define function get_classical_power_usage()
    #---------------------------------------------
    def get_classical_power_usage(self):

        if self.param is None:
            print("get_classical_power_usage: missing parameter param")
            return None, None
        param = self.param

        eco2AI_file = param['eco2ai_file_name']

        if eco2AI_file is None:
            return None, None

        # Read the CSV file
        try:
            df = pd.read_csv(eco2AI_file, sep =',')
        except Exception as e:
            print(f"Error reading the eco2AI file: {e}")
            return None, None

        # Access columns by name
        try:
            # Sum duration (convert seconds → hours)
            duration = df["duration(s)"].sum() / 3600.0
            param['duration'] = duration

            # Sum power consumption (already in kWh)
            classical_power_usage = df["power_consumption(kWh)"].sum()
            param['classical_power_usage'] = classical_power_usage

            print(
                f"\nClassical processing - Duration (h): {duration:.4f} "
                f"- Power consumption (kWh): {classical_power_usage:.4f}"
                )
        except Exception as e:
            print(f"Error accessing the duration and power consumption of the eco2AI file: {e}")
            return None, None
        
        return duration, classical_power_usage

    #-----------------------------------------------------------------------------------------------------------------------------------------
    # Define function step_3() - Execute using Qiskit Primitives
    #
    # We run the circuit on the target hardware, or we just read in 100k samples drawn from ``ibm_brisbane`` at an earlier time 
    # and collect samples for ground state energy estimation of the $N_2$ molecule with basis set `6-31g`.
    #
    # As shown in the tutorial [Improving energy estimation of a chemistry Hamiltonian with SQD]
    # (https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb), 
    # we can use `Fake_Sherbrooke`, a fake 127-qubit backend from qiskit_ibm_runtime to emulate a real device, 
    # and generate random samples drawn from the uniform distribution. This approach is however not recommended, 
    # see [A case against uniform sampling]
    # (https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).
    #-----------------------------------------------------------------------------------------------------------------------------------------
    def step_3(self):

        if self.param is None:
            print("step_3: missing parameter param")
            return None
        param = self.param

        if param['run_on_QPU']:
            #-------------------------------------------
            # Retrieve parameters from param dictionary
            #-------------------------------------------
            isa_circuit = self.param['isa_circuit']
            
            if isinstance(self.backend, AerSimulator):
                #------------------------------
                # Simulating with AerSimulator
                #------------------------------
                print("Simulating with AerSimulator")
                job = self.backend.run([isa_circuit], shots=param['nshots'])
                result = job.result()
                counts = result.get_counts(isa_circuit)
                bit_array = BitArray.from_counts(counts)
            
            else:
                job = param['job']

                if job is not None:
                    try:
                        result = job.result()
                    except Exception as e:
                        print(f"Error retrieving job result: {e}")
                else:
                    #----------------------------------------------------
                    # Running the quantum circuit on the target hardware
                    #----------------------------------------------------
                    print("Running the quantum circuit on the target hardware: ", self.backend.name)

                    # Migrate from backend.run to Qiskit Runtime primitives
                    # https://docs.quantum.ibm.com/migration-guides/qiskit-runtime
                    job = self.sampler.run([isa_circuit], shots=param['nshots'])
                    param['job'] = job
                    
                    param['job_id'] = job.job_id()
                    print("\njob id:", param['job_id'])

                    #-----------------------------------------------------------------------------
                    # Monitor job
                    # https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.providers.JobStatus
                    #-----------------------------------------------------------------------------
                    timeout = param['timeout']
                    poll_interval = param['poll_interval']

                    t0 = time.time()          # start time
                    t1 = t0                   # time when status is QUEUED

                    while True:
                        try:
                            status = job.status()
                        except Exception as e:
                            print(f"Error retrieving job status: {e}")
                            time.sleep(poll_interval)
                            continue

                        if status == "QUEUED" and t1 == t0:
                            t1 = time.time()
                            print(f"Waiting qpu time = {t1 - t0:.2f}, status = {status}")

                        elif status in ["VALIDATING", "RUNNING"]:
                            print(f"status = {status}")

                        elif status in ["CANCELLED", "DONE", "ERROR"]:
                            t2 = time.time()
                            print(f"Executing QPU time = {t2 - t1:.2f}, status = {status}")
                            break

                        if time.time() - t0 > timeout:
                            print("Job monitoring timed out.")
                            break

                        time.sleep(poll_interval)

                    #--------------------------------
                    # Wait until the job is complete
                    #--------------------------------
                    try:
                        result = job.result()
                    except Exception as e:
                        print(f"Error retrieving job result: {e}")
                        result = None

                if result is None:
                    print("\nThe job running the quantum circuit has failed")
                    #---------------------------------------------------------------------------------------------
                    # load a bit_array from a file or generate random samples drawn from the uniform distribution
                    #---------------------------------------------------------------------------------------------
                    bit_array = self.get_bit_array()

                else:
                    # Display Qiskit Runtime usage and power consumption
                    QPU_usage, QPU_power_consumption = self.get_QPU_usage()
                    
                    # Get results for the first (and only) PUB
                    pub_result = result[0]
                                
                    # Get counts from the result
                    # https://docs.quantum.ibm.com/migration-guides/qiskit-runtime-examples#2-get-counts-from-the-result
                    counts = pub_result.data.meas.get_counts()

                    #---------------------------
                    # Get bit_array from counts
                    #---------------------------
                    bit_array = BitArray.from_counts(counts)
                    
            #---------------------
            # Save bit array file
            #---------------------
            save_bit_array_file = param['save_bit_array_file']
            if save_bit_array_file != None:
                print("\nSaving bit array in file: ", save_bit_array_file)
                # https://numpy.org/doc/2.1/reference/generated/numpy.save.html
                np.save(save_bit_array_file, bit_array)
        
        else:
            #---------------------------------------------------------------------------------------------
            # load a bit_array from a file or generate random samples drawn from the uniform distribution
            #---------------------------------------------------------------------------------------------
            bit_array = self.get_bit_array()
        
        self.param['bit_array'] = bit_array
        
        return bit_array

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Define function post_process() - Post-process and return result to desired classical format
    #
    # The first iteration of the self-consistent configuration recovery procedure uses the raw samples after post-selection 
    # on symmetries as input to the diagonalization process to #obtain an estimate of the average orbital occupancies.
    #
    # Subsequent iterations use these occupancies to generate new configurations from raw samples that violate the symmetries 
    # (i.e., are incorrect). These configurations are collected and then subsampled to produce batches, which are subsequently 
    # used to project the Hamiltonian and compute a ground-state estimate using an eigenstate solver.
    #
    # We use the `diagonalize_fermionic_hamiltonian` function defined in fermion.py,
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py.
    #
    # The solver included in the SQD addon uses PySCF's implementation of selected CI, 
    # specifically pyscf.fci.selected_ci.kernel_fixed_space,
    # https://pyscf.org/pyscf_api_docs/pyscf.fci.html#pyscf.fci.selected_ci.kernel_fixed_space)
    #
    # The following parameters are user-tunable:
    # * `max_iterations`: Limit on the number of configuration recovery iterations.
    # * `num_batches`: The number of batches of configurations to subsample (i.e., the number of separate calls to the eigenstate solver).
    # * `samples_per_batch`: The number of unique configurations to include in each batch.# 
    # * `max_cycles`: The maximum number of Davidson cycles run by the eigenstate solver.
    #
    # * `occupancies_tol`: Numerical tolerance for convergence of the average orbital occupancies. If the maximum change in absolute value 
    # of the average occupancy of an orbital between iterations is smaller than this value, then the configuration recovery loop will exit, 
    # if the energy has also converged (see the ``energy_tol`` argument).
    #
    # * `energy_tol`: Numerical tolerance for convergence of the energy. If the change in energy between iterations is smaller than this value, 
    # then the configuration recovery loop will exit, if the occupancies have also converged (see the ``occupancies_tol`` argument).
    #-------------------------------------------------------------------------------------------------------------------------------------------
    def post_process(self):

        if self.param is None:
            print("post_process: missing parameter param")
            return None
        param = self.param
        
        #--------------------------------------
        # Get parameters from param dictionary
        #--------------------------------------
        bit_array = param['bit_array']

        #---------------------
        # Molecular integrals
        #---------------------
        hcore = param['hcore']                                         # The one-body tensor of the Hamiltonian
        eri = param['eri']                                             # The two-body tensor of the Hamiltonian
        nelec = param['nelec']                                         # The numbers of alpha and beta electrons.
        nuclear_repulsion_energy = param['nuclear_repulsion_energy']   # Nuclear repulsion energy
        num_orbitals = param['num_orbitals']                           # The number of spatial orbitals

        #-------------
        # SQD options
        #-------------
        energy_tol = param['energy_tol']                    # Numerical tolerance for convergence of the energy
        occupancies_tol = param['occupancies_tol']          # Numerical tolerance for convergence of the average orbital occupancies
        max_iterations = param['max_iterations']            # Limit on the number of configuration recovery iterations

        #---------------------------
        # Eigenstate solver options
        #---------------------------
        num_batches = param['num_batches']                   # The number of batches of configurations to subsample
        samples_per_batch = param['samples_per_batch']       # The number of bitstrings to include in each subsampled batch of bitstrings
        symmetrize_spin = param['symmetrize_spin']           # Whether to always merge spin-alpha and spin-beta CI strings into a single list
        carryover_threshold = param['carryover_threshold']   # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
        max_cycle = param['max_cycle']                       # The maximum number of Davidson cycles run by the eigenstate solver
        rng = np.random.default_rng(param['seed'])           # A seed for the pseudorandom number generator
        spin_sq = param['spin_sq']

        #-------------------------------------------------------------------------------------
        # Pass options to the built-in eigensolver. If you just want to use the defaults,
        # you can omit this step, in which case you would not specify the sci_solver argument
        # in the call to diagonalize_fermionic_hamiltonian below.
        # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
        #-------------------------------------------------------------------------------------
        sci_solver = partial(solve_sci_batch, spin_sq=spin_sq, max_cycle=max_cycle)

        # Initialize result_history list to capture intermediate results
        result_history = [] 

        #--------------------------
        # Define callback function
        #--------------------------
        def callback(results: list[SCIResult]): 
            result_history.append(results)
            iteration = len(result_history)
            print(f"Iteration {iteration}")
            for i, result in enumerate(results):
                print(f"\tSubsample {i}")
                print(f"\t\tEnergy: {result.energy + nuclear_repulsion_energy}")
                print(f"\t\tSubspace dimension: {np.prod(result.sci_state.amplitudes.shape)}")

        #---------------------------------------------------------------------------------------------------------------------------------------
        # Diagonalize fermionic hamiltonian
        # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
        # Output: List of (energy, sci_state, occupancies) triplets, where each triplet contains the result of the corresponding diagonalization
        #----------------------------------------------------------------------------------------------------------------------------------------
        result = diagonalize_fermionic_hamiltonian(
            hcore,                                   # The one-body tensor of the Hamiltonian
            eri,                                     # The two-body tensor of the Hamiltonian
            bit_array,                               # Array of sampled bitstrings
            samples_per_batch=samples_per_batch,     # The number of bitstrings to include in each subsampled batch of bitstrings
            norb=num_orbitals,                       # The number of spatial orbitals
            nelec=nelec,                             # The numbers of alpha and beta electrons
            num_batches=num_batches,                 # The number of batches of configurations to subsample
            energy_tol=energy_tol,                   # Numerical tolerance for convergence of the energy
            occupancies_tol=occupancies_tol,         # Numerical tolerance for convergence of the average orbital occupancies
            max_iterations=max_iterations,           # Limit on the number of configuration recovery iterations
            sci_solver=sci_solver,                   # Selected configuration interaction solver function
            symmetrize_spin=symmetrize_spin,         # Whether to always merge spin-alpha and spin-beta CI strings into a single list
            carryover_threshold=carryover_threshold, # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            callback=callback,                       # A callback function to be called after each configuration recovery
            seed=rng,                                # A seed for the pseudorandom number generator
        )

        self.param['result'] = result
        self.param['result_history'] = result_history

        #----------------------------------
        # Return result and result_history
        #----------------------------------
        return(result, result_history)

    #---------------------------------------------
    # Define function plot_energy_and_occupancy()
    #---------------------------------------------
    def plot_energy_and_occupancy(self):

        if self.param is None:
            print("plot_energy_and_occupancy: missing parameter param")
            return None
        param = self.param

        #--------------------------------------
        # Get parameters from param dictionary
        #--------------------------------------
        result = param['result']
        result_history = param['result_history']
        exact_energy = param['exact_energy']
        nuclear_repulsion_energy = param['nuclear_repulsion_energy']
        chem_accuracy = param['chem_accuracy'] # Chemical accuracy (+/- 1 milli-Hartree)

        #------------------------
        # Data for energies plot
        #------------------------
        x1 = range(len(result_history))
        min_e = [
            min(result, key=lambda res: res.energy).energy + nuclear_repulsion_energy
            for result in result_history
        ]
        e_diff = [abs(e - exact_energy) for e in min_e]
        yt1 = [1.0, 1e-1, 1e-2, 1e-3, 1e-4]

        #----------------------------------------
        # Data for avg spatial orbital occupancy
        #----------------------------------------
        y2 = np.sum(result.orbital_occupancies, axis=0)
        x2 = range(len(y2))
    
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    
        # Plot energies
        axs[0].plot(x1, e_diff, label="energy error", marker="o")
        axs[0].set_xticks(x1)
        axs[0].set_xticklabels(x1)
        axs[0].set_yticks(yt1)
        axs[0].set_yticklabels(yt1)
        axs[0].set_yscale("log")
        axs[0].set_ylim(min(yt1))
        axs[0].axhline(y=chem_accuracy, color="#BF5700", linestyle="--", label="chemical accuracy")
        axs[0].set_title("Approximated Ground State Energy Error vs SQD Iterations")
        axs[0].set_xlabel("Iteration Index", fontdict={"fontsize": 12})
        axs[0].set_ylabel("Energy Error (Ha)", fontdict={"fontsize": 12})
        axs[0].legend()

        #------------------------
        # Plot orbital occupancy
        #------------------------
        axs[1].bar(x2, y2, width=0.8)
        axs[1].set_xticks(x2)
        axs[1].set_xticklabels(x2)
        axs[1].set_title("Avg Occupancy per Spatial Orbital")
        axs[1].set_xlabel("Orbital Index", fontdict={"fontsize": 12})
        axs[1].set_ylabel("Avg Occupancy", fontdict={"fontsize": 12})

        param['SQD_energy'] = min_e[-1]
        param['Absolute_error'] = e_diff[-1]
        
        print(f"Exact energy (Ha): {exact_energy:.5f}")
        print(f"SQD energy (Ha): {min_e[-1]:.5f}")
        print(f"Absolute error (Ha): {e_diff[-1]:.5f}")

        #-----------
        # Save plot
        #-----------
        plt.tight_layout()
        plt.savefig("plot_energy_and_occupancy.png")

        param['plt_energy'] = plt

        #--------------
        # Display plot
        #--------------
        display(plt.show())

        #-------------------------
        # Stop the eco2AI tracker
        #-------------------------
        tracker = self.param['eco2ai_tracker']
        if tracker is not None:
            tracker.stop()

        # Return the figure
        return plt

    #---------------------------------------
    # Define function run() - Run all steps
    #---------------------------------------
    def run(self):

        if self.param is None:
            print("run: missing parameter param")
            return None
        param = self.param

        config_filename = param['config_filename']
        basis = param['basis']
        
        if config_filename is not None:
            print(f"\nFind an approximation to the molecule defined in configuration {config_filename} in the {basis} basis set")

        t0 = time.time()  # ⏱️ Start timing

        #--------
        # STEP 1
        #--------
        print("\nstep_1 - Perform a CCSD calculation")
        try:
            plt_circuit = self.step_1()
            t1 = time.time()
            print(f"✅ step_1 completed in {t1 - t0:.2f} seconds")
        except Exception as e:
            return(f"Error in step 1: {e}")

        #--------
        # STEP 2
        #--------
        print("\nstep_2 - Optimize the circuit for a target hardware")
        try:
            isa_circuit = self.step_2()
            t2 = time.time()
            print(f"✅ step_2 completed in {t2 - t1:.2f} seconds")
        except Exception as e:
            return(f"Error in step 2: {e}")

        #--------
        # STEP 3
        #--------
        print("\nstep_3 - Execute using Qiskit Primitives or generate random samples")
        try:
            bit_array = self.step_3()
            t3 = time.time()
            print(f"✅ step_3 completed in {t3 - t2:.2f} seconds")
        except Exception as e:
            return(f"Error in step 3: {e}")

        #--------------
        # POST PROCESS
        #--------------
        print("\npost_process - Self-consistent configuration recovery procedure")
        result, result_history = self.post_process()
        t4 = time.time()
        print(f"✅ post_process completed in {t4 - t3:.2f} seconds")

        #---------------------------
        # Plot energy and occupancy
        #---------------------------
        print("\nplot_energy_and_occupancy")
        try:
            plt_energy = self.plot_energy_and_occupancy()
            t5 = time.time()
            print(f"✅ plot_energy_and_occupancy completed in {t5 - t4:.2f} seconds")
        except Exception as e:
            return(f"Error in plot_energy_and_occupancy: {e}")

        print(f"\nSQD configuration saved into {config_filename} and simulation complete.")
        print(f"\nExact energy (Ha): {self.param['exact_energy']:.5f}")
        print(f"SQD energy (Ha): {self.param['SQD_energy']:.5f}")
        print(f"Absolute error (Ha): {self.param['Absolute_error']:.5f}")

        #------------------------------------------------------------------------------------
        # Print a rough estimate of the energy consumption of the quantum device
        # **Assumption**
        # A ballpark figure for a typical modern IBM-class superconducting quantum computer 
        # (including cryogenics + support, while idle or lightly used): ~ 15–25 kW. 
        # Source: [Green quantum computing, Capgemini, 8 May 2023]
        # (https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/).
        #------------------------------------------------------------------------------------
        QPU_usage, QPU_power_consumption = self.get_QPU_usage()

        #----------------------------------------------------------------------------------------
        # Read the CSV file and print the duration and power consumption of classical processing
        #----------------------------------------------------------------------------------------
        duration, classical_power_usage = self.get_classical_power_usage()
    
        return "simulation complete"