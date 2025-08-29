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
#------------------------------------------------------------------------------------------------------------------------

# Import common packages first
import psutil
import sys
import os
import time
import numpy as np
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

# Import qiskit classes
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
                 backend_name = None,                              # IBM cloud backend name
                 do_plot_gate_map = True,                          # Whether to plot the gate map 
                 load_bit_array_file = None,                       # If provided, function step_3 will load samples from this file
                 save_bit_array_file = None,                       # If provided, function step_3 will save samples into this file
                 n_ancillary_qubits = 0,                           # Number of ancillary qubits
                 run_on_QPU = False,                               # Whether to run the quantum circuit on the target hardware
                 nshots = 1000,                                    # Number of shots
                 opt_level = 1,                                    # Optimization level
                 poll_interval = 5,                                # Poll interval in seconds for job monitor
                 timeout = 600,                                    # Time out for job monitor
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
                 compute_exact_energy = True,                      # Whether to compute the exact energy with PySCF, cas.run().e_tot
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
        print("Backend name:", backend_name)
        print("do plot gate map: ", do_plot_gate_map)
        
        if load_bit_array_file != None:
            print("load_bit_array_file:", load_bit_array_file)
            
        if save_bit_array_file != None:
            print("save_bit_array_file:", save_bit_array_file)
            
        print("Run on QPU:", run_on_QPU)

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
        print("Compute exact energy: ", compute_exact_energy)

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
            "circuit": None,                                 # Quantum circuit returned by function step_1
            "n_qubits": 2,                                   # Number of qubits set by function step_1
            "isa_circuit": None,                             # ISA circuit returned by function step_2
            "bit_array": None,                               # Bit array returned by function step_3
            "result": None,                                  # result returned by post_process function
            "result_history": None,                          # result_history returned by post_process function
            "SQD_energy": None,                              # SQD energy returned by post_process function
            "Absolute_error": None,                          # Absolute error returned by post_process function
            "do_plot_gate_map": do_plot_gate_map,            # Whether to plot the gate map  
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
            "backend_name": backend_name,                    # IBM cloud backend name
            "spin_a_layout": [],                             # Spin a layout set by __init__
            "spin_b_layout": [],                             # Spin b layout set by __init__
            "n_ancillary_qubits": n_ancillary_qubits,        # Number of ancillary qubits
            "run_on_QPU": run_on_QPU,                        # Whether to run the quantum circuit on the target hardware
            "nshots": nshots,                                # Number of shots
            'opt_level':opt_level,                           # Optimization level
            "poll_interval": poll_interval,                  # Poll interval in seconds for job monitor
            "timeout": timeout,                              # Time out in seconds for gob monitor
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            "mol": pyscf.gto.Mole(),                         # Molecule
            "basis": basis,                                  # Basis set
            "spin": spin,                                    # Spin
            "charge": charge,                                # Charge
            "atom": atom,                                    # Atom configuration, PySCF Initializing a molecule
            "symmetry": symmetry,                            # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry 
            "n_frozen": n_frozen,                            # Number of frozen orbitals
            #---------------------
            # Molecular integrals
            #---------------------
            "num_orbitals": 2,                                # The number of spatial orbitals
            "hcore": 0.0,                                     # The one-body tensor of the Hamiltonian
            "nuclear_repulsion_energy": 0.0,                  # Nuclear repulsion energy
            "eri": 0.0,                                       # The two-body tensor of the Hamiltonian
            "nelec": 1,                                       # The numbers of alpha and beta electrons
            "compute_exact_energy": compute_exact_energy,     # Whether to compute the exact energy with PySCF, cas.run().e_tot
            "exact_energy": 0.0,                              # Exact energy
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

    #-------------------------------------------------------------------------------
    # Define function build_molecule
    # Define and initialize a molecule with pySCF
    # PySCF Initializing a molecule
    # https://pyscf.org/user/gto.html#initializing-a-molecule
    # https://pyscf.org/pyscf_api_docs/pyscf.gto.html#pyscf.gto.mole.MoleBase.build
    #-------------------------------------------------------------------------------
    def build_molecule(self):
        param = self.param
        mol = pyscf.gto.Mole()
        mol.build(
            atom = param['atom'],
            basis = param['basis'],
            spin = param['spin'],
            charge = param['charge'],
            symmetry = param['symmetry']
        )
        self.param['mol'] = mol
        return
    
    def setup_backend(self):
        #-------------------------------------
        # Define the function setup_backend() 
        #-------------------------------------
        #-------------------------------------------------------------------------------------------
        # Instantiate the service
        # Once the account is saved on disk, you can instantiate the service without any arguments:
        # https://docs.quantum.ibm.com/api/migration-guides/qiskit-runtime
        #-------------------------------------------------------------------------------------------
        try:
            service = QiskitRuntimeService()
        except:
            service = None
        
        backend_name = self.param['backend_name']
        n_qubits = self.param['n_qubits']
        
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
                    self.param['backend_name'] = self.backend.name

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
                        self.param['backend_name'] = self.backend.name
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
                self.param['backend_name'] = self.backend.name
                self.sampler = StatevectorSampler()
                print("\nUsing AerSimulator with method statevector and noise model from", backend_name)

            if not isinstance(self.backend, AerSimulator): 
                print(f"Backend name: {self.backend.name}\n"
                      f"Version: {self.backend.version}\n"
                      f"Number of qubits: {self.backend.num_qubits}")
        
        if isinstance(self.backend, AerSimulator):
            # Check that there is enough memory to perform a simulation with AerSimulator
            self.check_size()
            
        else:
            # Setup the zig-zag pattern for qubit interactions
            self.setup_zig_zag()

            # Plot gate map
            self.plot_gate_map()
            
            print(f"Backend name: {self.backend.name}\n"
                  f"Version: {self.backend.version}\n"
                  f"Number of qubits: {self.backend.num_qubits}"
                 )

    def setup_zig_zag(self):
        #--------------------------------------------------
        # Setup the zig-zag pattern for qubit interactions
        #--------------------------------------------------
        if isinstance(self.backend, AerSimulator):
            return

        if self.backend.num_qubits == 127:
            spin_a_layout = [0,14,18,19,20,33,39,40,41,53,60,61,62,72,81,82,83,92,102,103,104,111,122,123,124]
            spin_b_layout = [2,3,4,15,22,23,24,34,43,44,45,54,64,65,66,73,85,86,87,93,106,107,108,112,126]
        
        elif self.backend.num_qubits == 133:
            spin_a_layout = [2,3,4,16,23,24,25,35,44,45,46,55,65,66,67,74,86,87,88,94,107,108,109,113,128,127]
            spin_b_layout = [0,15,19,20,21,34,40,41,42,54,61,62,63,73,82,83,84,93,103,104,105,112,124,123,122,121]

        else:
            spin_a_layout=[]
            spin_b_layout=[]
        
        self.param['spin_a_layout'] = spin_a_layout
        self.param['spin_b_layout'] = spin_b_layout
        
        print("spin_a_layout: ", spin_a_layout)
        print("spin_b_layout: ", spin_b_layout)
        
        return

    def plot_gate_map(self):
        #--------------------------------------------------------------------------------------------------
        # Plot gate map for the backend with the spin-a layout shown in red, and the spin-b layout in blue
        # The 'Graphviz' library is required to use 'plot_coupling_map'.  
        # To install, follow the instructions at https://graphviz.org/download/
        # On Debian, Ubuntu, run the following command: 
        # sudo apt install graphviz
        #--------------------------------------------------------------------------------------------------
        param = self.param
        
        spin_a_layout = param['spin_a_layout']
        spin_b_layout = param['spin_b_layout']
        
        if param['do_plot_gate_map']:
            try:
                qubit_color = []
                for i in range(self.backend.num_qubits):
                    if i in spin_a_layout:
                        qubit_color.append("#FFA07A") # light red
                    elif i in spin_b_layout:
                        qubit_color.append("#87CEEB") # sky blue
                    else:
                        qubit_color.append("#888888") # gray
                
                print("\nGate map for :", self.backend.name, ", the spin-a layout is shown in red, and the spin-b layout in blue.\n")
                
                fig = plot_gate_map(self.backend, qubit_color=qubit_color, plot_directed=False)
                
                # Display the gate map
                display(fig)

                # Save the figure to a file
                fig.savefig("plot_gate_map.png", bbox_inches='tight')

                # Close the figure to free memory
                plt.close(fig)
            except:
                print("The 'Graphviz' library is required to use 'plot_coupling_map'")
                print(sys.exc_info())
            
        return

    #---------------------------------------------------------------------------------------------------------------------
    # Define the function check_size() which checks that there is enough memory to perform a simulation with AerSimulator
    #---------------------------------------------------------------------------------------------------------------------
    def check_size(self):

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

        if self.param == None:
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

        #---------------------------------
        # Compute exact energy
        # Default to converged SCF energy
        #---------------------------------
        param['exact_energy'] = scf
        
        if self.param['compute_exact_energy']:
            try:
                param['exact_energy'] = cas.run().e_tot
                print("exact energy: ", param['exact_energy']) 
            except:
                print("\nstep_1 - Failed to compute exact energy - Default to converged SCF energy")

        #----------------------------------------------------------------
        # Get CCSD calculation t2 amplitudes for initializing the ansatz
        # https://en.wikipedia.org/wiki/Coupled_cluster#Cluster_operator
        #----------------------------------------------------------------
        ccsd = pyscf.cc.CCSD(scf, frozen=[i for i in range(mol.nao_nr()) if i not in active_space]).run()
        t1 = ccsd.t1
        t2 = ccsd.t2
        print("\n")

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
        print("step_1 - Number of qubits:", n_qubits)

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

        # Draw the circuit using matplotlib
        fig = circuit.draw("mpl", fold=-1)

        # Display the figure
        display(fig)

        # Save the figure to a file
        fig.savefig("quantum_circuit.png", bbox_inches='tight')
        
        # Return the figure
        return fig

    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Define function step_2() - Optimize problem for quantum execution
    #
    # Next, we optimize the circuit for a target hardware. We generate a staged pass manager using the function generate\_preset\_pass\_manager,
    # https://docs.quantum.ibm.com/api/qiskit/transpiler_preset#generate_preset_pass_manager) from qiskit with the specified `backend` 
    # and `initial_layout`.
    #
    # We set the `pre_init` stage of the staged pass manager to `ffsim.qiskit.PRE_INIT`. It includes qiskit transpiler passes that decompose gates 
    # into orbital rotations and then merges the orbital rotations, resulting in fewer gates in the final circuit.
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    def step_2(self):

        if self.param == None:
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

        #------------------------
        # Get backend parameters 
        #------------------------
        spin_a_layout = param['spin_a_layout']
        spin_b_layout = param['spin_b_layout']

        k = param['n_qubits'] // 2
        
        OK = True
        if k <= len(spin_a_layout):
            spin_a_layout = spin_a_layout[:k]
        else:
            print("step_2 - provided spin_a_layout is too small")
            OK = False

        if k <= len(spin_b_layout):
            spin_b_layout = spin_b_layout[:k]
        else:
            print("step_2 - provided spin_b_layout is too small")
            OK = False

        if OK:
            initial_layout = spin_a_layout + spin_b_layout
            #--------------------------------
            # Generate a staged pass manager
            #--------------------------------
            print("Generating a staged pass manager with initial_layout = spin_a_layout + spin_b_layout")
            pass_manager = generate_preset_pass_manager(optimization_level=3, backend=self.backend, initial_layout=initial_layout)
        else:
            print("Generating a staged pass manager with no initial_layout")
            pass_manager = generate_preset_pass_manager(optimization_level=3, backend=self.backend)

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

        return isa_circuit

    #----------------------------------------------------------------------------------------------------
    # Define the function get_bit_array() which loads a bit_array from a file or generates random samples 
    # drawn from the uniform distribution
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb
    #-----------------------------------------------------------------------------------------------------
    def get_bit_array(self):

        param = self.param
        
        load_bit_array_file = param['load_bit_array_file']
        
        if load_bit_array_file != None and os.path.isfile(load_bit_array_file):
            print("Reading in samples from bit array file: ", load_bit_array_file)
            # https://numpy.org/doc/2.2/reference/generated/numpy.load.html
            bit_array = np.load(load_bit_array_file, allow_pickle=True).item()
            
        else:
            print("Generating 10k samples drawn from the uniform distribution")
            rng = np.random.default_rng(param['seed'])
            bit_array = generate_bit_array_uniform(10_000, param['num_orbitals'] * 2, rand_seed=rng)
            
        return bit_array

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

        if self.param == None:
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
                #----------------------------------------------------
                # Running the quantum circuit on the target hardware
                #----------------------------------------------------
                print("Running the quantum circuit on the target hardware: ", self.backend.name)

                # Migrate from backend.run to Qiskit Runtime primitives
                # https://docs.quantum.ibm.com/migration-guides/qiskit-runtime
                job = self.sampler.run([isa_circuit], shots=param['nshots'])
                print("\njob id:", job.job_id())

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

                if result is not None:
                    # Get results for the first (and only) PUB
                    pub_result = result[0]
                                
                    # Get counts from the result
                    # https://docs.quantum.ibm.com/migration-guides/qiskit-runtime-examples#2-get-counts-from-the-result
                    counts = pub_result.data.meas.get_counts()

                    #---------------------------
                    # Get bit_array from counts
                    #---------------------------
                    bit_array = BitArray.from_counts(counts)
                
                else:
                    #---------------------------------------------------------------------------------------------
                    # load a bit_array from a file or generate random samples drawn from the uniform distribution
                    #---------------------------------------------------------------------------------------------
                    bit_array = self.get_bit_array()
                    
            #---------------------
            # Save bit array file
            #---------------------
            save_bit_array_file = param['save_bit_array_file']
            if save_bit_array_file != None:
                print("Saving bit array in file: ", save_bit_array_file)
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

        if self.param == None:
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

        if self.param == None:
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

        #--------------
        # Display plot
        #--------------
        param['SQD_energy'] = min_e[-1]
        param['Absolute_error'] = e_diff[-1]
        
        print(f"Exact energy: {exact_energy:.5f} Ha")
        print(f"SQD energy: {min_e[-1]:.5f} Ha")
        print(f"Absolute error: {e_diff[-1]:.5f} Ha")
        
        plt.tight_layout()
        display(plt.show())

        #-----------
        # Save plot
        #-----------
        plt.savefig("plot_energy_and_occupancy.png")

        # Return the figure
        return plt

    #---------------------------------------
    # Define function run() - Run all steps
    #---------------------------------------
    def run(self):
        
        print("\nstep_1 - Perform a CCSD calculation")
        try:
            plt_circuit = self.step_1()
        except Exception as e:
            return(f"Error in step 1: {e}", None)
        
        print("\nstep_2 - Optimize the circuit for a target hardware")
        try:
            isa_circuit = self.step_2()
        except Exception as e:
            return(f"Error in step 2: {e}", None)
        
        print("\nstep_3 - Execute using Qiskit Primitives or generate random samples drawn from the uniform distribution ")
        try:
            bit_array = self.step_3()
        except Exception as e:
            return(f"Error in step 3: {e}", None)
        
        print("\npost_process - Self-consistent configuration recovery procedure")
        try:
            result, result_history = self.post_process()
        except Exception as e:
            return(f"Error in post_process: {e}", None)
        
        print("\nplot_energy_and_occupancy")
        try:
            plt = self.plot_energy_and_occupancy()
        except Exception as e:
            return(f"Error in plot_energy_and_occupancy: {e}", None)
        
        return result, result_history