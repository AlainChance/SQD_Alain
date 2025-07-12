# SQD_Alain
Find an approximation to the ground state of a chemistry Hamiltonian with Sample-based Quantum Diagonalization (SQD).

The Python file `SQD_Alain.py` and the Jupyter notebooks `SQD_N2.ipynb`, `SQD_H2O.ipynb`, `SQD_LiH.ipynb`, `SQD_CH2_Triplet.ipynb`, `SQD_CH2_Singlet.ipynb` are compatible with Python 3.13, Qiskit v2.1, Qiskit runtime version: 0.40 and Qiskit Runtime V2 primitives. 
|||
|-|-|
|**Author:** |Alain Chancé|
|**Date:** |July 1, 2025|
|**Version:** |**1.00**<br/>*Details see at the end of this notebook*|
|**Credit:**|
The Python file `SQD_Alain.py` and the Jupyter notebooks `SQD_N2.ipynb`, `SQD_H2O.ipynb`, `SQD_LiH.ipynb`, `SQD_CH2_Triplet.ipynb`, `SQD_CH2_Singlet.ipynb` are derived from the tutorial [Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb) which is distributed under the [Apache License 2.0](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/LICENSE.txt).
|**References:**|
[Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb)
[Sample-based Quantum Diagonalization (SQD), Quantum diagonalization algorithms](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview)
[Quantum diagonalization algorithms, IBM Quantum Platform, Quantum Learning](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms)
[Ieva Liepuoniute, Kirstin D. DoneyJavier Robledo MorenoJoshua A. JobWilliam S. FriendGavin O. Jones, Quantum-Centric Computational Study of Methylene Singlet and Triplet States, J. Chem. Theory Comput. 2025, 21, 10, 5062–5070](https://doi.org/10.1021/acs.jctc.5c00075)
[Antonio Mezzacapo: Quantum diagonalization methods for lattice models, IQuS - The InQubator For Quantum Simulation, Feb 19, 2025](https://www.youtube.com/watch?v=b1fhh71hY2g)
[J. Robledo-Moreno et al., "Chemistry Beyond Exact Solutions on a Quantum-Centric Supercomputer" (2024), arXiv:quant-ph/2405.05068](https://arxiv.org/abs/2405.05068)
[M. Motta et al., “Bridging physical intuition and hardware efficiency for correlated electronic states: the local unitary cluster Jastrow ansatz for electronic structure” (2023), Chem. Sci., 2023, 14, 11213](https://pubs.rsc.org/en/content/articlehtml/2023/sc/d3sc02516k)
[Quantum-Selected Configuration Interaction: classical diagonalization of Hamiltonians in subspaces selected by quantum computers, arXiv:2302.11320, quant-ph](https://doi.org/10.48550/arXiv.2302.11320)
[Introduction to Qiskit patterns](https://quantum.cloud.ibm.com/docs/en/guides/intro-to-patterns)
[Keeper L. Sharkey and Alain Chancé, Quantum Chemistry and Computing for the Curious: Illustrated with Python and Qiskit® code, Packt 2022](https://a.co/d/6YCwgPb)
[Companion Jupyter notebook for Chapter 5 of the book Quantum Chemistry and Computing for the Curious: Illustrated with Python and Qiskit® code]( https://github.com/AlainChance/Quantum-Chemistry-and-Computing-for-the-Curious/blob/main/Chapter_05_Variational_Quantum_Eigensolver_(VQE)_algorithm_V4.ipynb)
<br/>

## SQD process using Qiskit patterns
Using the [Qiskit patterns](https://quantum.cloud.ibm.com/docs/en/guides/intro-to-patterns) framework, the SQD process can be described in four steps, [arXiv:2405.05068
](https://arxiv.org/abs/2405.05068):

$$\begin{array}{|c|c|c|}
\hline
\text{Step} &\text{Purpose} &\text{Method}\\
\hline
\text{1} &\text{Map classical inputs to a quantum problem} &\text{Generate an ansatz for estimating the ground state}\\
\hline
\text{2} &\text{Optimize problem for quantum execution} &\text{Transpile the ansatz for the backend}\\
\hline
\text{3} &\text{Execute experiments using Qiskit Primitives} &\text{Draw samples from the ansatz using the Sampler primitive}\\
\hline
\text{4} &\text{Post-process and return results to desired classical format} &\text{Self-consistent configuration recovery loop}\\
\hline
\end{array}$$

The probabilistic self-consistent configuration recovery is an iterative procedure that runs in a loop which comprises the following steps:
- Post-process the full set of bitstring samples, using prior knowledge of particle number and the average orbital occupancy calculated on the most recent iteration.
- Probabilistically create batches of subsamples from recovered bitstrings.
- Project and diagonalize the molecular Hamiltonian over each sampled subspace.
- Save the minimum ground state energy found across all batches and update the average orbital occupancy.

The SQD workflow with self-consistent configuration recovery is depicted in the following diagram.

![SQD diagram](https://raw.githubusercontent.com/Qiskit/qiskit-addon-sqd/7fcec2a686bfe115560db20b840cf2b185ae06f6/docs/_static/images/sqd_diagram.png)

Source: [Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb)

### The Jupyter notebook [SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_N2.ipynb) finds an approximation to the ground state of the $N_2$ molecule
The `SQD class` is imported from [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) and instantiated as follows:
```
from SQD_Alain import SQD

SQD_N2 = SQD(
            #-------------
            # Run options
            #-------------
            backend_name = "ibm_brisbane",                            # IBM cloud backend name
            do_plot_gate_map = True,                                  # Whether to plot the gate map 
            load_bit_array_file = "N2_ibm_brisbane_bitarray.npy",     # If provided, function step_3 will load samples from this file
            save_bit_array_file = None,                               # If provided, function step_3 will save samples into this file
            n_ancillary_qubits = 0,                                   # Number of ancillary qubits
            run_on_QPU = False,                                       # Whether to run the quantum circuit on the target hardware
            nshots = 1000,                                            # Number of shots
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis = "6-31g",                                          # Basis set
            atom = [["N", (0.0, 0.0, 0.0)], ["N", (1.0, 0.0, 0.0)]],  # Atom configuration, PySCF Initializing a molecule
            symmetry = "Dooh",                                        # linear molecular symmetries, https://pyscf.org/user/gto.html#point-group-symmetry
            n_frozen = 2,                                             # Number of frozen orbitals
            compute_exact_energy = True,                              # Whether to compute the exact energy with PySCF, cas.run().e_tot
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            chem_accuracy = 1e-3,                                     # Chemical accuracy (+/- 1 milli-Hartree)
            energy_tol = 3e-5,                                        # Numerical tolerance for convergence of the energy
            occupancies_tol = 1e-3,                                   # Numerical tolerance for convergence of the average orbital occupancies
            max_iterations = 15,                                      # Limit on the number of configuration recovery iterations
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches = 5,                                          # The number of batches to subsample in each configuration recovery iteration
            samples_per_batch = 300,                                  # The number of bitstrings to include in each subsampled batch of bitstrings
            symmetrize_spin = True,                                   # Whether to always merge spin-alpha and spin-beta CI strings into a single list
            carryover_threshold = 1e-4,                               # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            max_cycle = 200,                                          # The maximum number of Davidson cycles run by the eigenstate solver
            seed = 24,                                                # A seed for the pseudorandom number generator
            spin_sq = 0.0)                                            # spin square
```
The last cell displays the following two plots:

![plot_energy_and_occupancy_N2](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_N2.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. The usual chemical accuracy is typically defined as `1 kcal/mol` $\approx$ `1.6 mHa`. However in the context of quantum chemistry, chemical accuracy is often approximated as `1 mHa`. The estimation of the energy can be improved by drawing more samples from the circuit, increasing the number of samples per batch and increasing the number of iterations.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

### The Jupyter notebook [SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_H2O.ipynb) finds an approximation to the ground state of the $H2O$ molecule
For $H2O$ in 6-31G, the FCI ground state energy is approximately `–76.20` Hartree. The experimental non-relativistic Born–Oppenheimer energy is about `–76.438` Hartree.

The `SQD class` is imported from [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) and instantiated as follows:
```
from SQD_Alain import SQD

SQD_H2O = SQD(
            #-------------
            # Run options
            #-------------
            backend_name = None,                                      # IBM cloud backend name
            do_plot_gate_map = True,                                  # Whether to plot the gate map 
            load_bit_array_file = None,                               # If provided, function step_3 will load samples from this file
            save_bit_array_file = None,                               # If provided, function step_3 will save samples into this file
            n_ancillary_qubits = 0,                                   # Number of ancillary qubits
            run_on_QPU = False,                                       # Whether to run the quantum circuit on the target hardware
            nshots = 1000,                                            # Number of shots
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis = "6-31g",                                          # Basis set
            atom = [["O", (0.0, 0.0, 0.0)], ["H", (0.0, 1.0, 0.0)], ["H", (0.0, 0.0, 1.0)]], # Atom configuration, PySCF Initializing a molecule
            symmetry = True,                                          # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry
            n_frozen = 0,                                             # Number of frozen orbitals
            compute_exact_energy = True,                              # Whether to compute the exact energy with PySCF, cas.run().e_tot
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            chem_accuracy = 1e-3,                                     # Chemical accuracy (+/- 1 milli-Hartree)
            energy_tol = 3e-5,                                        # Numerical tolerance for convergence of the energy
            occupancies_tol = 1e-3,                                   # Numerical tolerance for convergence of the average orbital occupancies
            max_iterations = 10,                                      # Limit on the number of configuration recovery iterations
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches = 5,                                          # The number of batches to subsample in each configuration recovery iteration
            samples_per_batch = 300,                                  # The number of bitstrings to include in each subsampled batch of bitstrings
            symmetrize_spin = True,                                   # Whether to always merge spin-alpha and spin-beta CI strings into a single list
            carryover_threshold = 1e-4,                               # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            max_cycle = 200,                                          # The maximum number of Davidson cycles run by the eigenstate solver
            seed = 24,                                                # A seed for the pseudorandom number generator
            spin_sq = 0.0)                                            # spin square
```
The last cell displays the following two plots:

![plot_energy_and_occupancy_H2O.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_H2O.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

### The Jupyter notebook [SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_LiH.ipynb) finds an approximation to the ground state of the  $LiH$ molecule
The `SQD class` is imported from [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) and instantiated as follows:
```
from SQD_Alain import SQD

SQD_LiH = SQD(
            #-------------
            # Run options
            #-------------
            backend_name = None,                                      # IBM cloud backend name
            do_plot_gate_map = True,                                  # Whether to plot the gate map 
            load_bit_array_file = None,                               # If provided, function step_3 will load samples from this file
            save_bit_array_file = None,                               # If provided, function step_3 will save samples into this file
            n_ancillary_qubits = 0,                                   # Number of ancillary qubits
            run_on_QPU = False,                                       # Whether to run the quantum circuit on the target hardware
            nshots = 1000,                                            # Number of shots
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis = "6-31g",                                          # Basis set
            atom = [["Li", (0.0, 0.0, 0.0)], ["H", (1.0, 0.0, 1.5474)]],  # Atom configuration, PySCF Initializing a molecule
            symmetry = True,                                          # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry
            n_frozen = 0,                                             # Number of frozen orbitals
            compute_exact_energy = True,                              # Whether to compute the exact energy with PySCF, cas.run().e_tot
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            chem_accuracy = 1e-3,                                     # Chemical accuracy (+/- 1 milli-Hartree)
            energy_tol = 3e-5,                                        # Numerical tolerance for convergence of the energy
            occupancies_tol = 1e-3,                                   # Numerical tolerance for convergence of the average orbital occupancies
            max_iterations = 2,                                       # Limit on the number of configuration recovery iterations
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches = 5,                                          # The number of batches to subsample in each configuration recovery iteration
            samples_per_batch = 300,                                  # The number of bitstrings to include in each subsampled batch of bitstrings
            symmetrize_spin = True,                                   # Whether to always merge spin-alpha and spin-beta CI strings into a single list
            carryover_threshold = 1e-4,                               # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            max_cycle = 200,                                          # The maximum number of Davidson cycles run by the eigenstate solver
            seed = 24,                                                # A seed for the pseudorandom number generator
            spin_sq = 0.0)                                            # spin square  
```
The last cell displays the following two plots:

![plot_energy_and_occupancy_LiH.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_LiH.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first two orbitals with high probability.

### The Jupyter notebook [SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CH2_Triplet.ipynb) finds an approximation to the ground state of the triplet $CH_2$ molecule
The `SQD class` is imported from [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) and instantiated as follows:
```
SQD_CH2_T = SQD(
            #-------------
            # Run options
            #-------------
            backend_name = None,                                      # IBM cloud backend name
            do_plot_gate_map = True,                                  # Whether to plot the gate map 
            load_bit_array_file = None,                               # If provided, function step_3 will load samples from this file
            save_bit_array_file = None,                               # If provided, function step_3 will save samples into this file
            n_ancillary_qubits = 0,                                   # Number of ancillary qubits
            run_on_QPU = False,                                       # Whether to run the quantum circuit on the target hardware
            nshots = 1000,                                            # Number of shots
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis = "6-31g",                                          # Basis set
            atom = [["C", (0.0, 0.0, 0.0)], ["H", (0.0, 0.0, 1.1160)], ["H", (0.9324, 0.0, -0.3987)]], # Atom configuration, PySCF Initializing a molecule
            spin = 2,                                                 # Triplet state (2S = 2 → S = 1)
            symmetry = True,                                          # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry
            n_frozen = 0,                                             # Number of frozen orbitals
            compute_exact_energy = True,                              # Whether to compute the exact energy with PySCF, cas.run().e_tot
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            chem_accuracy = 1e-3,                                     # Chemical accuracy (+/- 1 milli-Hartree)
            energy_tol = 3e-5,                                        # Numerical tolerance for convergence of the energy
            occupancies_tol = 1e-3,                                   # Numerical tolerance for convergence of the average orbital occupancies
            max_iterations = 10,                                      # Limit on the number of configuration recovery iterations
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches = 5,                                          # The number of batches to subsample in each configuration recovery iteration
            samples_per_batch = 300,                                  # The number of bitstrings to include in each subsampled batch of bitstrings
            symmetrize_spin = False,                                  # Whether to always merge spin-alpha and spin-beta CI strings into a single list
            carryover_threshold = 1e-4,                               # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            max_cycle = 200,                                          # The maximum number of Davidson cycles run by the eigenstate solver
            seed = 24,                                                # A seed for the pseudorandom number generator
            spin_sq = 0.0)                                            # spin square
```
* **Bond length C–H** ≈ 1.116 Å (first H) and ≈ 1.01 Å (second H), depending on the method.
* **H–C–H angle** ≈ 135°, reflecting the open-shell triplet nature.
* The geometry is asymmetric, with **one hydrogen atom out of the molecular plane**—this reflects the true triplet ground state ($^3B_1$) of $CH_2$.

The last cell displays the following two plots:

![plot_energy_and_occupancy_CH2_Triplet.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_CH2_Triplet.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

### The Jupyter notebook [SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CH2_Singlet.ipynb) finds an approximation to the ground state of the singlet $CH_2$ molecule
The `SQD class` is imported from [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) and instantiated as follows:
```
from SQD_Alain import SQD

SQD_CH2_S = SQD(
            #-------------
            # Run options
            #-------------
            backend_name = None,                                      # IBM cloud backend name
            do_plot_gate_map = True,                                  # Whether to plot the gate map 
            load_bit_array_file = None,                               # If provided, function step_3 will load samples from this file
            save_bit_array_file = None,                               # If provided, function step_3 will save samples into this file
            n_ancillary_qubits = 0,                                   # Number of ancillary qubits
            run_on_QPU = False,                                       # Whether to run the quantum circuit on the target hardware
            nshots = 1000,                                            # Number of shots
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis = "6-31g",                                          # Basis set
            atom = [["C", (0.0, 0.0, 0.0)], ["H", (0.0, 0.0, 1.0780)], ["H", (1.0250, 0.0, -0.3630)]], # Atom configuration, PySCF Initializing a molecule
            spin = 0,                                                 # Singlet state
            symmetry = True,                                          # Point group symmetry, https://pyscf.org/user/gto.html#point-group-symmetry
            n_frozen = 0,                                             # Number of frozen orbitals
            compute_exact_energy = True,                              # Whether to compute the exact energy with PySCF, cas.run().e_tot
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            chem_accuracy = 1e-3,                                     # Chemical accuracy (+/- 1 milli-Hartree)
            energy_tol = 3e-5,                                        # Numerical tolerance for convergence of the energy
            occupancies_tol = 1e-3,                                   # Numerical tolerance for convergence of the average orbital occupancies
            max_iterations = 10,                                      # Limit on the number of configuration recovery iterations
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches = 5,                                          # The number of batches to subsample in each configuration recovery iteration
            samples_per_batch = 300,                                  # The number of bitstrings to include in each subsampled batch of bitstrings
            symmetrize_spin = True,                                   # Whether to always merge spin-alpha and spin-beta CI strings into a single list
            carryover_threshold = 1e-4,                               # Threshold for carrying over bitstrings with large CI weight from one iteration of configuration recovery to the next
            max_cycle = 200,                                          # The maximum number of Davidson cycles run by the eigenstate solver
            seed = 24,                                                # A seed for the pseudorandom number generator
            spin_sq = 0.0)                                            # spin square
```
* The above geometry gives a bond length of $\~1.09$ Å and an H–C–H angle of $\~104°$, which is reasonable for bent $CH_2$.

The last cell displays the following two plots:

![plot_energy_and_occupancy_CH2_Singlet.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_CH2_Singlet.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first four orbitals with high probability.
