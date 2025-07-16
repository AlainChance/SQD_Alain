# SQD_Alain
Find an approximation to the ground state of a chemistry Hamiltonian with Sample-based Quantum Diagonalization (SQD).

The Python file `SQD_Alain.py` and the following list of Jupyter notebooks are compatible with Python 3.13, Qiskit v2.1, Qiskit runtime version: 0.40 and Qiskit Runtime V2 primitives:
* [SQD_Alain.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.ipynb)
* [SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CH2_Singlet.ipynb)
* [SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CH2_Triplet.ipynb)
* [SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CO2.ipynb)
* [SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_H2O.ipynb)
* [SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_LiH.ipynb)
* [SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_N2.ipynb)

|||
|-|-|
|**Author:** |Alain Chancé|
|**Date:** |July 1, 2025|
|**Version:** |**1.00**<br/>*Details see at the end of this notebook*|
|**License:** |[MIT License](https://github.com/AlainChance/SQD_Alain/blob/main/LICENSE)|
|**Credit:**|
The Python file `SQD_Alain.py` and the Jupyter notebooks listed above are derived from the tutorial [Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb) which is distributed under the [Apache License 2.0](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/LICENSE.txt).
|**References:**|
[Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb)
[Sample-based Quantum Diagonalization (SQD), Quantum diagonalization algorithms](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview)
[Quantum diagonalization algorithms, IBM Quantum Platform, Quantum Learning](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms)
[Antonio Mezzacapo: Quantum diagonalization methods for lattice models, IQuS - The InQubator For Quantum Simulation, Feb 19, 2025]( https://www.youtube.com/watch?v=b1fhh71hY2g)
[Javier Robledo-Moreno, Gavin Jones, Roberto Lo Nardo, Robert Davis, Lockheed Martin & IBM combine quantum computing with classical HPC in new research, IBM Quantum Research Blog, 22 May 2025](https://www.ibm.com/quantum/blog/lockheed-martin-sqd)
[Ieva Liepuoniute, Kirstin D. Doney, Javier Robledo Moreno, Joshua A. Job, William S. Friend, Gavin O. Jones, Quantum-Centric Computational Study of Methylene Singlet and Triplet States, J. Chem. Theory Comput. 2025, 21, 10, 5062–5070](https://doi.org/10.1021/acs.jctc.5c00075)
[Danil Kaliakin, Akhil Shajan, Fangchun Liang, Kenneth M. Merz Jr., Implicit Solvent Sample-Based Quantum Diagonalization, J. Phys. Chem. B 2025, 129, 23, 5788–5796](https://doi.org/10.1021/acs.jpcb.5c01030)
[Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb](https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb)
[J. Robledo-Moreno et al., "Chemistry Beyond Exact Solutions on a Quantum-Centric Supercomputer" (2024), arXiv:quant-ph/2405.05068](https://arxiv.org/abs/2405.05068)
[M. Motta et al., “Bridging physical intuition and hardware efficiency for correlated electronic states: the local unitary cluster Jastrow ansatz for electronic structure” (2023), Chem. Sci., 2023, 14, 11213](https://pubs.rsc.org/en/content/articlehtml/2023/sc/d3sc02516k)
[Quantum-Selected Configuration Interaction: classical diagonalization of Hamiltonians in subspaces selected by quantum computers, arXiv:2302.11320, quant-ph](https://doi.org/10.48550/arXiv.2302.11320)
[Introduction to Qiskit patterns](https://quantum.cloud.ibm.com/docs/en/guides/intro-to-patterns)
[Keeper L. Sharkey and Alain Chancé, Quantum Chemistry and Computing for the Curious: Illustrated with Python and Qiskit® code, Packt 2022](https://a.co/d/6YCwgPb)
[Companion Jupyter notebook for Chapter 5 of the book Quantum Chemistry and Computing for the Curious: Illustrated with Python and Qiskit® code]( https://github.com/AlainChance/Quantum-Chemistry-and-Computing-for-the-Curious/blob/main/Chapter_05_Variational_Quantum_Eigensolver_(VQE)_algorithm_V4.ipynb)
[Taha Selim, Ad van der Avoird, Gerrit C. Groenenboom, State-to-state rovibrational transition rates for CO2 in the bend mode in collisions with He atoms, J. Chem. Phys. 159, 164310 (2023)](https://doi.org/10.1063/5.0174787)
<br/>

# SQD process using Qiskit patterns
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

## Self-consistent configuration recovery procedure
The probabilistic self-consistent configuration recovery is an iterative procedure that runs in a loop which comprises the following steps:
- Post-process the full set of bitstring samples, using prior knowledge of particle number and the average orbital occupancy calculated on the most recent iteration.
- Probabilistically create batches of subsamples from recovered bitstrings.
- Project and diagonalize the molecular Hamiltonian over each sampled subspace.
- Save the minimum ground state energy found across all batches and update the average orbital occupancy.

## SQD workflow
The SQD workflow with self-consistent configuration recovery is depicted in the following diagram.

![SQD diagram](https://raw.githubusercontent.com/Qiskit/qiskit-addon-sqd/7fcec2a686bfe115560db20b840cf2b185ae06f6/docs/_static/images/sqd_diagram.png)

Source: [Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb)

## Step 1: Map classical inputs to a quantum problem
We perform a CCSD calculation. The [t1 and t2 amplitudes](https://en.wikipedia.org/wiki/Coupled_cluster#Cluster_operator) from this calculation will be used to initialize the parameters of the `LUCJ` ansatz circuit.

We use [ffsim](https://github.com/qiskit-community/ffsim) to create the ansatz circuit. ffsim is a software library for simulating fermionic quantum circuits that conserve particle number and the Z component of spin. This category includes many quantum circuits used for quantum chemistry simulations. By exploiting the symmetries and using specialized algorithms, ffsim can simulate these circuits much faster than a generic quantum circuit simulator.

Depending on the number of unpaired electrons, we use either:
* the spin-balanced variant of the unitary cluster Jastrow (UCJ) ansatz, ffsim class [UCJOpSpinBalanced](https://qiskit-community.github.io/ffsim/api/ffsim.html#ffsim.UCJOpSpinBalanced).
* the spin-unbalanced variant of the unitary cluster Jastrow (UCJ) ansatz, ffsim class [UCJOpSpinUnbalanced](https://qiskit-community.github.io/ffsim/api/ffsim.html#ffsim.UCJOpSpinUnbalanced).

## Step 2: Optimize problem for quantum execution
Next, we optimize the circuit for a target hardware. We generate a staged pass manager using the [generate\_preset\_pass\_manager](https://docs.quantum.ibm.com/api/qiskit/transpiler_preset#generate_preset_pass_manager) function from qiskit with the specified `backend` and `initial_layout`.

We set the `pre_init` stage of the staged pass manager to `ffsim.qiskit.PRE_INIT`. It includes qiskit transpiler passes that decompose gates into orbital rotations and then merges the orbital rotations, resulting in fewer gates in the final circuit.

## Step 3: Execute using Qiskit Primitives
As shown in the tutorial [Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb), we can use `Fake_Sherbrooke`, a fake 127-qubit backend from `qiskit_ibm_runtime` to emulate a real device, and generate random samples drawn from the uniform distribution. This approach is however not recommended, see [A case against uniform sampling](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).

## Step 4: Post-process and return result to desired classical format
The first iteration of the self-consistent configuration recovery procedure uses the raw samples, after post-selection on symmetries, as input to the diagonalization process to obtain an estimate of the average orbital occupancies.

Subsequent iterations use these occupancies to generate new configurations from raw samples that violate the symmetries (i.e., are incorrect). These configurations are collected and then subsampled to produce batches, which are subsequently used to project the Hamiltonian and compute a ground-state estimate using an eigenstate solver.

We use the `diagonalize_fermionic_hamiltonian` function defined in [qiskit-addon-sqd/qiskit\_addon\_sqd/fermion.py](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py).

The solver included in the SQD addon uses PySCF's implementation of selected CI, specifically [pyscf.fci.selected_ci.kernel_fixed_space](https://pyscf.org/pyscf_api_docs/pyscf.fci.html#pyscf.fci.selected_ci.kernel_fixed_space)

The following parameters are user-tunable:
* `max_iterations`: Limit on the number of configuration recovery iterations.
* `num_batches`: The number of batches of configurations to subsample (i.e., the number of separate calls to the eigenstate solver).
* `samples_per_batch`: The number of unique configurations to include in each batch.
* `max_cycles`: The maximum number of Davidson cycles run by the eigenstate solver.
* `occupancies_tol`: Numerical tolerance for convergence of the average orbital occupancies. If the maximum change in absolute value of the average occupancy of an orbital between iterations is smaller than this value, then the configuration recovery loop will exit, if the energy has also converged (see the ``energy_tol`` argument).
* `energy_tol`: Numerical tolerance for convergence of the energy. If the change in energy between iterations is smaller than this value, then the configuration recovery loop will exit, if the occupancies have also converged (see the ``occupancies_tol`` argument).

## Plot energy and occupancy
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. The usual chemical accuracy is typically defined as `1 kcal/mol` $\approx$ `1.6 mHa`. However in the context of quantum chemistry, chemical accuracy is often approximated as `1 mHa`. The estimation of the energy can be improved by drawing more samples from the circuit, increasing the number of samples per batch and increasing the number of iterations.

The second plot shows the average occupancy of each spatial orbital after the final iteration. 

## The Jupyter notebook [SQD_Alain.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.ipynb) is a basis for your own experiments.

## The Jupyter notebook [SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CH2_Singlet.ipynb) finds an approximation to the carbene singlet state $CH_2$ molecule in the `6-31G` basis set.
The first excited state $CH_2$ is called a carbene singlet state, where the two electrons in the carbon atom are paired and one orbital is empty. In a singlet state, the two electrons are paired with opposite spins, resulting in a total spin of zero. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

We import the `SQD class` from [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) and we create an instance with singlet $CH_2$ molecule definition, `6-31G` basis set, `max_iterations = 10`, `num_batches = 5`, `samples_per_batch = 300`, `energy_tol = 3e-5`.

The PySCF atom configuration gives a bond length of ≈ 1.09 Å and an $H–C–H$ angle of ≈ 104°, which is reasonable for bent $CH_2$.

The last cell displays the following two plots:

![plot_energy_and_occupancy_CH2_Singlet.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_CH2_Singlet.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first four orbitals with high probability.

## The Jupyter notebook [SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CH2_Triplet.ipynb) finds an approximation to the ground state of the triplet $CH_2$ molecule.
In its ground state, $CH_2$ is a diradical that adopts a triplet electronic structure configuration, where the carbon atom's outer shell contains two unpaired electrons which have parallel spins, leading to a total spin of one. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

We create an instance of the SQD class with triplet $CH_2$ molecule definition, `6-31G` basis set, `max_iterations = 10`, `num_batches = 5`, `samples_per_batch = 300`, `energy_tol = 3e-5`. 

In the PySCF list of options, we set `spin = 2`, the number of unpaired electrons of a molecule. Spin symmetrization will not be performed because the number of alpha and beta electrons is unequal. The [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) class automatically sets `symmetrize_spin = False` when the PySCF `spin` option is not zero.

The PySCF atom configuration gives the following:
* Bond length $C–H$ ≈ 1.116 Å (first H) and ≈ 1.01 Å (second H), depending on the method.
* $H–C–H$ angle ≈ 135°, reflecting the open-shell triplet nature.
* The geometry is asymmetric, with one hydrogen atom out of the molecular plane—this reflects the true triplet ground state ($^3B_1$) of $CH_2$.

This geometry is based on high-quality quantum chemistry references and is widely used in benchmarking.

The last cell displays the following two plots:

![plot_energy_and_occupancy_CH2_Triplet.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_CH2_Triplet.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. It is consistent with a triplet electronic structure configuration, where the carbon atom's outer shell contains two unpaired electrons which have parallel spins. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

## The Jupyter notebook [SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_CO2.ipynb) finds an approximation to the ground state of the singlet $CO_2$ molecule.
Carbon dioxide is a linear molecule with the carbon atom in the center and two oxygen atoms symmetrically placed along the z-axis.

* The bond length used here is **1.160 Å**, which is close to the experimental value (≈1.16 Å for $C=O$).
* The molecule is placed along the **z-axis** to maintain linearity and benefit from symmetry.
* `charge=0` and `spin=0` correspond to a neutral, singlet CO₂ molecule.

We create an instance of the SQD class with singlet $CO_2$ molecule definition, `sto3g` basis set, `max_iterations = 10`, `num_batches = 5`, `samples_per_batch = 300`, `energy_tol = 3e-5`.

The last cell displays the following two plots:

![plot_energy_and_occupancy_N2](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_CO2.png)

The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first eleven orbitals with high probability.

## The Jupyter notebook [SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_H2O.ipynb) finds an approximation to the ground state of the $H_2O$ molecule.
For $H_2O$ in 6-31G, the FCI ground state energy is approximately `–76.20` Hartree. The experimental non-relativistic Born–Oppenheimer energy is about `–76.438` Hartree.

The last cell displays the following two plots:

![plot_energy_and_occupancy_H2O.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_H2O.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

## The Jupyter notebook [SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_LiH.ipynb) finds an approximation to the ground state of the  $LiH$ molecule.
The last cell displays the following two plots:

![plot_energy_and_occupancy_LiH.png](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_LiH.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first two orbitals with high probability.

## The Jupyter notebook [SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_N2.ipynb) finds an approximation to the ground state of the $N_2$ molecule.
The last cell displays the following two plots:

![plot_energy_and_occupancy_N2](https://github.com/AlainChance/SQD_Alain/blob/main/plot_energy_and_occupancy_N2.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.
