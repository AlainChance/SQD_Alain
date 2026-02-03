# SQD_Alain
Demonstrating chemistry simulations with Sample-based Quantum Diagonalization (SQD).

|||
|-|-|
|**Author:** |Alain Chancé|
|**Date:** |July 1, 2025|
|**Version:** |**1.00**|
|**License:** |[MIT License](https://github.com/AlainChance/SQD_Alain/blob/main/LICENSE)|
<br/>

The Python files [SQD_Alain.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py), [SQD_Alain_Gradio.py](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain_Gradio.py) and the Jupyter notebooks in this repository `SQD_Alain` are compatible with Python 3.13, Qiskit v2.1, Qiskit runtime version: 0.40 and Qiskit Runtime V2 primitives.

# Credits
## Tutorial: "Improving energy estimation of a chemistry Hamiltonian with SQD"
The Jupyter notebooks in this repository `SQD_Alain` are derived from the tutorial [Improving energy estimation of a chemistry Hamiltonian with SQD](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb), which is distributed under the [Apache License 2.0](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/LICENSE.txt).

## Automated selection of "zig-zag" layout
The module `zigzag.py` is derived from the tutorial [Sample-based quantum diagonalization of a chemistry Hamiltonian](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization), specifically Step 2: *Optimize the problem for quantum hardware execution*.

```python
def create_lucj_zigzag_layout()
def get_zigzag_physical_layout()
initial_layout, _ = get_zigzag_physical_layout(num_orbitals, backend=backend)
```

# Additions by Alain Chancé
## SQD Gradio User Interface 
The module `SQD_Alain_Gradio.py` provides a simple web interface for configuring SQD simulations that leverages Gradio.

*“Gradio is an open-source Python package that allows you to quickly build a demo or web application for your machine learning model, API, or any arbitrary Python function. You can then share a link to your demo or web application in just a few seconds using Gradio's built-in sharing features. No JavaScript, CSS, or web hosting experience needed!”* [Quickstart](https://www.gradio.app/guides/quickstart).    

## Statistics about power consumption
The SQD class integrates the [eco2AI](https://github.com/sb-ai-lab/Eco2AI) tracking feature, a python library which accumulates statistics about power consumption and CO2 emission during running code. The Eco2AI is licensed under a [Apache licence 2.0](https://www.apache.org/licenses/LICENSE-2.0).

# Installation
## Requirements
Be sure you have the following installed:

* Qiskit SDK v2.1 or later, with visualization support (`pip install 'qiskit[visualization]'`)
* 'qiskit-aer' library (`pip install qiskit-aer`)
* Qiskit runtime 0.40 or later (`pip install qiskit-ibm-runtime`)
* Qiskit addon: sample-based quantum diagonalization (SQD) (`pip install qiskit_addon_sqd`)
* PySCF (`pip install pyscf`)
* ffsim (`pip install ffsim`)
* The 'Graphviz' library is required to use 'plot_coupling_map' (`sudo apt install graphviz`)
* Gradio (`pip install --upgrade gradio`)
* [eco2AI](https://github.com/sb-ai-lab/Eco2AI) is optional (`pip install eco2ai`)

## Clone the repository `SQD_Alain`
`git clone https://github.com/AlainChance/SQD_Alain`

# Setup your own SQD simulation
There are two options:
* Classical simulation by generating random samples drawn from a uniform distribution. *“While this approach may work for small problems, it tends to fail for larger and more practical problems.”* See [A Case Against Uniform Sampling](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).

* Simulation on a real QPU. The number of shots is highly system-and hardware-dependent. Check [Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb](https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb) for suggestions.

## Manual configuration
Use [SQD_manual.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_manual.ipynb) as template.
Keep [SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_H2O.ipynb) as reference.

### Import required dependencies
```python
from SQD_Alain import SQD
from SQD_simulation import run_simulation, get_QPU_power

# Import QiskitRuntimeService
from qiskit_ibm_runtime import QiskitRuntimeService
```

### Setup the configuration dictionary
See an example configuration here: [CH2_Triplet_ibm_fez_manual/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb)

### Classical simulation
* `run_on_QPU`: `False` to perform a classical simulation by loading a bit array file or by generating random samples drawn from the uniform distribution.
* `nsamples`: number of random samples (x1,000) to be drawn from a uniform distribution.
* `load_bit_array_file`: `None` or a real file name.
  * `None` → generate random samples drawn from the uniform distribution.
  * `bit array file name`, e.g. `N2_ibm_brisbane_bitarray.npy` → load a bit array file from a previous run.

### Simulation on a real QPU
* `run_on_QPU`: `True` to perform a simulation on a real QPU.
* `backend_name`: `None` or a real backend name.
* `nshots`: Number of shots.

### Run simulation
```python
result = run_simulation(config=config, do_print=True)
```

## Configuration with the SQD Gradio user interface
Use [SQD_Alain.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.ipynb) as template.

### Loading a configuration from a Json file
To load a configuration enter its name in the field `Json configuration filename` and then click on the button **Load configuration from a Json file**. It will update values for all components of the Gradio user interface.

The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none.

### Classical simulation
* Checkbox `run_on_QPU`: leave unticked to perform a classical simulation by loading a bit array file or by generating random samples drawn from the uniform distribution.
* Slider `Number of samples (x1,000)`: select the number of random samples (x1,000) to be drawn from a uniform distribution.
* Textbox `Load bit array file name (or 'None')`:
  * `None` → generate random samples drawn from the uniform distribution.
  * `bit array file name`, e.g. `N2_ibm_brisbane_bitarray.npy` → load a bit array file from a previous run.

See for exemple: ![CH2_Singlet_cc-pvdz/SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz/SQD_config_and_run.png)

### Simulation on a real QPU
* Checkbox `run_on_QPU`: ticked to perform a simulation on a real QPU.
* Dropdown listbox `Select IBM Backend Name (or 'None' for least busy)`: select `None` or a real backend name.
* Slider `Number of shots`: select number of shots.

See for exemple: ![CH2_Singlet_cc-pvdz_ibm_fez/SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz_ibm_fez/SQD_config_and_run.png)

### Run simulation
Click on the button **Save SQD configuration into a Json file and run simulation**.

The configuration is automatically saved into a Json configuration file, see for exemple: [CH2_Triplet_6-31g_ibm_fez
/SQD_CH2_Triplet.json](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.json).

### Status window
The status window displays error messages and the following information at the end of a successful simulation:
* SQD configuration saved into (Json configuration file name) and simulation complete.
* Exact energy (Ha)
* SQD energy (Ha)
* Absolute error (Ha)
* QPU usage (s):  (h): 
* Rough estimate for QPU power consumption (kWh)
* Classical processing - Duration (h) - Power consumption (kWh)

See for example: ![CH2_Triplet_6-31g_ibm_fez/CH2_Triplet_ibm_fez - Gradio status window.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/CH2_Triplet_ibm_fez%20-%20Gradio%20status%20window.png)

### Plot energy and occupancy window
* The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. The usual chemical accuracy is typically defined as `1 kcal/mol` $\approx$ `1.6 mHa`. However in the context of quantum chemistry, chemical accuracy is often approximated as `1 mHa`. The estimation of the energy can be improved by drawing more samples from the circuit, increasing the number of samples per batch and increasing the number of iterations.

* The second plot shows the average occupancy of each spatial orbital after the final iteration.

# Energetics Analysis
* **Assumption:**
  A ballpark estimate for a typical modern IBM-class superconducting quantum computer (including cryogenics and supporting infrastructure, while idle or lightly used) is approximately **15–25 kW**.
  Source: [*Green Quantum Computing*, Capgemini, 8 May 2023](https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/).

* The SQD class integrates the [eco2AI](https://github.com/sb-ai-lab/Eco2AI) tracking feature, a python library which accumulates statistics about power consumption and CO2 emission during running code. The Eco2AI is licensed under a [Apache licence 2.0](https://www.apache.org/licenses/LICENSE-2.0).

## Manual configuration

### Provide a ballpark figure for the QPU power consumption
Specify a value between `15` and `25` kW:

```python
power_QPU = get_QPU_power(power_QPU=25.0)  # Default value: 25 kW
```

### Define the eco2AI project in the configuration dictionary
See an example configuration here: [CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb)

```python
config = {
        ...
        #-------------------------------------
        # eco2AI Tracker options
        # https://github.com/sb-ai-lab/Eco2AI
        #-------------------------------------
        "do_eco2ai": True,                                      # Whether to track energy usage with eco2AI
        "project_name": "SQD_Alain",                            # Project name
        "experiment_description": "SQD_CH2_Triplet_ibm_fez",    # Experiment description
        "eco2ai_file_name": "SQD_CH2_Triplet_ibm_fez.csv",      # eco2AI file name
        #---------------------------------------------------------------------------------
        # Ballpark figure (kW) for the power consumption of the IBM cloud backend
        # "The power consumption of a quantum computer is about 15-25kW"
        # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
        #---------------------------------------------------------------------------------
        "power_QPU": power_QPU,                                   # Ballpark figure (kW)
        ...
```

## Configuration with the SQD Gradio user interface
### Provide a ballpark figure for the QPU power consumption
* Slider `QPU power consumption (kW)`: Select a value between `15` and `25` kW.

### Enter the eco2AI parameters
* Checkbox `eco2AI`: Tick to enable eco2AI tracking. The following text boxes will appear:

  * Project_name: Project name, e.g., `SQD_Alain`
  * Experiment_description: Experiment description, e.g., `SQD_CH4_Triplet`
  * eco2ai_file_name: eco2AI file name, e.g., `SQD_CH4_Triplet.ibm_fez.csv`

See for exemple: [CH2_Triplet_6-31g_ibm_fez/SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/SQD_config_and_run.png)

## Result
A rough estimate for the QPU power consumption is presented as follows::
```
A rough estimate for QPU power consumption is computed as power QPU (kW) * QPU usage (h)
power QPU (kW): 25
QPU usage (s): 2.00, (h): 0.0006
Rough estimate for QPU power consumption: (kWh) = 0.0139
```
The classical power consumption is shown as follows:
```
Classical processing - Duration (h): 0.0279 - Power consumption (kWh): 0.0014
```

# What is SQD?
“Sample-based quantum diagonalization (SQD) is a chemistry simulation technique that uses the quantum computer to extract a noisy distribution of possible electronic configurations of a molecule in the form of bitstrings. 

Then, it runs an iterative correction entirely on a classical computer to cluster these bitstrings around the correct electron configuration. 

This technique, for a large enough noisy quantum computer with low enough error rates, has the potential to perform more efficiently than classical approximation methods.”

Source: [Demonstrating a true realization of quantum-centric supercomputing](https://www.ibm.com/quantum/blog/supercomputing-24)

## SQD process using Qiskit patterns
Using the [Qiskit patterns](https://quantum.cloud.ibm.com/docs/en/guides/intro-to-patterns) framework, the SQD process can be described in four steps, Ref. [arXiv:2405.05068
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

Here $\tilde{\mathcal{X}}$ is a set of noisy samples obtained from execution on a real QPU which contain configurations represented as bitstrings.

The wave function is supported in a set of basis states $\mathcal{S} = \{|x\rangle \}$ whose size does not increase exponentially with the size of the problem. The quantum computing device produces samples of the members of $\mathcal{S}$ only. A quantum circuit prepares the state $|\Psi\rangle$ in the quantum device. The Jordan-Wigner encoding is used, hence members of the computational basis represent Fock states (electronic configurations/determinants). The quantum circuit is sampled in the computational basis, yielding the set of noisy configurations $\tilde{\mathcal{X}}$. The configurations are represented by bitstrings. The set $\tilde{\mathcal{X}}$ is then passed into the classical post-processing block, where the [self-consistent configuration recovery technique](https://arxiv.org/abs/2405.05068) is used.

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
There are two options:
* Classical simulation by generating random samples drawn from a uniform distribution. *“While this approach may work for small problems, it tends to fail for larger and more practical problems.”* See [A Case Against Uniform Sampling](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).
* Simulation on a real QPU. The number of shots is highly system-and hardware-dependent. Check [Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb](https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb) for suggestions.

## Step 4: Post-process and return result to desired classical format
The first iteration of the self-consistent configuration recovery procedure uses the raw samples, after post-selection on symmetries, as input to the diagonalization process to obtain an estimate of the average orbital occupancies.

Subsequent iterations use these occupancies to generate new configurations from raw samples that violate the symmetries (i.e., are incorrect). These configurations are collected and then subsampled to produce batches, which are subsequently used to project the Hamiltonian and compute a ground-state estimate using an eigenstate solver.

We use the `diagonalize_fermionic_hamiltonian` function defined in [qiskit-addon-sqd/qiskit\_addon\_sqd/fermion.py](https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py).

The solver included in the `SQD addon` uses PySCF's implementation of selected CI, specifically [pyscf.fci.selected_ci.kernel_fixed_space](https://pyscf.org/pyscf_api_docs/pyscf.fci.html#pyscf.fci.selected_ci.kernel_fixed_space)

The following parameters of the `SQD class` are user-tunable:
* `max_iterations`: Limit on the number of configuration recovery iterations.
* `num_batches`: The number of batches of configurations to subsample (i.e., the number of separate calls to the eigenstate solver).
* `samples_per_batch`: The number of unique configurations to include in each batch.
* `max_cycles`: The maximum number of Davidson cycles run by the eigenstate solver.
* `occupancies_tol`: Numerical tolerance for convergence of the average orbital occupancies. If the maximum change in absolute value of the average occupancy of an orbital between iterations is smaller than this value, then the configuration recovery loop will exit, if the energy has also converged (see the ``energy_tol`` argument).
* `energy_tol`: Numerical tolerance for convergence of the energy. If the change in energy between iterations is smaller than this value, then the configuration recovery loop will exit, if the occupancies have also converged (see the ``occupancies_tol`` argument).

# Ethylene $C_2H_4$ molecule
* [C2H4_sto-6g/SQD_C2H4.ipynb/SQD_C2H4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/C2H4/C2H4_sto-6g/SQD_C2H4.ipynb)
* [C2H4_sto-6g_ibm_fez/SQD_C2H4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/C2H4/C2H4_sto-6g_ibm_fez/SQD_C2H4.ipynb)
* [C2H4_sto-6g_ibm_fez_manual/SQD_C2H4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/C2H4/C2H4_sto-6g_ibm_fez_manual/SQD_C2H4.ipynb)

## [C2H4_sto-6g_ibm_fez/SQD_C2H4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/C2H4/C2H4_sto-6g_ibm_fez/SQD_C2H4.ipynb) finds an approximation to the ethylene $C_2H_4$ molecule in the `sto-6g` basis set
Ethylene is a hydrocarbon that is widely used in the chemical industry. Source: [How to simulate the local unitary cluster Jastrow (LUCJ) ansatz](https://qiskit-community.github.io/ffsim/dev/how-to-guides/simulate-lucj.html).

The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_C2H4.json` and then click on the button **Load configuration from a Json file**.

* `Basis set`: `sto-6g`
* `Job id` : `d4s2sj3her1c73bcs1pg`
* `Run on QPU`: ticked
* `Backend name`: `ibm_fez`
* `Number of shots`: `50000`
* `Atom configuration`: `[['C', (0, 0, 0.6695)], ['C', (0, 0, -0.6695)], ['H', (0, 0.9289, 1.2321)], ['H', (0, -0.9289, 1.2321)], ['H', (0, 0.9289, -1.2321)], ['H', (0, -0.9289, -1.2321)]]`
* `spin`: `0`
* `Symmetry`: `d2h`
* `Number of frozen orbitals`:  `0`
* `Max recovery iterations`: `15`
* `Batches`: `5`
* `Samples per batch`: `300`
* `Compute exact energy`: ticked
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_C2H4`
* `File name`: `SQD_C2H4.csv`
* `Json configuration file name`: `SQD_C2H4.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/C2H4/C2H4_sto-6g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/C2H4/C2H4_sto-6g_ibm_fez/plot_energy_and_occupancy.png)

The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. The usual chemical accuracy is typically defined as `1 kcal/mol` $\approx$ `1.6 mHa`. However in the context of quantum chemistry, chemical accuracy is often approximated as `1 mHa`. The estimation of the energy can be improved by drawing more samples from the circuit, increasing the number of samples per batch and increasing the number of iterations.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first eight orbitals with high probability.

# Carbene singlet state $CH_2$ molecule
* [CH2_Singlet_6-31g/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g/SQD_CH2_Singlet.ipynb)
* [CH2_Singlet_6-31g_ibm_fez/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez/SQD_CH2_Singlet.ipynb)
* [CH2_Singlet_6-31g_ibm_fez_manual/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez_manual/SQD_CH2_Singlet.ipynb)
* [CH2_Singlet_cc-pvdz/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz/SQD_CH2_Singlet.ipynb)
* [CH2_Singlet_cc-pvdz_ibm_fez/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz_ibm_fez/SQD_CH2_Singlet.ipynb)

The first excited state $CH_2$ is called a carbene singlet state, where the two electrons in the carbon atom are paired and one orbital is empty. In a singlet state, the two electrons are paired with opposite spins, resulting in a total spin of zero. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

[SQD-12] Ieva Liepuoniute, Kirstin D. Doney, Javier Robledo Moreno, Joshua A. Job, William S. Friend, Gavin O. Jones, Quantum-Centric Computational Study of Methylene Singlet and Triplet States, J. Chem. Theory Comput. 2025, 21, 10, 5062–5070, [https://doi.org/10.1021/acs.jctc.5c00075](https://doi.org/10.1021/acs.jctc.5c00075).

The PySCF atom configuration gives a bond length of ≈ 1.09 Å and an $H–C–H$ angle of ≈ 104°, which is reasonable for bent $CH_2$.

## [CH2_Singlet_6-31g_ibm_fez/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez/SQD_CH2_Singlet.ipynb) finds an approximation to the carbene singlet state $CH_2$ molecule in the `6-31G` basis set.
The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_CH2_Singlet.json` and then click on the button **Load configuration from a Json file**.

* `Json configuration file name`: `SQD_CH2_Singlet.json`
* `IBM Backend Name (or 'None' for least busy)`: `ibm_fez`
* `Number of shots`: `1000`
* `Atom configuration`: `[['C', (0.0, 0.0, 0.0)], ['H', (0.0, 0.0, 1.078)], ['H', (1.025, 0.0, -0.363)]]`
* `spin`: `0`
* `basis`: `6-31G`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `0`
* `max_iterations`: `10`.

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first four orbitals with high probability.

## [CH2_Singlet_cc-pvdz_ibm_fez/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz_ibm_fez/SQD_CH2_Singlet.ipynb) finds an approximation to the carbene singlet state $CH_2$ molecule in the `cc-pVDZ` basis set.
The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_CH2_Singlet.json` and then click on the button **Load configuration from a Json file**.

* `Basis Set`: `cc-pvdz`
* `Job id` : `None`  Update according to your context
* `Save bit array file name`: `CH2_Singlet_bitarray.npy`
* `Run on QPU`: ticked
* `IBM Backend Name (or 'None' for least busy)`: `ibm_fez`
* `Number of shots`: `50000`
* `Atom configuration`: `[['C', (0.0, 0.0, 0.0)], ['H', (0.0, 0.0, 1.078)], ['H', (1.025, 0.0, -0.363)]]` 
* `Spin`: `0`
* `Symmetry`: `Dooh`
* `Number of frozen orbitals`:  `0`
* `compute_exact_energy`: ticked,
* `Max recovery iterations`: `20`
* `num_batches`: `5`,
* `Samples per batch`: `1000`
* `QPU power consumption (kW)`: `25`
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_CH2_Singlet`
* `File name`: `SQD_CH2_Singlet.csv`
* `Json configuration file name`: `SQD_CH2_Singlet.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first four orbitals with high probability.

## SQD_CH2_Singlet - Analysis across five experiments
* **Assumption:**
  A ballpark estimate for a typical modern IBM-class superconducting quantum computer (including cryogenics and supporting infrastructure, while idle or lightly used) is approximately **15–25 kW**.
  Source: [*Green Quantum Computing*, Capgemini, 8 May 2023](https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/).
* Classical simulation by generating random samples drawn from a uniform distribution. *“While this approach may work for small problems, it tends to fail for larger and more practical problems.”* See [A Case Against Uniform Sampling](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).
* Simulation on a real QPU. The number of shots is highly system-and hardware-dependent. Check [Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb](https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb) for suggestions.

1. [CH2_Singlet_6-31g/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g/SQD_CH2_Singlet.ipynb)
2. [CH2_Singlet_6-31g_ibm_fez/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez/SQD_CH2_Singlet.ipynb)
3. [CH2_Singlet_6-31g_ibm_fez_manual/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_6-31g_ibm_fez_manual/SQD_CH2_Singlet.ipynb)
4. [CH2_Singlet_cc-pvdz/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz/SQD_CH2_Singlet.ipynb)
5. [CH2_Singlet_cc-pvdz_ibm_fez/SQD_CH2_Singlet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Singlet/CH2_Singlet_cc-pvdz_ibm_fez/SQD_CH2_Singlet.ipynb)

* `QPU power consumption (kW)`: `25`

$$\begin{array}{|c|c|c|c|c|c|}
\hline
\text{Metric} & \text{1} & \text{2} & \text{3} & \text{4} & \text{5}\\
\hline
\text{Basis Set} & \text{6-31G} & \text{6-31G} & \text{6-31G} & \text{cc-pVDZ} & \text{cc-pVDZ}\\
\hline
\text{Backend} & \verb|No QPU| & \verb|ibm_fez| & \verb|ibm_fez| & \verb|No QPU| & \verb|ibm_fez|\\
\hline
\text{Number of shots} & 0 & 1000 & 50000 & 0 & 50000\\
\hline
\text{Number of samples x1,000} & 10 & 0 & 0 & 5000 & 0\\
\hline
\text{Number of qubits} & 0 & 26 & 26 & 0 & 48\\
\hline
\text{Circuit depth (ISA)} & 0 & 353 & 355 & 0 & 635\\
\hline
\text{Exact energy (Ha)} & -38.94147 & -38.94147 & -38.94147 & -39.02174 & -39.02174\\
\hline
\text{SQD energy (Ha)} & -38.94144 & -38.94138 & -38.94141 & -39.02139  & -39.02149\\
\hline
\text{Chemical accuracy iterations} & 2 & 3 & 1 & 6 & 5\\
\hline
\text{Qiskit Runtime usage (s)} & 0.00 & 2.00 & 15.00 & 0.00 & 16.00\\
\hline
\text{IBM QPU power (kWh)} & 0.0000 & 0.0139 & 0.1042 & 0.0000 & 0.1111\\
\hline
\text{Classical processing duration (h)} & 0.0256 & 0.0276 & 0.0324 & 4.3108 & 1.6587\\
\hline
\text{Classical processing power (kWh)} & 0.0013 & 0.0014 & 0.0013 & 0.1223 & 0.0887\\
\hline
\end{array}$$

In the `6-31G` basis set, the first simulation consumes the least energy, while using roughly the same amount of energy for classical processing on the same laptop. The third simulation achieves chemical accuracy with the fewest iterations.

In the `cc-pVDZ` basis set, which is the most accurate, the fifth simulation achieves chemical accuracy with the fewest iterations while using roughly the same total amount of energy as the fourth simulation.

# Methylene triplet $CH_2$ molecule
* [CH2_Triplet_6-31g/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g/SQD_CH2_Triplet.ipynb)
* [CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.ipynb)
* [CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb)
* [CH2_Triplet_cc-pvdz/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_cc-pvdz/SQD_CH2_Triplet.ipynb)

In its ground state, $CH_2$ is a diradical that adopts a triplet electronic structure configuration, where the carbon atom's outer shell contains two unpaired electrons which have parallel spins, leading to a total spin of one. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

[SQD-12] Ieva Liepuoniute, Kirstin D. Doney, Javier Robledo Moreno, Joshua A. Job, William S. Friend, Gavin O. Jones, Quantum-Centric Computational Study of Methylene Singlet and Triplet States, J. Chem. Theory Comput. 2025, 21, 10, 5062–5070, [https://doi.org/10.1021/acs.jctc.5c00075](https://doi.org/10.1021/acs.jctc.5c00075).

The methylene triplet $^3B_1$ has:
- one doubly occupied σ‑type orbital  
- **two singly occupied orbitals** (the “diradical” pair)  
- the rest doubly occupied core/valence orbitals

The PySCF atom configuration gives the following:
* Bond length $C–H$ ≈ 1.116 Å (first H) and ≈ 1.01 Å (second H), depending on the method.
* $H–C–H$ angle ≈ 135°, reflecting the open-shell triplet nature.
* The geometry is asymmetric, with one hydrogen atom out of the molecular plane—this reflects the true triplet ground state ($^3B_1$) of $CH_2$.

This geometry is based on high-quality quantum chemistry references and is widely used in benchmarking.

We set `spin` to `2`, the number of unpaired electrons of a molecule. Spin symmetrization will not be performed because the number of alpha and beta electrons is unequal. The [SQD](https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain.py) class automatically sets `symmetrize_spin = False` when the PySCF `spin` option is not zero.

## [CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.ipynb) finds an approximation to the ground state of the triplet $CH_2$ molecule in the `6-31G` basis set.
The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_CH2_Triplet.json` and then click on the button **Load configuration from a Json file**.

* `Basis set`: `6-31G`
* `Job id` : `None`                                          Update according to your own context.
* `Run on QPU`: ticked
* `IBM Backend Name (or 'None' for least busy)`: `ibm_fez`
* `Number of shots`: `1000`
* `Atom configuration`: `[['C', (0.0, 0.0, 0.0)], ['H', (0.0, 0.0, 1.116)], ['H', (0.9324, 0.0, -0.3987)]]`
* `spin`: `2`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `0`
* `Max recovery iterations`: `10`
* `Batches`: `5`
* `Samples per batch`: `300`
* `Compute exact energy`: ticked
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_CH2_Triplet`
* `File name`: `SQD_CH2_Triplet.csv`
* `Json configuration file name`: `SQD_CH2_Triplet.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. It is consistent with a triplet electronic structure configuration, where the carbon atom's outer shell contains two unpaired electrons which have parallel spins. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

## [CH2_Triplet_cc-pvdz/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_cc-pvdz/SQD_CH2_Triplet.ipynb) finds an approximation to the ground state of the triplet $CH_2$ molecule in the `cc-pVDZ` basis set.
We freeze 1 doubly occupied core orbital (Carbon 1s). This is standard across SCI, CC, CAS, and most correlated methods.

The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_CH2_Triplet.json` and then click on the button **Load configuration from a Json file**.

* `Basis Set`: `cc-pvdz`
* `Save bit array file name`: `CH2_Triplet_bitarray`
* `Run on QPU`: unticked
* `Number of samples (×1,000)`: `10000`
* `Atom configuration`: `[['C', (0.0, 0.0, 0.0)], ['H', (0.0, 0.0, 1.116)], ['H', (0.9324, 0.0, -0.3987)]]`
* `spin`: `2`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `1`
* `Max recovery iterations`: `20`
* `Batches`: `5`
* `Samples per batch`: `1000`
* `Compute exact energy`: ticked
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_CH2_Triplet`
* `File name`: `SQD_CH2_Triplet.csv`
* `Json configuration file name`: `SQD_CH2_Triplet.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_cc-pvdz/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_cc-pvdz/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. It is consistent with a triplet electronic structure configuration, where the carbon atom's outer shell contains two unpaired electrons which have parallel spins. Source: [Lockheed Martin & IBM combine quantum computing with classical HPC in new research](https://www.ibm.com/quantum/blog/lockheed-martin-sqd).

## SQD_CH2_Triplet - Analysis across four experiments
* **Assumption:**
  A ballpark estimate for a typical modern IBM-class superconducting quantum computer (including cryogenics and supporting infrastructure, while idle or lightly used) is approximately **15–25 kW**.
  Source: [*Green Quantum Computing*, Capgemini, 8 May 2023](https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/).
* Classical simulation by generating random samples drawn from a uniform distribution. *“While this approach may work for small problems, it tends to fail for larger and more practical problems.”* See [A Case Against Uniform Sampling](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).
* Simulation on a real QPU. The number of shots is highly system-and hardware-dependent. Check [Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb](https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb) for suggestions.

1. [CH2_Triplet_6-31g/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g/SQD_CH2_Triplet.ipynb)
2. [CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez/SQD_CH2_Triplet.ipynb)
3. [CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_6-31g_ibm_fez_manual/SQD_CH2_Triplet.ipynb)
4. [CH2_Triplet_cc-pvdz/SQD_CH2_Triplet.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH2_Triplet/CH2_Triplet_cc-pvdz/SQD_CH2_Triplet.ipynb)

* `QPU power consumption (kW)`: `25`

$$\begin{array}{|c|c|c|c|c|}
\hline
\text{Metric} & \text{1} & \text{2} & \text{3} & \text{4}\\
\hline
\text{Basis Set} & \text{6-31G} & \text{6-31G} & \text{6-31G} & \text{cc-pVDZ}\\
\hline
\text{Backend} & \verb|No QPU| & \verb|ibm_fez| & \verb|ibm_fez| & \verb|No QPU|\\
\hline
\text{Number of shots} & 0 & 1,000 & 50,000 & 0\\
\hline
\text{Number of samples x1,000} & 10 & 0 & 0 & 10,000\\
\hline
\text{Number of qubits} & 0 & 26 & 26 & 0\\
\hline
\text{Circuit depth (ISA)} & 0 & 1515 & 1380 & 0\\
\hline
\text{Exact energy (Ha)} & -38.96994 & -38.96994 & -38.96994 & -39.03015\\
\hline
\text{SQD energy (Ha)} & -38.96979 & -38.96964 & -38.96986 & -39.02969\\
\hline
\text{Chemical accuracy iterations} & 5 & 6 & 4 & 5\\
\hline
\text{Qiskit Runtime usage (s)} & 0.00 & 2.00 & 17.00 & 0.00\\
\hline
\text{IBM QPU power (kWh)} & 0.0000 & 0.0139 & 0.1181 & 0.0000\\
\hline
\text{Classical processing duration (h)} & 0.0302 & 0.0260 & 0.0423 & 7.4716\\
\hline
\text{Classical processing power (kWh)} & 0.0015 & 0.0013 & 0.0016 & 0.1334\\
\hline
\end{array}$$

In the `6-31G` basis set, the first simulation consumes the least energy, while using roughly the same amount of energy for classical processing on the same laptop. The third simulation achieves chemical accuracy with the fewest iterations.

In the `cc-pVDZ` basis set, which is the most accurate, the only simulation consumes only energy for classical processing.

# Methane $CH_4$ molecule
* [CH4_6-31g/SQD_CH4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH4/CH4_6-31g/SQD_CH4.ipynb)
* [CH4_6-31g_ibm_fez/SQD_CH4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH4/CH4_6-31g_ibm_fez/SQD_CH4.ipynb)
* [CH4/CH4_6-31g_ibm_fez_manual](https://github.com/AlainChance/SQD_Alain/blob/main/CH4/CH4_6-31g_ibm_fez_manual/SQD_CH4.ipynb)

## [CH4_6-31g_ibm_fez/SQD_CH4.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CH4/CH4_6-31g_ibm_fez/SQD_CH4.ipynb) finds an approximation to the ground state of the $CH_4$ molecule at equilibrium in the `6-31g` basis set.
Methane is a chemical compound with the chemical formula $CH_4$, the simplest alkane, and the main constituent of natural gas.

The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_CH4.json` and then click on the button **Load configuration from a Json file**.

* `Basis set`: `6-31G`
* `Job id` : `None`                                          Update according to your own context.
* `Save bit array file name`: `CH4_ibm_fez_bitarray.npy`
* `Run on QPU`: ticked
* `Backend name`: `ibm_fez`
* `Number of shots`: `100000`
* `Atom configuration`: `[['C', (0.0, 0.0, 0.0)], ['H', (0.0, 0.0, 1.089)], ['H', (1.0267, 0.0, -0.363)], ['H', (-0.5133, -0.8892, -0.363)], ['H', (-0.5133, 0.8892, -0.363)]]`
* `spin`: `0`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `0`
* `Max recovery iterations`: `15`
* `Batches`: `5`
* `Samples per batch`: `300`
* `Compute exact energy`: ticked
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_CH4_ibm_fez`
* `File name`: `SQD_CH4_ibm_fez.csv`
* `Json configuration file name`: `SQD_CH4.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH4/CH4_6-31g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/CH4/CH4_6-31g_ibm_fez/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 
The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

# Carbone dioxide $CO_2$ molecule
* [CO2_sto-6g/SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g/SQD_CO2.ipynb)
* [CO2_sto-6g_ibm_fez/SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g_ibm_fez/SQD_CO2.ipynb)

Carbon dioxide is a linear molecule with the carbon atom in the center and two oxygen atoms symmetrically placed along the z-axis.

* The bond length used here is **1.160 Å**, which is close to the experimental value (≈ 1.16 Å for $C=O$).
* The molecule is placed along the **z-axis** to maintain linearity and benefit from symmetry.
* `charge=0` and `spin=0` correspond to a neutral, singlet $CO_2$ molecule.

## [CO2_sto-6g_ibm_fez/SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g_ibm_fez/SQD_CO2.ipynb) finds an approximation to the ground state of the $CO_2$ molecule at equilibrium in the `sto-6g` basis set.
The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_CO2.json` and then click on the button **Load configuration from a Json file**.

* `Basis set`: `sto-6g`
* `Job id` : `None`    Update according to your own context.
* `Run on QPU`: ticked
* `Backend name`: `ibm_fez`
* `Number of shots`: `75000`
* `Atom configuration`: `[['O', (0.0, 0.0, -1.16)], ['C', (0.0, 0.0, 0.0)], ['O', (0.0, 0.0, 1.16)]]`
* `spin`: `0`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `0`
* `Max recovery iterations`: `10`
* `Batches`: `5`
* `Samples per batch`: `300`
* `Compute exact energy`: ticked
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_CO2`
* `File name`: `SQD_CO2_ibm_fez.csv`
* `Json configuration file name`: `SQD_CO2_ibm_fez.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g_ibm_fez/plot_energy_and_occupancy.png)

The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first eleven orbitals with high probability.

## CO2 - Analysis across two experiments
* **Assumption:**
  A ballpark estimate for a typical modern IBM-class superconducting quantum computer (including cryogenics and supporting infrastructure, while idle or lightly used) is approximately **15–25 kW**.
  Source: [*Green Quantum Computing*, Capgemini, 8 May 2023](https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/).
* Classical simulation by generating random samples drawn from a uniform distribution. *“While this approach may work for small problems, it tends to fail for larger and more practical problems.”* See [A Case Against Uniform Sampling](https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview#32-a-case-against-uniform-sampling).
* Simulation on a real QPU. The number of shots is highly system-and hardware-dependent. Check [Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb](https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb) for suggestions.

1. [CO2_sto-6g/SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g/SQD_CO2.ipynb)  

2. [CO2_sto-6g_ibm_fez/SQD_CO2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/CO2/CO2_sto-6g_ibm_fez/SQD_CO2.ipynb)

* `QPU power consumption (kW)`: `25`

$$\begin{array}{|c|c|c|}
\hline
\text{Metric} & \text{1} & \text{2}\\
\hline
\text{Basis Set} & \text{sto-6g} & \text{sto-6g}\\
\hline
\text{Backend} & \texttt{No QPU} & \verb|ibm_fez|\\
\hline
\text{Number of shots} & 0 & 75000\\
\hline
\text{Number of samples x1,000} & 10 & 0\\
\hline
\text{Number of qubits} & 0 & 26\\
\hline
\text{Circuit depth (ISA)} & 0 & 371\\
\hline
\text{Exact energy (Ha)} & -187.0652 & -187.0652\\
\hline
\text{SQD energy (Ha)} & -187.0650 & -187.0650\\
\hline
\hline
\text{Chemical accuracy iterations} & 4 & 4\\
\hline
\text{Qiskit Runtime usage (s)} & 0.00 & 21.00\\
\hline
\text{IBM QPU power (kWh)} & 0.00 & 0.1458\\
\hline
\text{Classical processing duration (h)} & 0.2585 & 0.2825\\
\hline
\text{Classical processing power (kWh)} & 0.0141 & 0.0142\\
\hline
\end{array}$$

The first simulation consumes the least energy, while using roughly the same amount of energy for classical processing on the same laptop. The second simulation achieves chemical accuracy with the same number of iterations.

# Water $H_2O$ molecule
* [H2O_6-31g/SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/H2O/H2O_6-31g/SQD_H2O.ipynb)
* [H2O_6-31g_ibm_fez/SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/H2O/H2O_6-31g_ibm_fez/SQD_H2O.ipynb)
* [H2O_6-31g_ibm_fez_manual/SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/H2O/H2O_6-31g_ibm_fez_manual/SQD_H2O.ipynb)

## [H2O_6-31g_ibm_fez/SQD_H2O.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/H2O/H2O_6-31g_ibm_fez/SQD_H2O.ipynb) finds an approximation to the ground state of the $H_2O$ molecule.
For $H_2O$ in 6-31G, the FCI ground state energy is approximately `–76.20` Hartree. The experimental non-relativistic Born–Oppenheimer energy is about `–76.438` Hartree.

The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_H2O.json` and then click on the button **Load configuration from a Json file**.

* `Basis set`: `6-31G`
* `Job id` : `None`  Update according to your context
* `Run on QPU`: ticked
* `Backend name`: `ibm_fez`
* `Number of shots`: `10000`
* `Atom configuration`: `[['O', (0.0, 0.0, 0.0)], ['H', (0.0, 1.0, 0.0)], ['H', (0.0, 0.0, 1.0)]]`
* `spin`: `0`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `0`
* `max_iterations`: `10`
* `Json configuration file name`: `SQD_H2O.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/H2O/H2O_6-31g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/H2O/H2O_6-31g_ibm_fez/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

## Lithium hydride $LiH$ molecule
* [LiH_6-31g/SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g/SQD_LiH.ipynb)
* [LiH_6-31g_ibm_fez/SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g_ibm_fez/SQD_LiH.ipynb)
* [LiH_6-31g_ibm_fez_manual/SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g_ibm_fez_manual/SQD_LiH.ipynb)
* [LiH_6-31g_ibm_torino/SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g_ibm_torino/SQD_LiH.ipynb)

## [LiH_6-31g_ibm_fez/SQD_LiH.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g_ibm_fez/SQD_LiH.ipynb) finds an approximation to the ground state of the $LiH$ molecule.

The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_LiH.json` and then click on the button **Load configuration from a Json file**.

* `Json configuration file name`: `SQD_LiH.json`
* `IBM Backend Name (or 'None' for least busy)`: `ibm_fez`
* `Number of shots`: `1000`
* `Atom configuration`: `[['Li', (0.0, 0.0, 0.0)], ['H', (1.0, 0.0, 1.5474)]]`
* `spin`: `0`
* `basis`: `6-31G`
* `Symmetry`: `True`
* `Number of frozen orbitals`:  `0`
* `max_iterations`: `10`.

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy.png](https://github.com/AlainChance/SQD_Alain/blob/main/LiH/LiH_6-31g_ibm_fez/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`. 

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first two orbitals with high probability.

# Nitrogen gas = Dinitrogen $N_2$ molecule
* [N2_6-31g/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_6-31g/SQD_N2.ipynb)
* [N2_6-31g_ibm_fez/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_6-31g_ibm_fez/SQD_N2.ipynb)
* [N2_6-31g_ibm_fez_manual/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_6-31g_ibm_fez_manual/SQD_N2.ipynb)
* [N2_cc-pvdz_ibm_fez/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez/SQD_N2.ipynb)
* [N2_cc-pvdz_ibm_fez_manual/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez_manual/SQD_N2.ipynb)

The two nitrogen atoms are connected by a triple covalent bond (N≡N), one of the strongest bonds in chemistry. This makes $N_2$ chemically inert, hard to break apart (requires high energy), and useful as a protective atmosphere in labs and industry.

The reference energy from SCI calculation performed separately is: `exact_energy` = `-109.22690201485733`. 
[Sample-based quantum diagonalization of a chemistry Hamiltonian, https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization).

See [SQD-3] Javier Robledo-Moreno et al., [Chemistry Beyond the Scale of Exact Diagonalization on a Quantum-Centric Supercomputer](https://doi.org/10.48550/arXiv.2405.05068), arXiv:2405.05068 [quant-ph], 13 Jul 2025.

## [N2_6-31g_ibm_fez/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_6-31g_ibm_fez/SQD_N2.ipynb) finds an approximation to the ground state of the nitrogen gas $N_2$ molecule in the `6-31g` basis set.
The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_N2.json` and then click on the button **Load configuration from a Json file**.

* `Basis set`: `6-31G`
* `Job Id`:`None`                                           Update according to your own context
* `Save bit array file name`: `N2_ibm_fez_bitarray.npy`
* `Run on QPU`: ticked
* `IBM Backend Name (or 'None' for least busy)`: `ibm_fez`
* `Number of shots`: `10000`
* `Atom configuration`: `[['N', (0.0, 0.0, 0.0)], ['N', (1.0, 0.0, 0.0)]]` 
* `load_bit_array_file`: `None`
* `save_bit_array_file`: `None`
* `spin`: `0`
* `basis`: `6-31G`
* `Symmetry`: `Dooh`
* `Number of frozen orbitals`:  `2`
* `max_iterations`: `20`
* `Batches`: `5`
* `Samples per batch`: `300`
* `Compute exact energy`: ticked
* `eco2AI tracker`: `ticked`
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_N2_ibm_fez`
* `File name`: `SQD_N2.csv`
* `Json configuration file name`: `SQD_N2.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_6-31g_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_6-31g_ibm_fez/plot_energy_and_occupancy.png)
The first plot shows an estimation of the ground state energy within $\approx$ `1 milli-Hartree (mHa)`.

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first five orbitals with high probability.

## [N2_cc-pvdz_ibm_fez/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez/SQD_N2.ipynb) finds an approximation to the ground state of the nitrogen gas $N_2$ molecule in the `cc-pVDZ` basis set.
The field `Json configuration filename` is initialized with the name of first file ending with `.json` in the current directory or `None` if there is none. Check that it contains the value `SQD_N2.json` and then click on the button **Load configuration from a Json file**.

* `Basis Set`: `cc-pvdz`
* `Job id` : `None`  Update according to your context
* `Save bit array file name`: `N2_cc-pvdz_bitarray.npy`
* `Run on QPU`: ticked
* `IBM Backend Name (or 'None' for least busy)`: `ibm_fez`
* `Number of shots`: `60000`
* `Atom configuration`: `[['N', (0.0, 0.0, 0.0)], ['N', (1.0, 0.0, 0.0)]]` 
* `Spin`: `0`
* `Symmetry`: `Dooh`
* `Number of frozen orbitals`:  `2`
* `compute_exact_energy`: unticked,
* `Exact energy` = `-109.2269`
* `Max recovery iterations`: `20`
* `num_batches`: `5`
* `Samples per batch`: `1000`
* `QPU power consumption (kW)`: `25`
* `eco2AI tracker`: ticked
* `Project name`: `SQD_Alain`
* `Experiment description`: `SQD_N2_cc-pvdz`
* `File name`: `SQD_N2_cc-pvdz.csv`
* `Json configuration file name`: `SQD_N2.json`

Click on the button **Save SQD configuration into a Json file and run simulation**.

![SQD_config_and_run.png](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez/SQD_config_and_run.png)

![plot_energy_and_occupancy](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez/plot_energy_and_occupancy.png)

The first plot shows an estimate of the ground-state energy within approximately **1 milli-Hartree (mHa)**. Chemical accuracy is usually defined as **`1 kcal/mol` $\approx$ `1.6 mHa`**. However, in the context of quantum chemistry, chemical accuracy is often approximated as **`1 mHa`**. The energy estimate can be further improved by drawing more samples from the circuit, increasing the number of samples per batch, and increasing the number of iterations.

These results are better than those shown in the **“Visualize the results”** section of the tutorial *Sample-based Quantum Diagonalization of a Chemistry Hamiltonian* ([link](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization#visualize-the-results)). In that tutorial, the authors report:

> “The first plot shows that after a couple of iterations we estimate the ground-state energy within approximately **`41 mHa`** (chemical accuracy is typically accepted to be **`1 kcal/mol` $\approx$ `1.6 mHa`**). The energy can be improved by allowing more iterations of configuration recovery or by increasing the number of samples per batch.”

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first eight orbitals with high probability.

## [N2_cc-pvdz_ibm_fez_manual/SQD_N2.ipynb](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez_manual/SQD_N2.ipynb) finds an approximation to the ground state of the nitrogen gas $N_2$ molecule in the `cc-pVDZ` basis set.
Define the configuration dictionary as follows:

* `job_id`: `None`                                       Update according to your own context.
* `save_bit_array_file`: `N2_ibm_fez_bitarray.npy`
* `run_on_QPU`: `True`
* `nshots`: `50000`
* `do_eco2ai`: `True`
* `project_name`: `SQD_Alain`
* `experiment_description`: `N2_cc-pvdz`
* `eco2ai_file_name`: `SQD_N2_cc-pvdz.csv`
* `power_QPU`: `power_QPU`,
* `atom`: `[['N', (0.0, 0.0, 0.0)], ['N', (1.0, 0.0, 0.0)]]` 
* `spin`: `0`
* `basis`: `cc-pvdz`
* `symmetry`: `Dooh`
* `n_frozen`:  `2`
* `compute_exact_energy`: `False`,
* `exact_energy` = `-109.22690201485733`
* `max_iterations`: `20`
* `num_batches`: `5`
* `samples_per_batch`: `1000`

![plot_energy_and_occupancy](https://github.com/AlainChance/SQD_Alain/blob/main/N2/N2_cc-pvdz_ibm_fez_manual/plot_energy_and_occupancy.png)

The first plot shows an estimate of the ground-state energy within approximately **1 milli-Hartree (mHa)**. Chemical accuracy is usually defined as **`1 kcal/mol` $\approx$ `1.6 mHa`**. However, in the context of quantum chemistry, chemical accuracy is often approximated as **`1 mHa`**. The energy estimate can be further improved by drawing more samples from the circuit, increasing the number of samples per batch, and increasing the number of iterations.

These results are better than those shown in the **“Visualize the results”** section of the tutorial *Sample-based Quantum Diagonalization of a Chemistry Hamiltonian* ([link](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization#visualize-the-results)). In that tutorial, the authors report:

> “The first plot shows that after a couple of iterations we estimate the ground-state energy within approximately **`41 mHa`** (chemical accuracy is typically accepted to be **`1 kcal/mol` $\approx$ `1.6 mHa`**). The energy can be improved by allowing more iterations of configuration recovery or by increasing the number of samples per batch.”

The second plot shows the average occupancy of each spatial orbital after the final iteration. Both the spin-up and spin-down electrons occupy the first eight orbitals with high probability.

# References
## Sample-based quantum diagonalization (SQD)

[SQD-1] Sample-based quantum diagonalization (SQD) overview, https://quantum.cloud.ibm.com/docs/en/guides/qiskit-addons-sqd

[SQD-2] Osama M. Raisuddin, Haimeng Zhang, Mario Motta, Fabian M. Faulstich, From Promise to Practice: Benchmarking Quantum Chemistry on Quantum Hardware, arXiv:2512.01012 [quant-ph], 30 Nov 2025, https://doi.org/10.48550/arXiv.2512.01012

[SQD-3] Javier Robledo-Moreno et al., Chemistry Beyond the Scale of Exact Diagonalization on a Quantum-Centric Supercomputer, arXiv:2405.05068 [quant-ph], 13 Jul 2025, https://doi.org/10.48550/arXiv.2405.05068

[SQD-4] IBM Quantum Webinar Series: Quantum-Centric Supercomputing Algorithms and Use Cases, July 2025, https://www.youtube.com/watch?v=gkyH6X_yQyA

[SQD-5] Improving energy estimation of a chemistry Hamiltonian with SQD, https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/tutorials/01_chemistry_hamiltonian.ipynb

[SQD-6] Sample-based Quantum Diagonalization (SQD), Quantum diagonalization algorithms, https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview

[SQD-7] Quantum diagonalization algorithms, IBM Quantum Platform, Quantum Learning, https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms

[SQD-8] Sample-based quantum diagonalization of a chemistry Hamiltonian, https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization

[SQD-9] Samuele Piccinelli et al., Quantum chemistry with provable convergence via randomized sample-based quantum diagonalization, 4 Aug 2025, arXiv:2508.02578 [quant-ph], https://doi.org/10.48550/arXiv.2508.02578

[SQD-10] Antonio Mezzacapo: Quantum diagonalization methods for lattice models, IQuS - The InQubator For Quantum Simulation, Feb 19, 2025, https://www.youtube.com/watch?v=b1fhh71hY2g

[SQD-11] Javier Robledo-Moreno, Gavin Jones, Roberto Lo Nardo, Robert Davis, Lockheed Martin & IBM combine quantum computing with classical HPC in new research, IBM Quantum Research Blog, 22 May 2025, https://www.ibm.com/quantum/blog/lockheed-martin-sqd

[SQD-12] Ieva Liepuoniute, Kirstin D. Doney, Javier Robledo Moreno, Joshua A. Job, William S. Friend, Gavin O. Jones, Quantum-Centric Computational Study of Methylene Singlet and Triplet States, J. Chem. Theory Comput. 2025, 21, 10, 5062–5070, https://doi.org/10.1021/acs.jctc.5c00075

[SQD-13] Danil Kaliakin, Akhil Shajan, Fangchun Liang, Kenneth M. Merz Jr., Implicit Solvent Sample-Based Quantum Diagonalization,J. Phys. Chem. B 2025, 129, 23, 5788–5796, https://doi.org/10.1021/acs.jpcb.5c01030

[SQD-14] Deploy and Run the SQD IEF-PCM Function Template, qiskit-function-templates/chemistry/sqd_pcm/deploy_and_run.ipynb, https://github.com/qiskit-community/qiskit-function-templates/blob/main/chemistry/sqd_pcm/deploy_and_run.ipynb

[SQD-15] J. Robledo-Moreno et al., "Chemistry Beyond Exact Solutions on a Quantum-Centric Supercomputer" (2024), arXiv:quant-ph/2405.05068, https://arxiv.org/abs/2405.05068

[SQD-16] M. Motta et al., “Bridging physical intuition and hardware efficiency for correlated electronic states: the local unitary cluster Jastrow ansatz for electronic structure” (2023). Chem. Sci., 2023, 14, 11213, https://pubs.rsc.org/en/content/articlehtml/2023/sc/d3sc02516k

[SQD-17] Keita Kanno, Masaya Kohda, Ryosuke Imai, Sho Koh, Kosuke Mitarai, Wataru Mizukami, Yuya O. Nakagawa, Quantum-Selected Configuration Interaction: classical diagonalization of Hamiltonians in subspaces selected by quantum computers, arXiv:2302.11320 [quant-ph], https://doi.org/10.48550/arXiv.2302.11320

[SQD-18] Introduction to Qiskit patterns, https://quantum.cloud.ibm.com/docs/en/guides/intro-to-patterns

[SQD-19] Qiskit addon: sample-based quantum diagonalization (SQD), https://github.com/Qiskit/qiskit-addon-sqd

[SQD-20] Bounding the subspace dimension, https://github.com/Qiskit/qiskit-addon-sqd/blob/main/docs/how_tos/choose_subspace_dimension.ipynb

## Energetics of quantum computing

[EN-1] informatique quantique état de l’art, perspective et défis, Olivier Ezratty, SFGP, Paris, 5 novembre 2025, https://www.oezratty.net/Files/Conferences/Olivier%20Ezratty%20Informatique%20Quantique%20SFGP%20Nov2025.pdf 

[EN-2] Q2B25 Paris | Olivier Ezratty, Academic, Co Founder, Free Electron, EPITA, Quantum Energy Initiative, 
QC Ware, September 24-25 2025, https://www.youtube.com/watch?v=JVtm3pbesnA

[EN-3] Green quantum computing, Capgemini, 8 May 2023, https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/

## Quantum chemistry

[CHEM-1] Keeper L. Sharkey and Alain Chancé, Quantum Chemistry and Computing for the Curious: Illustrated with Python and Qiskit® code, Packt 2022, https://a.co/d/6YCwgPb

[CHEM-2] Companion Jupyter notebook for Chapter 5 of the book Quantum Chemistry and Computing for the Curious: Illustrated with Python and Qiskit® code, https://github.com/AlainChance/Quantum-Chemistry-and-Computing-for-the-Curious/blob/main/Chapter_05_Variational_Quantum_Eigensolver_(VQE)_algorithm_V4.ipynb

[CHEM-3] Taha Selim, Ad van der Avoird, Gerrit C. Groenenboom, State-to-state rovibrational transition rates for CO2 in the bend mode in collisions with He atoms, J. Chem. Phys. 159, 164310 (2023), https://doi.org/10.1063/5.0174787