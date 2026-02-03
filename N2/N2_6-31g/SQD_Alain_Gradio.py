#--------------------------------------------------------------------------------
# # Sample-based Quantum Diagonalization (SQD) by Alain Chancé

## MIT License

# Copyright (c) 2025 Alain Chancé

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#--------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# Install gradio
# Gradio is an open-source Python package that allows you to quickly build a demo or web application for your machine learning model, 
# API, or any arbitrary Python function. You can then share a link to your demo or web application in just a few seconds using 
# Gradio's built-in sharing features. No JavaScript, CSS, or web hosting experience needed! https://www.gradio.app/guides/quickstart
#-------------------------------------------------------------------------------------------------------------------------------------
import gradio as gr

from SQD_Alain import SQD

# Import QiskitRuntimeService
from qiskit_ibm_runtime import QiskitRuntimeService

import json
import os
import ast
import time

#----------------------------------------------------------------------------------------------------
# Define a function that converts atom coordinates from character strings to numeric values (floats)
#----------------------------------------------------------------------------------------------------
def convert_atom_coordinates(atom_list):
    converted = []
    for element, coords in atom_list:
        # Convert each coordinate to float if it's a string
        numeric_coords = tuple(float(c) if isinstance(c, str) else c for c in coords)
        converted.append([element, numeric_coords])
    return converted

#--------------------------------------------------------------------------
# Define a function that loads a configuration dictionary from a Json file
#--------------------------------------------------------------------------
def load_sqd_config(filename=None):
    try:
        with open(filename, "r") as f:
            config = json.load(f)

        if not isinstance(config, dict):
            raise TypeError("Loaded config is not a dictionary.")

        if "atom" in config:
            atom_raw = config["atom"]

            # If it's a string, evaluate it safely
            if isinstance(atom_raw, str):
                atom_parsed = ast.literal_eval(atom_raw)
            else:
                atom_parsed = atom_raw

            config["atom"] = convert_atom_coordinates(atom_parsed)

        return config

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except json.JSONDecodeError:
        print(f"Error: File '{filename}' contains invalid JSON.")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return None

#--------------------------------------------------------------------------
# Define a function that loads a configuration dictionary from a Json file
# and returns updated values for each UI component
#--------------------------------------------------------------------------
def load_and_reset_config(config_filename):
    
    config = load_sqd_config(filename=config_filename)
    
    if config is None:
        return [gr.update(value="Invalid configuration.")] + [None] * 23

    # Extract values from config
    backend = config.get("backend_name", "None")
    job_id = config.get("job_id", "None")                    # job_id of a previously run job
    do_plot = config.get("do_plot_gate_map", True)
    load_bit = config.get("load_bit_array_file", "None")
    save_bit = config.get("save_bit_array_file", "None")
    run_qpu = config.get("run_on_QPU", False)
    nshots = config.get("nshots", 1000)                      # Number of shots
    nsamples = config.get("nsamples", 10000)/1000            # Number of samples (×1,000) to be drawn from the uniform distribution
    #-------------------------------------
    # eco2AI Tracker options
    # https://github.com/sb-ai-lab/Eco2AI
    #-------------------------------------
    do_eco2ai = config.get("do_eco2ai", True)                                   # Whether to track energy usage with eco2AI
    project_name = config.get("project_name", "SQD_Alain")                      # Project name
    experiment_desc = config.get("experiment_description", "SQD_experiment")    # Experiment description
    eco2ai_file_name = config.get("eco2ai_file_name", "SQD.csv")                # eco2AI file name
    #---------------------------------------------------------------------------------
    # Ballpark figure (kW) for the power consumption of the IBM cloud backend
    # "The power consumption of a quantum computer is about 15-25kW"
    # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
    #---------------------------------------------------------------------------------
    power_QPU = config.get("power_QPU", 25)                                     # QPU power consumption (kW)
    #---------------------------------------------------------
    # PySCF options
    # https://pyscf.org/user/gto.html#initializing-a-molecule
    #---------------------------------------------------------
    basis = config.get("basis", "6-31g")
    atom = str(config.get("atom", []))
    spin = config.get("spin", 0)
    symmetry = str(config.get("symmetry", "True"))
    frozen = config.get("n_frozen", 0)
    #---------------------
    # Molecular integrals
    #---------------------
    compute_exact_energy = config.get("compute_exact_energy", True)
    exact_energy = config.get("exact_energy", -109.23)
    #----------------------------------------------------------------------------------
    # SQD options
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
    #----------------------------------------------------------------------------------
    max_iter = config.get("max_iterations", 10)
    #-----------------------------------------------------------------------------------
    # Eigenstate solver options
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
    #-----------------------------------------------------------------------------------
    num_batches = config.get("num_batches", 5)
    samples_per_batch = config.get("samples_per_batch", 300)
    
    # Return updates for each component
    return (
        gr.update(value="Configuration loaded."),  # status
        gr.update(value=backend),
        gr.update(value=job_id),                   # job_id of a previously run job
        gr.update(value=do_plot),
        gr.update(value=load_bit),
        gr.update(value=save_bit),
        gr.update(value=run_qpu),
        gr.update(value=nshots),                   # Number of shots
        gr.update(value=nsamples),                 # Number of samples (×1,000) to be drawn from the uniform distribution
        gr.update(value=do_eco2ai),
        gr.update(value=project_name),
        gr.update(value=experiment_desc),
        gr.update(value=eco2ai_file_name),
        gr.update(value=power_QPU),                # QPU power consumption (kW)
        gr.update(value=basis),
        gr.update(value=atom),
        gr.update(value=spin),
        gr.update(value=symmetry),
        gr.update(value=frozen),
        gr.update(value=compute_exact_energy),
        gr.update(value=exact_energy),
        gr.update(value=max_iter),
        gr.update(value=num_batches),
        gr.update(value=samples_per_batch)
    )

#-------------------------------------------------------------------------------
# Define a function that returns "None" and a list of real operational backends
#-------------------------------------------------------------------------------
def list_backends(n_qubits=22):
    try:
        service = QiskitRuntimeService()
    except:
        service = None

    r_backends = ["None"]
    
    if service is not None:
        try:
            backends = service.backends(None, min_num_qubits=n_qubits, simulator=False, operational=True)
        except Exception as e:
            backends = []
        
        r_backends = ["None"] + [f"{backend.name}" for backend in backends]

    return r_backends

#--------------------------------------------------------------------------------------------------
# Define a function that returns the name of first file ending with .json in the current directory 
# or None if there is none
#--------------------------------------------------------------------------------------------------
def list_json_files():
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]

    return None if json_files is None else json_files[0]

#------------------------------------------------------
# Define a function that process inputs from Gradio UI
#------------------------------------------------------
def process_inputs(
    #---------------------------------
    # Parameters in config dictionary
    #---------------------------------
    # Run options
    #-------------
    config_filename,
    backend_name,
    job_id,                    # job_id of a previously run job
    do_plot_gate_map,
    load_bit_array_file,
    save_bit_array_file,
    run_on_QPU,
    nshots,                    # Number of shots
    nsamples,                  # Number of samples (×1,000) to be drawn from the uniform distribution
    #-------------------------------------
    # eco2AI Tracker options
    # https://github.com/sb-ai-lab/Eco2AI
    #-------------------------------------
    do_eco2ai,                 # Whether to track energy usage with eco2AI
    project_name,              # Project name
    experiment_description,    # Experiment description
    eco2ai_file_name,          # eco2AI file name
    #---------------------------------------------------------------------------------
    # Ballpark figure (kW) for the power consumption of the IBM cloud backend
    # "The power consumption of a quantum computer is about 15-25kW"
    # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
    #---------------------------------------------------------------------------------
    power_QPU,                  # QPU power consumption (kW)
    #---------------------------------------------------------
    # PySCF options
    # https://pyscf.org/user/gto.html#initializing-a-molecule
    #---------------------------------------------------------
    basis,
    atom,
    spin,
    symmetry,
    n_frozen,
    #---------------------
    # Molecular integrals
    #---------------------
    compute_exact_energy,
    exact_energy,
    #----------------------------------------------------------------------------------
    # SQD options
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
    #----------------------------------------------------------------------------------
    max_iterations,
    #-----------------------------------------------------------------------------------
    # Eigenstate solver options
    # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
    #-----------------------------------------------------------------------------------
    num_batches,
    samples_per_batch,
    #-------------------------------------
    # Parameters not in config dictionary
    #-------------------------------------
    run_simulation
): 
    config = {
        #-------------
        # Run options
        #-------------
        "config_filename": None if config_filename == "None" else config_filename,
        "backend_name": None if backend_name == "None" else backend_name,
        "job_id": None if job_id in ["None", ""] else job_id,                                      # job_id of a previously run job
        "do_plot_gate_map": do_plot_gate_map,
        "load_bit_array_file": None if load_bit_array_file in ["None", ""] else load_bit_array_file,
        "save_bit_array_file": None if save_bit_array_file in ["None", ""] else save_bit_array_file,
        "n_ancillary_qubits": 0,
        "run_on_QPU": run_on_QPU,
        "nshots": nshots,                                    # Number of shots
        "nsamples": nsamples*1000,                           # Number of samples to be drawn from the uniform distribution
        #-------------------------------------
        # eco2AI Tracker options
        # https://github.com/sb-ai-lab/Eco2AI
        #-------------------------------------
        "do_eco2ai": do_eco2ai,                              # Whether to track energy usage with eco2AI
        "project_name": project_name,                        # Project name
        "experiment_description": experiment_description,    # Experiment description
        "eco2ai_file_name": eco2ai_file_name,                # eco2AI file name
        #---------------------------------------------------------------------------------
        # Ballpark figure (kW) for the power consumption of the IBM cloud backend
        # "The power consumption of a quantum computer is about 15-25kW"
        # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
        #---------------------------------------------------------------------------------
        "power_QPU": power_QPU,                              # QPU power consumption (kW)
        #---------------------------------------------------------
        # PySCF options
        # https://pyscf.org/user/gto.html#initializing-a-molecule
        #---------------------------------------------------------
        "basis": basis,
        "atom": atom,
        "spin": spin,
        "symmetry": symmetry,
        "n_frozen": n_frozen,
        #---------------------
        # Molecular integrals
        #---------------------
        "compute_exact_energy": compute_exact_energy,
        "exact_energy": exact_energy,
        #----------------------------------------------------------------------------------
        # SQD options
        # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
        #----------------------------------------------------------------------------------
        "chem_accuracy": 1e-3,
        "energy_tol": 3e-5,
        "occupancies_tol": 1e-3,
        "max_iterations": max_iterations,
        #-----------------------------------------------------------------------------------
        # Eigenstate solver options
        # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
        #-----------------------------------------------------------------------------------
        "num_batches": num_batches,
        "samples_per_batch": samples_per_batch,
        "symmetrize_spin": True,
        "carryover_threshold": 1e-4,
        "max_cycle": 200,
        "seed": 24,
        "spin_sq": 0.0
    }
    
    #-------------------------------------------------
    # Write configuration dictionary into a Json file
    #-------------------------------------------------
    try:
        with open(config_filename, "w") as f:
            json.dump(config, f, indent=4)
    except Exception as e:
        return(f"Error saving configuration: {e}", None)
    
    #--------------------------------------------------------------------------------
    # Convert the atom we get from Gradio UI from a character string to a tuple 
    # and convert atom coordinates from character strings to numeric values (floats). 
    #--------------------------------------------------------------------------------
    # If atom is a string, evaluate it safely
    if isinstance(atom, str):
        atom_parsed = ast.literal_eval(atom)
    else:
        atom_parsed = atom

    if atom_parsed is not None:
        config["atom"] = convert_atom_coordinates(atom_parsed)
    else:
        return(f"Invalid atom configuration {atom}.", None)

    #----------------------------------------------------------------------------------------
    # The symmetry we get from Gradio UI are expected to be "True", "Dooh", "Coov", or "D2h"
    # https://pyscf.org/user/gto.html#point-group-symmetry
    # Convert "True" and "False" to corresponding booleans
    #-----------------------------------------------------------------------------------------
    if isinstance(symmetry, str):
        symmetry = symmetry.strip()
        if symmetry == "True":
            config["symmetry"] = True
        elif symmetry == "False":
            config["symmetry"] = False

    if run_simulation:
        #-----------------------------------------------------------------------------
        # Create an instance of the SQD_Alain class from the configuration dictionary
        #-----------------------------------------------------------------------------
        try:
            My_SQD = SQD(**config)
        except Exception as e:
            return(f"Error creating SQD instance: {e}", None)

        #---------------
        # Run all steps
        #---------------
        try:
            result = My_SQD.run()
        except Exception as e:
            return(f"Error running SQD instance: {e}", None)
    
        #---------------------------
        # Energy and occupancy plot
        #---------------------------
        plt_energy = My_SQD.param['plt_energy']
        try:
            gr_plt_energy = gr.Plot(plt_energy)
        except:
            gr_plt_energy = None

        text = f"SQD configuration saved into {config_filename} and {result}."
        text += f"\nExact energy (Ha): {My_SQD.param['exact_energy']:.5f}"
        text += f"\nSQD energy (Ha): {My_SQD.param['SQD_energy']:.5f}"
        text += f"\nAbsolute error (Ha): {My_SQD.param['Absolute_error']:.5f}"

        #----------------------------------------------------
        # Print Qiskit Runtime usage in seconds and in hours
        #----------------------------------------------------
        QPU_usage = My_SQD.param['QPU_usage']
        if QPU_usage is not None:
            text += f"\nQPU usage (s): {QPU_usage:.2f}, (h): {QPU_usage/3600.0:.4f}"

        #------------------------------------------------------------------------------------
        # Print a rough estimate of the energy consumption of the quantum device
        # **Assumption**
        # A ballpark figure for a typical modern IBM-class superconducting quantum computer 
        # (including cryogenics + support, while idle or lightly used): ~ 15–25 kW. 
        # Source: [Green quantum computing, Capgemini, 8 May 2023]
        # (https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/).
        #------------------------------------------------------------------------------------
        QPU_power_consumption = My_SQD.param['QPU_power_consumption']
        if QPU_power_consumption is not None:
            text += f"\nRough estimate for QPU power consumption (kWh): {QPU_power_consumption:.4f}"

        #------------------------------------------------------------------
        # Print the duration and power consumption of classical processing
        #-------------------------------------------------------------------
        duration = My_SQD.param['duration']
        classical_power_usage = My_SQD.param['classical_power_usage']
        if duration is not None and classical_power_usage is not None:
            text += f"\nClassical processing - Duration (h): {duration:.4f} "
            text += f"- Power consumption (kWh): {classical_power_usage:.4f}"

        #--------
        # Return
        #--------
        return(text, gr_plt_energy)      # Energy and occupancy plot
  
    else:
        return(f"SQD configuration saved into {config_filename}.", None)

#-----------
# Gradio UI
#-----------
# Get list of backends
backend_options = list_backends()
    
with gr.Blocks() as demo:
    gr.Markdown("## Configure And Run SQD Simulation")

    with gr.Row():
        basis = gr.Textbox(label="Basis Set", value="6-31g")

        job_id = gr.Textbox(label="Job id or 'None'", value='None') # job_id of a previously run job
        
        load_bit_array_file = gr.Textbox(
            label="Load bit array file name or 'None'",
            value='None'
        )
        save_bit_array_file = gr.Textbox(
            label="Save bit array file name or 'None'",
            value='None'
        )

    with gr.Row():
        run_on_QPU = gr.Checkbox(label="Run on QPU", value=True)

        backend_name = gr.Dropdown(
        label="IBM Backend Name (or 'None' for least busy)",
        choices=backend_options,
        value="None",
        visible=True
        )

        nshots = gr.Slider(1000, 200000, step=1000, label="Number of shots", value=1000, visible=True)
        do_plot_gate_map = gr.Checkbox(label="Plot Gate Map", value=True, visible=True)

        run_on_QPU.change(
            fn=lambda checked: [gr.update(visible=checked)] * 3,
            inputs=run_on_QPU,
            outputs=[backend_name, nshots, do_plot_gate_map]
        )

        nsamples = gr.Slider(10, 10000, step=10, label="Number of samples (×1,000)", value=10, visible=False) # Hidden by default
        
        run_on_QPU.change(
            fn=lambda checked: gr.update(visible=not checked),
            inputs=run_on_QPU,
            outputs=nsamples
        )

    with gr.Row():
        atom_json = gr.Textbox(
            label="Atom Configuration (Python list format)",
            value='[["C", (0.0, 0.0, 0.0)], ["H", (0.0, 0.0, 1.1160)], ["H", (0.9324, 0.0, -0.3987)]]'
        )
        symmetry = gr.Textbox(
            label="Point group symmetry 'True', 'Dooh', 'Coov', or'D2h'",
            value="True")

    with gr.Row():
        spin = gr.Slider(0, 4, step=1, label="Spin", value=2)
        n_frozen = gr.Slider(0, 10, step=1, label="Number of Frozen Orbitals", value=0)
        max_iterations = gr.Slider(10, 20, step=1, label="Max recovery iterations", value=10)
        num_batches = gr.Slider(1, 10, step=1, label="Batches", value=5)
        samples_per_batch = gr.Slider(300, 1000, step=100, label="Samples per batch", value=300)

    with gr.Row():
        compute_exact_energy = gr.Checkbox(label="Compute Exact Energy", value=True)
        
        exact_energy = gr.Slider(
            minimum=-250.0,
            maximum=-1.17,
            value=-109.2269,
            step=0.0001, # float step → float output label="Exact Energy"
            label="Exact Energy",
            visible=False
        )
        
        compute_exact_energy.change(
            fn=lambda checked: gr.update(visible=not checked),
            inputs=compute_exact_energy,
            outputs=exact_energy
        )

    with gr.Row():
        power_QPU = gr.Slider(15, 25, step=1, label="QPU power consumption (kW)", value=25, visible=True)
        
        do_eco2ai = gr.Checkbox(label="eco2AI tracker", value=True)

        project_name = gr.Textbox(label="Project name", value="SQD_Alain", visible=True)
        experiment_description = gr.Textbox(label="Experiment description", value="SQD_experiment", visible=True)
        eco2ai_file_name = gr.Textbox(label="eco2AI file name", value="SQD.csv", visible=True)

        do_eco2ai.change(
        fn=lambda checked: [gr.update(visible=checked)] * 3,
        inputs=do_eco2ai,
        outputs=[project_name, experiment_description, eco2ai_file_name]
    )

    with gr.Row():
        config_filename = gr.Textbox(label="Json configuration file name", value=list_json_files())
        run_simulation = gr.Checkbox(label="Run simulation", value=True)

    load_btn = gr.Button("Load SQD configuration from a Json file")
    submit_btn = gr.Button("Save SQD configuration into a Json file and run simulation")
    output = gr.Textbox(label="Status", lines=8, max_lines=20)
    energy_plot = gr.Plot(label="Energy and Occupancy Plot")

    load_btn.click(
        fn=load_and_reset_config,
        inputs=[config_filename],
        outputs=[
            output,                # status textbox
            #-------------
            # Run options
            #-------------
            backend_name,
            job_id,                # job_id of a previously run job
            do_plot_gate_map,
            load_bit_array_file,
            save_bit_array_file,
            run_on_QPU,
            nshots,                 # Number of shots
            nsamples,               # Number of samples (×1,000) to be drawn from the uniform distribution
            #-------------------------------------
            # eco2AI Tracker options
            # https://github.com/sb-ai-lab/Eco2AI
            #-------------------------------------
            do_eco2ai,
            project_name,
            experiment_description,
            eco2ai_file_name,
            #---------------------------------------------------------------------------------
            # Ballpark figure (kW) for the power consumption of the IBM cloud backend
            # "The power consumption of a quantum computer is about 15-25kW"
            # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
            #---------------------------------------------------------------------------------
            power_QPU,                  # QPU power consumption (kW)
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis,
            atom_json,
            spin,
            symmetry,
            n_frozen,
            #---------------------
            # Molecular integrals
            #---------------------
            compute_exact_energy,
            exact_energy,
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            max_iterations,
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches,
            samples_per_batch,
        ]
    )
    
    submit_btn.click(
        fn=process_inputs,
        inputs=[
            # Parameters in config dictionary
            config_filename,
            #-------------
            # Run options
            #-------------
            backend_name,
            job_id,                # job_id of a previously run job
            do_plot_gate_map,
            load_bit_array_file,
            save_bit_array_file,
            run_on_QPU,
            nshots,                 # Number of shots
            nsamples,               # Number of samples (×1,000) to be drawn from the uniform distribution
            #-------------------------------------
            # eco2AI Tracker options
            # https://github.com/sb-ai-lab/Eco2AI
            #-------------------------------------
            do_eco2ai,
            project_name,
            experiment_description,
            eco2ai_file_name,
            #---------------------------------------------------------------------------------
            # Ballpark figure (kW) for the power consumption of the IBM cloud backend
            # "The power consumption of a quantum computer is about 15-25kW"
            # https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/
            #---------------------------------------------------------------------------------
            power_QPU,               # QPU power consumption (kW)
            #---------------------------------------------------------
            # PySCF options
            # https://pyscf.org/user/gto.html#initializing-a-molecule
            #---------------------------------------------------------
            basis,
            atom_json,
            spin,
            symmetry,
            n_frozen,
            #---------------------
            # Molecular integrals
            #---------------------
            compute_exact_energy,
            exact_energy,
            #----------------------------------------------------------------------------------
            # SQD options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #----------------------------------------------------------------------------------
            max_iterations,
            #-----------------------------------------------------------------------------------
            # Eigenstate solver options
            # https://github.com/Qiskit/qiskit-addon-sqd/blob/main/qiskit_addon_sqd/fermion.py
            #-----------------------------------------------------------------------------------
            num_batches,
            samples_per_batch,
            # Parameters not in config dictionary
            run_simulation
        ],
        outputs=[output, energy_plot]
    )

#-------------------------
# Launch Gradio interface
#-------------------------
demo.launch()