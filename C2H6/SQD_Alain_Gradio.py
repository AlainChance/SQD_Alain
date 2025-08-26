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

import json
import ast

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

#------------------------------------------------------
# Define a function that process inputs from Gradio UI
#------------------------------------------------------
def process_inputs(
    # Parameters in config dictionary
    backend_name,
    do_plot_gate_map,
    load_bit_array_file,
    save_bit_array_file,
    run_on_QPU,
    basis,
    atom,
    spin,
    symmetry,
    n_frozen,
    compute_exact_energy,
    max_iterations,
    # Parameters not in config dictionary
    load_config_file,
    config_filename,
    run_simulation,
    molecule_name
): 
    config = {
        "backend_name": None if backend_name == "None" else backend_name,
        "do_plot_gate_map": do_plot_gate_map,
        "load_bit_array_file": None if load_bit_array_file == "None" else load_bit_array_file,
        "save_bit_array_file": None if save_bit_array_file == "None" else save_bit_array_file,
        "n_ancillary_qubits": 0,
        "run_on_QPU": run_on_QPU,
        "nshots": 1000,
        "basis": basis,
        "atom": atom,
        "spin": spin,
        "symmetry": symmetry,
        "n_frozen": n_frozen,
        "compute_exact_energy": compute_exact_energy,
        "chem_accuracy": 1e-3,
        "energy_tol": 3e-5,
        "occupancies_tol": 1e-3,
        "max_iterations": max_iterations,
        "num_batches": 5,
        "samples_per_batch": 300,
        "symmetrize_spin": True,
        "carryover_threshold": 1e-4,
        "max_cycle": 200,
        "seed": 24,
        "spin_sq": 0.0
    }

    if load_config_file:
        #--------------------------------------------------
        # Load a configuration dictionary from a Json file
        #--------------------------------------------------       
        config = load_sqd_config(filename = config_filename)
        
        if config is None:
            return(
                f"Invalid configuration.",
                None
            )

        atom = config["atom"]

    else:
        #-------------------------------------------------
        # Write configuration dictionary into a Json file
        #-------------------------------------------------
        try:
            with open(config_filename, "w") as f:
                json.dump(config, f, indent=4)
        except Exception as e:
            return(
                f"Error saving configuration: {e}",
                None
            )
    
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
        return(
            f"Invalid atom configuration {atom}.",
            None
        )

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
            return f"Error creating SQD instance: {e}"

        #---------------
        # Run all steps
        #---------------
        print(f"\nFind an approximation to the {molecule_name} molecule in the {basis} basis set")
        
        print("\nstep_1 - Perform a CCSD calculation")
        plt_circuit = My_SQD.step_1()
        
        print("\nstep_2 - Optimize the circuit for a target hardware")
        isa_circuit = My_SQD.step_2()
        
        print("\nstep_3 - Execute using Qiskit Primitives or generate random samples drawn from the uniform distribution ")
        bit_array = My_SQD.step_3()
        
        print("\npost_process - Self-consistent configuration recovery procedure")
        result, result_history = My_SQD.post_process()
        
        print("\nplot_energy_and_occupancy")
        plt_energy = My_SQD.plot_energy_and_occupancy()

    #--------
    # Return
    #--------
    if run_simulation:
        #---------------------------
        # Energy and occupancy plot
        #---------------------------
        try:
            gr_plt_energy = gr.Plot(plt_energy)
        except:
            gr_plt_energy = None
        
        return (
            f"SQD configuration saved into {config_filename} and simulation complete.",
            gr_plt_energy      # Energy and occupancy plot
        )
    else:
        return (
            f"SQD configuration saved into {config_filename}.",
            None
        )
        
#-----------
# Gradio UI
#-----------
with gr.Blocks() as demo:
    gr.Markdown("## Configure And Run SQD Simulation")

    with gr.Row():
        config_filename = gr.Textbox(label="Json configuration file name", value="SQD_CH2_Triplet.json")
        load_config_file = gr.Checkbox(label="Load configuration file", value=False)
        run_simulation = gr.Checkbox(label="Run simulation", value=True)
        molecule_name = gr.Textbox(
            label="Molecule name",
            value="CH2_Triplet")

    with gr.Row():
        do_plot_gate_map = gr.Checkbox(label="Plot Gate Map", value=True)
        load_bit_array_file = gr.Textbox(
            label="Load bit array file name or 'None'",
            value='None'
        )
        save_bit_array_file = gr.Textbox(
            label="Save bit array file name or 'None'",
            value='None'
        )
        compute_exact_energy = gr.Checkbox(label="Compute Exact Energy", value=True)

    with gr.Row():
        run_on_QPU = gr.Checkbox(label="Run on QPU", value=False)
        backend_name = gr.Textbox(label="IBM Backend Name (or 'None')", value="None")
        basis = gr.Textbox(label="Basis Set", value="6-31g")

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
        max_iterations = gr.Slider(10, 20, step=1, label="Limit on the number of configuration recovery iterations", value=10)
    
    submit_btn = gr.Button("Save SQD configuration into a Json file and run simulation")
    output = gr.Textbox(label="Status")
    energy_plot = gr.Plot(label="Energy and Occupancy Plot")

    submit_btn.click(
        fn=process_inputs,
        inputs=[
            # Parameters in config dictionary
            backend_name,
            do_plot_gate_map,
            load_bit_array_file,
            save_bit_array_file,
            run_on_QPU,
            basis,
            atom_json,
            spin,
            symmetry,
            n_frozen,
            compute_exact_energy,
            max_iterations,
            # Parameters not in config dictionary
            load_config_file,
            config_filename,
            run_simulation,
            molecule_name
        ],
        outputs=[output, energy_plot]
    )

#-------------------------
# Launch Gradio interface
#-------------------------
demo.launch()