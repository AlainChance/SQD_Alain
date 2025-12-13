#------------------------------
# Import required dependencies
#------------------------------
import json
import ast
import time

from SQD_Alain import SQD

#----------------------------------------------------------------------------------------------------
# Define a function that converts atom coordinates from character strings to numeric values (floats)
# This function is copied from the SQD_Alain_gradio.py file
# https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain_Gradio.py
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
# This function is copied from the SQD_Alain_gradio.py file
# https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain_Gradio.py
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

#------------------------------------------------------------------------
# Define a function that runs a SQD simulation
# This function is derived from the SQD_Alain_gradio.py file
# https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain_Gradio.py
#------------------------------------------------------------------------
def run_simulation(config=None, config_filename=None, do_print=True):

    #-------------------------------------------------
    # Write configuration dictionary into a Json file
    #-------------------------------------------------
    if config_filename is not None:
        try:
            with open(config_filename, "w") as f:
                json.dump(config, f, indent=4)
        except Exception as e:
            return(f"Error saving configuration: {e}", None)

    #-----------------------------------------------------------
    # Retrieve basis and atom from the configuration dictionary
    #-----------------------------------------------------------
    basis = config['basis']
    atom = config['atom']
    symmetry = config['symmetry']

    #---------------------
    # Print configuration
    #---------------------
    if do_print:
        print(json.dumps(config, indent=4))

    #--------------------------------------------------------------------------------
    # Convert the atom we get from a character string to a tuple 
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

    #------------------------------------------------------------------------
    # The symmetry we get is expected to be "True", "Dooh", "Coov", or "D2h"
    # https://pyscf.org/user/gto.html#point-group-symmetry
    # Convert "True" and "False" to corresponding booleans
    #------------------------------------------------------------------------
    if isinstance(symmetry, str):
        symmetry = symmetry.strip()
        if symmetry == "True":
            config["symmetry"] = True
        elif symmetry == "False":
            config["symmetry"] = False
    
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
    if config_filename is not None:
        print(f"\nFind an approximation to the molecule defined in configuration {config_filename} in the {basis} basis set")

    t0 = time.time()  # ⏱️ Start timing

    # STEP 1
    print("\nstep_1 - Perform a CCSD calculation")
    try:
        plt_circuit = My_SQD.step_1()
        t1 = time.time()
        print(f"✅ step_1 completed in {t1 - t0:.2f} seconds")
    except Exception as e:
        return(f"Error in step 1: {e}", None)

    # STEP 2
    print("\nstep_2 - Optimize the circuit for a target hardware")
    try:
        isa_circuit = My_SQD.step_2()
        t2 = time.time()
        print(f"✅ step_2 completed in {t2 - t1:.2f} seconds")
    except Exception as e:
        return(f"Error in step 2: {e}", None)

    # STEP 3
    print("\nstep_3 - Execute using Qiskit Primitives or generate random samples")
    try:
        bit_array = My_SQD.step_3()
        t3 = time.time()
        print(f"✅ step_3 completed in {t3 - t2:.2f} seconds")
    except Exception as e:
        return(f"Error in step 3: {e}", None)

    # POST PROCESS
    print("\npost_process - Self-consistent configuration recovery procedure")
    result, result_history = My_SQD.post_process()
    t4 = time.time()
    print(f"✅ post_process completed in {t4 - t3:.2f} seconds")

    # PLOT
    print("\nplot_energy_and_occupancy")
    try:
        plt_energy = My_SQD.plot_energy_and_occupancy()
        t5 = time.time()
        print(f"✅ plot_energy_and_occupancy completed in {t5 - t4:.2f} seconds")
    except Exception as e:
        return(f"Error in plot_energy_and_occupancy: {e}", None)

    text = f"SQD configuration saved into {config_filename} and simulation complete."
    text += f"\nExact energy: {My_SQD.param['exact_energy']:.5f} Ha"
    text += f"\nSQD energy: {My_SQD.param['SQD_energy']:.5f} Ha"
    text += f"\nAbsolute error: {My_SQD.param['Absolute_error']:.5f} Ha"

    print(text)
    
    return config