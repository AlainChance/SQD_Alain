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

#-----------------------------------------------------------------------------------
# Define function get_QPU_power() 
# Prompt for a ballpark figure for the QPU power consumption
# A ballpark figure for a typical modern IBM-class superconducting quantum computer 
# (including cryogenics + support, while idle or lightly used): ~ 15–25 kW. 
# Source: Green quantum computing, Capgemini, 8 May 2023,
# https://www.capgemini.com/insights/expert-perspectives/green-quantum-computing/.
#-----------------------------------------------------------------------------------
def get_QPU_power(power_QPU=25.0):

    power_QPU = float(power_QPU)
    
    if power_QPU >= 15.0 and power_QPU <= 25.0:
        return power_QPU

    text = "Enter a rough estimate between 15 and 25 (kW) for QPU power consumption"

    while True:
        try:
            power_QPU = float(input(text))
            print("You entered:", power_QPU)
            power_QPU = min(max(15.0, power_QPU), 25.0)
            break
        except ValueError:
            print("Please enter a rough estimate between 15 and 25 (kW)")

    print(f"QPU Power consumption (kW):{power_QPU}")

    return power_QPU

#------------------------------------------------------------------------
# Define a function that runs a SQD simulation
# This function is derived from the SQD_Alain_gradio.py file
# https://github.com/AlainChance/SQD_Alain/blob/main/SQD_Alain_Gradio.py
#------------------------------------------------------------------------
def run_simulation(config=None, do_print=True):

    if config is None:
        return None

    #----------------------------------------------------------------------------
    # Retrieve configuration file name (.json) from the configuration dictionary
    #----------------------------------------------------------------------------
    try:
        config_filename = config['config_filename']
    except Exception as e:
        return(f"Error retrieving config_filename: {e}")

    #-------------------------------------------------
    # Write configuration dictionary into a Json file
    #-------------------------------------------------
    if config_filename is not None:
        try:
            with open(config_filename, "w") as f:
                json.dump(config, f, indent=4)
        except Exception as e:
            return(f"Error saving configuration: {e}")

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
        return(f"Invalid atom configuration {atom}.")

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
        return(f"Error creating SQD instance: {e}")

    #---------------
    # Run all steps
    #---------------
    result = My_SQD.run()
    
    return result