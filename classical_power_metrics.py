#--------------------------------------------------------------------------------
# Classical_power_metrics module by Alain Chancé

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

#------------------------------
# Import required dependencies
#------------------------------
import os
import pandas as pd
import re

#-------------------------------------------------------------------------
# Define function get_eco2AI_file which prompts for a eco2AI filename 
#-------------------------------------------------------------------------
def get_eco2AI_file(eco2AI_file=None):
    
    if eco2AI_file is None:
        prompt = True
    else:
        if os.path.isfile(eco2AI_file):
            prompt = False
        else:
            print("Missing eco2AI file: ", eco2AI_file)
            prompt = True

    if prompt:
        while True:
            eco2AI_file = input("Enter eco2AI filename (.csv): ").strip()
            if eco2AI_file.lower().endswith((".csv")):
                print("You entered: ", eco2AI_file)
                if os.path.isfile(eco2AI_file):
                    break
                else:
                    print("Missing eco2AI file: ", eco2AI_file)    
            else:
                print("Invalid file type. Please enter a .csv file.")
        
    return eco2AI_file

#---------------------------------------------------------------------------------
# Define function get_classical_power_usage()
#
# Input:
# - eco2ai_file: eco2ai file name or None (default)
# 
# Prompt for a eco2AI filename (.csv) and retrieve duration and power consumption
#---------------------------------------------------------------------------------
def get_classical_power_usage(eco2AI_file=None):

    eco2AI_file = get_eco2AI_file(eco2AI_file)

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

        # Sum power consumption (already in kWh)
        power_usage = df["power_consumption(kWh)"].sum()

        print(
            f"\nClassical processing - Duration (h): {duration:.4f} "
            f"- Power consumption (kWh): {power_usage:.4f}"
            )
    except Exception as e:
        print(f"Error accessing the duration and power consumption of the eco2AI file: {e}")
        return None, None
        
    return duration, power_usage