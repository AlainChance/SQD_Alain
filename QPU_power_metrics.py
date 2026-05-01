#--------------------------------------------------------------------------------
# QPU_power_metrics module by Alain Chancé

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
import re

from qiskit_ibm_runtime import QiskitRuntimeService

#---------------------------------------------------------------
# Define function get_valid_string() which prompts for a string 
# (letters and digits only, no spaces or special characters)
#---------------------------------------------------------------
def get_valid_string(text="input text"):
    while True:
        user_input = input(text)

        # Regex: only allows A-Z, a-z, 0-9
        if re.fullmatch(r"[A-Za-z0-9]+", user_input):
            print("Valid input job_id:", user_input)
            return user_input
        else:
            print("❌ Invalid input. Please try again.")

#------------------------------------------------------------
# Define function get_QPU_usage() 
# Prompt for a job, get job_id and fetch metrics from server
#------------------------------------------------------------
def get_QPU_usage(job_id=None):

    if job_id is None:
        job_id = get_valid_string(
            text="Enter job_id as a string (letters and digits only, no spaces or special characters): "
        )

    try:
        service = QiskitRuntimeService()
    except Exception as e:
        print(f"Error creating an instance of QiskitRuntimeService(): {e}")
        return None

    try:
        job = service.job(job_id)
    except Exception as e:
        print(f"Error getting job: {e}")
        return None

    try:
        metrics = job.metrics()           # Fetch metrics from server
        usage = metrics.get("usage", {})
        QPU_usage = usage.get("quantum_seconds")
        print(f"\nQiskit Runtime usage: {QPU_usage:.2f} seconds\n")
    except Exception as e:
        print(f"Error fetching job metrics: {e}")
        return None
        
    return QPU_usage

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

#-------------------------------------------------------------------------
# Define function get_QPU_power_usage()
#
# Input:
# - job_id: job id of a job or None (default)
# - power_QPU: a ballpark figure for the QPU power consumption
# 
# Prompt for a job, get job_id, retrieve QPU usage, get power_QPU
# and print a rough estimate of the power consumption on the IBM QPU (kWh)
#--------------------------------------------------------------------------
def get_QPU_power_usage(job_id=None, power_QPU=25.0):

    # Get QPU usage (s)
    QPU_usage = get_QPU_usage(job_id)

    if QPU_usage is None:
        return None, None
    print(f"\nQiskit Runtime usage: {QPU_usage:.2f} seconds")

    try:
        print("\nA rough estimate of the energy consumption is computed as power_QPU (kW) * QPU_usage (s) * 1/3600 (h/s)")
        print(f"power_QPU (kW): {power_QPU}")
        print(f"QPU_usage/3600 (h): {QPU_usage/3600:.4f}")
        QPU_power_consumption = power_QPU*QPU_usage/3600
        print(f"Rough estimate of the power consumption on the IBM QPU (kWh) = {QPU_power_consumption:.4f}")
    except Exception as e:
        print(f"Error computing a rough estimate of the power consumption on the IBM QPU : {e}")
        return QPU_usage, None

    return QPU_usage, QPU_power_consumption