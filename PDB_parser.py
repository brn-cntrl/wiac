"""
Temporary script for testing PDB sonification strategies. The final data parsing will be 
implemented in Unity (C# or C++) making the OSC and tkinter libraries eventually unnecessary.
Parameters passed here to the OSC parser are equivalent to the @hv_param variables used in 
Pure Data and the Heavy compiler.

This version sonifies only the data pertaining to R Groups. Future iterations will make use 
of data from protein backbone data as well. 

Brian Cantrell, Worldbuilding Media Lab. Sept, 2020.
"""
from pythonosc import osc_message_builder
from pythonosc import osc_bundle_builder
from pythonosc import udp_client
import os
import time

#################################################################################
# OSC SETUP #####################################################################
port = 5005
ip = "127.0.0.1"
client = udp_client.SimpleUDPClient(ip, port)
BFactorMsg = 0
hydroMsg = 0
# TODO: test to see if boolean values work with PD
newAsymMsg = False
newEntityMsg = False
newAAMsg = False

#################################################################################
#GLOBAL VARIABLES ###############################################################
# Column indices (WARNING: THESE MIGHT CHANGE ACCORDING TO FILE TYPE)
BFactorCol = 14 #Column for temperature (B Factor)
atomCol = 3 # Column for atom name
AACol = 5 # Column for amino acid label
chainCol = 8 # Column for the chain number
asymCol = 6 # Column for the asymmetry chain letter
entityCol = 7 # Column for the entity number

# Lists to hold the relevant data
RGroupAAs = [] # List to hold R Group amino acids (used to index hydroVals)
RGroupChains = [] # List to hold R Group chain numbers (corresponding to amino acid)
RGroupBFactors = [] # List to hold R Group B factors
RGroupAsyms = [] # List to hold R Group asym values
RGroupEntities = [] # List to hold R Group entity values

# Variables to hold temp values
newAsymVal = ' '
oldAsymVal = ' '
newEntityVal = 0
oldEntityVal = 20000 # Set oldEntityVal to something high to prevent false positives
newAAVal = ' '
oldAAVal = ' '
newChainVal = 0
oldChainVal = 20000 # Set high to prevent false positives

# Booleans for "new AA," "new asym," and "new entity" messages
newAA = True
newAsym = True
newEntity = True
newChain = True
oldAA = False
oldAsym = False
oldEntity = False
oldChain = False

PDBFile = " " # Hold the selected file
speed = .3 # Global speed for the OSC engine
 # Create an iterator to use for timing
iterator = 0

# Dictionary of hydrophobicity values
hydroVals = { 
        "PHE": 100,
        "ILE":  99,
        "TRP":  97,
        "LEU":  97,
        "VAL":  76,
        "MET":  74,
        "TYR":  63,
        "CYS":  49,
        "ALA":  41,
        "THR":  13,
        "HCS":   0,
        "HIS":   8,
        "SER":  -5,
        "GLN": -10,
        "ARG": -14,
        "LYS": -23,
        "ASN": -28,
        "GLU": -31,
        "PRO": -46,
        "ASP": -55 
    }

#################################################################################
# Open a file
def openFile(data):
    file = open(data, 'r')

    # make global variables visible
    
    global RGroupAAs
    global RGroupBFactors
    global RGroupAsyms
    global RGroupEntities
    global RGroupChains
    
    # Clear lists just in case...
    RGroupAAs.clear()
    RGroupBFactors.clear()
    RGroupAsyms.clear()
    RGroupEntities.clear()
    RGroupChains.clear()
    
    if file is not None:
        # Split the file up by line
        theFile = file.readlines()
         # Loop over file and select for ATOM entries
        for line in theFile:
            entry = line.split()
            # Select for only those lines listing ATOM coordinates
            if entry[0] == "ATOM":
                # Ignore the lines with backbone entries
                atomName = str(entry[atomCol])
                atomName.strip()
                
                if atomName == "CA" or atomName == "C" or atomName == "N" or atomName == "O" or atomName == "CB":
                    pass
                # Push R Group amino acid labels and B Factor values to lists
                else:
                    #print(atomName)
                    RGroupAAs.append(entry[AACol])
                    RGroupBFactors.append(entry[BFactorCol])
                    RGroupAsyms.append(entry[asymCol])
                    RGroupEntities.append(entry[entityCol])
                    RGroupChains.append(entry[chainCol])
    
    else:
        println("file is none")
    
#################################################################################
# Start the OSC engine
def runOSC():
    global iterator
    global newAsymVal
    global oldAsymVal
    global newEntityVal
    global oldEntityVal
    global newChainVal
    global oldChainVal
    global newAAVal
    global oldAAVal

    global newAsym
    global newEntity
    global newChain
    global newAA

    global BFactorMsg
    global hydroMsg

    # Set B Factor message
    BFactorMsg = RGroupBFactors[iterator]
    client.send_message("/BFactor", float(BFactorMsg))
    # Retrieve hydrophobicity from table by indexing the amino acids list
    aminoAcid = RGroupAAs[iterator]
    # Set hydrophobicity message
    hydroMsg = hydroVals[aminoAcid]
    # Set new asym message
    # Check if new asym is same as old asym
    newAsymVal = RGroupAsyms[iterator]
    newEntityVal = RGroupEntities[iterator]
    #newAA = RGroupAAs[iterator]
    newChainVal = RGroupChains[iterator]

    # New asym logic
    if newAsymVal != oldAsymVal:
        newAsym = True
        client.send_message("/newAsym", newAsym)
        oldAsymVal = newAsymVal
        newAsym = False
    else:
        pass
    
    # New entity logic
    if newEntityVal != oldEntityVal:
        newEntity = True
        client.send_message("/newEntity", newEntity)
        oldEntityVal = newEntityVal
        newEntity = False
    else:
        pass
    
    # Redundancy with newChain
    '''
    # New AA logic
    if newAAVal != oldAAVal:
        newAA = True
        client.send_message("/newAA", newAA)
        oldAAVal = newAAVal 
        newAA = False
    else:
        pass
    '''

    # New chain # logic
    if newChainVal != oldChainVal:
        newChain = True
        client.send_message("/newChain", newChain)
        oldChainVal = newChainVal
        newChain = False
    else:
        pass
    
    # Set BFactor and Hydrophibicity messages
    print(BFactorMsg)
    
    client.send_message("/hydrophobicity", hydroMsg)

    # Increment iterator    
    iterator += 1
    # Take modulus to loop back to beginning of lists
    iterator %= len(RGroupAAs)
    # Sleep to control speed
    time.sleep(speed)

#################################################################################
# Main Loop
def main():
    pdbFile = input("\nPlease enter the path/name of the .cif or .pdb file and press 'Enter': \n\n")
    openFile(pdbFile)
    # TODO: Implement better quit code
    print("Press ctrl+c to quit.")
    while True:
        runOSC()
      
if __name__ == '__main__': main()
