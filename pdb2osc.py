"""
Temporary script for testing PDB sonification strategies. The final data parsing will be 
implemented in Unity (C# or C++) making the OSC and tkinter libraries eventually unnecessary.
Parameters passed here to the OSC parser are equivalent to the @hv_param variables used in 
Pure Data and the Heavy compiler.

This version sonifies only the data pertaining to R Groups. Future iterations may make use 
of data from protein backbone data as well. 

NOTE: When porting to C# for Unity, all references to OSC should be scrubbed. All variables
and the importing of data should be made internal to C# and Unity with no calls over UDP.

Brian Cantrell, Worldbuilding Media Lab. Sept, 2020.
Updated Oct. 2020.
"""
from pythonosc import osc_message_builder
from pythonosc import osc_bundle_builder
from pythonosc import udp_client
import os
import time
import json

#################################################################################
# OSC SETUP #####################################################################
port = 5005
ip = "127.0.0.1"
client = udp_client.SimpleUDPClient(ip, port)
BFactorMsg = 0
hydroMsg = 0

# TODO: test to see if boolean values work with PD
# TODO: Get rid of global variables and make script OO with get methods 

newAsymMsg = False
newEntityMsg = False
newAAMsg = False
structTypeMsg = 5 # Initialize a value higher than the one used for message (0,1,2)
categoryMsg = 8 # Initialize a value higher than the one used for message (0,1,2,3,4)
iteratorMsg = 0

#################################################################################
#GLOBAL VARIABLES ###############################################################
# Column indices (WARNING: THESE MIGHT CHANGE ACCORDING TO FILE TYPE)
BFactorCol = 14 #Column for temperature (B Factor)
atomCol = 3 # Column for atom name
AACol = 5 # Column for amino acid label
chainCol = 8 # Column for the chain number
asymCol = 6 # Column for the asymmetry chain letter
entityCol = 7 # Column for the entity number
beginHelix = 5 # Column number for the start of the helix (gives amino acid number)
endHelix = 9 # Column number for the end of the helix (gives amino acid number)
beginSheet = 4 # Column number for the beginning AA of the sheet
endSheet = 8 # Column number for the ending AA of the sheet

# Lists to hold the relevant data
RGroupAAs = [] # List to hold R Group amino acids (used to index hydroVals)
RGroupChains = [] # List to hold R Group chain numbers (corresponding to amino acid)
RGroupBFactors = [] # List to hold R Group B factors
RGroupAsyms = [] # List to hold R Group asym values
RGroupEntities = [] # List to hold R Group entity values
helices = [] # List to hold helices for structType
sheets = [] # List to hold sheets for structType
loops = [] # List to hold loops for structType

# Variables to hold temp values
newAsymVal = ' '
oldAsymVal = ' '
newEntityVal = 0
oldEntityVal = 20000 # Set oldEntityVal to something high to prevent false positives
newAAVal = ' '
oldAAVal = ' '
newChainVal = 0
oldChainVal = 20000 # Set high to prevent false positives
newHydroVal = 0
oldHydroVal = 2000

# Booleans for "new AA," "new asym," and "new entity" messages
newAA = True
newAsym = True
newEntity = True
newChain = True
oldAA = False
oldAsym = False
oldEntity = False
oldChain = False
newHydro = True
oldHydro = False
# Boolean for determining whether the next entries are sheet structures
isSheet = False

PDBFile = " " # Hold the selected file
speed = .4 # Global speed for the OSC engine
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
        "GLY":   0,
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
    # Dictionary to Index Amino Acid Categories
aaCategories = {
        "PHE": 0, # "aromatic"
        "ILE": 1, # "aliphatic"
        "TRP": 0, # "aromatic"
        "LEU": 1, # "aliphatic"
        "VAL": 1, # "aliphatic"
        "MET": 2, # "polar"
        "TYR": 0, # "aromatic"
        "CYS": 2, # "polar"
        "ALA": 1, # "aliphatic"
        "THR": 2, # "polar"
        "HCS": 4, # "unique"
        "GLY": 4, # "unique"
        "HIS": 3, # "charged"
        "SER": 2, # "polar"
        "GLN": 2, # "polar"
        "ARG": 3, # "charged"
        "LYS": 3, # "charged"
        "ASN": 2, # "polar"
        "GLU": 3, # "charged"
        "PRO": 4, # "unique"
        "ASP": 3, # "charged" 
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
    global helices
    global beginHelix
    global endHelix

    global isSheet
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

            # Here is where we determine whether the current amino acid is part of a helix, loop, or sheet
            # Check for helix entry first
            if entry[0] == "HELX_P":
                firstHAA = int(entry[beginHelix])
                lastHAA = int(entry[endHelix])
                # Push the first and last AA number and every number between into helices list
                for i in range(firstHAA, lastHAA+1):
                    helices.append(i)

            # Fill sheets list
            if entry[0] == "_struct_sheet_range.end_auth_seq_id":
                # Assign boolean to get sheet values until "#" is reached.
                isSheet = True
            
            if isSheet == True and entry[0] != "_struct_sheet_range.end_auth_seq_id" and entry[0] != "#":
            
                firstSAA = int(entry[beginSheet])
                lastSAA = int(entry[endSheet])
                for i in range(firstSAA, lastSAA+1):
                    sheets.append(i)
            elif entry[0] == "#":
                isSheet = False

            
            # Move on to atomic coordinates for rest of data
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
        print("file is none")
    
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
    global newHydroVal
    global oldHydroVal

    global newAsym
    global newEntity
    global newChain
    global newAA
    global newHydro

    global BFactorMsg
    global hydroMsg
    global structTypeMsg
    global categoryMsg

    # Set B Factor message
    BFactorMsg = RGroupBFactors[iterator]
    client.send_message("/BFactor", float(BFactorMsg))
    
    # Retrieve hydrophobicity from table by indexing the amino acids list
    aminoAcid = RGroupAAs[iterator]
    
    # Set hydrophobicity message DO NOT DELETE
    hydroMsg = hydroVals[aminoAcid]
    client.send_message("/hydrophobicity", hydroMsg)
    # Set new asym message
    # Check if new asym is same as old asym
    newAsymVal = RGroupAsyms[iterator]
    newEntityVal = RGroupEntities[iterator]
    
    #newHydroVal = hydroVals[aminoAcid] DO NOT DELETE
    #newAA = RGroupAAs[iterator]
    
    newChainVal = RGroupChains[iterator]
    
    """ DO NOT DELETE
    # New hydrophobicity logic
    if newHydroVal != oldHydroVal:
        newHydro = True
        client.send_message("/hydrophobicity", hydroMsg)
        oldHydroVal = newHydroVal
        newHydro = False
    else:
        pass
    """

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

    # New chain number logic. 
    # In this conditional, we check to see if we are in a new chain and get its structure type
    # Then send both over OSC

    if newChainVal != oldChainVal:
        newChain = True
        client.send_message("/newChain", newChain)    
        
        # Retrieve category type and assign to category message
        categoryMsg = aaCategories[aminoAcid]
        client.send_message("/category", categoryMsg)
        
        # Get structure type of the chain (0 = helix, 1 = sheet, 2 = loop) and set struct message.
        ncv = int(newChainVal)
        if (ncv in helices):
            structTypeMsg = 0
        elif (ncv in sheets):
            structTypeMsg = 1
        else:
            structTypeMsg = 2
        client.send_message("/structType", structTypeMsg)
        
        # Set new chain to old chain to check for next new chain.
        oldChainVal = newChainVal
        newChain = False
    else:
        pass
    
    # Set Hydrophibicity messages
    print("Category: ", categoryMsg, " ", "Structure Type: ", structTypeMsg, " ", "Asym: ", 
            newAsymVal, " ", "Entity: ", newEntityVal, " ", "Chain: ", newChainVal, "\n")
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
    
    """
    # Debug sheets and helices lists
    print("Sheets: \n")
    for i in sheets:
        print(i)
    print("Helices: \n")
    for i in helices:
        print(i)
    """
  
    # TODO: Implement better quit code
    print("Press ctrl+c to quit.")
    while True:
        runOSC()
      
if __name__ == '__main__': main()
