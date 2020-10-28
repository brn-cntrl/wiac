"""
Script for writing PDB data to JSON for use in C#/Unity. 
NOTE: .cif or .pdb files must first be made plain text in order for the script to open them.

This version creates a JSON file in which you will find the data entries for amino acid
category, structure type (helix/sheet/loop), chain number, asym letter, and entity number.
NOTE: All entries are for side chains only. Backbone structures are already being handled 
in the visual design.

IMPORTANT: The C# script should handle timing of the sonification as well as the triggering 
of events based on boolean values (e.g. when a new chain number, asym letter, or entity 
number is dectected).

Brian Cantrell, Worldbuilding Media Lab. Oct, 2020.
"""

import os
import json

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
siteNumCol = 1 # Column number for the site number

# Pyhon lists for processing data
RGroupAAs = [] # List to hold R Group amino acids (used to index hydroVals)
RGroupChains = [] # List to hold R Group chain numbers (corresponding to amino acid)
RGroupBFactors = [] # List to hold R Group B factors
RGroupAsyms = [] # List to hold R Group asym values
RGroupEntities = [] # List to hold R Group entity values
categories = []
structureType = []
helicesList = [] # List to hold helices for structType
sheetsList = [] # List to hold sheets for structType
loopsList = [] # List to hold loops for structType
siteList = [] # List to hold site numbers

# JSON lists to hold the relevant data
jsonData = {}
"""
jsonData['categories'] = []
jsonData['AAs'] = [] # List to hold R Group amino acids (used to index hydroVals)
jsonData['chainNums'] = [] # List to hold R Group chain numbers (corresponding to amino acid)
jsonData['asyms'] = [] # List to hold R Group asym values
#jsonData['entities'] = [] # List to hold R Group entity values
#jsonData['helices'] = [] # List to hold helices for structType
#jsonData['sheets'] = [] # List to hold sheets for structType
jsonData['loops'] = [] # List to hold loops for structType
"""

# Sentry bool for struct type logic
isSheet = False

# String to hold the PDB file.
PDBFile = " " 

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
    global jsonData
    global beginHelix
    global endHelix
    global isSheet
    global RGroupAAs
    global RGroupChains 
    global RGroupBFactors 
    global RGroupAsyms 
    global RGroupEntities 
    global helicesList
    global sheetsList 
    global loopsList

    # Clear lists just in case...
    RGroupAAs.clear()
    RGroupAsyms.clear()
    RGroupEntities.clear()
    RGroupChains.clear()
    categories.clear()
    helicesList.clear()
    sheetsList.clear()
    loopsList.clear()
    
    # Open the .cif file NOTE: file must first be made plain text.
    if file is not None:
        # Split the file up by line
        theFile = file.readlines()
         # Loop over file 
        for line in theFile:
            # Split each line on white space
            entry = line.split()

            ###############################################################################################
            # Determine whether the current amino acid is part of a helix, loop, or sheet
            # Check for helix entry first and fill helices list
            if entry[0] == "HELX_P":
                firstHAA = int(entry[beginHelix])
                lastHAA = int(entry[endHelix])
                # Push the first and last AA number and every number between into helices list
                for i in range(firstHAA, lastHAA+1):
                    helicesList.append(i)
                    
            # Fill sheets list
            # First check for last header
            if entry[0] == "_struct_sheet_range.end_auth_seq_id":
                # Assign boolean to get sheet values until "#" is reached.
                # This is necessary because the "#" tag marks the beginning and end of section
                # TODO: Check to see if sheets are marked out by headers in the first column in other formats.
                isSheet = True
            
            if isSheet == True and entry[0] != "_struct_sheet_range.end_auth_seq_id" and entry[0] != "#":
                firstSAA = int(entry[beginSheet])
                lastSAA = int(entry[endSheet])
                for i in range(firstSAA, lastSAA+1):
                    sheetsList.append(i)
            
            # NOTE: "loops" list comprises any atom not in either sheets or helices. This will be handled 
            # when looping through the ATOM segment of the data and assigning JSON entries.

            # Check for "#" to set boolean back to false (probably unnecessary but extra security)
            if entry[0] == "#":
                isSheet = False

            ###############################################################################################
            # ATOMs
            # Select for only those lines listing ATOM coordinates
            if entry[0] == "ATOM":
                # Ignore the lines with backbone entries
                atomName = str(entry[atomCol])
                atomName.strip()

                if atomName == "CA" or atomName == "C" or atomName == "N" or atomName == "O" or atomName == "CB":
                    pass
                # Push R Group amino acid labels and B Factor values to lists
                else:
                    # Data values
                    RGroupAAs.append(entry[AACol])
                    RGroupAsyms.append(entry[asymCol])
                    RGroupEntities.append(entry[entityCol])
                    RGroupChains.append(entry[chainCol])
                    
                    # Get the site number to identify atom
                    siteList.append(entry[siteNumCol]) 
    
    else:
        print("file is none")
    
#################################################################################
def pushJSON():
    global jsonData
    jsonData["entries"] = []
    JSsiteNum = " "
    JSasym = " "
    JScat = " "
    JSstruct = " "
    JSentity = " "

    for i in range(len(siteList)):
        # Set values for JSON entries
        
        # Set site number 
        JSsiteNum = siteList[i]

        # Set the structure type
        curChain = int(RGroupChains[i])
        if curChain in helicesList:
            JSstruct = 0 # 0 for helix
        elif RGroupChains[i] in sheetsList:
            JSstruct = 1 # 1 for sheets
        else:
            JSstruct = 2 # 2 for loops

        # Set the category
        curAA = RGroupAAs[i]
        JScat = aaCategories[curAA]

        # Set the asym value
        JSasym = RGroupAsyms[i]

        # Set the entity
        JSentity = RGroupEntities[i]

        # Append data to JSON
        jsonData["entries"].append({
            "site": JSsiteNum,
            "asym": JSasym,
            "entity": JSentity,
            "structType": JSstruct,
            "category": JScat
        })
    
    # Open text file and dump json data to file
    with open('data.txt', 'w') as outputfile:
        json.dump(jsonData, outputfile)

#################################################################################      

# Main Loop
def main():
    pdbFile = input("\nPlease enter the path/name of the .cif or .pdb file and press 'Enter': \n\n")
    openFile(pdbFile)
    pushJSON()
    print("JSON file created.")
    
if __name__ == '__main__': main()
