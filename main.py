#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 14:02:19 2023

@author: sara
"""
from math import sqrt
import freesasa
import numpy as np
import shutil
import os
import subprocess
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.pose.carbohydrates import glycosylate_pose, glycosylate_pose_by_file, idealize_last_n_glycans_in_pose
from pyrosetta.rosetta.core.pose import pose_from_saccharide_sequence
from rosetta.protocols.analysis import *
from rosetta.core.select.residue_selector import *
from rosetta.protocols.carbohydrates import *
from rosetta.protocols.calc_taskop_movers import *
from rosetta.core.simple_metrics.metrics import *
from rosetta.core.simple_metrics.composite_metrics import *
from rosetta.core.simple_metrics.per_residue_metrics import *
from pyrosetta.rosetta.core.select.residue_selector import ResiduePropertySelector
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select import get_residues_from_subset

'''
BIOPYTHON citation: 
    Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. 
    Bioinformatics 19: 2308–2310
'''


# double = input('Would you like to include double mutation sites {Y/N}? (*Warning: depending on protein length, this can add significant runtime: ')
# chain = input('Please input which protein chain (if any) to use for prediction: ')
# size = int(input('Enter the number of residues to be within each final combination (integer): '))
# number = int(input('Enter the number of combinations to be produced (integer): '))
# priority_sites_in = input('Enter the residues to be glycosylated (separated by a comma): ')
# uncovered_sites_in = input('Enter the residues to be left uncovered (separated by a comma): ')
# model_glycans = input('Would you like to model glycans {Y/N}?: ')



# if len(priority_sites_in)>0:
#     priority_sites = priority_sites_in.split(',')
    
# else:
#     priority_sites = []


# if len(uncovered_sites_in)>0:
#     uncovered_sites = uncovered_sites_in.split(',')
# else:
#     uncovered_sites = []
    

    
    

#priority_sites = []


#check that priority_sites and uncovered_sites have no overlaps
#if set(priority_sites).intersection(set(uncovered_sites))!= set():
#    print('Error! Please re-input residues to be glycosylated and uncovered, and ensure no overlaps.')


def surface_residues(pdb, min_SASA=30):
    '''
    DESCRIPTION: this function takes a protein input and outputs a 
    list of all the surface residues if  their SASA > min_SASA, set at 
    default to 2.5 Å
    **based on: http://pymolwiki.org/index.php/FindSurfaceResidues
    and modified to use fresasa module to find SASA of residues
    
    ** freesasa is used for SASA calcuation and cited here:
        Simon Mitternacht (2016) FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Research 5:189. (doi: 10.12688/f1000research.7931.1)
        
    Parameters
    ----------
    pdb : FILE
        pdb file of protein.
    
    Returns
    -------
    surface_residues: [List] of surface residues in the protein.

    '''
    surface_residues = []
    
    prot = freesasa.Structure(pdb)
    
    #calculate SASA
    areas = freesasa.calc(prot)
   
    #from areas get residue Areas
    resAreas = areas.residueAreas()
    
    #loop through residues areas and check if SASA is above threshold
    #area dict stores SASA of surface residues:
    areaDict = {}
    for chain in resAreas:
        for residue, residueArea in resAreas[chain].items():
            if residueArea.total > min_SASA:
                resID = residue + ' ' + chain
                surface_residues.append(resID)
                areaDict[resID] = residueArea.total   
    return surface_residues, areaDict


def residue_info(surface_residues, pdb, destination):
    
    '''
    DESCRIPTION: this function generates dictionaries with information
    about a list of input residues, including a dictionary of the 
    type of each residue and the coordinates
    made to reduce the number of times a pdb file needs to be parsed
    
    Parameters
    ----------
    surface_residues: LIST
        list of strings: 'resnum chain'
        surface residues of interest
        
    pdb : FILE
        pdb file of protein.
    
    Returns
    -------
    surface_dict: dictionary of {'resnum chain': residue type}
    surface_coords: dictionary of {'resnum chain': coordinates}
    map_dict: nested dictionary of {chain: {resnum: index}}
    jwalk_map: dictionary that maps old pdb file numbering to new pdb file numbering 
                {chain: {original pdb num: jwalk pbd id}}
    
    '''
    
    
    parser = PDBParser()
    res_dict = {}
    res_coords = {}
    full_prot = parser.get_structure('protein', pdb)
    # need to go from structure --> model
    prot = full_prot[0]
    
    for residue in prot.get_residues():
        resID = residue.get_id()
        chain = residue.get_full_id()[2]
        #check that residue ID does not contain any letters (no insertion code)
        if resID[0] == ' ': #make sure that HET field is empty (amino acids only)
            if resID[2] == ' ': 
                resnum = resID[1]
                resname = str(resnum) + ' ' + chain
                restype = prot[chain][resnum].get_resname()
                CA_coords = prot[chain][resnum]['CA'].get_coord()
                #store res type
                res_dict[resname] = restype
                #store coordinates
                res_coords[resname] = list(CA_coords)
            else: #specially reference amino acids with letters in their name ("insertion code")
                resnum = resID[1]
                insertion = resID[2]
                resname = str(resnum) + insertion + ' ' + chain
                restype = prot[chain][(' ', resnum, insertion)].get_resname()
                CA_coords = prot[chain][(' ', resnum, insertion)]['CA'].get_coord()
                #store res type
                res_dict[resname] = restype
                #store coordinates
                res_coords[resname] = list(CA_coords)
    
    #also generate map of residue names to numerical sequence, by chain
    map_dict = {}
    for chain in prot.get_list():
        #get all residues of chain
        cur_resdict = {}
        indx = 1 #start index counting AT ONE to match NetNGlyc
        for res in chain.get_list():
            if res.get_id()[0] == ' ': #only use amino acids 
                res_key = str(res.get_id()[1])
                if res.get_id()[2] != ' ': #will only add insertion code if residue has it
                    res_key += res.get_id()[2]
                #if name in surface_residues: #only add to dict if surface residue **TOOK THIS OUT BC OF VALUE ERROR
                cur_resdict[res_key] = indx
            indx += 1
            map_dict[chain.get_id()] = cur_resdict
    
    #create new PDB file with corrected numbering for jwalk, and store map to revert back to original numbering
    n = PDBIO()
    n.set_structure(full_prot)
    
    jwalk_map = {}
    for chain in prot.get_list():
        cur_resdict = {}
        i = 1
        cur_chain = chain.get_id()
        for res in chain.get_list():
            original_name = '' #store original name of residue
            org_id = res.get_id()
            if org_id[0] != ' ':
                original_name+=org_id[0]
            original_name += str(org_id[1])
            if org_id[2] != ' ':
                original_name+=org_id[2]
            try:
                new_id = (org_id[0], i, ' ')
                res.id = new_id
            #in case residue with this index already exists, temporarily change
            #index of next residue   
            except ValueError: 
                next_res = chain[(org_id[0], i, ' ')]
                next_res.id = (org_id[0], i+10000+ord(org_id[0]), ' ') #set new ID to one that is not already occupied
                res.id = new_id
            cur_resdict[original_name] = new_id[1]
            i+=1
        jwalk_map[cur_chain] = cur_resdict
    n.save(destination+'jwalk_protein.pdb')     
    return res_dict, res_coords, map_dict, jwalk_map

  
#consecutive_residues based on two functions in original script + combined  
def consecutive_residues(resMap):
    """
    DESCRIPTION:
        This function finds three or more consecutive residues and 
        adds them into a single separate list together
	
	PARAMETERS:
        resMap,      dictionary of residues by chain
        
    RETURNS:
		[list]  of lists of consecutive numbers in the input residue numbers
                or empty list if no consecutive residue numbers found
    """
    final_list = []
    #iterate through each chain in resMap
    for chain in resMap:
        name_list = [key for key in resMap[chain].keys()]
        cur_chain = [] #keep track of consecutive residues in each chain
        for residue, indx in resMap[chain].items(): #loop through each residue in chain
        #check if both consective values are in the index values
            consec_list = []
            if (indx+1 in resMap[chain].values()) and (indx+2 in resMap[chain].values()):
                i  = name_list.index(residue) #get index of first residue in set of 3
                consec_list.append(residue + ' ' + chain)
                consec_list.append(name_list[i+1] + ' ' + chain)
                consec_list.append(name_list[i+2] + ' ' + chain)
            if len(consec_list) ==3:
                final_list.append(consec_list)
  
    return final_list 


#consensus_sequences is from original script
def consensus_sequences():
    """
    DESCRIPTION:
        This function creates possible combinations of N-glycosylation 
        consensus sequences that may be found in an amino acid sequence.

	RETURNS:
		[list] of possible consensus sequences for N-linked glycosylation
    """
    
    consensus_sequences = []
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    for aa in d.keys():
        for end in ['SER','THR']:
            sequence = ['ASN']
            sequence.append(aa)
            sequence.append(end)
            consensus_sequences.append(sequence)
    return consensus_sequences

def consensus_sequences_similarAA():
    """
    DESCRIPTION:
        This function creates possible combinations of N-glycosylation 
        consensus sequences that may be found in an amino acid sequence.

	RETURNS:
		[list] of possible consensus sequences for N-linked glycosylation
    """
    
    consensus_sequences = []
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    for aa in d.keys():
        for end in ['SER','THR','ASN','GLN','GLY','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']:
            for beginning in ['SER','THR','ASN','GLN','GLY','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']:
                sequence = []
                sequence.append(beginning)
                sequence.append(aa)
                sequence.append(end)
                consensus_sequences.append(sequence)
    return consensus_sequences


#identify_native_site from original script
def identify_native_site(target, query, res_Dict):
    """
    DESCRIPTION:
        This function identify native glycosylation sites in EpitopeCA.txt by 
        identifying whether the consensus sequence lies in its sequence.
	
	PARAMETERS:
        target,       list of consecutive residue numbers from EpitopeCA.txt
        query,        list of possible consensus sequence combinations
        res_Dict,     dictionary of {residue number: residue name} from CA file.
	
	RETURNS:
		nested [list] of residue numbers of native glycosylation sites.
    """
    native_sites = [] 
    for segment in target:
        if segment[0] in res_Dict.keys() and segment[1] in res_Dict.keys() and segment[2] in res_Dict.keys():
            resOne = res_Dict[segment[0]]
            resTwo = res_Dict[segment[1]]
            resThree = res_Dict[segment[2]]
            if [resOne, resTwo, resThree] in query:
                native_sites.append(segment)
    return native_sites


#id_single_mutation from original script
def id_single_mutation(target, query, res_Dict):
    """
    DESCRIPTION:
        This function identify novel glycosylation sites after single mutation
        in EpitopeCA.txt by identifying whether 2 of three amino acids in 
        the consensus sequence lies in its sequence.
	
	PARAMETERS:
        target,       list of consecutive residue numbers from EpitopeCA.txt
        query,        list of possible consensus sequence combinations
        res_Dict,     dictionary of {residue number: residue name} from CA file.
	
	RETURNS:
		nested [list] of residue numbers of glycosylation sites that can be singly mutated.
    """
    mutant_sites = [] 
    for segment in target:
        if segment[0] in res_Dict.keys() and segment[1] in res_Dict.keys() and segment[2] in res_Dict.keys():
            resOne = res_Dict[segment[0]]
            resTwo = res_Dict[segment[1]]
            resThree = res_Dict[segment[2]]
            list1 = [resOne, resTwo, resThree]
            if list1 not in query: #ensure only single mutation sites are added, do not double count native sites
                for subquery in query:
                    if subquery[:2]==list1[:2] or subquery[1:]==list1[1:]:
                        mutant_sites.append(segment)
                        break
    return mutant_sites


#id_single_mutation
def id_double_mutation(target, query, res_Dict):
    """
    DESCRIPTION:
        This function identify novel glycosylation sites after single mutation
        in EpitopeCA.txt by identifying whether 2 of three amino acids in 
        the consensus sequence lies in its sequence.
	
	PARAMETERS:
        target,       list of consecutive residue numbers from EpitopeCA.txt
        query,        list of possible consensus sequence combinations
        res_Dict,     dictionary of {residue number: residue name} from CA file.
	
	RETURNS:
		nested [list] of residue numbers of glycosylation sites that can be singly mutated.
    """
    mutant_sites = [] 
    for segment in target:
        resOne = res_Dict.get(segment[0])
        resTwo = res_Dict.get(segment[1])
        resThree = res_Dict.get(segment[2])
        list1 = [resOne, resTwo, resThree]
        if list1 not in query:
            for subquery in query:
                if subquery[0]==list1[0] or subquery[1]==list1[1] or subquery[2] == list1[2]:
                    mutant_sites.append(segment)
                    break
    return mutant_sites
    
def id_double_mutation_similarAA(target, query, query_nat, res_Dict):
    """
    DESCRIPTION:
        This function identify novel glycosylation sites after single mutation
        in EpitopeCA.txt by identifying whether 2 of three amino acids in 
        the consensus sequence lies in its sequence.
	
	PARAMETERS:
        target,       list of consecutive residue numbers from EpitopeCA.txt
        query,        list of possible consensus sequence combinations
        res_Dict,     dictionary of {residue number: residue name} from CA file.
	
	RETURNS:
		nested [list] of residue numbers of glycosylation sites that can be singly mutated.
    """
    mutant_sites = [] 
    for segment in target:
        resOne = res_Dict.get(segment[0])
        resTwo = res_Dict.get(segment[1])
        resThree = res_Dict.get(segment[2])
        list1 = [resOne, resTwo, resThree]
        if list1 not in query_nat:
            for subquery in query:
                if subquery[0]==list1[0] and subquery[1]==list1[1] and subquery[2] == list1[2]:
                    mutant_sites.append(segment)
                    break
    return mutant_sites


#dist from original script
def dist(a, b):
    """
    DESCRIPTION:
        This function calculates distance between two points (a, b)
	
	PARAMETERS:
        a,      point 1 (x1, y1, z1)
        b,      point 2 (x2, y2, z2)	
        
	RETURNS:
        Distance between two points as float.
    """
    
    d = [float(a[0]) - float(b[0]), float(a[1]) - float(b[1]), float(a[2]) - float(b[2])]
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) 


#distance_n_atoms from original script
def distance_n_atoms(cords):
    """
    DESCRIPTION:
        This function calculates distance between a set of points. 
	
	PARAMETERS:
        cords,      list of alpha carbon coordinates from possible glycosylation sites
        
	RETURNS:
        Coord_matrix,   nested [dictionary] of distances between coordinates
    """

    
    D = {}
    
    for n_atom1, cords1 in cords.items():
        D[n_atom1] = {}
        for n_atom2, cords2 in cords.items():
            D[n_atom1][n_atom2] = dist(cords1, cords2) 
    
    return D



def allSASD(reslist, jwalk_pdbfile, Jwalk_path, jwalk_map):
    
    '''
    This program uses Jwalk, cited here:
        
    Sinnott et al., Combining Information from Crosslinks and Monolinks in the Modeling
    of Protein Structures, Sturcture (2020), https://doi.org/10/1016/j.str.2020.05.012
    
    The Importance of Non-accessible Crosslinks and Solvent Accessible Surface Distance in Modeling Proteins with Restraints From Crosslinking Mass Spectrometry.
    J Bullock, J Schwab, K Thalassinos, M Topf. Mol Cell Proteomics. 15, 2491–2500, 2016


    Purpose: This function uses Jwalk to find the SAS distances between a list
    of residues 
    
    Parameters
    ----------
    reslist: list of surface residues to calculate SASD for
    
    jwalk_pdbfile : FILE
        name of pdb file which contains the residues in JWALK renumbered format.
    
    Jwalk_path: STRING
        string of the path location to Jwalk download
    
    jwalk_map: dictionary that maps new pdb file numbering to original pdb file numbering 
                {chain: {jwalk pdb num: original num}}
                
    Returns
    -----------
    SASD_dict: dict of SASD for each residue pair in residue list
            key: (res1, res2) - strings
            value: SASD - float 
        
    '''

    aa_file = Jwalk_path + '/aa_inputs.txt'
    
    #first, translate list of residues to Jwalk pdb file numbering scheme
    newlist = []
    #store backwards dict {jwalk res: org res}
    backwards = {}
    for res in reslist:
        resnum = res.split()[0]
        chain = res.split()[1]
        newres = str(jwalk_map[chain][resnum]) + ' ' + chain
        newlist.append(newres)
        backwards[newres] = res
    #fix reslist so that numbering is jwalk compatible
    reslist = newlist

    visited = {} #keep track of visited pairs to avoid calculating dist twice
    #reformat input information for Jwalk as txt file
    with open(aa_file, 'w') as inputfile: 
         line_indx = 1    #keep track of which line residue pair will be calculated at
         for res1 in reslist:
             #create second loop to calculate distances
             for res2 in reslist:
                 if ((backwards[res1], backwards[res2]) not in visited.values()) and ((backwards[res2], backwards[res1]) not in visited.values()) and (res1!=res2):
                     #extract residue number and chain from inputs
                     res_num1 = res1.split()[0]
                     chain1 = res1.split()[1]
                
                     res_num2 = res2.split()[0]
                     chain2 = res2.split()[1]
                    
                     inputfile.write(res_num1 + '|' + chain1 + '|' + res_num2 + '|' +chain2 + '|' + '\n')
                     #mark that this distance pair has been visited
                     visited[line_indx] = (backwards[res1],backwards[res2])
                     line_indx += 1
    
    #copy pdb file into cd of Jwalk
    jwalkpdb = Jwalk_path + '/pdbfile.pdb'
    shutil.copyfile(jwalk_pdbfile, jwalkpdb)

    #access terminal to run Jwalk
    os.chdir(Jwalk_path)
    os.system('python Jwalk.v2.1.py -i pdbfile.pdb -xl_list aa_inputs.txt -max_dist 50')           
    
    #extract distance information from outputted txt file
    dist_output = Jwalk_path + '/Jwalk_results/pdbfile_crosslink_list.txt'
    
    #store distances for pairs of residues
    SASD_dict = {}
            
    with open(dist_output, 'r') as distance_info:
        #check that there is an SASD:
        all_lines = distance_info.readlines()
        distance_lines = all_lines[1:]
        for line in distance_lines:
            #check if residue information at current line matches with the input  
            #(means SASD was calculated)
            
            line_data = line.split()
            cur_res1 = line_data[2].split('-')[1]+ ' ' + line_data[2].split('-')[2]
            cur_res2 = line_data[3].split('-')[1]+ ' ' + line_data[3].split('-')[2]
            SASD= line.split()[4]
            dist = float(SASD)
            #translate each res from original to Jwalk numbering
            new_res1 = backwards[cur_res1]
            new_res2 = backwards[cur_res2]
            SASD_dict[(new_res1, new_res2)] = dist

    #default dist to 100 for pairs where no SASD was found
    for res_pair in visited.values():
        if (res_pair[1], res_pair[0]) not in SASD_dict and (res_pair[0], res_pair[1]) not in SASD_dict:
            SASD_dict.setdefault(res_pair, 100.0)
    return SASD_dict


def calc_overlap(sphere_rad, distances, priority_sites):
    """
    DESCRIPTION:
        This function identifies and creates clusters 
	
	PARAMETERS:
        sphere_rad,     input radius of glycan size
        distances,      distance dictionary between coordinate pairs
        priority_sites   list of residues that must be in each combo
	RETURNS:
        all_clusters, Dictionary of {site number: non overlapping sites}
    """
    print('START', flush = True)
    all_sites = list(distances.keys())
    #create nested list starting with first site
    combinations = [[all_sites[0]]]
    #use visited to keep track of which combinations are already included 
    visited = set(map(frozenset, combinations))
    i = 0
    while i < len(combinations):
        cur_combination = combinations[i]
        #only check against residues that aren't already in current combination
        for site1 in all_sites:
            if site1 not in cur_combination:
                distance_to_cur_combination = {site2: distances[site1][site2] for site2 in cur_combination}
                if any(distance <= sphere_rad for distance in distance_to_cur_combination.values()):
                    clash_set = {residue for residue in cur_combination if distances[site1][residue] <= sphere_rad}
                    #create new list that does not include clashing residues, plus add in new site
                    new_cluster = [res for res in cur_combination if res not in clash_set] + [site1]
                    #append new combination to ongoing list of combinations
                    #first check that cluster is not already in list of combinations, or that longer version of cluster is not in combinations
                    if frozenset(new_cluster) not in visited:
                        combinations.append(new_cluster)
                        visited.add(frozenset(new_cluster))
            # no residues in cluster cause a clash
                else:
                    cur_combination.append(site1)
                    #make sure this combination is not already in the queue
                    if frozenset(cur_combination) in visited:
                        cur_combination.pop()
                    else:
                        visited.add(frozenset(cur_combination))
        i = i+1
        if i%1000==0:
            print(i, flush=True)
            print(len(combinations))

    #create set of priority residues to eliminate clusters that do not contain them
    priority_set = set(priority_sites)
    if priority_sites == [' K']:
        return combinations
    else:
        all_combos = {frozenset(c) for c in combinations if frozenset(c).intersection(priority_set) == priority_set}
        #print('Difference in lengths?: ')
        #print(combinations)
        #print(len(combinations), flush=True)
        #for c in combinations:
        #    if frozenset(c) not in all_combos and frozenset(c).intersection(priority_set)==priority_set:
        #        all_combos.append(set(c))
        final_combos = [list(combo) for combo in all_combos]
        return final_combos

def reduce_clusters(sphere_rad, distances, combinations, priority_sites):
    """
    DESCRIPTION:
        This function reduces size of clusters until they have no remaining clashes under new clash radius
	
	PARAMETERS:
        sphere_rad,     input radius of glycan size
        distances,      distance dictionary between coordinate pairs
        combinations,   nested list of combinations for clusters created from previous overlap radius
        priority_sites   list of residues that must be in each combo
        
	RETURNS:
        new_clusters, list of new combinations that have no clashes based on starting combination and sphere_rad
    """
    #sort to ensure consistent results:
    for combo in combinations:
        combo.sort()
    combinations.sort(key=len, reverse=True)

    final_clusters = []
    visited = [] #use visited to keep track of every combination that has been iterated through
    for combination in combinations:
        if (len(combinations))%1000==0:
            print(len(combinations))
        visited.append(set(combination))
        clashes=False
        for res1 in combination:
            for res2 in combination:
                if res1!=res2 and distances[res1][res2]<=sphere_rad:
                    clashes = True
                    #remove clashing residues and append to list of combinations
                    new_combo1 = [res for res in combination if res!=res1]
                    new_combo2 = [res for res in combination if res!=res2]
                    if set(new_combo1) not in visited:
                        combinations.append(new_combo1)
                        visited.append(set(new_combo1))
                    if set(new_combo2) not in visited:
                        combinations.append(new_combo2)
                        visited.append(set(new_combo2))
                    break
            if clashes:
                break
        if not clashes:
            if set(combination) not in final_clusters:
                final_clusters.append(set(combination))
    #create priority set to make sure all priority residues are in clusters
    if priority_sites == [' K']:
        return [list(cluster) for cluster in final_clusters]
    else:
        priority_set = set(priority_sites)
        return [list(cluster) for cluster in final_clusters if cluster.intersection(priority_set)==priority_set]


def write_fasta(pdb, destination):
    '''
    DESCRIPTION: take pdb file and write FASTA sequence manually to a txt file
    (ensures that all residues that appear in PDB match those submitted to
     NetNGlyc)
    
    Parameters:
        pdb: pdb file of protein
        destination: path where user wants to send txt file of FASTA seq, should 
        end in "/"
    
    Returns:
        None
        
    '''
    #residue dictionary to convert 3 letter code to one letter code
    AACodeDict = {'ALA': 'A', 'ARG': 'R', 'ASN':'N', 'ASP':'D', 'CYS': 'C', 
                    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS':'H', 'HYP': 'O', 
                    'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 
                    'PRO': 'P', 'GLP': 'U', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                    'TYR': 'Y', 'VAL': 'V'}
    #create txt file for output
    name = pdb.split('.')[0]
    prot = name.split('/')[-1]
    filename = destination + prot + '_fasta.txt'
    fasta = open(filename, 'w')
    
    #open pdb file and loop through by chain
    parser = PDBParser()
    #protein model from pdb file
    prot = parser.get_structure('prot', pdb)[0]
    for chain in prot.get_list():
        #check that chain has amino acids; else skip
        check = False
        for n in chain.get_list():
            name = n.get_resname()
            if name in AACodeDict.keys():
                check = True
        if check:
            fasta.write('>Chain' + '_' + chain.get_id() +'\n')
            for residue in prot[chain.get_id()].get_list():
                if residue.get_resname() in AACodeDict:
                    letter_code = AACodeDict[residue.get_resname()]
                    fasta.write(letter_code) #add AA code to fasta file
        fasta.write('\n') #new line for next chain
    fasta.close()
    return None


def mutate_fasta(res_to_mutate, allRes, map_dict, mutant_sites, original_fasta, destination):
    '''
    DESCRIPTION: take pdb file and write FASTA sequence manually to a txt file
    with residue to mutate changed to appropriate residue to ensure consensus sequence
    
    Parameters:
        res_to_mutate: string ('number chain' e.g. '20 A') of residue to mutate (based on PDB numbering NOT mapped integer numbering)
            just provides first position residue in set of 3, may not actually be the residue that will
            change (could be third residue)
        allRes: dict of surface residues (number and chain) and corresponding amino acid values
        map_dict : nested dict
            nested dictionary of {chain: {resnum: index}} 
        destination: path where user wants to send txt file of FASTA seq, should end in '/'
        mutant_sites: list of list of sites that can be consensus sequence with mutation
        original_fasta: txt file of fasta seq for original protein sequence
    Returns:
        None
        
    '''
    #residue dictionary to convert 3 letter code to one letter code
    AACodeDict = {'ALA': 'A', 'ARG': 'R', 'ASN':'N', 'ASP':'D', 'CYS': 'C', 
                    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS':'H', 'HYP': 'O', 
                    'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 
                    'PRO': 'P', 'GLP': 'U', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                    'TYR': 'Y', 'VAL': 'V'}
    
    resnum = res_to_mutate.split()[0] #get number of residue to mutate
    reschain = res_to_mutate.split()[1] #get chain of residue to mutate
    mutation = ''
    #get sequence of 3 residues which mutation will occur in
    seq = [site for site in mutant_sites if site[0] == res_to_mutate] 
    #determine residue to change current residue to for consensus sequence
    if allRes[res_to_mutate] == 'ASN': #if first residue is already N, we need to mutate last residue to 'T'
        mutation = 'T'
        to_change = [seq[0][2]] #store residue name of AA actually being mutated
    #determine if second residue is S or T, in this case change first res to N
    elif allRes[seq[0][2]] == 'THR' or allRes[seq[0][2]] == 'SER':
        mutation = 'N'
        to_change = [res_to_mutate] #residue that will change is same as input
    #if first residue and second residue both need to change
    else:
        mutation1 = 'N'
        mutation2 = 'T'
        to_change = [res_to_mutate, seq[0][2]]
    #only one change needs to be made
    if len(to_change) == 1:
        #get index of residue to change from map dict
        chain = to_change[0].split()[1]
        num = to_change[0].split()[0]
        indx = map_dict[chain][num]
        
        
        #create txt file for output
        filename = destination + 'fasta_' + resnum + '_' + reschain + '.txt'
        newfile = open(filename, 'w')
        og_file = open(original_fasta, 'r')
        newfile.write('>Chain' + ' ' + chain +'\n')
        #keep track of if reached chain of mutation yet
        j = False
        for line in og_file.readlines():
            check =  line.split('_')
            if len(check)>1:
                if check[1].split('\n')[0] == chain:
                    j = True
            elif j:
                #starting at index 1, copy all letters from seq into new fasta, except 
                #letter at index of mutation
                i=1
                for c in line:
                    if i!= indx:
                        newfile.write(c)
                    else:
                        newfile.write(mutation)
                    i=i+1
                break #break so loop does not repeat once fasta is written for chain of interest
    
        newfile.close()
        og_file.close()

    #need to make mutation at both sites
    elif len(to_change) == 2:
        #get index of residue to change from map dict
        to_change1 = to_change[0]
        to_change2 = to_change[1]
        chain = to_change1.split()[1] #same chain for both res, use first one
        num1 = to_change1.split()[0]
        indx1 = map_dict[chain][num1]
        num2 = to_change2.split()[0]
        indx2 = map_dict[chain][num2]
        
        #create txt file for output
        filename = destination + 'fasta_' + resnum + '_' + reschain + '.txt'
        newfile = open(filename, 'w')
        og_file = open(original_fasta, 'r')
        newfile.write('>Chain' + ' ' + chain +'\n')
        #keep track of if reached chain of mutation yet
        j = False
        for line in og_file.readlines():
            check =  line.split('_')
            if len(check)>1:
                if check[1].split('\n')[0] == chain:
                    j = True
            elif j:
                #starting at index 1, copy all letters from seq into new fasta, except 
                #letter at index of mutation
                i=1
                for c in line:
                    if i!= indx1 and i!= indx2:
                        newfile.write(c)
                    else:
                        if i == indx1:
                            newfile.write(mutation1)
                        elif i == indx2:
                            newfile.write(mutation2)
                    
                    i=i+1
                break #break so loop does not repeat once fasta is written for chain of interest
    
        newfile.close()
        og_file.close()
        
    return None
  
  
def run_netNglyc_chain(fasta, map_dict, netnglyc_loc, threshold = 6):
    '''
    

    Parameters
    ----------
    fasta : file
        input file of protein fasta sequence.
    map_dict : nested dict
        nested dictionary of {chain: {resnum: index}}
    threshold : integer, optional
        number of neural nets that must agree for site to count as 'likely' glycosylated. The default is 6.

    Returns
    -------
    predicted_sites: list
        list of strings ('resnum chain') predicted by netNglyc to be glycosylated.

    '''
    #obtain output by running fasta file in NetNGlyc
    #file path assumes starting in HyperImmunISE folder, append name of fasta file to input
    command = 'cd '+netnglyc_loc+' && ./netNglyc ' + fasta 
    output = subprocess.run(command, shell=True, stdout=subprocess.PIPE) #CHANGE LATER Shell = TRUE

    results = output.stdout.decode().split('\n')
    #obtain name of file sent to NetNGlyc to identify lines with scoring output
    seqname = fasta.split('.')[0]
    print(seqname)
    
    if len(seqname.split('_'))>2:
        #get current chain from filename
        chain = seqname.split("_")[2]
        print(chain)
        #create reverse dictionary to map numerical index --> pdb resnum for translating netNglyc results
        rev_dict = {value: key for (key, value) in map_dict[chain].items()}
        #use list to store results of netNglyc scoring for each residue
        predicted_sites = []
        #loop through each line of output
    
        for i in results:
            terms=i.split()
            if len(terms) > 0 and terms[0].split('_')[0] =='Chain':
                #terms[0] is seqname
                #terms[1] is amino acid number of N in consensus sequence
                #terms[2] is sequence 
                #terms[3] is potential
                #terms[4] is neural net agreement 
                cur_index = int(terms[1])
                potential = float(terms[3])
                score = int((terms[4].split('/')[0])[1]) #first integer of score /9
                #check that score is greater than threshold, default = 6
                if (score>= threshold) and (potential >0.5):
                    resid = str(rev_dict[cur_index]) + ' ' + chain #residue number and chain
                    predicted_sites.append(resid)
        
    return predicted_sites

           
def run_netNglyc_all(fasta, map_dict, netnglyc_loc, threshold = 6):
    '''
    r

    Parameters
    ----------
    fasta : file
        input file of protein fasta sequence.
    map_dict : nested dict
        nested dictionary of {chain: {resnum: index}}
    threshold : integer, optional
        number of neural nets that must agree for site to count as 'likely' glycosylated. The default is 6.

    Returns
    -------
    predicted_sites: list
        list of strings ('resnum chain') predicted by netNglyc to be glycosylated.

    '''
    #obtain output by running fasta file in NetNGlyc
    #file path assumes starting in HyperImmunISE folder, append name of fasta file to input
    command = 'cd '+netnglyc_loc+' && ./netNglyc ' + fasta 
    output = subprocess.run(command, shell=True, stdout=subprocess.PIPE) #CHANGE LATER Shell = TRUE

    results = output.stdout.decode().split('\n')
    #obtain name of file sent to NetNGlyc to identify lines with scoring output
    predicted_sites = []
    
    #loop through results
    cur_chain = ''
    rev_dict = {}
    for i in results:
        terms=i.split()
        #terms[0] is seqname
        #terms[1] is amino acid number of N in consensus sequence
        #terms[2] is sequence 
        #terms[3] is potential
        #terms[4] is neural net agreement
        
        if len(terms) > 0 and terms[0].split('_')[0] =='Chain':
            #check if on same chain or still on current chain
            #same chain as previously
            if terms[0].split('_')[1] == cur_chain:
                 
                cur_index = int(terms[1])
                potential = float(terms[3])
                score = int((terms[4].split('/')[0])[1]) #first integer of score /9
                #check that score is greater than threshold, default = 6
                if (score>= threshold) and (potential>0.5):
                    resid = str(rev_dict[cur_index]) + ' ' + cur_chain #residue number and chain
                    predicted_sites.append(resid)

           
            #new chain
            else: 
                cur_chain = terms[0].split('_')[1]
                #create reverse dictionary to map numerical index --> pdb resnum for translating netNglyc results
                rev_dict = {value: key for (key, value) in map_dict[cur_chain].items()}
                cur_index = int(terms[1])
                potential = float(terms[3])
                score = int((terms[4].split('/')[0])[1]) #first integer of score /9
                #check that score is greater than threshold, default = 6
                if (score>= threshold) and (potential >0.5):
                    resid = str(rev_dict[cur_index]) + ' ' + cur_chain #residue number and chain
                    predicted_sites.append(resid)
        
    return predicted_sites


def surface_quant(pdb, glycopath, residues, all_surface_res,destination):
    """
    Parameters
    ----------
    res_dict : dictionary of {'resnum chain': residue type} 
       
    pdb : string (FULL file path)
        path to pdb file of protein
    glycopath : string
        path to glycopath download on computer
    residues : list
        list of residues within surface_dict to examine
    all_surface_res: list
        list of residues determined by freeSASA to be surface residues, use to calculate coverage parameter

    Returns
    -------
    coverage : INTEGER
        integer of sum of number of glycan atoms for surface residues
       
        ** This program uses GLYCO, cited here:
            Lee M, Reveiz M, Rawi R, Kwong PD, Chuang GY. GLYCO: 
                a tool to quantify glycan shielding of glycosylated proteins. 
                Bioinformatics. 2021 Nov 23;38(4):1152–4. 
                doi: 10.1093/bioinformatics/btab791. Epub ahead of print. 
                PMID: 34864901; PMCID: PMC8796370.
    """
    #write surface residues into text file readable by GLYCO
    os.chdir(glycopath)
    #rerun residue info to get correct residue type for each residue number, since mutated for each protein
    res_dict, sc, md, jd = residue_info(all_surface_res, pdb,destination)

    with open('surfaces.txt', 'w') as outfile:
        for cur_res in residues:
            temp = cur_res.split()
            outfile.write(res_dict[cur_res] + ' ' + temp[1] + ' ' + temp[0] + "\n")
    
    
    #copy pdb file into cd of glyco
    
    glycopdb = glycopath + '/pdbfile.pdb'
    #only copy pdb file if it isn't already in path
    if not os.path.exists('glycopdb'):
        shutil.copyfile(pdb, glycopdb)
    
    
    #check if name of results folder is taken, if so iterate up an integer
    os.system('python glyco.py -pdb pdbfile.pdb -cutoff 20 -module sub -glycan BMA,AMA,BGL,NAG,MAN -residue surfaces.txt -out_folder res')
    
    #read output file to obtain result
    ans_file = glycopath + '/res/pdbfile_sub_glysum.txt'
    with open(ans_file, 'r') as sumfile:
        int_string = sumfile.readlines()[0]
        coverage = float(int_string)
    #delete file and results folder to clear for next run (if function will be run again)
    os.remove(ans_file)
    os.rmdir(glycopath + '/res')

    return coverage


def coverage_rank(distance_n_atoms, clusters, i, radius = 17.5):
    '''
    PURPOSE: refine possible cluster combinations based on predicted
    coverage -- only high coverage combinations remain 

    Parameters
    ----------
    distance_n_atoms : dictionary of {site: {residue: distance}} using linear
    distance calculation (all surface residues)
       
    clusters : dict
        Dictionary of {site number: non overlapping sites}
    
    radius: int
        distance cutoff of nearby residues considered to be covered by glycan
    
    i: int
        final number of combinations needed
    
    Returns
    -------
    refined_clusters: dict
    Dictionary of {site number: non overlapping sites} (i.e. combos/clusters)
    with lower-coverage combinations eliminated
    '''
    #store max rank to take top percentile of coverage combinations

    refined_clusters = []
      
    ranking_dict = {}
    list_dict = {} #map string combos to list combos
    
    n_residues = len(distance_n_atoms.keys()) #store total number of surface residues
    #calculate approximate glycan coverage for each cluster
    for cluster in clusters:
        covered_residues = []
        #check coverage of each residue in the combination
        print(cluster)
        for res in cluster:
            for res2, dist in distance_n_atoms[res].items():
                if dist < radius and res2 not in covered_residues:
                    covered_residues.append(res2)
        rank = len(covered_residues)/n_residues
        cluster.sort() #sort for consistent results
        ranking_dict[str(cluster)] = rank
        list_dict[str(cluster)] = cluster

    #sort dictionary by value (ranking, greatest to least)
    ordered_dict = sorted(ranking_dict.items(), key = lambda x:x[1], reverse=True)

    #add clusters to final list while index is below required number of clusters
    indx = 0
    for tup in ordered_dict:
        if int(indx) < int(i):
            refined_clusters.append(tup[0])
        else: 
            break
        indx+=1
        
    final_clusters = []
    for clust in refined_clusters:
        final_clusters.append(list_dict[clust])
    print()
    return final_clusters
   
 
def iterate_clusters(distance_n_atoms, clusters, siteDistances, priority_sites, i, j, initial_r):
    '''
    PURPOSE: iteratively increase radius of overlap to reduce size of combinations in clusters
    and check coverage to reduce number of combinations

    Parameters
    ----------
    distance_n_atoms : dictionary of {site: {residue: distance}} using linear
    distance calculation (all surface residues)
    
    priority_sites : list
        list of N residue of user designated priority sites that must be glycosylated
        
    clusters : dict
        Dictionary of {site number: non overlapping sites}
    
    
    i: int
        final number of clusters (j different combinations)
    j: int
        final size of each cluster (each cluster constitutes i or less residues)
        
    initial_r: int
        first radius of overlap to determine overlapping residues in clusters
    
    
    siteDistances: dict
      distance dictionary between coordinate pairs
    
    Returns
    next_clusters: dict
        Dictionary of {site number: non overlapping sites} with j site numbers, combinations of size i
   
    '''

    #check if another iteration is needed by checking i and j
    print(i)
    print(j)
    cur_i = len(clusters)
    cur_j = 0
    print(cur_i)
    #obtain maximum length of clusters
    for val in clusters:
        if len(val) > cur_j: #check if lenth of cluster is greater than current max length
            cur_j = len(val)
            
    print(cur_j)
    #case 1: clusters too large AND too many clusters
    if  int(cur_i) > int(i) and int(cur_j) > int(j):
        print('i and j too large',flush=True)
        print('num clusters: ' + str(cur_i) + ' len clust: ' + str(cur_j))
        #increase overlap radius
        new_r = initial_r+1
        print('new radius: ' + str(new_r))
        new_clusters_long = reduce_clusters(new_r, siteDistances, clusters, priority_sites)
        #reduce number of clusters
        new_clusters = refine_clusters(new_clusters_long)
        print('new clusters, refined clusters')
        print(len(new_clusters_long), len(new_clusters))

        next_clusters = iterate_clusters(distance_n_atoms, new_clusters,
                                        siteDistances, priority_sites, i, j, new_r)
        
    elif int(cur_j) > int(j): 
        #obtain new overlaps with increased overlap radius size
        print('too many residues in cluster:',flush=True)
        print('num clusters: ' + str(cur_i) + ' len clust: ' + str(cur_j))
        new_r = initial_r+1
        print('new radius: ' + str(new_r))
        new_clusters = reduce_clusters(new_r, siteDistances, clusters, priority_sites)
        next_clusters = iterate_clusters(distance_n_atoms, new_clusters,
                                        siteDistances, priority_sites, i, j, new_r)
        
    #case 3: clusters correct size, too many clusters:
    elif int(cur_i) > int(i):
        print('too many clusters, filtering by coverage: ',flush=True)
        print('num clusters: ' + str(cur_i) + ' len clust: ' + str(cur_j))
        next_clusters = coverage_rank(distance_n_atoms, clusters, i)
        
    else:
        next_clusters = clusters
    
    return next_clusters


def refine_clusters(cluster_list, max_len=0):
    '''
    

    Parameters
    ----------
    cluster_list : nested LIST
        Nested list of all combinations
    max_len : integer, optional
        Minimum length of a combination. The default is 0.

    Returns
    -------
    new_list : LIST
        Nested list of combinations with minimum length of max_len.

    PURPOSE: reduce number of combinations in cluster_list so that only clusters of 
    longer length are included
    
    '''
    new_list = []
    #if max_len at default value, adjust to actual max length
    print(cluster_list)
    if max_len ==0:
        max_len = max([len(cluster) for cluster in cluster_list])
    #store each residue that was included in original clusters
    initial_res = set(np.concatenate(cluster_list))
    for cluster in cluster_list:
        if len(cluster)<max_len:
            pass
        else:
            new_list.append(cluster)
    final_residues= set(np.concatenate(new_list))
    #check that every residue in final residues and initial residues match, else do recursion
    if initial_res!=final_residues and max_len>1:
        new_list=refine_clusters(cluster_list, max_len-1)
    
    return new_list

def add_glycans(pdb, combo_list, glycan, native_sites, destination, model_glycans):
    '''
    

    Parameters
    ----------
    pdb : file
        pdb file of protein.
    combo_list : nested list
        list of list of site combinations for adding glycans.
    glycan: string
        indicates glycan to be added
    native_sites:
        list of native glycosylation sites identified previously
    model_glycans:
        STRING 'Y' or 'N'
        only if 'Y' will have rosetta do glycan sampling of conformations
        
    Returns
    -------
    updated_combos list
        list of lists of combos to account for any residues removed due to disulfide bonds.
    
    Purpose:
        generates glycosylated pdb files for each combination
        
    '''
    # code is largely from : https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/13.02-Glycan-Modeling-and-Design.ipynb
    # and https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/13.01-Glycan-Trees-Selectors-and-Movers.ipynb
    
    options = '''
    -edensity::score_symm_complex false
    -cryst::crystal_refine
    -include_sugars
    -auto_detect_glycan_connections
    -alternate_3_letter_codes pdb_sugar
    -maintain_links
    -write_pdb_link_records
    -write_glycan_pdb_codes
    -ignore_unrecognized_res
    -ignore_zero_occupancy false
    -load_PDB_components false
    -pdb_comments
    -scorefile_format json
    -skip_connect_info
    -use_input_sc
    -jd2:delete_old_poses
    -other_pose_to_scorefile
    '''
    
    init(" ".join(options.split('\n')))
    
    first_n_sites = [site[0] for site in native_sites]
    
    base_protein = pose_from_pdb(pdb)
    
    #check for any disulfide bonded sites, skip these:
    p = ResiduePropertySelector()
    prop = ResidueProperty(63)
    p.set_property(prop)
    p.apply(base_protein)
    x = list(get_residues_from_subset(p.apply(base_protein))) #residues to SKIP and remove from final combo

    digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'] #for checking if insertion code present
    #loop through each combination of sites
    updated_combos = [] #store any updated combos
    for combo in combo_list:
        temp_pose = base_protein.clone()
        #loop through each site in combo
        rosetta_combo = [] #store combo with rosetta numbering
        #store list of residues to remove
        res_to_remove = []
        for site in combo:
            chain = site.split()[1]
            num = site.split()[0]
            #check if insertion code is present!
            if num[-1] not in digits: #insertion code is present
                pose_num = base_protein.pdb_info().pdb2pose(chain, int(num[0:-1]), str(num)[-1])
            else: #insertion code not present
                pose_num = base_protein.pdb_info().pdb2pose(chain, int(num))
            #check if disulfide bond -- if so, move on to next site
            if pose_num in x:
                res_to_remove.append(site)
                continue
            
            #check if native site; if not, mutate sequence
            if site not in first_n_sites:
                n_res = ResidueIndexSelector()
                n_res.set_index(pose_num)
                #check type of third res to mutate to either S or T
                sequon = CreateSequenceMotifMover(n_res)
                third_res = base_protein.residue(pose_num+2)
                if third_res.name()=='SER':
                    sequon.set_motif("N[-]S")
                else: #if third res is THR or anything else
                    sequon.set_motif("N[-]T")
                sequon.apply(temp_pose)
            rosetta_combo.append(pose_num)
            
        updated_combo = [s for s in combo if s not in res_to_remove]
        updated_combos.append(updated_combo)
        
        #add glycans
        glycans = SimpleGlycosylateMover()
        glycans.set_glycosylation(glycan)
        #loop through sites again, add glycans to each site
        for site in rosetta_combo:
            n_res = ResidueIndexSelector()
            n_res.set_index(site)
            glycans.set_residue_selector(n_res)
            glycans.apply(temp_pose)
            
        #glycan sampling IF user input
        if model_glycans == 'Y':
            GTM = GlycanTreeModeler()
            GTM.apply(temp_pose)
            
        #output new protein pdb file
        pdb_name = destination + 'Output_' + str(combo_list.index(combo))+'.pdb'
        temp_pose.dump_pdb(pdb_name)

    return updated_combos



# RUNNING THE SCRIPT
def process(netnglyc_loc, jwalk_loc, chain, size, number, priority_sites, uncovered_sites, model_glycans,pdb_list,path,destination):
    for i in range(len(pdb_list)):
        netnglyc_location=netnglyc_loc[i]
        jwalk_location=jwalk_loc[i]
        download_directory=path[i]
        prot_name = pdb_list[i].split('.')[0] #name of protein of interest
        residues, areaDict = surface_residues(download_directory)
        #print(residues, flush=True)
        #priority_sites[i]=[priority_sites[i]]
        #uncovered_sites[i]=[uncovered_sites[i]]
        allRes, allCoords, resMap, jwalk_map = residue_info(residues, download_directory, destination[i])
        # SEARCHING SURFACE RESIDUES FOR GLYCOSYLATION SEQUENCE
        query=consensus_sequences()
        query_similarAA=consensus_sequences_similarAA()
        targets = consecutive_residues(resMap)
        #reduce targets to only sequons starting with a surface residue and not having a LYS residue
        s_targets = [target for target in targets if target[0] in residues]
        surface_sites = identify_native_site(s_targets, query, allRes)

        
        # IDENTIFY POSSIBLE MUTATION SITES 
        single_mutant_sites = id_single_mutation(s_targets, query, allRes)

        #add double mutant sites into consideration if double indicated:
        double_mutant_sites = id_double_mutation_similarAA(s_targets, query_similarAA, query, allRes)
        for site in double_mutant_sites:
            if site not in single_mutant_sites:
                single_mutant_sites.append(site)
           
        mutant_sites=single_mutant_sites
        #linear distance matrix
        #use for calculating coverage
        sCoord = {k:v for (k,v) in allCoords.items() if k.split()[1]==chain[i]}
        d_linear = distance_n_atoms(sCoord)
        #also create a linear distance matrix only for surface residues (or, user designated glycosylation site):
        sCoord_surf = {k:v for (k,v) in allCoords.items() if k.split()[1]==chain[i] and (k in residues or k in priority_sites[i])}
        d_linear_surf = distance_n_atoms(sCoord_surf)

        
        #COMMENT BEGIN

        #NETNGLYC ANALYSIS
        #first, write original fasta file with no mutations; obtain netNglyc output
        write_fasta(download_directory,destination[i])
        fasta_file = destination[i] + prot_name + '_fasta.txt' #naming according to write_fasta function
        initial_sites_all = run_netNglyc_all(fasta_file, resMap, netnglyc_location) 


        #obtain initial sites for specified chain, only if surface residue or user specified site
        initial_sites = [site for site in initial_sites_all if site.split()[1]==chain[i] and (site in residues or site in priority_sites[i])]

        native_sequons = [trip for trip in surface_sites if trip[0] in initial_sites] #CHECK that native site was also passed by NetNGlyc
        #compile all residues into one list
        
        all_native_seq = [site for seq in native_sequons for site in seq]

        
        # RESNUM OF ASN RESIDUE FOR NEW SITE, only for CHAIN specified and not overlapping with native glycosylation sites

        first_res = [site[0] for site in mutant_sites if site[0].split()[1]==chain[i] and site[0] not in all_native_seq and site[1] not in all_native_seq and site[2] not in all_native_seq]

        #print(d_linear)
        #knock out mutant sites that are within too close a radius of native sites and region to remain uncovered:
        rad = 10
        #print(priority_sites[i])
        if len(priority_sites[i])<2:
            avoid_sites = initial_sites + uncovered_sites[i]
        else:
            avoid_sites = initial_sites + uncovered_sites[i] + priority_sites[i]
        temp_sites = [site1 for site1 in first_res if all(d_linear[site1][site2] > rad for site2 in avoid_sites)]
        #for site1 in first_res:
        #    no_clashes = True
        #    for site2 in avoid_sites:
        #        #if site clashes with site to avoid
        #        if d_linear[site1][site2]<=rad:
        #            no_clashes = False
        #    if no_clashes:
        #        temp_sites.append(site1)
        #update temp sites
        first_res = temp_sites
        #create cur_sites for adding new sites from mutations
        if len(priority_sites[i])<2:
            cur_sites = initial_sites
        else:
            cur_sites = initial_sites+priority_sites[i]

        #loop through each mutant site, write fasta file with mutation and feed into netNglyc
        for res in first_res:
            mutate_fasta(res, allRes, resMap, mutant_sites, fasta_file, destination[i])
            resnum = res.split()[0]
            reschain = res.split()[1]
            mutate_file = destination[i] + 'fasta_' + resnum + '_' + reschain + '.txt'
            new_sites = run_netNglyc_chain(mutate_file, resMap, netnglyc_location)
            add_list = [site for site in new_sites if site not in cur_sites and (site in residues or priority_sites[i])]
            #add unique sites to final list
            if len(add_list) > 0:
                cur_sites.append(add_list[0])
            if len(add_list) > 1:
                for add in add_list:
                    cur_sites.append(add)
        #Jwalk distances:
        dmat = allSASD(cur_sites, destination[i]+'jwalk_protein.pdb', jwalk_location, jwalk_map)
        siteDistances = {}
        for r1 in cur_sites:
            siteDistances[r1] = {}
            for r2 in cur_sites:
                if r1!=r2:
                    if ((r1, r2)) in dmat.keys():
                        tempkey = (r1, r2)
                        siteDistances[r1][r2] = dmat[tempkey]
                    elif ((r2, r1)) in dmat.keys():
                        tempkey = (r2, r1)
                        siteDistances[r1][r2] = dmat[tempkey]
                    else:
                        print((r1, r2))
     


        #sites available for coverage are all surface sites not in uncovered_sites
        covered_sites = [site for site in residues if site not in uncovered_sites[i]]
        print(siteDistances,flush=True)
        # OVERLAPPING and not overlapping SITES FOR EACH SITE NUMBER (DICT)
        first_overlaps = calc_overlap(17.5, siteDistances, priority_sites[i]) #to calc overlap, use siteDistances
        print(first_overlaps)
        #reduce length of combinations
        next_overlaps = refine_clusters(first_overlaps)
        
            
        #to obtain more accurate estimation of nearby residues
        
        print(number[i])
        print(size[i])
        iterate_cl = iterate_clusters(d_linear_surf, next_overlaps, siteDistances, priority_sites[i], number[i], size[i], 17.5)
        print('Modeling Glycans')
        #glycosylate top 10 files
        final_combos = add_glycans(download_directory, iterate_cl, 'fucosylated_full', surface_sites, destination[i], model_glycans[i])
        print('Running Initial Ranking',flush=True)
        #run glyco to score combinations
        output_dict = {}
        j = 0
        n_residues = len(d_linear_surf.keys())
        #for cluster in final_combos:
            #get name of pdb file
        #    cur_pdb = destination[i] + 'Output_'+ str(j)+'.pdb'
            #first calculate coverage over surface residues:
        #    coverage1 = surface_quant(cur_pdb, destination[i]+'GLYCO-main', covered_sites, residues,destination[i])
            #next calculate coverage over residues to remain uncovered:
        #    coverage2 = surface_quant(cur_pdb, destination[i]+'GLYCO-main', uncovered_sites[i], residues,destination[i])
        #    output_dict[j] = [coverage1, coverage2, str(cluster)]
        #    j+=1
            
            
        for cluster in final_combos:
            covered_residues = []
            cur_pdb = destination[i] + 'Output_'+ str(j)+'.pdb'
            for res in cluster:
                for res2, dist in d_linear_surf[res].items():
                    if dist < 15 and res2 not in covered_residues:
                        covered_residues.append(res2)
            rank = len(covered_residues)/n_residues
            cluster.sort()
            output_dict[j] = [rank, str(cluster)]
            j+=1
            
        print('Outputting Constructs')
        #create sets to avoid assigning different ranks to equivalent scores
        s1 = set([v[0] for k,v in output_dict.items()])
        #s2 = set([int(v[1]) for k,v in output_dict.items()])
        #final ranking metric is sum of position for both rankings
        sort1 = sorted(list(s1),reverse=True)
        #sort2 = sorted(list(s2)) #for coverage on uncovered residues, want lowest coverage to rank highest

        final_ranks = {}

        for combo, scores in output_dict.items():
            r1 = sort1.index(scores[0])
            #r2 = sort2.index(scores[1])
            #total = r1+r2
            final_ranks[combo] = r1

        sorted_outputs = [k for k,v in sorted(final_ranks.items(), key=lambda x: x[1])]
        j = 1
        for output in sorted_outputs:
            print(str(j)+ ':')
            print('Output_'+str(output)+ ": " + str(output_dict[output][1]))
            print('Coverage 1:')
            print(output_dict[output][0])
            j+=1
        j = 1
        with open("output.txt", 'a') as file:
            for output in sorted_outputs:
                file.write(str(j)+ ':'+'\n')
                file.write('Output_'+str(output)+ ": " + str(output_dict[output][1])+'\n')
                file.write('Coverage:')
                file.write(str(output_dict[output][0])+'\n')
                j+=1
            file.close()



if __name__ == "__main__":
	arguments = sys.argv[1:]
	argnum = 0
	netnglyc_loc = []
	jwalk_loc = []
	chain = []
	size = []
	number = []
	priority_sites = []
	uncovered_sites = []
	model_glycans = []
	pdb_list = []
	path = []
	destination = []
	while True:
		if argnum >= len(arguments):
			break
		if arguments[argnum] == "-netnglyc_loc":
			while arguments[argnum+1] != "-jwalk_loc":
				argnum += 1
				netnglyc_loc.append(arguments[argnum])
		if arguments[argnum] == "-jwalk_loc":
			while arguments[argnum+1] != "-chain":
				argnum += 1
				jwalk_loc.append(arguments[argnum])
		if arguments[argnum] == "-chain":
			while arguments[argnum+1] != "-size":
				argnum += 1
				chain.append(arguments[argnum])
		if arguments[argnum] == "-size":
			while arguments[argnum+1] != "-number":
				argnum += 1
				size.append(arguments[argnum])
		if arguments[argnum] == "-number":
			while arguments[argnum+1] != "-priority_sites":
				argnum += 1
				number.append(arguments[argnum])
		if arguments[argnum] == "-priority_sites":
			while arguments[argnum+1] != "-uncovered_sites":
				argnum += 1
				priority_sites.append(arguments[argnum])
		if arguments[argnum] == "-uncovered_sites":
			while arguments[argnum+1] != "-model_glycans":
				argnum += 1
				uncovered_sites.append(arguments[argnum])
		if arguments[argnum] == "-model_glycans":
			argnum += 1
			model_glycans.append(arguments[argnum])
		if arguments[argnum] == "-pdb_list":
			while arguments[argnum+1] != "-path":
				argnum += 1
				pdb_list.append(arguments[argnum])
		if arguments[argnum] == "-path":
			while arguments[argnum+1] != "-destination":
				argnum += 1
				path.append(arguments[argnum])
		if arguments[argnum] == "-destination":
			argnum += 1
			destination.append(arguments[argnum])
			break
#		else:
#			fastq_files.append(arguments[argnum])
#			f = os.path.expanduser(fastq_files[-1])			
#			if not os.path.isfile(f):
#				raise Exception('The provided file {0} does not exist'.format(f))
#			fastq_files[-1] = os.path.abspath(f)			
		argnum += 1	
	print(netnglyc_loc, flush = True)
	print(jwalk_loc, flush = True)
	print(chain, flush = True)
	print(size, flush = True)
	print(model_glycans, flush = True)
	print(pdb_list, flush = True)
	print(path, flush = True)
	print(destination, flush = True)
	for place in range(0,len(priority_sites)):
		priority_sites[place]=priority_sites[place].split(',')
		temp_priority_sites = []
		for site in priority_sites[place]:
			site=str(site)
			if len(site.split('-'))==2:
				s1 = site.split('-')[0]
				start = int(s1.split(' ')[0])
				s2 = site.split('-')[1]
				end = int(s2.split(' ')[0])
				i = start
				while i<=end:
					next_site = str(i) + ' ' + chain[place]
					temp_priority_sites.append(next_site)
					i=i+1
			else:
				temp_priority_sites.append(str(site)+ ' ' + chain[place])
		priority_sites[place] = temp_priority_sites
	print(priority_sites, flush = True)
	for place in range(0,len(uncovered_sites)):
		uncovered_sites[place]=uncovered_sites[place].split(',')
		temp_uncovered_sites = []
		for site in uncovered_sites[place]:
			site=str(site)
			if len(site.split('-'))==2:
				s1 = site.split('-')[0]
				start = int(s1.split(' ')[0])
				s2 = site.split('-')[1]
				end = int(s2.split(' ')[0])
				i = start
				while i<=end:
					next_site = str(i) + ' ' + chain[place]
					temp_uncovered_sites.append(next_site)
					i=i+1
			else:
				temp_uncovered_sites.append(str(site)+ ' ' + chain[place])
		uncovered_sites[place] = temp_uncovered_sites
	print(uncovered_sites, flush = True)
	if len(netnglyc_loc)<1:
		raise Exception("Error: Please provide path to the NetNGlyc folder")
	if len(netnglyc_loc)<1:
		raise Exception("Error: Please provide path to the XLM-Tools folder")
	if len(chain)<1:
		raise Exception("Error: Please input chain to glycosylate")
	if len(size)<1:
		raise Exception("Error: Please input the number of targeted glycosylation sites")
	if len(number)<1:
		raise Exception("Error: Please input the number of designs to output")
	if len(model_glycans)<1:
		raise Exception("Error: Please select if you do (Y) or do not (N) want to sample glycans to increase PDB accuracy")
	if len(pdb_list)<1:
		raise Exception("Error: Please input list of PDB files to glycosylate")
	if len(path)<1:
		raise Exception("Error: Please input path of PDB files to glycosylate")
	if len(destination)<1:
		raise Exception("Error: Please specify output desination")
	process(netnglyc_loc, jwalk_loc, chain, size, number, priority_sites, uncovered_sites, model_glycans,pdb_list,path,destination)
    