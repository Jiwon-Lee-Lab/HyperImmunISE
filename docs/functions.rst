Functions
=========

surface_residues(pdb, min_SASA)
-------------------------------

This function takes a protein input and outputs a 
list of all the surface residues if  their SASA > min_SASA, set at 
default to 2.5 Å

based on: http://pymolwiki.org/index.php/FindSurfaceResidues
and modified to use freesasa module to find SASA of residues

freesasa is used for SASA calculation and cited here:

Simon Mitternacht (2016) FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Research 5:189. (doi: 10.12688/f1000research.7931.1)

Parameters:

pdb : FILE
   pdb file of protein.

Returns:

surface_residues: [List] of surface residues in the protein.

residue_info(surface_residues, pdb, destination)
------------------------------------------------

This function generates dictionaries with information
about a list of input residues, including a dictionary of the 
type of each residue and the coordinates
made to reduce the number of times a pdb file needs to be parsed

Parameters:

surface_residues: LIST
   list of strings: 'resnum chain'
   surface residues of interest

pdb : FILE
   pdb file of protein.

Returns:

surface_dict: dictionary of {'resnum chain': residue type}

surface_coords: dictionary of {'resnum chain': coordinates}

map_dict: nested dictionary of {chain: {resnum: index}}

jwalk_map: dictionary that maps old pdb file numbering to new pdb file numbering::

   {chain: {original pdb num: jwalk pdb id}}

consecutive_residues(resMap)
----------------------------

This function finds three or more consecutive residues and 
adds them into a single separate list together

Parameters:

resMap, dictionary of residues by chain

Returns:

[list] of lists of consecutive numbers in the input residue numbers
or empty list if no consecutive residue numbers found

consensus_sequences_similarAA()
-------------------------------

This function creates possible combinations of N-glycosylation 
consensus sequences that may be found in an amino acid sequence.

Returns:

[list] of possible consensus sequences for N-linked glycosylation

identify_native_site(target, query, res_Dict)
---------------------------------------------

This function identify native glycosylation sites in EpitopeCA.txt by 
identifying whether the consensus sequence lies in its sequence.

Parameters:

target, list of consecutive residue numbers from EpitopeCA.txt

query, list of possible consensus sequence combinations

res_Dict, dictionary of {residue number: residue name} from CA file.

Returns:

nested [list] of residue numbers of native glycosylation sites.

id_double_mutation_similarAA(target, query, query_nat, res_Dict)
----------------------------------------------------------------

This function identify novel glycosylation sites two mutations
in EpitopeCA.txt by identifying whether 1 of three amino acids in 
the consensus sequence lies in its sequence.

Parameters:

target, list of consecutive residue numbers from EpitopeCA.txt

query,list of possible consensus sequence combinations

res_Dict, dictionary of {residue number: residue name} from CA file.

Returns:

nested [list] of residue numbers of glycosylation sites that can be singly mutated.

dist(a, b)
----------

This function calculates distance between two points (a, b)

Parameters:

a,  point 1 (x1, y1, z1)
b,  point 2 (x2, y2, z2)   

Returns:

Distance between two points as float.

distance_n_atoms(cords)
-----------------------

This function calculates distance between a set of points. 

Parameters:
cords, list of alpha carbon coordinates from possible glycosylation sites

Returns:
Coord_matrix, nested [dictionary] of distances between coordinates

allSASD(reslist, jwalk_pdbfile, Jwalk_path, jwalk_map)
------------------------------------------------------

This program uses Jwalk, cited here:

Sinnott et al., Combining Information from Crosslinks and Monolinks in the Modeling
of Protein Structures, Structure (2020), https://doi.org/10/1016/j.str.2020.05.012

The Importance of Non-accessible Crosslinks and Solvent Accessible Surface Distance in Modeling Proteins with Restraints From Crosslinking Mass Spectrometry.
J Bullock, J Schwab, K Thalassinos, M Topf. Mol Cell Proteomics. 15, 2491–2500, 2016

Purpose: This function uses Jwalk to find the SAS distances between a list
of residues 

Parameters:

reslist: list of surface residues to calculate SASD for

jwalk_pdbfile : FILE
   name of pdb file which contains the residues in JWALK renumbered format.

Jwalk_path: STRING
   string of the path location to Jwalk download

jwalk_map: dictionary that maps new pdb file numbering to original pdb file numbering 
   {chain: {jwalk pdb num: original num}}

Returns:

SASD_dict: dict of SASD for each residue pair in residue list
   key: (res1, res2) - strings
   value: SASD - float 

calc_overlap(sphere_rad, distances, priority_sites)
---------------------------------------------------

This function identifies and creates clusters 

Parameters:

sphere_rad, input radius of glycan size

distances, distance dictionary between coordinate pairs

priority_sites, list of residues that must be in each combo

Returns:

all_clusters, Dictionary of {site number: non overlapping sites}

reduce_clusters(sphere_rad, distances, combinations, priority_sites)
--------------------------------------------------------------------

This function reduces size of clusters until they have no remaining clashes under new clash radius

Parameters:

sphere_rad, input radius of glycan size

distances, distance dictionary between coordinate pairs

combinations, nested list of combinations for clusters created from previous overlap radius

priority_sites, list of residues that must be in each combo

Returns:

new_clusters, list of new combinations that have no clashes based on starting combination and sphere_rad

write_fasta(pdb, destination)
-----------------------------

Takes pdb file and write FASTA sequence manually to a txt file
(ensures that all residues that appear in PDB match those submitted to
NetNGlyc)

Parameters:

pdb: pdb file of protein
destination: path where user wants to send txt file of FASTA seq, should 
end in "/"

Returns:

None

mutate_fasta(res_to_mutate, allRes, map_dict, mutant_sites, original_fasta, destination)
----------------------------------------------------------------------------------------

Takes pdb file and write FASTA sequence manually to a txt file
with residue to mutate changed to appropriate residue to ensure consensus sequence

Parameters:

res_to_mutate: string ('number chain' e.g. '20 A') of residue to mutate (based on PDB numbering NOT mapped integer numbering)
   Provides the first position residue in set of 3, may not actually be the residue that will change (could be third residue)

allRes: dict of surface residues (number and chain) and corresponding amino acid values

map_dict: nested dict::

   nested dictionary of {chain: {resnum: index}} 

destination: path where user wants to send txt file of FASTA seq, should end in '/'

mutant_sites: list of list of sites that can be consensus sequence with mutation

original_fasta: txt file of fasta seq for original protein sequence

Returns:

None

run_netNglyc_all(fasta, map_dict, threshold = 6, netnglyc_loc)
--------------------------------------------------------------

Evaluates likelihood of glycosylation using NetNGlyc

Parameters:

fasta : file
   input file of protein fasta sequence.
map_dict : nested dict
   nested dictionary of {chain: {resnum: index}}
threshold : integer, optional
   number of neural nets that must agree for site to count as 'likely' glycosylated. The default is 6.
netnglyc_loc : string
   path to netnglyc folder

Returns:

predicted_sites: list
   list of strings ('resnum chain') predicted by netNglyc to be glycosylated.


run_netNglyc_chain(fasta, map_dict, threshold = 6,netnglyc_loc)
---------------------------------------------------------------

Evaluates likelihood of glycosylation using NetNGlyc

Parameters:

fasta : file
   input file of protein fasta sequence.
map_dict : nested dict
   nested dictionary of {chain: {resnum: index}}
threshold : integer, optional
   number of neural nets that must agree for site to count as 'likely' glycosylated. The default is 6.
netnglyc_loc : string
   path to netnglyc folder

Returns:

predicted_sites: list
   list of strings ('resnum chain') predicted by netNglyc to be glycosylated.

coverage_rank(distance_n_atoms, clusters, i, radius = 17.5)
-----------------------------------------------------------

Refine possible cluster combinations based on predicted
coverage -- only high coverage combinations remain 

Parameters:

distance_n_atoms : dictionary of {site: {residue: distance}} using linear
   distance calculation (all surface residues)
   
clusters : dict
   Dictionary of {site number: non overlapping sites}

radius: int
   distance cutoff of nearby residues considered to be covered by glycan

i: int
   final number of combinations needed

Returns:

refined_clusters: dict
   Dictionary of {site number: non overlapping sites} (i.e. combos/clusters) with lower-coverage combinations eliminated

iterate_clusters(distance_n_atoms, clusters, siteDistances, priority_sites, i, j, initial_r)
--------------------------------------------------------------------------------------------

Iteratively increase radius of overlap to reduce size of combinations in clusters
and check coverage to reduce number of combinations

Parameters:

distance_n_atoms : dictionary of {site: {residue: distance}} using linear distance calculation (all surface residues)

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

Returns:

next_clusters: dict
   Dictionary of {site number: non overlapping sites} with j site numbers, combinations of size i

refine_clusters(cluster_list, max_len=0)
----------------------------------------

Reduce number of combinations in cluster_list so that only clusters of 
longer length are included

Parameters:

cluster_list: nested LIST
   Nested list of all combinations
max_len: integer, optional
   Minimum length of a combination. The default is 0.

Returns:

new_list: LIST
   Nested list of combinations with minimum length of max_len.

add_glycans(pdb, combo_list, glycan, native_sites, destination, model_glycans)
------------------------------------------------------------------------------

Generates glycosylated pdb files for each combination

Parameters:

pdb : file
   pdb file of protein.
combo_list : nested list
   list of list of site combinations for adding glycans.
glycan: string
   indicates glycan to be added
native_sites:
   list of native glycosylation sites identified previously
model_glycans:
   STRING 'Y' or 'N', only if 'Y' will have rosetta do glycan sampling of conformations

Returns:

updated_combos list
   list of lists of combos to account for any residues removed due to disulfide bonds.




