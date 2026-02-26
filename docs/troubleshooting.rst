Troubleshooting and FAQs
=======================

I am getting an "open: can't stat file" error
---------------------------------------------

This is a known bug when running NetNGlyc on modern architecture, please refer to the following guide: 
https://squanderingti.me/blog/2020/10/28/extreme-debugging.html

I am not seeing the scores of the designed constructs
-----------------------------------------------------

By default the scores are output in console. If these jobs are being run on a computing cluster, ensure that the cluster keeps record of the console
output (typically in a .out file). These reports can be long due to outputs from various simulations but searching for "Output" should land you in the right place.

Are there restrictions on the use of the software?
--------------------------------------------------

Currently, this software is available without a license to any user. It may be distributed and edited freely, however, 
please cite this paper if you plan on using this in a manuscript.

Does HyperImmunISE support multi-chain analysis?
------------------------------------------------

The current implementation of HyperImmunISE does not support generating structures with multiple chains.

Why should I model glycans using pyRosetta?
------------------------------------------

We highly recommend allowing pyRosetta to model glycans for the designs. This increases the computational time significantly
however the outputted PDB will contain a more accurate glycan positioning. This is particularly critical for evaluation using Amber
and scoring may differ if this option is not selected.

Can I select what types of glycans are added to my designs?
-----------------------------------------------------------

Currently there is no implementation to allow for different glycans. However, manually editing the "add_glycans()" function and 
the corresponding library of glycans in pyRosetta can achieve this if needed. 

What value should I select for defining the desired number of glycans?
----------------------------------------------------------------------

This value is intended to vary from case to case. The algorithm will require a minimum of 10 angstroms of spacing between glycans and thus the
density of glycans is inherently limited by this. In our experience, molecules of ~25-50kDa can be adequately covered by between 10-20 glycans.

There is not a PDB structure available for my target, can I still use this program?
-----------------------------------------------------------------------------------

Unfortunately, HyperImmunISE requires a PDB to model glycans. Tools including AlphaFold (https://deepmind.google/technologies/alphafold/) and 
Rosetta (https://www.rosettacommons.org/software) can predict these structures and can be used to generate a PDB structure.

Can I run HyperImmunISE locally?
--------------------------------

HyperImmunISE can be run on a local Linux system, however, due to the computational power required for this analysis, this is not recommended.

Can I run Amber simulations without a GPU?
------------------------------------------

Yes, Amber simulations can be run without a GPU, however analysis will be slower. This can be accomplished by running the "sander" function instead of "pmemd".

