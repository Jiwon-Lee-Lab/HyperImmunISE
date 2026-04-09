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

Why does glycan modelling fail with "Unrecognized sugar 3-letter code 'fuc'"?
-----------------------------------------------------------------------------

This usually indicates a mismatch between the HyperImmunISE glycan preset and the installed PyRosetta carbohydrate database.
HyperImmunISE now requires the preset ``fucosylated_full`` to exist directly in your local PyRosetta carbohydrate database before glycan modelling starts.

To register it manually, add the following entry to ``pyrosetta/database/chemical/carbohydrates/common_glycans/common_names.txt``::

  fucosylated_full  fucosylated_full.iupac

Then create ``pyrosetta/database/chemical/carbohydrates/common_glycans/fucosylated_full.iupac`` with this content::

  b-D-GlcpNAc-(1->2)-a-D-Manp-(1->3)-[a-D-GlcpNAc-(1->2)-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-[a-L-Fucp-(1->6)]-b-D-GlcpNAc-

If either the registration line or the ``.iupac`` file is missing, HyperImmunISE will stop before running PyRosetta.

What value should I select for defining the desired number of glycans?
----------------------------------------------------------------------

This value is intended to vary from case to case. The algorithm will require a minimum of 10 Å of spacing between glycans and thus the
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

Why can ``Output_0.pdb`` or the selected glycosylation-site combination change between runs?
---------------------------------------------------------------------------------------------

There are two main reasons this can happen:

1. PyRosetta can introduce stochastic behavior. If glycan modeling or related packing steps are
   allowed to sample different trajectories, repeated runs can produce different structural
   outputs even when the same candidate sites are provided.
2. Earlier versions of the site-selection pipeline did not define a deterministic tie-break rule
   when multiple candidate combinations received the same coverage score. In that situation, the
   order of equally scored combinations could vary, and a different combination could be passed
   downstream to PyRosetta.

If you see output variation, it is reasonable to run the same input more than once to confirm
whether the top result is stable. If repeated runs converge on the same site combination and
similar structures, that is a useful sign that the optimization is behaving consistently.

