# HyperImmunISE

HyperImmunISE is a computational framework for designing glycan-mediated masking strategies to modulate protein immunogenicity.  
The platform integrates structural modelling, accessibility analysis, and rule-based site selection to support rational construct design.

It is intended for researchers in antibody engineering, vaccine design, and computational immunology who require reproducible, automatable workflows.

Full documentation, advanced configuration, and dependency details are available at:  
https://hyperimmunise.readthedocs.io/en/latest/

## Quick Start

1. Create the conda environment

conda env create -f hyperimmunise.yml
conda activate hyperimmunise

2. Install external dependencies

Download and install the following tools and make note of their absolute paths:

- NetNGlyc (https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/)
- Jwalk from XLM-Tools (https://github.com/Topf-Lab/XLM-Tools)

3. Run the example

python main.py \
-netnglyc_loc "ABSOLUTE DIRECTORY OF netNglyc-1.0" \
-jwalk_loc "ABSOLUTE DIRECTORY OF XLM-Tools/" \
-chain "K" \
-size 1 \
-number 1 \
-priority_sites "94" \
-uncovered_sites "99-105,221-228" \
-model_glycans "N" \
-pdb_list "HG_HAhead" \
-path "ABSOLUTE DIRECTORY OF HyperImmunISE/example/test_HAhead.pdb" \
-destination "ABSOLUTE DIRECTORY OF target directory"

4. Expected output

The run will generate:

- `Output_0.pdb` — the designed glycosylated protein structure  
- `*_coverage_rank_searchTrajectory_scores.csv` — a log of intermediate glycosylation-site combinations evaluated during iterative pruning

*_coverage_rank_children_scores.csv records all intermediate glycosylation-site combinations evaluated during iterative pruning. Each row corresponds to one candidate combo scored by coverage_rank() and reports the iteration number, predicted surface coverage (fraction and count within the current radius cutoff), and the exact set of sites in that combo.


5. etc
Rerunning jobs after forced termination Sometimes glycan modeling takes longer than expected and the job may need to be stopped. When rerunning, it’s safer to remove the Jwalk_results folder first. This is optional, but helps avoid potential issues.