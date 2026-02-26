Usage
=====

Installing HyperImmunISE and Dependencies
---------------------------------------

To use HyperImmunISE, first pull from git: https://github.com/Jiwon-Lee-Lab/HyperImmunISE.git

.. code-block:: console

  $ git clone https://github.com/Jiwon-Lee-Lab/HyperImmunISE.git

HyperImmunISE requires Python 3 (https://www.python.org/) and the following dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_kmp_llvm
  - appdirs=1.4.4=pyhd3eb1b0_0
  - biopython=1.78=py310h7f8727e_0
  - blas=1.0=mkl
  - bzip2=1.0.8=h7b6447c_0
  - ca-certificates=2023.08.22=h06a4308_0
  - certifi=2023.11.17=py310h06a4308_0
  - cudatoolkit=11.8.0=h6a678d5_0
  - icu=73.1=h6a678d5_0
  - intel-openmp=2021.4.0=h06a4308_3561
  - ld_impl_linux-64=2.38=h1181459_1
  - libboost=1.82.0=h109eef0_2
  - libboost-python=1.82.0=py310hcb52e73_6
  - libffi=3.3=he6710b0_2
  - libgcc-ng=13.2.0=h807b86a_3
  - libstdcxx-ng=13.2.0=h7e041cc_3
  - libuuid=1.0.3=h7f8727e_2
  - llvm-openmp=14.0.6=h9e868ea_0
  - lz4-c=1.9.4=h6a678d5_0
  - mako=1.2.3=py310h06a4308_0
  - markupsafe=2.1.1=py310h7f8727e_0
  - mkl=2021.4.0=h06a4308_640
  - mkl-service=2.4.0=py310h7f8727e_0
  - mkl_fft=1.3.1=py310hd6ae3a3_0
  - mkl_random=1.2.2=py310h00e6091_0
  - ncurses=6.3=h7f8727e_2
  - numpy=1.24.3=py310hd5efca6_0
  - numpy-base=1.24.3=py310h8e6c178_0
  - openssl=1.1.1w=h7f8727e_0
  - pip=21.2.4=py310h06a4308_0
  - platformdirs=3.10.0=py310h06a4308_0
  - pycuda=2023.1=py310hd456308_0
  - pyrosetta=2023.33+release.9c16e13=py310_0
  - python=3.10.4=h12debd9_0
  - python_abi=3.10=2_cp310
  - pytools=2023.1.1=pyhd8ed1ab_0
  - readline=8.1.2=h7f8727e_1
  - setuptools=61.2.0=py310h06a4308_0
  - six=1.16.0=pyhd3eb1b0_1
  - sqlite=3.38.3=hc218d9a_0
  - tk=8.6.11=h1ccaba5_1
  - typing_extensions=4.7.1=py310h06a4308_0
  - tzdata=2022a=hda174b7_0
  - wheel=0.37.1=pyhd3eb1b0_0
  - xz=5.4.2=h5eee18b_0
  - zlib=1.2.13=h5eee18b_0
  - zstd=1.5.5=hc292b87_0
  - pip:
      - freesasa==2.1.0

The provided "hyperimmunise.yml" can also be used to directly import the corresponding packages.
It is recommended that these are installed into a conda environment.

Jwalk, a software used to calculate solvent accessible surface distances 
(Bullock, Joshua Matthew Allen, et al. Molecular & Cellular Proteomics, 2016)
can be downloaded as part of XLM-Tools (https://github.com/Topf-Lab/XLM-Tools/).

NetNGlyc, a software to predict if a glycosylation sequon will contain a glycan
(Gupta R, Brunak S. Pac Symp Biocomput., 2002) can be downloaded 
(https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/. If the primary link is unavailable, use https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netNglyc&version=1.0_copy&packageversion=1.0d&platform=Linux). 
This software may need to be edited to work on modern operating systems, a guide
to a hotfix to allow for this can be found here:
https://squanderingti.me/blog/2020/10/28/extreme-debugging.html

Amber, a molecular dynamic simulation tool (R. Salomon-Ferrer, D.A. Case, R.C. Walker. WIREs 
Comput. Mol. Sci., 2013) can be downloaded from https://ambermd.org/.
Amber simulation can be edited based on user needs, the default framework contained here
assumes that Amber is installed with gpu compatibility and a gpu is available for simulations.

Necessary System Architecture
-----------------------------

This script is intended to be run on a computing cluster and is currently written to accommodate 
a CentOS Linux architecture with both CPU and GPU nodes. Other architectures may be compatible
but have not been tested with this program. It is highly recommended that this software be run 
on a high-performance computing cluster, as the process is memory intensive. 

Designing Constructs with HyperImmunISE
---------------------------------------

Constructs are designed by running the 'main.py' script with the following flags:

-netnglyc_loc <file path to the NetNGlyc installation>

-jwalk_loc <file path to the XLM-Tools folder>

-chain <chain in PDB file to be used for analysis, contained in quotations>

-size <integer with the targeted number of glycans to add>

-number <integer with the number of designs to create>

-priority_sites <residue numbers of site which must contain a glycan>

-uncovered_sites <residues number of sites which should remain unblocked by glycans>"99-105,221-228" 

-model_glycans <whether Rosetta should model glycans, "Y" or "N"> 

-pdb_list <name of PDB file, no extension (.pdb) necessary>

-path <path to the PDB construct>

-destination <path to the output folder>

Residue numbers for the "priority_sites" and "uncovered_sites" can be entered as a range (e.g., "99-105")
and/or as a sequence (e.g., "99-105,221-228,330"). Numbering is based off of the PDB numbering.

Output designs will be written in the destination folder as .pdb files and scoring evaluations will be output
in text in the console.

An example input is provided below:

.. code-block:: console

  $ python main.py -netnglyc_loc "/NetNGlyc/netNglyc-1.0" -jwalk_loc "/XLM-Tools-master" -chain "K" -size 10 -number 10 -priority_sites "" -uncovered_sites "99-105,221-228" -model_glycans "Y" -pdb_list "mono_head" -path "/mono_head.pdb" -destination "/"

Evaluating HyperImmunISE Constructs in Amber
--------------------------------------------

With the installation of Amber, navigate to the installed directory:

.. code-block:: console

  $ source amber.sh

Made folder called Simulations:

.. code-block:: console

  $ mkdir Simulations/
  $ cd Simulations/

Copy your desired .pdb structure to this folder.

Navigate to the AmberTools/src/ folder and run cpptraj:

.. code-block:: console

  $ cpptraj
  $ cd Simulations/
  $ parm <path_to_file>/<filename>.pdb
  $ prepareforleap crdset MyCrd name Final out <path_to_file>/leap.<filename>.in leapunitname mol pdbout <path_to_file>/<filename>.cpptraj.pdb nowat noh keepaltloc highestocc
  $ quit

Edit the resulting .in file to have the following:

At the top:

source leaprc.protein.ff19SB

source leaprc.GLYCAM_06j-1

source leaprc.water.opc3

loadAmberParams frcmod.ionslm_126_opc

and at the bottom:

solvateOct mol OPC3BOX 8.0

addIons mol Cl- 0

addIons mol Na+ 0

saveAmberParm mol <filename>.prmtop <filename>.inpcrd

quit

Navigate to the Amber bin folder and then run tleap and parmed:

.. code-block:: console

  $ tleap -f <path_to_file>/leap.<filename>.in
  $ parmed
  $ parm <filename>.prmtop
  $ HMassRepartition
  $ outparm <filename>_hmr.prmtop
  $ quit

Move the resulting .prmtop, inpcrd, and _hrm.prmtop files to the Simulations folder and navigate there. 

Make the following files:

File - minimize.in

Minimize

&cntrl

imin=1,

ntx=1,

ntxo=1,

irest=0,

maxcyc=5000,

ncyc=2500,

ntpr=250,

ntwx=0,

cut=8.0,

/

File - heat.in

Heat

&cntrl

imin=0,

ntx=1,

ntxo=1,

irest=0,

nstlim=11000,

dt=0.004,

ntf=2,

ntc=2,

tempi=200.0,

temp0=310.0,

ntpr=100,

ntwx=100,

ntwr=1000,

cut=8.0,

ntb=1,

ntp=0,

ntt=2,

vrand=1000

ig=-1,

/

&wt type='TEMP0', istep1=0, istep2=10000, value1=200.0,

value2=310.0 /

&wt type='TEMP0', istep1=10001, istep2=11000, value1=310.0,

value2=310.0 /

&wt type='END' 

File - equilibrate.in

Equilibrate w/ HMR

&cntrl

imin=0,

ntx=1,

ntxo=1,

irest=0,

nstlim=25000000,

dt=0.004,

ntf=2,

ntc=2,

tempi=310.0,

temp0=310.0,

ntpr=100000,

ntwx=50000,

ntwr=1000000,

cut=8.0,

ntb=2,

ntp=1,

ntt=2,

vrand=10000,

ig=-1,

/

The following commands are best run on a computing cluster:

.. code-block:: console

  $ module load cuda
  $ cd <path_to_amber>
  $ source amber.sh
  $ cd Simulations
  $ ../bin/sander -O -i minimize.in -o minimize.out -p <filename>_hmr.prmtop -c <filename>.inpcrd -r min.rst7
  $ ../bin/pmemd.cuda -O -i heat.in -o heat.out -p <filename>_hmr.prmtop -c min.rst7 -r heat.rst7 -x heat.nc
  $ ../bin/pmemd.cuda -O -i equilibrate.in -o equilibrate.out -p <filename>_hmr.prmtop -c heat.rst7 -r equilibrate.rst7 -x equilibrate.nc

Then load the conda environment which contains a mdtraj environment (install with pip or conda if you don't have it) and run the following:

.. code-block:: python

  import mdtraj as mdtraj
  import os
  os.chdir('Simulations')
  mtraj = mdtraj.load("equilibrate.nc", top="<filename>_hmr.prmtop")
  mtraj=simul.remove_solvent()
  mtraj.save_pdb("<filename>_nosolvent.PDB")
  sr = mdtraj.shrake_rupley(mtraj, probe_radius=1, mode='residue')
  avg_sr = np.mean(sr,0)
  np.savetxt("avg_sr.csv",avg_sr,delimiter=",")
  
The corresponding degree of blocking will be shown in the "avg_sr.csv" file and the numbering will follow the PDB numbering on your input design.


