# IndelPolymorphismsModulateDistalRNAProteinInteractions
This repository includes all code related to the manuscript "Indel polymorphisms modulate distal RNA-protein interactions through RNA secondary structure changes" by Carlos Owusu-Ansah, Elan Shatoff, and Ralf Bundschuh

## Manuscript Figures

- **[`generate_manuscript_figures.ipynb`](https://colab.research.google.com/github/bundschuhlab/IndelPolymorphismsModulateDistalRNAProteinInteractions/blob/main/generate_manuscript_figures.ipynb)**: This notebook generates all the figures used in the manuscript. It processes data from various directories to visualize the effects of indels on protein binding affinity. The output of this notebook is saved in the **`figures/`** folder. You can use the following link to open the notebook in google colab and run the code that generates the figures.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bundschuhlab/IndelPolymorphismsModulateDistalRNAProteinInteractions/blob/main/generate_manuscript_figures.ipynb)

## Simulation Data Folders

- **`randomData/`**: Contains data files related to synthetic transcripts of different lengths, including 150 and 250 nucleotide transcripts. Files are named with the format `DEL_<footprint>_<iterations>_<omit_region>_<seq_length>_dev.dat`

- **`dataHuman/`**: Contains data on human transcripts selected based on their proximity (<200 nucleotides) to known RNA-binding sites. Subdirectories include:
  - `likeRandomNoIndels/`: synthetic indels | arbitrary binding sites
  - `likeRandom/`: synthetic indels near (aprox 25 nts) real indels | arbitrary binding sites
  - `varySites/`: natural indels | arbitrary binding sites
  - `badSites/`: natural indels | HuR binding sites translated by 300 nucleotides.
  - `SameSites/`: natural indels | HuR binding sites

- **`dataHumanNoBS/`**: Contains human transcript data where the transcripts were not selected based on proximity to known binding sites. Subdirectories mirror those in `dataHuman/`.



## Code Folder

The **`code/`** folder contains the code used to generate the simulation data. It consists of four subdirectories:

- **`cppFiles/`**: Contains C++ source code that directly call the Vienna Package to fold RNA transcript fragment and compute various metrics
- **`executables/`**: Contains compiled versions of the C++ files.
- **`pythonFiles/`**: Contains Python scripts that call executables from `executables/` as subroutines.
- **`notebooks/`**: Parses output from `pythonFiles/` to prepare data for generating figures.

### Helper Scripts

- **`pythonFiles/combineHuRandIndels.py`**: Identifies and stores indels proximal (within 200 nts) to HuR binding locations.
- **`pythonFiles/GetSequence.py`**: Extracts sequences containing a specified indel.
- **`pythonFiles/join_indels.py`**: Merges multiple indel-related data files into a single file.
- **`pythonFiles/splitHuRplusIndels.py`**: Splits HuR binding location data and indel data into separate files per chromosome.

### Analysis Scripts

- **`pythonFiles/run_indels.py`**: Simulates indels and binding sites in random transcript fragments to compute $\Delta\Delta G$.
  - Calls the executable `DEL.x`, compiled from `indels.cpp`.
- **`pythonFiles/treatLikeRandom.py`**: Simulates synthetic indels and arbitrary binding sites in human sequences.
  - Calls the executable `indel_human`, compiled from `indel_human.cpp`.
- **`pythonFiles/varySites.py`**: Simulates real indels and arbitrary binding sites.
  - Calls the executable `varySites`, compiled from `varySites.cpp`.
- **`pythonFiles/GetKdandBestSites.py`**: Computes HuR dissociation constant and binding probabilities for sequences containing an HuR binding site.
  - Uses in vitro RNA-binding experiment data (e.g., RNAcompete).
  - Calls the executable `findConstant`, compiled from `reduced_elan.cpp`.
- **`pythonFiles/addDDG.py`**: Computes $\Delta\Delta G$ (change in binding free energy associated with HuR binding at its preferred site).
  - Calls the executable `findEnergy`, compiled from `reduced_indels.cpp`.

