# runListCompare: maximum likelihood sequence comparison, corrected for recombination

`runListCompare.py` and associated scripts provide a Python wrapper for generating maximum likelihood phylogenies from a list of fasta consensus sequence files obtained from mapping to the same reference. The script enables large numbers of samples to be initially handled in parallel and clustered with similar sequences based on a SNP threshold before calculating maximum likelihood trees for each cluster using either [PhyML](http://www.atgc-montpellier.fr/phyml/) or [IQTree](http://www.iqtree.org). Correction for recombination is done with [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML).



## Requirements

- [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)



## Installation

1. Download and decompress the [latest release](https://github.com/davideyre/runListCompare/releases/latest), and `cd` into it
2. Install dependencies using Conda:
   - `conda env create -f conda.yml`
3. Activate and test installation
   - `conda activate runlistcompare`
   - `pytest` (takes ~2mins)



### Python 2.7 version

1. Download and decompress the [release](https://github.com/davideyre/runListCompare/archive/0.2.7.tar.gz), and `cd` into it
2. `conda env create -f conda_python2.yml`
3. `conda activate runlistcompare2`



## Usage

```
python runListCompare.py tests/data/test.ini
```

Here `test.ini` is an ini file containing the required settings. It is advisable to run the above command to test that things are working with the included demo data. Input sequences are listed in a tab separated format, an example is provided in `tests/data/test.seqlist.txt` – the first column can be up to 8 characters in length and is used for the tip labels of the final trees.



## Output files

- `align_snps.fa`, `align_positions.txt`, `align-compare.txt` are the variable sites, their position in the reference genome and the raw pairwise snp difference between samples
- the `cluster` folder contains variable site alignments for each cluster of related samples
- the `cluster_ml` folder contains the output maximum likelihood phylogenies, e.g. `cluster_1_phyml_tree_scaled.tree` is the phyML generated tree scaled to have branch lengths in SNPs and the `cluster_1_cf_scaled.tree` is the ClonalFrameML corrected tree scaled to have branch lengths in SNPs
- the `recomb_corr` folder contains alignment of variable sites with recombination removed for use as input with other software, e.g. BEAST
- the `ML_distances.txt` and `CF_distances.txt` files contain the pairwise distances obtained from the maximum likelihood and ClonalFrameML phylogenies

---

David Eyre & Bede Constantinides  
david.eyre@bdi.ox.ac.uk   
17 April 2019