# haplot
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11001623-blue)](https://doi.org/10.5281/zenodo.11001623)

A set of easy-to-use plot functions can be used for genomic visual analysis

## Functions
- [x] manhattan and qq plot
- [x] LD heatmap plot
- [x] haplotype geographic distribution plot
- [x] haplotype network plot
- [x] gene structure plot
- [x] haplotype boxplot
- [x] nucleotide diversity and Fst plot
- [x] gene structure with haplotype plot

## Dependencies
haplot support python >= 3.10, with following packages installed:
- numpy >= 1.22.3
- pandas >= 1.4.2
- matplotlib >= 3.6.2
- scipy >= 1.9.3
- geopandas >= 0.13.2
- networkx >= 3.1

## Installation
```bash
git clone https://github.com/swu1019lab/haplot.git
cd haplot
# recommended to install the package after installing the dependencies with pip mirror
python setup.py install
```