
# Alnclst

Alnclst is a command line utility (and python module) for clustering pre-aligned nucleotide sequences.
Two alignment strategies are presented.


## Threshold based

The threshold based clustering algorithm is inspired by [UCLUST](http://drive5.com/usearch/manual/uclust_algo.html).
It differs in that

* Sequences must be pre-aligned for distance evalutation (normalized Hamming):
    * this is a disadvantage if you don't have or need an alignment, especially if it would be expensive to compute one;
    * however, it can be an advantage if you have a trusted alignment, from which you would like more informed clustering.
* Iterative [recentering](http://drive5.com/usearch/manual/recenter.html) is supported out of the box, resulting in better clusters, particularly for finding OTUs or cluster representatives (see the `-r|recenterrings` flag).
* It is possible to specify a minimum cluster size, which is desireable for some applications (see [Matsen, Small, et al. 2014](http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1003493)).
  Specifying `-m|--min-per-cluster` setting from the CLI makes Alnclst merge small clusters to nearest neighboring clusters (defined by centroid distance) until no clusters are left smaller than this specified size.

This method is available via the `alnclst.py threshold` subcommand (add `-h` for further help).


## KMeans clustering

Alnclst (as of `v0.1.0`) also supports K-means clustering.
Aside from being able to specify K -- as one might expect -- it also supports specifying a number of clustering batches to run via the `-b|--batches` option.
The one with the smallest average distance from sequences to their cluster centers is taken as the definitive clustering.

This method is available via `alnclst.py kmeans` (add `-h` for further help).


## Installation

This program is written in python, so you will need an installation of python on your system.
Currently, this has only been tested in python 2.7.

You will also need the `scipy`, `numpy`, and `biopython` libraries, as well as 'setuptools'.
On an Ubuntu system, these can all be installed by running

    sudo apt-get install python-pip python-scipy python-numpy python-biopython

You can alternatively, after ensuring you have [pip installed](http://www.pip-installer.org/en/latest/installing.html), install the other python packages using `pip` and perhaps follow any package/OS specific instructions as available on the web.


## Issues / TODO

* Currently, Alnclst does best when input sequences overlap well.
* The current implementation uses consensus sequences as centers (which helps a bit with the above), but it would be nice to allow for centroids as well.
* Enable `--min-per-cluster` for K-means.

