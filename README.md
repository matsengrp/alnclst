
# Alnclst

Alnclst is a command line utility (and python module) for clustering pre-aligned nucleotide sequences.
The method implemented herein is inspired by [UCLUST](http://drive5.com/usearch/manual/uclust_algo.html), but differs in that

* Sequences must be pre-aligned
* Distance measures can take advantage of a trusted alignment
* Iterative [recentering](http://drive5.com/usearch/manual/recenter.html) is supported out of the box


## Pre-alignment

One of the advantages of UCLUST is that its distance measure does not require pre-alignment of sequences.
This lets us cluster very large and possibly messy sequence sets, particularly useful in metagenomics.
However, in cases where we already have a sequence alignment we trust, it would be nice to use this to inform distance inferences.


## Iterative recentering

As [described on the Drive5 webiste](http://drive5.com/usearch/manual/recenter.html), recentering is a useful strategy for obtaining better clustering results for certain applications.
We have baked recentering into the core of Alnclst, so that no extra work is required from the user to recenter their clusters.


## Merge small clusters

For some applications (see Matsen, Small, et al. 2014, in press), it desirable not to have clusters with a very small number of sequences.
Alnclst makes it possible to merge smaller clusters to their nearest neighboring clusters until no more clusters below a certain size are left.
See the 


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


