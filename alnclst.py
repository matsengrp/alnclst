#!/usr/bin/env python
"""Alnclst: A simple UCLUST-like clustering algorithm for pre-aligned sequences with built in
recentering."""

import itertools as it
from numpy import mean
from Bio import SeqIO, Align, SeqRecord
from Bio.Align import AlignInfo
import argparse
import random
import csv

__version__ = "0.1.1"

__verbose__ = False


def hamming_dist(seq1, seq2):
    """Normalized hamming distance that ignores deletion characters. """
    diffs = 0
    length = 0
    for x, y in zip(str(seq1), str(seq2)):
        if x == '-' or y == '-':
            continue
        elif x != y:
            diffs += 1
        length += 1
    try:
        return float(diffs) / length
    except:
        return 0.5


class Clustering(object):
    def __init__(self, seqrecords, threshold, consensus_threshold):
        """Create a new clustering, using a UCLUST-esque greedy cluster creation mechanism."""
        # iterate through each cluster and find a representative
        self.clusters = list()
        self.threshold = threshold
        self.consensus_threshold = consensus_threshold
        for record in seqrecords:
            try:
                min_dist, clst = min((c.distance(record), c) for c in self.clusters)
            except ValueError:
                # If there aren't any clusters yet, we need to create a new one; this is a somewhat hackish
                # way of doihng so
                min_dist = threshold + 1
            if min_dist < threshold:
                clst.add(record)
            else:
                new_cluster = Cluster(record, consensus_threshold)
                self.clusters.append(new_cluster)

    def __repr__(self):
        return "Clustering(size %s):[\n%s\n]" % (len(self.clusters),
                "\n".join("    " + str(c) for c in self.sorted_clusters()))

    def merge_small_clusters(self, min_per_cluster):
        """Merges the smallest clusters with nearest neighbours until no cluster is smaller than min_per_cluster"""
        while any(c.size() < min_per_cluster for c in self.clusters):
            _, smallest = min((c.size(), c) for c in self.clusters)
            _, closest = min((c.distance(smallest.centroid), c) for c in self.clusters if c != smallest)
            closest.merge(smallest)
            self.clusters.remove(smallest)

    def mapping_iterator(self):
        """For spitting out to csv, essentially."""
        for i, cluster in enumerate(self.sorted_clusters()):
            for record in cluster.members:
                yield (i, record.name, cluster.distance(record))

    def sorted_clusters(self):
        """Generator for iterating over clusters in order of greatest to smallest"""
        return (c for _, c in sorted((-c.size(), c) for c in self.clusters))

    def recenter_iterator(self):
        """This iterates over cluster members of all clusters. First it yields each of the cluster centroids,
        in order of sorted_clusters. Next it yields each of the next most central members from each cluster,
        again in order of sorted clusters. It does this till there are no more cluster members."""
        cluster_iterators = (c.recenter_iterator() for c in self.sorted_clusters())
        return (r for r in it.chain(*it.izip_longest(*cluster_iterators)) if r)

    def recenter(self, n):
        """Does the requested number of recenterings."""
        clustering = self
        for i in xrange(n):
            clustering = Clustering(clustering.recenter_iterator(), threshold=clustering.threshold,
                    consensus_threshold=clustering.consensus_threshold)
        return clustering

    def write(self, handle):
        out_writer = csv.writer(handle)
        out_writer.writerow(('cluster_id', 'sequence', 'distance'))
        for row in self.mapping_iterator():
            out_writer.writerow(row)
        handle.close()


class KMeansClsutering(Clustering):
    """This class does a very simple KMeans clustering using some of the machinery of the other two classes.
    Should really refactor things so there is cleaner code inheritance and reuse."""
    def __init__(self, seqrecords, k, consensus_threshold, max_iters=50):
        print "Running KMeans for K=", k
        self.seqrecords = [sr for sr in seqrecords]
        self.clusters = [Cluster(sr, consensus_threshold) for sr in random.sample(self.seqrecords, k)]
        last_clustering = []
        itercount = 0
        while last_clustering != self.sorted_membernames() and itercount < max_iters:
            if __verbose__:
                print "On loop", itercount
                print self
                print ""
            itercount += 1
            last_clustering = self.sorted_membernames()
            self.reseed_empties()
            for c in self.clusters:
                c.recenter()
                c.members = []
            for sr in self.seqrecords:
                min_dist, clst = min((c.distance(sr), c) for c in self.clusters)
                clst.add(sr)
        print "KMeans completed in", itercount, "iterations"

    def average_distance(self):
        result = mean([c.distance(m) for c in self.clusters for m in c.members])
        info_str = "Average dist: %s\n" % result
        if __verbose__:
            print info_str
        return result

    def reseed_empties(self):
        empty_clusters = filter(lambda c: c.size() == 0, self.clusters)
        if empty_clusters:
            if __verbose__:
                print "Doing a reseed\n"
            dist_sort = sorted((-c.distance(m), m, c) for c in self.clusters for m in c.members)
            for (_, distal_m, distal_c), empty_c in it.izip(dist_sort, empty_clusters):
                distal_c.members.remove(distal_m)
                empty_c.add(distal_m)
                
    def sorted_membernames(self):
        return [sorted([m.name for m in c.members]) for c in self.sorted_clusters()]


class Cluster(object):
    """This class represents a specific cluster. Initally, gets stamped out just with a representative
    sequences and empty sequence list. Sequences have to be added..."""
    def __init__(self, centroid, consensus_threshold):
        self.centroid = centroid
        self.consensus_threshold = consensus_threshold
        self.members = list()
        if self.centroid.name != 'consensus':
            self.members.append(centroid)

    def __repr__(self):
        return "Cluster(size: %s members: %s centroid: %s)" % (self.size(),
                sorted([m.name for m in self.members]),
                self.centroid.seq)

    def distance(self, record):
        """Simple hamming distance between rep/center and seq arg."""
        return hamming_dist(self.centroid.seq, record.seq)

    def add(self, record):
        """Add a member to this cluster"""
        if record.name != 'consensus':
            self.members.append(record)

    def merge(self, cluster):
        """Merge two clusters"""
        self.members += cluster.members

    def average_distance(self, record):
        "Average distance to other members of the cluster."
        return mean([hamming_dist(r.seq, record.seq) for r in self.members])

    def recenter(self):
        """Assign a new centroid to the cluster"""
        self.centroid = self.consensus()
        return self.centroid

    def size(self):
        "Number of members..."
        return len(self.members)

    def consensus(self):
        aln = Align.MultipleSeqAlignment(self.members)
        info = AlignInfo.SummaryInfo(aln)
        seq = info.dumb_consensus(threshold=self.consensus_threshold, ambiguous='N')
        return SeqRecord.SeqRecord(seq, name='consensus', id='consensus')

    def recenter_iterator(self):
        """This method iterates over the sequences in the cluster in order of distance from centroid,
        starting from the centroid."""
        # Note that this is where the cluster recenter is happening. Not sure if this is really the best place
        # for it. I guess it doesn't REALLY matter too much as long as it happens before this iterator is
        # called...
        self.recenter()
        yield self.centroid
        for _, record in sorted((self.distance(r), r) for r in self.members):
            yield record


def threshold_handler(args):
    # Grab inputs, sort by ungapped length
    seqrecords = SeqIO.parse(args.alignment, 'fasta')
    seqrecords = (x for _, x in sorted((-len(sr.seq.ungap('-')), sr) for sr in seqrecords))

    # All of the work: cluster, recenter, and merge small clusters as necessary
    clustering = Clustering(seqrecords, args.threshold, args.consensus_threshold)
    print "Initial clustering complete. Starting recentering..."
    clustering = clustering.recenter(args.recenterings)
    print "Recenterings complete."
    if args.min_per_cluster:
        clustering.merge_small_clusters(args.min_per_cluster)
        print "Finished Merging small clusters"

    # Write output
    clustering.write(args.output)


def kmeans_handler(args):
    if args.seed:
        random.seed(args.seed)

    def clustering(args):
        seqrecords = SeqIO.parse(args.alignment, 'fasta')
        return KMeansClsutering(seqrecords, args.k, args.consensus_threshold, args.max_iters)

    if args.batches:
        print "Starting new batch"
        _, clusts = min((c.average_distance(), c) for c in (clustering(args) for i in xrange(args.batches)))
    else:
        clusts = clustering(args)

    clusts.write(args.output)


def add_common_args(parser):
    parser.add_argument('alignment', help="In fasta format, sequenced to be clustered")
    parser.add_argument('output', type=argparse.FileType('w'),
        help="Output csv file")
    parser.add_argument('-c', '--consensus-threshold', type=float, default=0.25,
        help="Consensus threshold for computing cluster consensus sequences")
    parser.add_argument('-v', '--verbose', action="store_true",
        help="Print out information during clustering")


def get_args():
    parser = argparse.ArgumentParser(description=__doc__ + " (version: %s)" % __version__)

    subparsers = parser.add_subparsers(title="clustering method",
        help="For further help, try `alnclst.py <subcommand> -h`")
    threshold_parser = subparsers.add_parser("threshold",
        help="Cluster based on distance thresholds for cluster membership decision; similar to UCLUST")
    kmeans_parser = subparsers.add_parser("kmeans",
        help="""KMeans clustering; specify number of clusters a priori and minimize average within
        cluster distnaces""")

    add_common_args(kmeans_parser)
    kmeans_parser.add_argument('-k', type=int, required=True, help="Number of clusters")
    kmeans_parser.add_argument('-M', '--max-iters', type=int, default=100,
        help="Maximum number of iterations.")
    kmeans_parser.add_argument('-b', '--batches', type=int,
        help="""Perform this many kmeans, and take the one with the smallest average_distance from cluster
        centers to corresponding members.""")
    kmeans_parser.add_argument('-s', '--seed',
        help="Specify random seed for initial set of clusters in NIPALS algorihtm")
    kmeans_parser.set_defaults(func=kmeans_handler)

    add_common_args(threshold_parser)
    threshold_parser.add_argument('-t', '--threshold', type=float, required=True,
        help="Roughly speaking, the maximum cluster 'radius', assuming -m is not specified")
    threshold_parser.add_argument('-r', '--recenterings', type=int, default=0,
        help="Number of recentering steps to perform")
    threshold_parser.add_argument('-m', '--min-per-cluster', type=int,
        help="""Minimum number of sequences allowed in clusters. Clusters smaller than this will be merged
        with the nearest neighbouring cluster until all clusters are at least this big.""")
    threshold_parser.set_defaults(func=threshold_handler)

    return parser.parse_args()


def main(args):
    global __verbose__
    __verbose__ = args.verbose
    args.func(args)


if __name__ == '__main__':
    main(get_args())


