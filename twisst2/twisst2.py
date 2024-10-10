#!/usr/bin/env python

import random
import argparse
import gzip
import sys
import itertools
import numpy as np
from twisst2.TopologySummary import *
from sticcs import sticcs
import cyvcf2

try: import tskit
except ImportError:
    pass

def comboGen(groups, max_iterations):
    n_combos = np.prod([len(t) for t in groups])
    if n_combos <= max_iterations:
        for combo in itertools.product(*groups):
            yield combo
    else:
        for i in range(max_iterations):
            yield tuple(random.choice(group) for group in groups)


def contains_interval(intervals, interval):
    return (intervals[:,0] <= interval[0]) & (interval[1] <= intervals[:,1])


def get_unique_intervals(intervals, verbose=True):
    
    intervals.sort(axis=1)
    intervals = intervals[intervals[:,1].argsort()]
    intervals = intervals[intervals[:,0].argsort()]
    
    intervals = np.unique(intervals, axis=0)
    
    starts = intervals[:,0]
    starts = np.unique(starts)
    n_starts = starts.shape[0]
    
    ends = intervals[:,1]
    ends.sort()
    ends = np.unique(ends)
    n_ends = ends.shape[0]
    
    current_start = starts[0]
    next_start_idx = 1
    next_end_idx = 0
    next_start = starts[next_start_idx]
    next_end = ends[next_end_idx]
    
    output = []
    
    while True:
        if next_start and next_start <= next_end:
            if next_start == current_start:
                #if we get here, the previous end was exactly 1 before the next start, so we jump to the next start
                next_start_idx += 1
                next_start = starts[next_start_idx] if next_start_idx < n_starts else None
            else:
                #we have to cut before this next start and start again
                output.append((current_start, next_start-1))
                current_start = next_start
                next_start_idx += 1
                next_start = starts[next_start_idx] if next_start_idx < n_starts else None
        else:
            #otherwise we use the end and start at next position
            output.append((current_start, next_end))
            current_start = next_end + 1
            next_end_idx += 1
            if next_end_idx == n_ends: break
            next_end = ends[next_end_idx]
    
    #retain only those that are contained within the starting intervals
    output = [interval for interval in output if np.any(contains_interval(intervals, interval))]
    
    return np.array(output)


#def split_weights(weights, intervals, new_intervals):
    
    #n_new_intervals = len(new_intervals)
    
    #n_topos = weights["weights"].shape[1]
    
    #_weights_ = np.zeros((n_new_intervals, n_topos), dtype=int)
    
    #k=0
    #max_k = len(intervals)
    #for j in range(n_new_intervals):
        #if intervals[k][0] <= new_intervals[j][0] and intervals[k][1] >= new_intervals[j][1]:
            ##k contains j
            #_weights_[j,:] += weights["weights"][k]
        
        #if intervals[k][1] == new_intervals[j][1]:
            ##both are ending, so advance k
            #k += 1
            #if k == max_k: break
    
    #weights_split = {"topos": weights["topos"], "weights": _weights_, "weights_norm":np.apply_along_axis(lambda x: x/x.sum(), 1, _weights_)}
    
    #return weights_split



def show_tskit_tree(tree, node_labels=None): print(tree.draw(format="unicode", node_labels=node_labels))

def is_bifurcating(tree):
    for node in tree.nodes():
        if len(tree.children(node)) > 2:
            return False
    return True


def list_all_topologies_tskit(n, include_polytomies = False):
    assert n <= 7, "There are 135135 rooted topologies with 8 tips (660032 including bifurcations). I doubt you want to list all of those."
    if include_polytomies: return list(tskit.all_trees(n))
    return [tree for tree in tskit.all_trees(n) if is_bifurcating(tree)]

def all_possible_patterns(n, lower=2, upper=None):
    if upper: assert upper <= n
    else: upper=n
    for k in range(lower,upper+1):
        for indices in list(itertools.combinations(range(n), k)):
            pattern = np.zeros(n, dtype=int)
            pattern[(indices,)] = 1
            yield pattern


def contains_array(arrayA, arrayB):
    assert arrayA.shape[1] == arrayB.shape[1]
    if arrayA.shape[0] < arrayB.shape[0]: return False
    for b in arrayB:
        contained = False
        for a in arrayA:
            if np.all(a == b): contained = True
        if not contained:
            return False
    return True

def all_possible_clusters(n, lower=2, upper=None):
    if upper: assert upper <= n
    else: upper=n-1
    ploidies = np.array([1]*n)
    #make first set of cluster pairs
    clusters = []
    size = 0
    while True:
        size += 1
        found_one = False
        for new_cluster in itertools.combinations(all_possible_patterns(n, lower, upper), size):
            fail=False
            new_cluster = np.array(new_cluster)
            for i in range(1,size):
                for j in range(i):
                    if not sticcs.passes_derived_gamete_test(new_cluster[i],new_cluster[j], ploidies):
                        fail = True
                        break
                if fail:
                    break
            
            if not fail:
                found_one = True
                #check if an old cluster is contained within the new one
                clusters = [cluster for cluster in clusters if not contains_array(new_cluster, cluster)] + [new_cluster]
        
        if not found_one:
            return clusters

def list_all_topologies(n):
    assert n <= 7, "There are 135135 rooted topologies with 8 tips (660032 including bifurcations). I doubt you want to list all of those."
    clusters = all_possible_clusters(n)
    topos = [sticcs.patterns_to_tree(cluster) for cluster in clusters]
    return topos

def make_topoDict(n, unrooted=False):
    topos = []
    topoIDs = []
    for topo in list_all_topologies(n):
        topoSummary = TopologySummary(topo)
        ID = topoSummary.get_topology_ID(list(range(n)), unrooted=unrooted)
        if ID not in topoIDs: #this is because multiple topos are the same when unrooted, so we just use the first
            topoIDs.append(ID)
            topos.append(topo)
    
    return {"topos":topos, "topoIDs":topoIDs}

def make_topoDict_tskit(n, include_polytomies=False):
    topos = list_all_topologies_tskit(n, include_polytomies=include_polytomies)
    ranks = [t.rank() for t in topos]
    return {"topos":topos, "ranks":ranks}


class Topocounts:
    def __init__(self, topos, counts, totals=None, intervals=None, label_dict=None):
        assert counts.shape[1] == len(topos)
        if intervals is not None:
            assert intervals.shape[1] == 2
            assert counts.shape[0] == intervals.shape[0]
        self.topos = topos
        self.counts = counts
        if totals is not None:
            assert np.all(totals >= self.counts.sum(axis=1))
            self.totals = totals
        else: self.totals = counts.sum(axis=1)
        self.intervals = intervals
        self.label_dict = label_dict
    
    def write(self, outfile, include_topologies=True, include_header=True):
        nTopos = len(self.topos)
        if include_topologies:
            for x in range(nTopos): outfile.write("#topo" + str(x+1) + " " + self.topos[x].as_newick(node_labels=self.label_dict) + "\n") 
        if include_header:
            outfile.write("\t".join(["topo" + str(x+1) for x in range(nTopos)]) + "\tother\n")
        #write counts
        output_array = np.column_stack([self.counts, self.totals - self.counts.sum(axis=1)]).astype(str)
        outfile.write("\n".join(["\t".join(row) for row in output_array]) + "\n")
    
    def write_window_data(self, outfile, include_header=True):
        if include_header: outfile.write("chrom\tstart\tend\n")
        outfile.write("\n".join(["\t".join(["chr1", str(self.intervals[i,0]), str(self.intervals[i,1])]) for i in range(len(self.totals))]) + "\n")

def stack_topocounts(topocounts_list):
    #get all unique intervals
    unique_intervals = get_unique_intervals(np.row_stack([tc.intervals for tc in topocounts_list]))
    
    n_unique_intervals = len(unique_intervals)
    
    _counts_ = np.zeros((n_unique_intervals, topocounts_list[0].counts.shape[1]), dtype=int)
    _totals_ = np.zeros(n_unique_intervals, dtype=int)
    
    #for each iteration we check each of its intervals
    for i in range(len(topocounts_list)):
        print(".", end="", file=sys.stderr, flush=True)
        k=0
        intervals = topocounts_list[i].intervals
        max_k = intervals.shape[0]
        for j in range(n_unique_intervals):
            #for each of the unique intervals - if nested within the iteration's intervals, add to the stack
            if intervals[k,0] <= unique_intervals[j,0] and intervals[k,1] >= unique_intervals[j,1]:
                #k contains j
                _counts_[j,:] += topocounts_list[i].counts[k]
                _totals_[j] += topocounts_list[i].totals[k]
            
            if intervals[k,1] == unique_intervals[j,1]:
                #both are ending, so advance k to the next interval in this iteration
                k += 1
                if k == max_k: break
    
    return Topocounts(topocounts_list[0].topos, _counts_, _totals_, unique_intervals, topocounts_list[i].label_dict)


def get_topocounts_tskit(ts, groups=None, group_names=None, topoDict=None, include_polytomies=False):
    
    if groups is None:
        #use populations from ts
        assert list(ts.populations()) != [], "Either specify groups or provide a treesequence with embedded population data."
        if group_names is None: group_names = [str(pop.id) for pop in ts.populations()]
        groups = [[s for s in ts.samples() if str(ts.get_population(s)) == t] for t in taxonNames]
    
    ngroups = len(groups)
    label_dict = dict(zip(range(ngroups), group_names if group_names else ("group"+str(i) for i in range(1, ngroups+1))))
    
    if not topoDict:
        topoDict = make_topoDict_tskit(ngroups, include_polytomies=include_polytomies)
    
    topos = topoDict["topos"]
    ranks = topoDict["ranks"]
    
    intervals = np.array([tree.interval for tree in ts.trees()], dtype= int if ts.discrete_genome else float)
    
    #intervals[:,0] +=1 #convert back to 1-based - I decided to keep this out of the function, because it should be optional to adjust that
    
    counts = np.zeros((ts.num_trees, len(ranks)), dtype=int)
    
    totals = np.zeros(ts.num_trees, dtype=int)
    
    counter_generator = ts.count_topologies(groups)
    
    key = tuple(range(ngroups)) #we only want to count subtree topologies with a tip for each of the groups
    
    for i, counter in enumerate(counter_generator):
        counts[i] = [counter[key][rank] for rank in ranks]
        totals[i] = sum(counter[key].values())
    
    return Topocounts(topos, counts, totals, intervals, label_dict)


def get_topocounts(trees, n_trees, groups, max_iterations, simplify=True, group_names=None, topoDict=None, unrooted=False):
    
    ngroups = len(groups)
    label_dict = dict(zip(range(ngroups), group_names if group_names else ("group"+str(i) for i in range(1, ngroups+1))))
    
    if not topoDict:
        topoDict = make_topoDict(ngroups, unrooted)
    
    topos = topoDict["topos"]
    topoIDs = topoDict["topoIDs"]
    
    counts = np.zeros((n_trees, len(topoIDs)), dtype=int)
    
    totals = np.zeros(n_trees, dtype=int)
    
    intervals = []
    
    leafGroupDict = makeGroupDict(groups) if simplify else None
    
    for i,tree in enumerate(trees):
        topoSummary = TopologySummary(tree, leafGroupDict)
        counts_dict = topoSummary.get_topology_counts(groups, max_iterations, unrooted)
        counts[i] = [counts_dict[ID] for ID in topoIDs]
        totals[i] = sum(counts_dict.values())
        intervals.append(tree.interval)
    
    intervals = np.array(intervals, dtype=float)
    
    return Topocounts(topos, counts, totals, intervals, label_dict)


def get_topocounts_stacking_sticcs(der_counts, positions, ploidies, groups, max_iterations, group_names=None,
                                   unrooted=False, second_chances=False, multi_pass=True, chrom_start=None, chrom_len=None, silent=True):
    
    comboGenerator = comboGen(groups, max_iterations = max_iterations)
    
    topocounts_iterations = []
    
    for iteration,combo in enumerate(comboGenerator):
        
        print(f"\nInferring treesequence for iteration {iteration}.", file=sys.stderr, flush=True)
        
        der_counts_sub = der_counts[:, combo]
        
        #ploidies
        ploidies_sub = ploidies[(combo,)]
        
        #Find non-missing genoypes for these individuals
        no_missing = np.all(der_counts_sub >= 0, axis=1)
        
        site_sum = der_counts_sub.sum(axis=1)
        variable = (1 < site_sum) & (site_sum < sum(ploidies_sub))
        
        usable_sites = np.where(no_missing & variable)
        
        #der counts and positions
        der_counts_sub = der_counts_sub[usable_sites]
        
        positions_sub = positions[usable_sites]
        
        patterns, matches, n_matches  = sticcs.get_patterns_and_matches(der_counts_sub)
        
        clusters = sticcs.get_clusters(patterns, matches, positions_sub, ploidies=ploidies_sub, second_chances=second_chances,
                                       seq_start=chrom_start, seq_len=chrom_len, silent=silent)
        
        trees = sticcs.infer_trees(patterns, ploidies_sub, clusters, multi_pass = multi_pass, silent=silent)
        
        print(f"\nCounting topologies for iteration {iteration}.", file=sys.stderr, flush=True)
        
        topocounts_iterations.append(get_topocounts(trees, len(clusters), groups = make_numeric_groups(ploidies_sub),
                                                    group_names=group_names, max_iterations=1024, unrooted=unrooted))
    
    print(f"\nStacking", file=sys.stderr)
    
    return stack_topocounts(topocounts_iterations)

#def merge_identical_weights

    ##print("\nMerging adjacent identical intervals\n", file=sys.stderr)
    ##unique_intervals_merged = [unique_intervals[0][:]]
    ##stacks_merged = [stacks[0]]
    ##onePercent = int(np.ceil(n_unique_intervals/100))
    
    ##for j in range(1, n_unique_intervals):
        
        ##if j % onePercent == 0: print(".", end="", file=sys.stderr, flush=True)
        
        ##if np.all(stacks_merged[-1] == stacks[j]):
            ##unique_intervals_merged[-1][1] = unique_intervals[j][-1]
        ##else:
            ##unique_intervals_merged.append(unique_intervals[j][:])
            ##stacks_merged.append(stacks[j])

def makeGroupDict(groups, names=None):
    groupDict = {}
    for x in range(len(groups)):
        for y in groups[x]: groupDict[y] = x if not names else names[x]
    return groupDict


def make_numeric_groups(sizes):
    partitions = np.cumsum(sizes)[:-1]
    return np.split(np.arange(sum(sizes)), partitions)


def main():
    ### parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_vcf", help="Input VCF file with DC field", action = "store", required=True)
    parser.add_argument("-o", "--out_prefix", help="Output file prefix", action = "store", required=True)
    
    parser.add_argument("--unrooted", help="Unroot topologies (results in fewer topologies)", action = "store_true")
    
    parser.add_argument("--ploidy", help="Sample ploidy if all the same. Use --ploidy_file if samples differ.", action = "store", type=int)
    parser.add_argument("--ploidy_file", help="File with samples names and ploidy as columns", action = "store")
    
    parser.add_argument("--max_iterations", help="Maximum iterations over combintions of individuals", action = "store", type=int, required=True)
    
    parser.add_argument("--allow_second_chances", help="Consider SNPs that are separated by incompatible SNPs", action='store_true')
    parser.add_argument("--single_pass", help="Single pass when building trees (only relevant for ploidy > 1, but not recommended)", action='store_true')
    
    #parser.add_argument("--inputTopos", help="Input file for user-defined topologies (optional)", action = "store", required = False)
    parser.add_argument("--output_topos", help="Output file for topologies used", action = "store", required = False)
    parser.add_argument("--group_names", help="Name for each group (separated by spaces)", action='store', nargs="+", required = True)
    parser.add_argument("--groups", help="Sample IDs for each individual (separated by commas), for each group (separated by spaces)", action='store', nargs="+")
    parser.add_argument("--groups_file", help="Optional file with a column for sample ID and group", action = "store", required = False)
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    
    args = parser.parse_args()
    
    group_names = args.group_names
    
    ngroups = len(group_names)
    
    assert ngroups >= 3, "Please specify at least three groups."
    
    if args.groups:
        assert len(args.groups) == ngroups, "Number of groups does not much number of group names"
        groups = [g.split(",") for g in args.groups]
    
    elif args.groups_file:
        groups = [[]]*ngroups
        with open(args.groups_file, "rt") as gf: groupDict = dict([ln.split() for ln in gf.readlines()])
        for sample in groupDict.keys():
            try: groups[group_names.index(groupDict[sample])].append(sample)
            except: pass
    
    assert min([len(g) for g in groups]) >= 1, "Please specify at least one sample ID per group."
    
    sampleIDs = [ID for group in groups for ID in group]
    assert len(sampleIDs) == len(set(sampleIDs)), "Each sample should only be in one group."
    
    #VCF file
    vcf = cyvcf2.VCF(args.input_vcf, samples=sampleIDs)
    
    for ID in sampleIDs: assert ID in vcf.samples, f"ID {ID} not found in vcf file header line."
    
    #the dac vcf reader does not give us the dac values in the right order, so we need to record the sample indices
    sample_indices = [vcf.samples.index(ID) for ID in sampleIDs]
    
    ploidies, ploidyDict = sticcs.parsePloidyArgs(args, sampleIDs)
    
    ploidies = np.array(ploidies)
    
    dac_generator = sticcs.parse_vcf_with_DC_field(vcf)
    
    chromLenDict = dict(zip(vcf.seqnames, vcf.seqlens))
    
    #get all topologies
    #if args.inputTopos:
    
    label_dict = dict(zip(range(ngroups), group_names))
    
    topoDict = make_topoDict(ngroups, unrooted=args.unrooted)
    
    if args.output_topos:
        with open(args.output_topos, "wt") as tf:
            tf.write("\n".join([t.as_newick(node_labels=label_dict) for t in topoDict["topos"]]) + "\n")
    
    sys.stderr.write("\n".join([t.as_newick(node_labels=label_dict) for t in topoDict["topos"]]) + "\n")
    
    groups_numeric = make_numeric_groups([len(g) for g in groups])
    
    print(f"\nReading first chromosome...", file=sys.stderr)
    
    for chrom, positions, der_counts in dac_generator:
        
        assert positions.max() <= chromLenDict[chrom], f"\tSNP at position {positions.max()} exceeds chromosome length {chromLenDict[chrom]} for {chrom}."
        
        print(f"\nAnalysing {chrom}. {positions.shape[0]} usable SNPs found.", file=sys.stderr)
        
        #correct sample order so that groups are together
        der_counts = der_counts[:,sample_indices]
        
        topocounts_stacked = get_topocounts_stacking_sticcs(der_counts, positions, ploidies=ploidies, groups=groups_numeric,
                                                            group_names=group_names, max_iterations=args.max_iterations,
                                                            unrooted=args.unrooted, multi_pass=not args.single_pass,
                                                            second_chances = args.allow_second_chances,
                                                            chrom_start=1, chrom_len=chromLenDict[chrom],
                                                            silent= not args.verbose)
        
        #could potnentially add step merging identical intervals here
        
        with gzip.open(args.out_prefix + "." + chrom + ".topocounts.tsv.gz", "wt") as outfile:
            topocounts_stacked.write(outfile)
        
        with gzip.open(args.out_prefix + "." + chrom + ".windowdata.tsv.gz", "wt") as outfile:
            topocounts_stacked.write_window_data(outfile)
    
    print(f"\nEnd of file reached. Looks like my work is done.", file=sys.stderr)

if __name__ == '__main__':
    main()
