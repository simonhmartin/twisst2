import numpy as np
from collections import defaultdict
from collections import deque
import itertools, random

class NodeChain(deque):
    def __init__(self, nodeList, dists=None):
        super(NodeChain, self).__init__(nodeList)
        if dists is None: self.dists = None
        else:
            assert len(dists) == len(self)-1, "incorrect number of iternode distances"
            self.dists = deque(dists)
        self._set_ = None
    
    def addNode(self, name, dist=0):
        self.append(name)
        if self.dists is not None: self.dists.append(dist)
    
    def addNodeLeft(self, name, dist=0):
        self.appendleft(name)
        if self.dists is not None: self.dists.appendleft(dist)
    
    def addNodeChain(self, chainToAdd, joinDist=0):
        self.extend(chainToAdd)
        if self.dists is not None:
            assert chainToAdd.dists is not None, "Cannot add a chain without distances to one with distances"
            self.dists.append(joinDist)
            self.dists.extend(chainToAdd.dists)
    
    def addNodeChainLeft(self, chainToAdd, joinDist=0):
        self.extendleft(chainToAdd)
        if self.dists is not None:
            assert chainToAdd.dists is not None, "Cannot add a chain without distances to one with distances"
            self.dists.appendleft(joinDist)
            self.dists.extendleft(chainToAdd.dists)
    
    def chopLeft(self):
        self.popleft()
        if self.dists is not None: self.dists.popleft()
    
    def chop(self):
        self.pop()
        if self.dists is not None: self.dists.pop()
    
    def fuseLeft(self, chainToFuse):
        new = NodeChain(self, self.dists)
        assert new[0] == chainToFuse[0], "No common nodes"
        i = 1
        while new[1] == chainToFuse[i]:
            new.chopLeft()
            i += 1
        m = len(chainToFuse)
        while i < m:
            new.addNodeLeft(chainToFuse[i], chainToFuse.dists[i-1] if self.dists is not None else None)
            i += 1
        return new
    
    def simplifyToEnds(self, newDist=None):
        if self.dists is not None:
            if not newDist: newDist = sum(self.dists)
            self.dists.clear()
        leftNode = self.popleft()
        rightNode = self.pop()
        self.clear()
        self.append(leftNode)
        self.append(rightNode)
        if self.dists is not None:
            self.dists.append(newDist)
    
    def setSet(self):
        self._set_ = set(self)


def getChainsToLeaves(tree, node=None, simplifyDict = None):
    if node is None: node = tree.root
    children = tree.children(node) #this syntax for tskit trees
    #children = tree.node_children[node]
    if len(children) == 0:
        #if it has no children is is a child
        #if it's in the simplifyDict or there is not simplifyDict
        #just record a weight for the node and return is as a new 1-node chain
        if simplifyDict is None or node in simplifyDict:
            chain = NodeChain([node])
            setattr(chain, "weight", 1)
            return [chain]
        else:
            return []
    #otherwise get chains for all children
    childrenChains = [getChainsToLeaves(tree, child, simplifyDict) for child in children]
    #now we have the chains from all children, we need to add the current node
    for childChains in childrenChains:
        for chain in childChains: chain.addNodeLeft(node)
    
    #if collapsing, check groups for each node
    if simplifyDict:
        nodeGroupsAll = np.array([simplifyDict[chain[-1]] for childChains in childrenChains for chain in childChains])
        nodeGroups = list(set(nodeGroupsAll))
        nGroups = len(nodeGroups)
        
        if (nGroups == 1 and len(nodeGroupsAll) > 1) or (nGroups == 2 and len(nodeGroupsAll) > 2):
            #all chains end in a leaf from one or two groups, so we can simplify.
            #first list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
            #Start by getting index of each chain for each group
            indices = [(nodeGroupsAll == group).nonzero()[0] for group in nodeGroups]
            #the new weight for each chain we keep will be the total node weight of all from each group 
            newWeights = [sum([chains[i].weight for i in idx]) for idx in indices]
            #now reduce to just a chain for each group 
            chains = [chains[idx[0]] for idx in indices]
            for j in range(nGroups):
                chains[j].simplifyToEnds()
                chains[j].weight = newWeights[j]
        
        #if we couldn't simply collapse completely, we might still be able to merge down a side branch
        #Side branches are child chains ending in a single leaf
        #If there is a lower level child branch that is itself a side branch, we can merge to it
        elif (len(childrenChains) == 2 and
            ((len(childrenChains[0]) == 1 and len(childrenChains[1]) > 1) or
            (len(childrenChains[1]) == 1 and len(childrenChains[0]) > 1))):
            chains,sideChain = (childrenChains[1],childrenChains[0][0]) if len(childrenChains[0]) == 1 else (childrenChains[0],childrenChains[1][0])
            #now check if any main chain is suitable (should be length 3, and the only one that is such. and have correct group
            targets = (np.array([len(chain) for chain in chains]) == 3).nonzero()[0]
            if len(targets) == 1 and simplifyDict[chains[targets[0]][-1]] == simplifyDict[sideChain[-1]]:
                #we have found a suitable internal chain to merge to
                targetChain = chains[targets[0]]
                newWeight = targetChain.weight + sideChain.weight
                targetChain.simplifyToEnds()
                targetChain.weight = newWeight
            else:
                #if we didn't find a suitable match, just add side chain
                chains.append(sideChain)
        else:
            #if there was no side chain, just list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
    #otherwise we are not collapsing, so just list all chains
    else:
        chains = [chain for childChains in childrenChains for chain in childChains]
    #now we have the chains from all children, we need to add the current node
    
    return chains


def comboGen(groups, max_iterations):
    n_combos = np.prod([len(t) for t in groups])
    if n_combos <= max_iterations:
        for combo in itertools.product(*groups):
            yield combo
    else:
        for i in range(max_iterations):
            yield tuple(random.choice(group) for group in groups)

def pairsDisjoint(pairOfPairs):
    if pairOfPairs[0][0] in pairOfPairs[1] or pairOfPairs[0][1] in pairOfPairs[1]: return False
    return True


pairs_generic = dict([(n, tuple(itertools.combinations(range(n),2))) for n in range(4,8)])

pairPairs_generic = dict([(n, [pairPair for pairPair in itertools.combinations(pairs_generic[n],2) if pairsDisjoint(pairPair)],) for n in range(4,8)])


class TopologySummary:
    #Summarises a tree as a set of NodeChain objects
    #Provides a convenient way to summarise the topology as a unique ID
    #If group data is provided, tree is only summarised in terms of chains
    # between nodes from distinct groups, and redundant chains (e.g. monophyletic clades)
    # are simplified and can be weighted accordingly with recorded leaf weights
    
    def __init__(self, tree, simplifyDict=None):
        self.root = tree.root
        self.simplifyDict = simplifyDict
        chains = getChainsToLeaves(tree, simplifyDict=self.simplifyDict)
        self.leafWeights = {}
        self.chains = {}
        #record each root-leaf chain with th
        for chain in chains:
            self.chains[(chain[0], chain[-1])] = chain
            self.leafWeights[chain[-1]] = chain.weight
        
        self.leavesRetained = set(self.leafWeights.keys())
        #now make chains for all pairs of leaves (or all pairs in separate groups if defined)
        if self.simplifyDict:
            leafPairs = [pair for pair in itertools.combinations(self.leavesRetained,2) if self.simplifyDict[pair[0]] != self.simplifyDict[pair[1]]]
        else:
            leafPairs = [pair for pair in itertools.combinations(self.leavesRetained,2)]
        
        for l0,l1 in leafPairs:
            self.chains[(l0,l1)] = self.chains[(tree.root,l0)].fuseLeft(self.chains[(tree.root,l1)])
        
        #add a reversed entry for each chain, and set the set for each one
        for pair in list(self.chains.keys()):
            self.chains[pair].setSet()
            self.chains[pair[::-1]] = self.chains[pair]
        
        #add weight 1 for root (do this now only because we don't want to include the root in the leavesRetained set)
        self.leafWeights[tree.root] = 1
    
    #def get_topology_ID(self, leaves, unrooted=False):
        #if leaves is None: leaves = sorted(self.leavesRetained)
        #if not unrooted: leaves = list(leaves) + [self.root]
        #n = len(leaves)
        #return tuple([self.chains[(leaves[pairs[0][0]],
                                   #leaves[pairs[0][1]],)]._set_.isdisjoint(self.chains[(leaves[pairs[1][0]],
                                                                                        #leaves[pairs[1][1]],)]._set_)
                                      #for pairs in pairPairs_generic[n]])
    
    def get_quartet_ID(self, quartet):
        if self.chains[(quartet[0],quartet[1],)]._set_.isdisjoint(self.chains[(quartet[2],quartet[3],)]._set_):
            return 0
        elif self.chains[(quartet[0],quartet[2],)]._set_.isdisjoint(self.chains[(quartet[1],quartet[3],)]._set_):
            return 1
        elif  self.chains[(quartet[0],quartet[3],)]._set_.isdisjoint(self.chains[(quartet[1],quartet[2],)]._set_):
            return 2
        return 3
    
    #def get_all_quartet_IDs(self, leaves=None, unrooted=False):
        #if leaves is None: leaves = sorted(self.leavesRetained)
        #if not unrooted: leaves = list(leaves) + [self.root]
        #return [self.get_quartet_ID(quartet) for quartet in itertools.combinations(leaves, 4)]
    
    def get_all_quartet_IDs(self, leaves=None, unrooted=False):
        if leaves is None: leaves = sorted(self.leavesRetained)
        if unrooted:
            return [self.get_quartet_ID(quartet) for quartet in itertools.combinations(leaves, 4)]
        else:
            return [self.get_quartet_ID(trio + (self.root,)) for trio in itertools.combinations(leaves, 3)]
    
    def get_topology_ID(self, leaves=None, unrooted=False):
        return tuple(self.get_all_quartet_IDs(leaves, unrooted))
    
    def get_topology_counts(self, groups, max_iterations, unrooted=False):
        _groups = [[leaf for leaf in group if leaf in self.leavesRetained] for group in groups]
        
        if self.simplifyDict is not None:
            nCombos = np.prod([len(g) for g in _groups])
            assert nCombos < max_iterations, "Tree simplification must be turned off when considering only a subset of combinations."
        
        combos = comboGen(_groups, max_iterations)
        counts = defaultdict(int)
        for combo in combos:
            comboWeight = np.prod([self.leafWeights[leaf] for leaf in combo])
            ID = self.get_topology_ID(combo, unrooted)
            counts[ID] += comboWeight
        
        return counts

def get_quartet_dist(tree1, tree2, unrooted=False, approximation_subset_size=None):
    topoSummary1 = TopologySummary(tree1)
    topoSummary2 = TopologySummary(tree2)
    
    if approximation_subset_size is None:
        quartetIDs1 = np.array(topoSummary1.get_all_quartet_IDs(unrooted=unrooted))
        quartetIDs2 = np.array(topoSummary2.get_all_quartet_IDs(unrooted=unrooted))
    else:
        #do approximate distance with random sets of quartets
        quartetIDs1 = np.zeros(approximation_subset_size, dtype=int)
        quartetIDs2 = np.zeros(approximation_subset_size, dtype=int)
        
        for i in range(approximation_subset_size):
            if unrooted:
                quartet = random.sample(topoSummary1.leavesRetained, 4)
                quartetIDs1[i] = topoSummary1.get_quartet_ID(quartet)
                quartetIDs2[i] = topoSummary2.get_quartet_ID(quartet)
            
            else:
                trio = random.sample(topoSummary1.leavesRetained, 3)
                quartetIDs1[i] = topoSummary1.get_quartet_ID(trio + [tree1.root])
                quartetIDs2[i] = topoSummary2.get_quartet_ID(trio + [tree2.root])
    
    dif = quartetIDs1 - quartetIDs2
    return np.mean(dif != 0)


def get_min_quartet_dist(tree1, tree2, inds, max_itr=10, unrooted=False):
    topoSummary1 = TopologySummary(tree1)
    topoSummary2 = TopologySummary(tree2)
    
    #quartet IDs for tree 1 are unchanging
    quartetIDs1 = np.array(topoSummary1.get_all_quartet_IDs(unrooted=unrooted))
    
    #for tree 2 we try with different permutations for each individual
    new_inds = inds[:]
    
    for itr in range(max_itr):
        
        for i in range(len(inds)):
            
            ind_orderings = list(itertools.permutations(inds[i]))
            
            dists = []
            
            for ind_ordering in ind_orderings:
                current_inds = new_inds[:]
                current_inds[i] = ind_ordering
                quartetIDs2 = np.array(topoSummary2.get_all_quartet_IDs(leaves=[i for ind in current_inds for i in ind], unrooted=unrooted))
                dif = quartetIDs1 - quartetIDs2
                dists.append(np.mean(dif != 0))
            
            new_inds[i] = ind_orderings[np.argmin(dists)]
        
        if itr > 0 and new_inds == previous_new_inds: break
        else:
            previous_new_inds = new_inds[:]
    
    #get final dist
    quartetIDs2 = np.array(topoSummary2.get_all_quartet_IDs(leaves=[i for ind in new_inds for i in ind], unrooted=unrooted))
    dif = quartetIDs1 - quartetIDs2
    
    return np.mean(dif != 0)

