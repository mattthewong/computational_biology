# Tree distance
# Ran Libeskind-Hadas
# December 1, 2012

import random

BITS = 16  # Number of bits used in Robinson-Foulds fingerprints

# This program assumes unrooted trees, represented in Newick format.
# A tree is an ordered pair (A, B) where the comma represents a root node
# and A and B are the left and right subtrees, respectively.  We assume that in the
# main tree, A is a leaf (we use 0 in the test tree examples below).  Since the tree is
# not actually rooted (there should be no internal nodes of degree 2), we can
# imagine short-cutting the root of the tree and replacing its two edges by a single
# edge:
#
#     *                                                 0
#    / \     becomes...          alternatively...       |
#   0   @               0 -- @                          @
#      / \                  / \                        / \


# The program comprises three parts.  The first part computes the Robinson-Foulds
# distance between a pair of trees.  The second part does aggregate clustering of tree.
# The third part computes the consensus partitions and specificity, but DOES NOT actually
# compute the consensus trees!  

# ----- Test Data -----

# Test Trees from the file MoretTrees.pptx

# Moret 1, tree distance 1

test1a = (0, (1, (2, (3, 4))))
test1b = (0, (1, (3, (2, 4))))

# Moret 2, tree distance 3

test2a = (0, ((1, (2, 3)), (7, (6, (4, 5)))))
test2b = (0, ((2, (1, 3)), (6, (4, (5, 7)))))

# ----- End Test Data -----

# Part 1:  Robinson-Foulds distance

def labelLeaves(leafList, leafLabelDict):
    ''' Takes a list of leaf numbers and an empty dictionary as input and
        populates the dictionary (as a side effect) with randomly generated fingerprints
        whose binary representations have length specified by the global variable BITS.
        The leaf label dictionary has keys that are leaf names and labels that are
        the randomly generated fingerprints.'''
    for leaf in leafList:
        label = random.randint(0, 2**BITS-1)
        leafLabelDict[leaf] = label

def makeLeafList(Tree):
    ''' Takes a tree as input and returns the list of leaves in that tree.'''
    if type(Tree) != tuple:
        return [Tree]
    else:                   
        leftTree = Tree[0]
        rightTree = Tree[1]
        return makeLeafList(leftTree) + makeLeafList(rightTree)

def labelInternal(Tree, internalLabels, leafLabelDict):
    ''' Takes a tree, an initially empty dictionary for the internal node labels, and
        the dictionary of leaf labels as input and populates the internal node labels
        dictionary as a side effect.  The internal node label dictionary has
        keys that are binary fingerprints.  The values are actually not important,
        but we have the tree as the value.'''
    
    if type(Tree) != tuple:# It's a leaf, so just return its label
        return leafLabelDict[Tree]
    else:
        leftTree = Tree[0]
        rightTree = Tree[1]
        leftLabel = labelInternal(leftTree, internalLabels, leafLabelDict)
        rightLabel = labelInternal(rightTree, internalLabels, leafLabelDict)
        treeLabel = leftLabel ^ rightLabel
        internalLabels[treeLabel] = Tree
        return treeLabel

def distance(Tree1, Tree2):
    ''' Computes the unweighted Robinson-Foulds distance between two trees.'''
    leaves = makeLeafList(Tree1)  # the two trees have the same leaves
    leafLabelDict = {}
    labelLeaves(leaves, leafLabelDict)
    labels1 = {}
    labelInternal(Tree1, labels1, leafLabelDict)
    labels2 = {}
    labelInternal(Tree2, labels2, leafLabelDict)
    counter = 0
    for label in labels1:
        if not label in labels2:
            counter += 1
    return counter

# Part 2:  Clustering

def pairwiseDistances(TreeList, D):
    ''' Takes a list of trees and an empty dictionary D as input and populates the
        dictionary with the pairwise distances for all pairs of trees in the given list.
        This dictionary is used to improve the run-time performance of the clustering
        algorithm since it allows us to avoid recomputing distances.'''
    for t1 in TreeList:
        for t2 in TreeList:
            d = distance(t1, t2)
            D[(t1, t2)] = d
            D[(t2, t1)] = d
        
def cluster(TreeList, minClusters, maxClusters):
    ''' Partition TreeList into clusters.  Gives solutions for all cluster sizes
        between minClusters (>=1) and maxClusters.
        
        Uses agglomeration clustering.  Distance between partitions defined by
        closestPair function.

        Returns a list called Stats whose elements are of the form (numClusters, clusterData)
        where numClusters is the number of clusters in this partition and clusterData
        is itself a list of tuples, one for cluster, of the form (cluster size, diameter, specificity).
        where specificity is a measure of how resolved
        the strict consensus tree is.  It is defined as #edges in strict consensus tree/#edges
        in original (binary) tree'''
    D = {}
    pairwiseDistances(TreeList, D)
    Partition = [ [Tree] for Tree in TreeList ]  # initial singleton partitioning
    Stats = []
    while True:
        currentNumClusters = len(Partition)
        if currentNumClusters <= minClusters: break
        print("Clusters remaining: "), currentNumClusters
        if minClusters <= currentNumClusters <= maxClusters:
            clusterData = [(len(cluster), diameter(cluster, D), specificity(cluster, 1.0)) for \
                           cluster in Partition]
            print(clusterData)
            Stats.append((currentNumClusters, clusterData))
        cluster1, cluster2 = closestPair(Partition, D)
        Partition.remove(cluster1)
        Partition.remove(cluster2)
        Partition.append(cluster1 + cluster2)
    return Stats	

def diameter(TreeList, D):
    ''' Computes the diameter of the set of points (trees) in TreeList with respect
        to the pairwise distances between trees that have been computed and passed in as
        dictionary D.'''
    max = 0
    for tree1 in TreeList:
        for tree2 in TreeList:
            	d = D[(tree1, tree2)]
            	if d > max: max = d
    return max

def clusterDistance(Cluster1, Cluster2, D):
    ''' return the Agg2 distance between Clusters (see paper by Stockham et al)'''
    d = 0.0
    for tree1 in Cluster1:
        for tree2 in Cluster2:
            d += D[(tree1, tree2)]
    return d/(len(Cluster1)*len(Cluster2))

def closestPair(Partition, D):
    ''' Takes a partition of trees and a distance dictionary as input and returns the two sets
        in that partition that are nearest. '''
    minDistance = float('inf') # Infinity
    for cluster1 in Partition:
        for cluster2 in Partition:
            if cluster1 != cluster2:
                d = clusterDistance(cluster1, cluster2, D)
                if d < minDistance:
                    minDistance = d
                    C1 = cluster1
                    C2 = cluster2
    return C1, C2

# Part 3:  Consensus trees

def listOfPartitions(tree):
    ''' Takes a tree as input and returns a list of all partitions induced by removal
        of edges in that tree.  However, for each partition, only the one without the
        root vertex of the original tree is presented.  The other is implicit.'''
    if type(tree) != tuple:
        return [(tree,)]
    elif len(tree) == 2:
        P0 = makeLeafList(tree[0])
        P0.sort()
        P1 = makeLeafList(tree[1])
        P1.sort()
        P0 = tuple(P0)
        P1 = tuple(P1)
        return [P0, P1] + listOfPartitions(tree[0]) + listOfPartitions(tree[1])      

def partitions(TreeList):
    ''' Takes a list of trees as input and builds a partitionDict dictionary where
        keys are ordered pairs (A, B) where A and B
        are tuples of leaves induced by removal of an edge from a tree in TreeList
        and tuples are sorted with the smallest leaf number always in A.
        The value is the number of times this key (partition) is found in TreeList.'''
    partitionDict = {}
    numLeaves = len(makeLeafList(TreeList[0]))
    for tree in TreeList:
        lop = listOfPartitions(tree)
        nontrivialPartitions = [p for p in lop if numLeaves-1 > len(p) > 1] # throw out leaves
        for partition in nontrivialPartitions:
            if partition in partitionDict:
                partitionDict[partition] += 1
            else: partitionDict[partition] = 1
    return partitionDict

def consensus(TreeList, threshold):
    ''' Given a list of trees in TreeList and a threshold value (e.g. 0.5 for majority consensus,
        1.0 for strict consensus, or values in between), returns the list of consensus partitions for
        that TreeList.'''
    numTrees = len(TreeList)
    leafList = makeLeafList(TreeList[0]) # all trees have same set of leaves
    numLeaves = len(leafList)
    partitionDict = partitions(TreeList)
    consensusPartitions = []
    for partition in partitionDict:
        if 1.0 * partitionDict[partition]/numTrees >= threshold:
            consensusPartitions.append(partition)
    return consensusPartitions

def specificity(cluster, threshold):
    numEdgesStrictConsensusTree = len(consensus(cluster, threshold))
    numOriginalEdges = len(makeLeafList(cluster[0])) - 3 # all trees have same number of leaves
    return 1.0 * numEdgesStrictConsensusTree / numOriginalEdges

    





        
        
        
    
