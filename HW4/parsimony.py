# Phylogenetic Trees via Maximum Parsimony
# Ran Libeskind-Hadas
# April 2016

# This program uses the Maximum Parsimony Method and Nearest Neighbor
# Interchange to find a local optimum maximum parsimony phylogenetic tree
# with branch lengths.

# The input comprises a list of tip names, a list of characters, and
# a dictionary that associates each tip with a character.  The output is
# a newick tree with branch lengths, corrected using Jukes-Cantor
# correction.

import random
import math
import grasses_seq as gs
Infinity = float('inf')


# --------------------------------------------------------------------------
# Start test inputs

t1 = ["Linde", "Atwood", "West"]
m1 = {"Linde": "AAA" , "Atwood": "AAT", "West": "AAG"}

t2 = ["Groody", "Froody", "Snoody", "Wumph", "Glumph", "Humph"]
m2 = {"Groody": "AACA", "Froody": "AACG", "Snoody": "AAGC", "Wumph": "CCAA", "Glumph": "CTAA", "Humph": "AAAA"}

t3 = ["Alien", "Zombie", "Human"]
m3 = {"Alien": "GA", "Zombie": "GT", "Human": "CT"}

t4 = ["Western", "Eastern", "Northern", "Southern"]
m4 = {"Western": "TTAA", "Eastern": "AAAA", "Northern": "AACT", "Southern": "TACG"}

t5 = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
m5 = {"1":"AAAT", "2":"AAAG", "3":"AAAC", "4":"ATAG", "5":"ATAC", "6":"ATAT", "7":"GGCC", "8":"GGTA", "9":"GGTC", "10":"GGCT", "11":"AAAA", "12":"AAAC"}

Groodyspecies = ["Groody", "Froody", "Snoody", "Wumph", "Glumph", "Humph"]
Groodymapping = {"Groody": "AACC", "Froody": "AACG", "Snoody": "AAGC", "Wumph": "CCAA", "Glumph": "CTAA", "Humph": "AAAA"}       
      
# End test inputs
# --------------------------------------------------------------------------

# This is the best single position function.
# It takes a tree, a character (A, T, C, or G) for the root, a position,
# a tip mapping, and memo dictionary as input and returns the best cost
# of a maximum parsimony solution for the tree.
def bestSinglePosition(tree, character, position, tipMapping, memo):
    root = tree[0]
    leftTree = tree[1]
    rightTree = tree[2]
    if leftTree == ():  # Tip
        if character == tipMapping[root][position]:
            return 0
        else:
            return Infinity
    else:
        if (tree, character, position) in memo: return memo[(tree, character, position)]
        bestSolution = Infinity
        for leftChar in ["A", "T", "C", "G"]:
            for rightChar in ["A", "T", "C", "G"]:
                bestLeft = bestSinglePosition(leftTree, leftChar, position, tipMapping, memo)
                bestRight= bestSinglePosition(rightTree, rightChar, position, tipMapping, memo)
                cost = bestLeft + bestRight
                if leftChar != character: cost += 1
                if rightChar != character: cost += 1
                if cost < bestSolution:
                    bestSolution = cost
        memo[(tree, character, position)] = bestSolution
    return bestSolution 

# This function calls best (above) four times, once for each possible
# symbol at the root node...
# For a given tree, position, and tipMapping, returns the
# solution of minimum cost for this tree and this position.
# The solution is of the form (cost, tree) where tree
# has all internal nodes changed to the character in the optimal
# solution.
def bestChar(tree, position, tipMapping):
    memo = {}
    bestSolution = Infinity
    for root in ["A", "T", "C", "G"]:
        option = bestSinglePosition(tree, root, position, tipMapping, memo)
        if option < bestSolution:
            bestSolution = option
    return bestSolution
        
# This function scores a tree by looping over each possible position in
# the tip sequence data, finding the best solution for each position,
# merging those separate trees into one tree, and then returns a tuple
# of the form (treeCost, bestTree) where treeCost is the total parsimony
# cost of that tree and bestTree is that tree with the full ancestral sequences
# at each internal node.
def maxParsimony(tree, tipMapping):
    keys = list(tipMapping.keys())
    positions = len(tipMapping[keys[0]])
    score = 0 
    for position in range(positions):
        score += bestChar(tree, position, tipMapping)
    return score

# -----------------------------------------------------------------------     
# NNI Heuristic

# Generate a single very lopsided RLR tree with the given set of tips.
# This can be used as the initial tree in our algorithm.
def oneTree(tipList):
    if len(tipList) == 1: return (tipList[0], (), ())
    else: return ("Anc", (tipList[0], (), ()), oneTree(tipList[1:]))

# Generate an initial tree by first shuffling the tipList and then calling
# the oneTree function above to build a tree with that shuffled tipList.
def randomTree(tipList):
    random.shuffle(tipList)
    return oneTree(tipList)

# The main NNI heuristic
def NNIheuristic(tipMapping, sampleSize):
    tipList = list(tipMapping.keys())
    bestTree = randomTree(tipList)
    bestScore = maxParsimony(bestTree,tipMapping)
    #print "Initial tree score ", bestScore
    while True:
        neighbors = allNNIs(bestTree)
        if sampleSize > len(neighbors): sampleSize = len(neighbors)
        sampleNeighbors = random.sample(neighbors, sampleSize)
        if len(sampleNeighbors) == 0: break
        scoredNeighbors = []
        for neighbor in sampleNeighbors:
            nScore  = maxParsimony(neighbor, tipMapping)
            scoredNeighbors.append((nScore, neighbor))
        scoredNeighbors.sort()
        if scoredNeighbors[0][0] >= bestScore: break
        bestScore = scoredNeighbors[0][0]
        #print " Best score so far ", bestScore
        bestTree = scoredNeighbors[0][1]
    #print "Best tree ", bestTree
    newickTree = RLRtoNewick(bestTree)
    #print "Final score: ", bestScore
    #print "Final newick tree: "
    return newickTree, bestScore

# Takes a root-left-right tree as input and returns all root-left-right
# trees obtained making the tree unrooted, re-rooting the tree at each
# internal branch, and computing the NNIs for that re-rooting.
def allNNIs(RLRtree):
    newickVersion = RLRtoNewick(RLRtree)
    allReroots = genAllTrees(newickVersion)
    newickNNIs = []
    RLRNNIs = []
    for tree in allReroots: 
        newickNNIs.extend(NNIs(tree))
    for tree in newickNNIs:
        RLRNNIs = list(map(lambda tree: newickToRLR(tree), newickNNIs))
    return RLRNNIs

# Takes a newick tree (A, B) as input.  Implicitly, that tree can be
# viewed as rooted on the branch between A and B.  This function returns
# a list of all possible re-rootings of this same tree, where the root
# moves to the left of the original root.
def genLeft(newickTree):
    left = newickTree[0]
    right = newickTree[1]
    if type(left) != tuple: return []  # We can't push forward left!
    else:
        leftleft = left[0]
        leftright = left[1]
        tree1 = (leftleft, (leftright, right))
        tree2 = (leftright, (leftleft, right))
        return [tree1, tree2] + genLeft(tree1) + genLeft(tree2)

# Takes a newick tree (A, B) as input.  Implicitly, that tree can be
# viewed as rooted on the branch between A and B.  This function returns
# a list of all possible re-rootings of this same tree, where the root
# moves to the right of the original root.
def genRight(newickTree):
    left = newickTree[0]
    right = newickTree[1]
    if type(right) != tuple: return [] # We can't push forward right!
    else:
        rightleft = right[0]
        rightright = right[1]
        tree1 = ((left, rightleft), rightright)
        tree2 = ((left, rightright), rightleft)
        return [tree1, tree2] + genRight(tree1) + genRight(tree2)

# Takes a newick tree (A, B) as input and returns a list of all possible
# re-rootings of this tree.
def genAllTrees(newickTree):
    return genLeft(newickTree) + genRight(newickTree)

# Take a newick tree ((A, B). (C, D)) as input and returns a list of the two
# NNI neighbors of that tree.
def NNIs(newickTree):
    candidates = filter(lambda tree: type(tree[0]) == tuple and \
                            type(tree[1]) == tuple, genAllTrees(newickTree))
    NNIlist = []
    for tree in candidates:
        A = tree[0][0]
        B = tree[0][1]
        C = tree[1][0]
        D = tree[1][1]
        NNIlist.append( ((A, C), (B, D)) )
        NNIlist.append( ((A, D), (C, B)) )
    return NNIlist

# Convert a root-left-right tree to plain newick format.
def RLRtoNewick(RLRTree):
    root = RLRTree[0]
    left = RLRTree[1]
    right = RLRTree[2]
    if len(left) == 0: return root
    else: return (RLRtoNewick(left), RLRtoNewick(right))

# Convert a plain newick tree to root-left-right format where the internal
# nodes are all named "Anc" for "Ancestor".
def newickToRLR(newickTree):
    if type(newickTree) != tuple: return (newickTree, (), ())
    else: 
        return ("Anc", newickToRLR(newickTree[0]), newickToRLR(newickTree[1]))
    

            
    
             
    
    


    
        
            
