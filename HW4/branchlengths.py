from __future__ import division
import numpy
import matplotlib.pyplot as plt
import ultrametric1 as ultra
import consistency
# Defines the possible characters at each site of a sequence
# With no argument, returns all of them
# With an argument consisting of one of them, returns all the rest (potential mutations)
def characters(existing=""):
    chars = ["A","T","C","G"]
    if existing in chars:
        chars.remove(existing)
    return chars

# Creates a random sequence using the characters defined in characters()
# Takes as input the sequence length n
# Returns a string of characters of length n
def randomSeq(n):
    seqList = [numpy.random.choice(characters()) for _ in range(n)]
    return "".join(seqList)

# Takes as input a sequence string and an expected number of changes mu per site
# Returns a new, changed sequence string (without modifying the original)
# The number of mutations is drawn from a Poisson distrubution with mean mu
def evolveSeq(seq, mu):
    seqList = [evolveSite(x, mu) for x in seq]
    return "".join(seqList)

# Takes as input a single character ch and an expected number of changes mu
# Draws a number of changes from a Poisson distribution with mean mu
# For each change, modifies the character to one of the other possibilities
# Returns the new, changed character
def evolveSite(ch, mu):
    nchanges = numpy.random.poisson(mu)
    for i in range(nchanges):
        ch = numpy.random.choice(characters(ch))
    return ch
        
# Takes two sequences and returns the percent difference between them, i.e.
# the number of changes divided by total sequence length
# Assumes they are of equal length
def seqDiff(seq1, seq2):
    if len(seq1) != len(seq2):
        print ("Sequences are not of the same length")
        return
    else:
        return sum(numpy.array(list(seq1))!=numpy.array(list(seq2)))/len(seq1)

# Takes the expected number of changes per site
# Returns the expected number of observed differences per site under the J-C model
def jukesCantorModel(mu):
    return 3/4 *(1 - numpy.exp(-4/3 * mu))

# Takes the number of observed differences per site
# Returns an estimate of the number of changes per site under the J-C model
def jukesCantorCorrection(observedDiff):
    return -3/4 * numpy.log(1 - 4/3 * observedDiff)

# Simulates sequence evolution over time: 
# plots the number of observed differences from start to end
# vs. the number of expected changes over that time period
def plotSeqChanges(maxChanges=3, stepChanges=0.01, seqLength=1000):
    mus = numpy.arange(0,maxChanges,stepChanges)
    s1 = randomSeq(seqLength)
    diffs = [seqDiff(s1,evolveSeq(s1,x)) for x in mus]
    predicted = [jukesCantorModel(x) for x in mus]
    plt.plot(mus,diffs,'ro',mus,predicted,'b-')
    plt.xlabel("Expected number of changes per site")
    plt.ylabel("Observed number of changes per site")
    plt.show()

# Example trees
sampleTree1 = [(0,""),[(9,""),[(1,"A"),[],[]],[(1,"B"),[],[]]],[(10,"C"),[],[]]]
sampleTree2 = [(0,""),[(1,""),[(1,"A"),[],[]],[(10,"B"),[],[]]],[(1,""),[(1,"C"),[],[]],[(10,"D"),[],[]]]]

# Simulates sequence evolution on a tree
# Takes as input a tree in [(BranchLength, Root name), Left tree, Right tree] format
# as well as a starting sequence and an expected number of changes per site per unit time
# Returns a dictionary with the tip labels as keys and simulated sequences as values
def simSeqOnTree(tree, seq, mu):
    newSeq = evolveSeq(seq, mu*tree[0][0])
    if len(tree[1])==0:
        return { tree[0][1] : newSeq }
    seqDict = simSeqOnTree(tree[1], newSeq, mu)
    seqDict.update(simSeqOnTree(tree[2], newSeq, mu))
    return seqDict
    
# Calculates the raw distance matrix for a set of sequences
# Takes as input a list of tip labels, and a dictionary mapping tip labels to sequences
# Returns a raw distance matrix (with rows and columns ordered by the tip label list)
def rawDistanceMatrix(tipList, seqDict):
    distances = []
    for i in range(len(tipList)):
        distances.append([])
        for j in range(len(tipList)):
            distances[i].append(seqDiff(seqDict[tipList[i]], seqDict[tipList[j]]))
    return distances

# Applies the Jukes Cantor correction to every element of a distance matrix
def correctDistanceMatrix(distMatrix):
    corrMatrix = jukesCantorCorrection(numpy.matrix(distMatrix))
    return corrMatrix.tolist()

# Takes a tree in [(BranchLength, Root name), Left tree, Right tree] format
# along with expected rate of change mu (per site, per unit time)
# and converts branch lengths (number of changes per site) to time 
# assuming a molecular clock with constant rate of change mu 
def branchLengthsAsTime(tree,mu):
    if len(tree[1])==0:
        return tree
    leftTree = branchLengthsAsTime(tree[1],mu)
    rightTree = branchLengthsAsTime(tree[2],mu)
    return [ tree[0]/mu, leftTree, rightTree ]

# Takes a list of species, a distance matrix, and a rate of change mu
# Returns an ultrametric tree in Newick format, with branch lengths in units of time
def buildUltrametricChronogram(tipList, distMatrix, mu):
    uDistMatrix, maxChange = ultra.ultrametrify(distMatrix)
    tree = ultra.buildTree(tipList, uDistMatrix)
    treeBLTime = branchLengthsAsTime(tree, mu)
    newick = ultra.toNewick(treeBLTime)
    return newick

# Simulates sequence evolution on the given tree
# Then calculates distance matrices (raw and J-C corrected)
# Then builds ultrametric trees for each and returns both
def simulateAndRecoverTrees(tree, mu, seqLength):
    seq = randomSeq(seqLength)
    seqDict = simSeqOnTree(tree, seq, mu)
    tipList = sorted(seqDict.keys())
    rawDist = rawDistanceMatrix(tipList, seqDict)
    corrDist = correctDistanceMatrix(rawDist)
    return buildUltrametricChronogram(tipList, rawDist, mu), buildUltrametricChronogram(tipList, corrDist, mu)
