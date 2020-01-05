# ultrametric.py
# Ran Libeskind-Hadas
# Last revised: March 2016
# 
# Program to build ultrametric trees using the recursive algorithm described
# by Gusfield (Algorithms on Trees, Strings, and Sequences), p.450
# This program also implements the "correction" algorithm that converts a
# non-ultrametric matrix to an ultrametric one using the minimum "downward"
# correction.

# ---------------------------------------------------------------------------

# START TEST DATA

L0 = ["A", "B", "C"]
M0 = [ [0, 3, 3],
       [3, 0, 2],
       [3, 2, 0]]

# Species list and distance matrix from Gusfield p. 450 (for testing purposes)
L1 = [ "A", "B", "C", "D", "E"]
M1 = [ [0, 8, 8, 5, 3],
       [8, 0, 3, 8, 8],
       [8, 3, 0, 8, 8],
       [5, 8, 8, 0, 5],
       [3, 8, 8, 5, 0]]

# Non-ultrametric data from homework
L2non =["Orang", "Chimp", "Neander", "San", "Yoruba", "Han", "HumRef"]
M2non =[[0.000000, 0.170199, 0.158900, 0.159732, 0.159246, 0.159845, 0.159257],
        [0.170199, 0.000000, 0.090661, 0.092981, 0.092281, 0.092203, 0.091658],
        [0.158900, 0.090661, 0.000000, 0.013081, 0.012099, 0.012282, 0.012528],
        [0.159732, 0.092981, 0.013081, 0.000000, 0.005936, 0.005996, 0.005876],
        [0.159246, 0.092281, 0.012099, 0.005936, 0.000000, 0.002357, 0.002478],
        [0.159845, 0.092203, 0.012282, 0.005996, 0.002357, 0.000000, 0.002175],
        [0.159257, 0.091658, 0.012528, 0.005876, 0.002478, 0.002175, 0.000000]]


L3non = ["Bovine", "Mouse", "Gibbon", "Orang", "Gorilla", "Chimp", "Human"]
M3non = [ [0.0000, 1.6866, 1.7198, 1.6606, 1.5243, 1.6043, 1.5905],
       [1.6866, 0.0000, 1.5232, 1.4841, 1.4465, 1.4389, 1.4629],
       [1.7198, 1.5232, 0.0000, 0.7115, 0.5958, 0.6179, 0.5583],
       [1.6606, 1.4842, 0.7115, 0.0000, 0.4631, 0.5061, 0.4710],
       [1.5243, 1.4465, 0.5958, 0.4631, 0.0000, 0.3484, 0.3083],
       [1.6043, 1.4389, 0.6179, 0.5061, 0.3484, 0.0000, 0.2692],
       [1.5905, 1.4629, 0.5583, 0.4710, 0.3083, 0.2692, 0.0000]]

# END TEST DATA

# ---------------------------------------------------------------------------

# START OF CODE FOR THE CASE THAT THE MATRIX IS ULTRAMETRIC.
# THIS PART OF HTE PROGRAM DETERMINES IF A MATRIX IS ULTRAMETRIC AND
# CONSTRUCTS THE ULTRAMETRIC TREE.

# WE USE THE FOLLOWING TREE FORMAT:  A TREE IS EITHER A TRIPLET
# [ROOT, LEFT, RIGHT] WHERE ROOT IS A DISTANCE AND LEFT AND RIGHT ARE 
# THEMSELVES SUBTREES OR THE TREE IS EMPTY [].  LATER, WE CONVERT OUR
# TREES INTO NEWICK FORMAT WITH BRANCH LENGTHS.

# Determine whether or not a tree is ultrametric
def ultrametric(Matrix):
    ''' returns True if the matrix is ultrametric and False otherwise. '''
    size = len(Matrix)
    for i in range(size):
        for j in range(size):
            for k in range(size):
                if i != j and i != k and j != k:
                    dists = [Matrix[i][j], Matrix[i][k], Matrix[j][k]]
                    maxdist = max(dists)
                    if dists.count(maxdist) != 2: return False
    return True

# Returns the index of the row of the distance matrix containing a 
# maximum element
def findMaxRowIndex(Matrix):
    return 0 # The maximum is always found in every row!

# Takes a list of indices (each between 0 and N-1), a list of N species,
# and a NxN distance matrix as input and returns BOTH a sublist and a 
# submatrix.  The sublist is the subset of the species list at the given 
# indices and the submatrix is the submatrix of distances between these 
# species.
def subMatrix(indexList, speciesList, Matrix):
    size = len(Matrix)
    newList = []
    newMatrix = []
    for row in indexList:  # build new species list and all rows of new matrix
        newList.append(speciesList[row])
        newRow = []
        for col in range(size):
            if col in indexList: newRow.append(Matrix[row][col])
        newMatrix.append(newRow)
    return newList, newMatrix
                
# Build the ultrametric tree for the given species list and distance matrix. 
# Returns an ultrametric tree.
def buildTree(speciesList, Matrix):
    size = len(speciesList) # this is N, the # of species and dim of Matrix
    if size == 1:
        return [speciesList[0], [], []]
    else:
        # The maxRowIndex is the index of the species furthest from 
        # all other species.
        maxRowIndex = findMaxRowIndex(Matrix) 

        # maxRow is the row of the matrix containing this maximum.
        maxRow = Matrix[maxRowIndex] 

        # Next, we wish to partition the remaining species into equivalence
        # classes with respect to distance from our maxRowIndex species.
        # We do this using a dictionary called distanceGroups where
        # the key is a distance and the value is the list of all indices 
        # (species) at this given distance.
        # Now we build this dictionary...
        distanceGroups = {}
        for col in range(size):
            dist = Matrix[maxRowIndex][col]
            if dist != 0:
                if not dist in list(distanceGroups.keys()): 
                    distanceGroups[dist] = [col]
                else: distanceGroups[dist] = distanceGroups[dist] + [col]

        # Next, we use recursion to build a list called subTreesWithDist
        # that contains tuples of the form (dist, tree) where tree is the 
        # subtree for all species in the same distance equivalence class and 
        # dist is that distance.
        subTreesWithDist = []
        for dist in list(distanceGroups.keys()):
            indexList = distanceGroups[dist]
            subS, subM = subMatrix(indexList, speciesList, Matrix)
            tree = buildTree(subS, subM)
            subTreesWithDist.append((dist,tree))
        
        # Next, we get all of the distances in descending order 
        distances = list(distanceGroups.keys())
        distances.sort()
        distances.reverse()
        
        # Similarly, we sort the list of subtrees by that distance
        subTreesWithDist.sort()
        subTreesWithDist.reverse()
        return assemble(distances, subTreesWithDist, speciesList[maxRowIndex])

# Assemble the subtrees into a single phylogenetic tree.  
def assemble(distances, subTrees, leftmostTip):
    if len(distances) == 1:
        return [distances[0], [leftmostTip, [], []], subTrees[0][1]]
    else:
        leftSubTree = subTrees[0][1] # remove distance annotation
        rightSubTree = assemble(distances[1:], subTrees[1:], leftmostTip)
        return [distances[0], leftSubTree, rightSubTree]

# END OF CODE FOR THE CASE THAT TREE IS ULTRAMETRIC        

# ---------------------------------------------------------------------------

# START OF MST-BASED CORRECTION ALGORITHM THAT TAKES A NON-ULTRAMETRIC MATRIX
# AND "CORRECTS" IT INTO AN ULTRAMETRIC ONE.

# Given a distance matrix, returns another matrix of the same size representing
# a MST for that Matrix.  This function uses Prim's algorithm.
def MST(Matrix):
    size = len(Matrix)
    MSTMatrix = []
    for row in range(size):
        MSTMatrix.append([0]*size)
    visitedList = [0]
    for iteration in range(size-1):
        fromV, toV = nearest(Matrix, visitedList)
        MSTMatrix[fromV][toV] = Matrix[fromV][toV]
        MSTMatrix[toV][fromV] = Matrix[toV][fromV]
        visitedList.append(toV)
    return MSTMatrix

# Helper function for Prim's algorithm.        
# Returns an ordered pair (fromV, toV) where from is in visitedList
# and to is not in visitedList and the edge Matrix[fromV][toV] is the
# minimum of all such edge costs.
def nearest(Matrix, visitedList):
    size = len(Matrix)
    bestCost = float('inf')
    for fromV in visitedList:
        for toV in range(size):
            if toV not in visitedList:
                if Matrix[fromV][toV] < bestCost:
                    bestCost = Matrix[fromV][toV]
                    bestFromV = fromV
                    bestToV = toV
    return bestFromV, bestToV

# Given the adjacency Matrix for the MST, a start vertex, an end vertex, and
# a list of visited vertices (initially empty), returns the maximum weight
# of an edge on the unique path from the start vertex to the end vertex.
def largestEdge(Matrix, start, end, visited):
    visited += [start]
    if start == end: return 0  # no edge on this path - so we're done!
    # If we've wandered down a path and gotten stuck (this path is "blocked")
    # return infinity to indicate that this path didn't find the end vertex
    if blocked(Matrix, start, visited):
        return float('inf')
    else:
        size = len(Matrix)
        result = None
        # The candidates list will contain the maximum edge weight on a 
        # path from start to end AND some infinities that were contributed 
        # by dead-end paths.  Ultimately, we'll filter out the infinities 
        # by asking for min(candidates) which will simply leave us with the 
        # actual maximum edge weight.
        candidates = []
        for neighbor in range(size):
            if Matrix[start][neighbor] > 0 and neighbor not in visited:
                maxSubPath = largestEdge(Matrix, neighbor, end,  visited)
                candidates = candidates + [max(Matrix[start][neighbor],
                                               maxSubPath)]
        return min(candidates)
                
# Given a MST adjacency matrix, the start vertex and a list of
# visited vertices, returns False if the start vertex has unvisited neighbors
# (our search path is _not_ blocked) and True otherwise (this path is blocked).
def blocked(Matrix, start, visited):
    size = len(Matrix)
    for neighbor in range(size):
        if Matrix[start][neighbor] > 0 and neighbor not in visited:
            return False
    return True

# Given a distance matrix, return the ultrametrified matrix and the max
# change of any entry during the conversion (max change in range from 0 to 100)
def ultrametrify(Matrix):
    MSTMatrix = MST(Matrix)
    size = len(Matrix)
    newMatrix = [[0]*size for rows in range(size)]
    maxChange = 0.0
    for row in range(size):
        for col in range(size):
            maxedge = largestEdge(MSTMatrix, row, col, [])
            newMatrix[row][col] = maxedge
        if Matrix[row][col] != 0:
	        change = 100.0*(Matrix[row][col]-maxedge)/Matrix[row][col]
	        if change > maxChange:
                 maxChange = change
        
    return newMatrix, maxChange


# END OF ULTRAMETRIC CORRECTION STEP
# ---------------------------------------------------------------------------

# START OF MAIN PROGRAM ALONG WITH CONVERSION OF TREE TO NEWICK FORMAT

# This program represents trees as [] or [Root, Left, Right] where Root is
# a distance and Left and Right are subtrees.  The next three functions
# are used to convert from this format to newick.

# Takes a tree in [Root, Left, Right] format and returns a newick cladogram
# (no branch lengths) as a tuple.
def toNewick(Tree):
    if len(Tree[1]) == 0: return Tree[0]
    else:
        return (toNewick(Tree[1]), toNewick(Tree[2]))

# Takes a tree in [Root, Left, Right] format as input and returns a
# newick tree with branch lengths.  This function first calls 
# computeBranchLengths to construct a new tree in the form 
# [(BranchLength, Root), Left, Right].  Then, this function calls 
# toNewickBranchLengthHelper on that annotated tree to construct
# the newick version.
def toNewickBranchLengths(Tree):
    Tannotated = computeBranchLengths(Tree[0], Tree)
    result = toNewickBranchLengthsHelper(Tannotated)
    return result.rstrip(":0.0")
    
# Takes a tree in [(BranchLength, Root), Left, Right] format (as computed
# by the helper function computeBranchLengths below) and returns
# a newick tree with branchlengths as a string.  This newick tree has a
# trailing :0 at the end which can be sliced off later if desired.
def toNewickBranchLengthsHelper(annotatedTree):
    branchLength = str(annotatedTree[0][0])
    rootName = str(annotatedTree[0][1])
    if len(annotatedTree[2]) == 0 : return rootName + ":" + branchLength
    else:
        return "("+ toNewickBranchLengthsHelper(annotatedTree[1]) + "," + \
            toNewickBranchLengthsHelper(annotatedTree[2]) + ")" + ":" + \
            branchLength

# Take the root of the tree and the tree in [Root, Left, Right] format
# and returns a new tree of the form [(BranchLength, Root), Left, Right]
# where BranchLength is the length of the branch entering the root.
# This is used as a helper function for toNewickBranchLengths.
def computeBranchLengths(parentValue, Tree):
    # branch length and node
    if len(Tree[1]) == 0: return [(0.5*parentValue, Tree[0]), [], []] 

    else:
        return [ ( 0.5*(parentValue-Tree[0]), Tree[0] ), 
                    computeBranchLengths(Tree[0], Tree[1]), 
                    computeBranchLengths(Tree[0], Tree[2])]

# This is the main program that takes a species list and a distance matrix
# as input.  The program tests if the matrix is ultrametric.  If so, it 
# returns a newick tree with branch lengths.  If not, the distance matrix
# is first "corrected" and then returns the newick tree for the corrected
# matrix.                                   
def main(speciesList, Matrix):
    if ultrametric(Matrix):
        print("Matrix is ultrametric")
        T = buildTree(speciesList, Matrix)
        print(toNewickBranchLengths(T))
    else:
        print("Matrix is not ultrametric")
        newMatrix, maxChange = ultrametrify(Matrix)
        print(" New Matrix is: ")
        print(newMatrix)
        print("Maximum change: ", maxChange, "%")
        T = buildTree(speciesList, newMatrix)
        print(toNewickBranchLengths(T))
        
# END OF MAIN PROGRAM
# ---------------------------------------------------------------------------


