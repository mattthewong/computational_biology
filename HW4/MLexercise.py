# Phylogenetic Trees via Maximum Likelihood
# MCB 118b class exercise
# April 2016

import math
# --------------------------------------------------------------------------
# Start test inputs

tips = {"Linde": "AAA" , "Atwood": "AAA", "West": "AAT"}
tree1 = ( (0, "Anc"),
              ( (1.0, "Linde"), (), ()),
              ( (1.0, "Anc"),
                ((1.0, "Atwood"), (), ()),
                ((1.0, "West"), (), ())
              )
        )

tree2 = ( (0, "Anc"),
              ( (0.1, "Linde"), (), ()),
              ( (0.1, "Anc"),
                ((0.1, "Atwood"), (), ()),
                ((1.0, "West"), (), ())
              )
        )

tree3 = ( (0, "Anc"),
              ( (0.01, "Linde"), (), ()),
              ( (0.01, "Anc"),
                ((0.01, "Atwood"), (), ()),
                ((1.0, "West"), (), ())
              )
        )

tree4 = ( (0, "Anc"),
              ( (0.01, "Linde"), (), ()),
              ( (0.01, "Anc"),
                ((0.01, "Atwood"), (), ()),
                ((10.0, "West"), (), ())
              )
        )

# End test inputs
# --------------------------------------------------------------------------


def jukesCantor(branchLength):
    return 3.0/4.0*(1 - math.exp(-4.0/3.0 * branchLength))

# This function takes a root-left-right tree as input where
# root is of the form (branchLength, Name), a character for the root,
# the position/site index, tipMapping dictionary that associates a species
# name to its sequence.  It returns the likelihood for that input.
# This implementation is not memoized.
def singlePosition(tree, rootChar, position, tipMapping):
    root = tree[0]
    branchLength = tree[0][0]
    rootName = tree[0][1]
    leftTree = tree[1]
    rightTree = tree[2]
    if len(leftTree) == 0:  # Tip
        if rootChar == tipMapping[rootName][position]:
            return 1.0 
        else:
            return 0.0
    else:
        score = 0.0
        leftBranchLength = leftTree[0][0]
        rightBranchLength = rightTree[0][0]
        
        for leftChar in ["A", "T", "C", "G"]:
            for rightChar in ["A", "T", "C", "G"]:
                leftScore = singlePosition(leftTree, leftChar, position, tipMapping)
                probObservedChange = jukesCantor(leftBranchLength)
                if rootChar == leftChar:
                    leftScore *= (1.0-probObservedChange)
                else:
                    leftScore *= (1.0/3.0)*probObservedChange
            

                rightScore = singlePosition(rightTree, rightChar, position, tipMapping)
                probObservedChange = jukesCantor(rightBranchLength)
                if rootChar == rightChar:
                    rightScore *=  (1.0-probObservedChange)
                else:
                    rightScore *= (1.0/3.0)*probObservedChange
                score += rightScore*leftScore
                    
        return score

def allRootChars(tree, position, tipMapping):
    score = 0.0
    for rootChar in ["A", "T", "C", "G"]:
        score += 0.25 * singlePosition(tree, rootChar, position, tipMapping)
    return score
        
def ML(tree, tipMapping):
    keys = list(tipMapping.keys())
    positions = len(tipMapping[keys[0]])
    score = 1.0
    for position in range(positions):
        score *= allRootChars(tree,position,tipMapping)
    return score

            
    
             
    
    


    
        
            
