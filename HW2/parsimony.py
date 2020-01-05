# -*- coding: utf-8 -*-

#Created on Fri Apr  8 16:16:44
#Matt Wong 
#this program determines maximum parsimony for positions multiple sequences in a tree, 
#entire trees, determining NNI's and performing heuristic searches.
import ultrametric1
import branchlengths
import random
def bestSinglePosition(tree, character, position, tipMapping, memo):
#returns the most parsimonious score (least cost) using memoization when passing a given character and position.
    if (tree, character, position) in memo: #if in memo, return cost.
        return memo[(tree, character, position)]
    
    elif tree[1] == ():#if tree itself is a leaf, return infinity.
        if tipMapping[tree[0]][position] != character:
            return float('inf')

        else:
            return 0
    else:
        baselist = ['A','T','C','G']
        leftscorelist = []
        rightscorelist = []
        for letter1 in baselist: #looping through potential letters for children on left side of tree.
            leftscore = bestSinglePosition(tree[1],letter1, position, tipMapping, memo) #recurse with memoization
            if letter1 != character:#if not equal, add a penalty of 1.
                leftscore += 1
            leftscorelist.append(leftscore)
               
        for letter2 in baselist: #looping through potential letters for children on right side of tree.
            rightscore = bestSinglePosition(tree[2],letter2, position, tipMapping, memo) #recurse with memoization
            if letter2 != character:#if not equal, add a penalty of 1.
                rightscore += 1
            rightscorelist.append(rightscore)    
        
        leftmin = min(leftscorelist)#access minimum cost of leftside
        rightmin = min(rightscorelist)#access minimum cost of rightside
        cost = leftmin + rightmin  
        memo[(tree, character, position)] = cost
        return cost #return cost.
        
            

def maxParsimony(tree, tipMapping):
#returns the most parsimoneous cost for a tree, evaluating all positions of the nodes.
    sumcost = []
    for sequence in list(tipMapping.values()):
        
        for index in range(len(sequence)):
            
            letter = list(sequence)[index]
            mycost = bestSinglePosition(tree,letter,index,tipMapping,{})
            
            sumcost.append(mycost)
    myset = set(sumcost)
    cost = sum(list(myset)[0:])     
    return cost 
    
def GenAllTrees(tree):
#generates all rerootings in newick format from the provided tree
    newtreelist = []
    if type(tree[0]) == tuple:
        return [tree]
    
    else:
        tree1 = (tree[0][0],(tree[0][1],tree[1]))
        tree2 = (tree[0][1],(tree[0][0],tree[1]))
        tree3 = (tree[1][0],(tree[1][1],tree[0]))
        tree4 = (tree[1][1],(tree[1][0],tree[0]))
        newtreelist.extend(getOneDirection(tree1))
        newtreelist.extend(getOneDirection(tree2))
        newtreelist.extend(getOneDirection(tree3))
        newtreelist.extend(getOneDirection(tree4))
    return newtreelist

def getOneDirection(tree):
#forces the recursion to go in one direction
    if type(tree[0]) != tuple:
        return [tree]
    else:
        leftTree =(tree[0][0],(tree[0][1],tree[1]))
        rightTree = (tree[0][1], (tree[0][0],tree[1]))
        return [tree] + getOneDirection(leftTree) + getOneDirection(rightTree)
        
        
def NewicktoRLR(newicktree):
#changes the newick tree provided to a tree with root left right format.
    if type(newicktree) != tuple:
        return (newicktree,(),())
    else:
        return ('Anc',(NewicktoRLR(newicktree[0])),(NewicktoRLR(newicktree[1])))
        
        
def AllNNIs(tree):
#generates all NNI's in a list using the NNI newick tree helper function with a regular tree as input.
    newick = ultrametric1.toNewick(tree)
    NNIlist = NNIs(newick)
    RLRlist = []
    for NNI in NNIlist:
        newtree = NewicktoRLR(NNI)
        RLRlist.append(newtree)
    return RLRlist
        
        
def NNIs(newicktree):
#generates a list of nearst neighbor interchange trees utilizing genalltrees function with a newicktree as input.
    Rerootings = GenAllTrees(newicktree)
    NNIlist = []
    for tree in Rerootings:
       
        if type(tree[0]) == tuple and type(tree[1]) == tuple:
            NNI1 = ((tree[0][0],tree[1][0]),(tree[1][1],tree[0][1]))
            NNI2 = ((tree[0][0],tree[1][1]),(tree[0][1],tree[1][0]))
            NNIlist.append(NNI1)
            NNIlist.append(NNI2)
    return NNIlist
        
def initmatrix(tipMapping):
#initializes a matrix of size by size with zeros.Needed to generate a MST.
    matrix = []
    for sequence in list(tipMapping.values()):    
        tempList = [] #create rows
        for sequence1 in list(tipMapping.values()):    
            hammingdistance = branchlengths.seqDiff(sequence,sequence1) 
            tempList.append(hammingdistance)
        matrix.append(tempList)
    return matrix
      
def NNIheuristic(tipMapping, sampleSize):
#returns the optimal NNI and score associated when randomly sampling with a provided sample size.
    #mymatrix = initmatrix(tipMapping) #initialize matrix
    initial = initmatrix(tipMapping)
    ultramatrix = ultrametric1.ultrametrify(initial)[0]#ultrametrify matrix
    ultratree = ultrametric1.buildTree(list(tipMapping.keys()),ultramatrix)#use ultrametrix to make tree
    
    NNIset = [] #initialize NNI list
    scoreset = [] #keep track of scores
    
    newick = ultrametric1.toNewick(ultratree)
    RLR = NewicktoRLR(newick)
    maxpars = maxParsimony(RLR,tipMapping) #calculate max pars for the tree
    scoreset.append(maxpars) #add to scoreset
    NNIset.append(RLR) #add tree to NNIset (count it)
    NNIlist = AllNNIs(RLR) #generate NNI list
    subsetNNI = random.sample(set(NNIlist),sampleSize) #subset list with random.sample & samplesize.
    for i in range(len(list(subsetNNI))): #loop through each NNI
        score = maxParsimony(subsetNNI[i],tipMapping) #get score
        scoreset.append(score)#add score
        NNIset.append(subsetNNI[i])#add tree
    optimalscore = min(scoreset) #get min of scoreset
    indexofscore = scoreset.index(optimalscore) #locate the index
    optimalNNI = NNIset[indexofscore] #use the index to get the optimal tree
    return (optimalNNI,optimalscore) #return both

    
    
    
       
       

      
       
       
      
   
           
           
   
            
            
        
            
        
             
        
        
            
                    
                        
                     
                         
                    
                        
                
              
                
                    
                        
           
               
        
                        
                        
                    
                    
                   
            
                
                
                    
            
            
                    
            
            