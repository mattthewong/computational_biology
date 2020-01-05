# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 21:24:00 2016

@author: Matt Wong
This program generates random sequences, allows for evolution of sequences 
with mu probability values, comparisons between sequences, and can simulate changes
over time with specific corrections such as Jukes-Cantor, purely observable changes, as 
well as observed vs. expected changes with parameters of sequence length, max # of changes, and stepwise
values.
"""
from __future__ import division
import numpy as np
from random import *
import math
import matplotlib.pyplot as plt
import ultrametric1 as ultra


def randomSeq(n):
#takes as input an integer representing the length of the sequence, 
#and returns a random string of As, Cs, Gs and Ts.
    sequence = ''
    bases = ["A","T","C","G"] #initialize a list of base options
    for i in range(n): #for every index in the range of the sequence
        index  = np.random.choice([0,1,2,3]) #randomly assign an index that is same length of bases.
        sequence += bases[index] #add that base.
        
    return sequence
    
def evolveSeq(seq,mu):
#returns a new sequence based on the expected number of changes per site 
#and the provided sequence.
    newsequence = ''
    
    bases = ["A","T","C","G"]
    
    for letter in range(len(seq)):#for index in the length of the sequence
        changes = np.random.poisson(lam = mu)#number of changes from a random selection on poisson dist.
        lettertochangeto = list(seq)[letter] #initialize a letter to change
        for x in range(changes):#for the total number of changes, select a letter other than the sequence[i] of that letter.
            if seq[letter] == 'A':
                index = np.random.choice([1,2,3])
                newletter = bases[index]
                lettertochangeto = newletter
                
            if seq[letter] == 'T':
                index = np.random.choice([0,2,3])
                newletter = bases[index]
                lettertochangeto = newletter
                
            if seq[letter] == 'C':
                index = np.random.choice([0,1,3])
                newletter = bases[index]
                lettertochangeto = newletter
                
            else:
                index = np.random.choice([0,1,2])
                newletter = bases[index] 
                lettertochangeto = newletter     
        newsequence += lettertochangeto #finally add that letter to the sequence
    return newsequence
        

def seqDiff(seq1,seq2):
#takes as input a starting sequence (a string of As, Cs, Gs and Ts) and an
# expected number of changes per site. For each site, it generates an
# actual number of changes as a random draw from a Poisson distribution
# with mean mu. 
    if len(seq1) != len(seq2): #if two sequences are not same length, cannot compare.
        return 'Sequences are not of the same length'
    else:
        seq1 = list(seq1) #intialize the sequence strings into lists.
        seq2 = list(seq2)
        differences = 0
        for index in range(len(seq1)): #loop through the indices of the sequence.
            if np.array_equal(seq1[index],seq2[index]) == False: #if the bases are not the same
                differences += 1 #increment the counter.
            else:
                differences = differences
        obsdiff = differences / len(seq1) #calculate the observed difference ratio.
        return obsdiff #return the observed difference ratio.

def jukesCantorModel(mu):
#takes as input an expected number of changes per site, 
#and returns the predicted number of observable changes 
#per site under the Jukes-Cantor model.
    observedchance = (3/4)*(1-(np.exp(-4/3 * mu))) #calculated observed chance.
    return observedchance

def plotSeqChanges(maxChanges, stepChanges, seqLength):
#makes a plot of observed vs. expected changes per site, 
#and compares this to the Jukes-Cantor prediction.takes arguments which
# specify the maximum number of expected changes on your plot (i.e., 
#the maximum value for the x axis), the step size between consecutive x values,
# and the length of the starting sequence. 
    changes = np.arange(0,maxChanges,stepChanges) #generates array of changes.
    changelist = [] 
    predicted = []
    observed = []
    for value in changes: #for each value in the change, append to list.
        changelist.append(value)
        
    for i in range(len(changelist)): #for each value in changeList
        change = changelist[i]
        mypredictedval = jukesCantorModel(change) #transform with Jukes Cantor
        predicted.append(mypredictedval) #add to list
        
    myseq = randomSeq(seqLength)
    for mu in changelist: #evolving the sequence and using purely seqDiff.
        evolved = evolveSeq(myseq,mu)
        difference = seqDiff(myseq,evolved)
        observed.append(difference)
        
    plt.plot(changelist,observed,'ro',changelist,predicted,'b-')#plot syntax.
    plt.xlabel("Observed Changes")
    plt.ylabel("Expected Changes")
    plt.show()

#PART TWO: Part 2 (50 points). Ultrametric trees on simulated data.

def simSeqOnTree(tree,seq,mu):
#takes as input a tree with branch lengths represented as time, 
#a starting sequence, and an expected number of changes per unit time. 
#The tree is represented in [(BranchLength, RootName), LeftTree, RightTree] 
#format, which is the same format returned by your computeBranchLengths function
# from last week. It then recursively simulates sequence evolution along each branch
# of the tree using evolveSeq, producing a sequence for each node in the tree.
# It returns a dictionary representing just the tips of the tree: each key is the 
#label for a tip, and the value is the simulated sequence for that tip.

    seqDict = {}
    simSeqOnTreeFunct(tree,seq,mu,seqDict) #call wrapper
    return seqDict

def simSeqOnTreeFunct(tree,seq,mu,dic):
#wrapper function for simSeqOnTree, taking in a provided dictionary as empty 
#to avoid recursion issues.
    if tree[1] == []: #if tree is a leaf
        myseq = evolveSeq(seq,mu) #first evolve the sequence
        dic[tree[0][1]] = myseq #the stored species and sequence added to dict.
    else:
        simSeqOnTreeFunct(tree[1],seq,mu,dic) #recurseive calls to either side of tree.
        simSeqOnTreeFunct(tree[2],seq,mu,dic)

def rawDistanceMatrix(tipList, seqDict):
# It constructs a raw (uncorrected) distance matrix using seqDiff to compare 
#each pair of species in tipList.
    matrix = [] #initialize matrix
    for member in tipList: #loop through the species
        myseq = seqDict[member]
        row = []
        for other in tipList:
            otherseq = seqDict[other]
            difference = seqDiff(myseq,otherseq) #calculates difference
            row.append(difference)  #appends difference
            
        matrix.append(row)
    return matrix #returns uncorrected matrix.

def correctDistanceMatrix(distMatrix):
#takes as input a raw (uncorrected) distance matrix (a list of lists) and creates
# a new matrix with the Jukes-Cantor correction applied to each entry. 
    expectedmatrix = []
    for index in distMatrix:
        row = []
        for mymu in index: #looping though each value in matrix
            newlambda = (-(3/4))*math.log1p((-(4/3)*mymu)) #defining lambda
            row.append(newlambda) #appending new corrected value
        expectedmatrix.append(row)
    return expectedmatrix #returning new matrix with JC applied.
    
def branchLengthsAsTime(tree,mu):
#It converts the original branch lengths, represented as average number of changes
# per site (i.e. as produced by a distance matrix), into units of time. 
    if tree[1] == []: #if it's a leaf, return the tree
        return tree
    else:
        timevalue = tree[0] / mu #timevalue calculation
        return [timevalue, branchLengthsAsTime(tree[1],mu),branchLengthsAsTime(tree[2],mu)] #recusion to both sides.

def simulateAndRecoverTrees(tree, mu, seqLength):
#takes as input a tree with branch lengths represented in units of time, 
#a rate of change (per site, per unit time), and a sequence length. The input
# tree should be represented in [(BranchLength, RootName), Left, Right] format.
# The function generates a random sequence of the specified length using randomSeq,
# evolves sequences along that tree using simSeqOnTree, takes the resulting set of 
#sequences and calculates a distance matrix using rawDistanceMatrix and then corrects
# it using correctDistanceMatrix. Finally, it creates ultrametric trees from both the
# raw and corrected distance matrices, converts the branch lengths back to units of
# time using branchLengthsAsTime, and returns the resulting trees as Newick trees.

    randomseq = randomSeq(seqLength) #first generate random sequence
    newdict = simSeqOnTree(tree,randomseq,mu) #create new dictionary
    species = list(newdict.keys()) #define species
    rawmatrix = rawDistanceMatrix(species,newdict) #generate raw matrix
    corrected = correctDistanceMatrix(rawmatrix) #generate corrected matrix
    
    ultrarawmatrix = ultra.ultrametrify(rawmatrix)[0] #ultrametrify raw matrix
    ultracorrectedmatrix = ultra.ultrametrify(corrected)[0] #ultrametrify corrected matrix
    
    ultrarawtree = ultra.buildTree(species,ultrarawmatrix) #build each tree
    ultracorrectedtree =ultra.buildTree(species,ultracorrectedmatrix)
    
    rawastime = branchLengthsAsTime(ultrarawtree,mu) #turn each tree into branched lengths as time
    correctedastime = branchLengthsAsTime(ultracorrectedtree,mu)
    
    
    return ultra.toNewickBranchLengths(rawastime), ultra.toNewickBranchLengths(correctedastime) #return them in Newick formatt.


    
    
        
        
    
        
    
            
    
        
    
    


    
    
    
    
        
    
        


    

            
        
    
    
    
    
