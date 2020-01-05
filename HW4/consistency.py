# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:36:42 2016

@author: mattthewong
"""
import treedistance
import branchlengths as bl
import ultrametric1 as ultra
import parsimony as pars
import matplotlib.pyplot as plt
def genspecies(tree):
#generates  all of the species of the tree in a list
    if type(tree) != tuple:
        return [tree]
    else:
        return genspecies(tree[0]) + genspecies(tree[1])


    
def leafinTree(tree,leaf):
#returns a boolean depending on whether or not the leaf is in the tree or not
     if leaf in genspecies(tree):
         return True
     else:
         return False

def makespecial(tree,leaf):
#helper function which takes a "special leaf" and inserts it as an outer branch of an unrooted tree.
    if tree[0] == leaf:
        return tree
    elif tree[1] == leaf:
        return (tree[1],tree[0])
    else:
         right = tree[0]
         left = tree[1]
         if leafinTree(right,leaf):
             return makespecial((right[0],(right[1],left)),leaf)
         else:
             return makespecial((left[0],(left[1],right)),leaf)
         
def samespecies(tree1,tree2):
#utilizes the genspecies function to make sure the two trees have the same species.
 counter = 0
 specieslist1 = genspecies(tree1)
 specieslist2 = genspecies(tree2)
 for species in specieslist1:
     if species in specieslist2:
         counter = counter
     else:
         counter +=1
 if counter == 0:
     return True
 else:
     return False
     
     
def sameUnrooted(tree1, tree2):
 #takes two Newick cladograms (without branch lengths) and returns true
 # if their unrooted topology is the same, and false otherwise.
 specieslist1 = genspecies(tree1)
 specieslist2 = genspecies(tree2)
 if samespecies(tree1,tree2):
     for index in range(len(specieslist1)):
         specialspecies = specieslist1[index]
         if specialspecies in specieslist2:
             tree1 = makespecial(tree1,specialspecies)
             tree2 = makespecial(tree2,specialspecies)
             if tree1[0] == specialspecies and tree2[0] == specialspecies:
                 if treedistance.distance(tree1,tree2) == 0:
                     return True
                 else:
                     return False
 return False
         
#PART 2
         
#def genspecies2(branchedtree):
#    if type(branchedtree) != list:
#        return [branchedtree]
#    else:
#        print(branchedtree)
#        return genspecies2(branchedtree[0]) + genspecies2(branchedtree[1]) 

def calcUltrametricConsistency(tree, mu, seqLengthList, jukes, reps):
#takes as input a tree with branch lengths represented in units of time, mu, a rate of
 #change (per site, per unit time), and a list of sequence lengths. It also takes a True/False
 #argument indicating whether the Jukes-Cantor correction should be applied. Finally, the reps
# argument indicates how many replicates will be run, i.e., how many times you will simulate
# and calculate the tree for each sequence length you have indicated. As in Assignment 2, the
# input tree should be represented in [(BranchLength, RootName), Left, Right] format. The output
# of this function is a list of numbers, the same length as seqLengthList, which indicate the number
# of times (out of reps attempts) that the ultrametric algorithm returned the true starting tree.
    accuracylist = []
    removedtree = removebranchlengths(tree)
    for length in seqLengthList:
        counter = 0
        for rep in range(reps):
            
            if jukes == True:
                newick = bl.simulateAndRecoverTrees(tree,mu,length)[1]
            
                outcome = sameUnrooted(removedtree,newick)
                if outcome == True:
                    counter += 1
                else:
                    counter = counter
            else:
                newick = bl.simulateAndRecoverTrees(tree,mu,length)[0]
                outcome = sameUnrooted(removedtree,newick)
                if outcome == True:
                    counter += 1
                else:
                    counter = counter
        accuracylist.append(counter)
    return accuracylist
                
def plotConsistency(seqList,consistencyList):
    
    plt.plot(seqList,consistencyList)
    plt.xlabel("Sequence Length")
    plt.ylabel("Number of Correct Trees out of 75")
    plt.show()
                
def removebranchlengths(tree):
#removes the branch lengths of a normal tree and returns the newick tree.
    if tree[1] == []:
        return tree[0][1]
    else:
        
        return (removebranchlengths(tree[1]),removebranchlengths(tree[2]))             
                
        
            
def calcParsimonyConsistency(tree, mu, seqLengthList, sampleSize, reps):
#takes as input a tree with branch lengths represented in units of time, mu, 
#a rate of change (per site, per unit time), and a list of sequence lengths. 
#It also takes an argument indicating the sample size, for the number of nearest
# neighbors to sample during the heuristic search. Finally, the reps argument indicates
# how many replicates will be run, i.e., how many times you will simulate and calculate
#the tree for each sequence length you have indicated. As in Assignment 2, the input tree
# should be represented in [(BranchLength, RootName), Left, Right] format. The output of 
#this function is a list of numbers, the same length as seqLengthList, which indicate the
# number of times (out of reps attempts) that the maximum parsimony algorithm returned 
#the true starting tree.
    accuracylist = []
    removedtree = removebranchlengths(tree)
    for length in seqLengthList:
        counter = 0
        sequence = bl.randomSeq(length)
        for rep in range(reps):
            mydict = bl.simSeqOnTree(tree,sequence,mu)
            best = pars.NNIheuristic(mydict,sampleSize)
            outcome = sameUnrooted(removedtree,best)
            if outcome == True:
                counter += 1
            else:
                counter = counter
        accuracylist.append(counter)
    return accuracylist
             












