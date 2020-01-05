# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 16:09:12 2016

@author: mattthewong
"""
#This program defines fingerprints for species of a tree as a binary sequence specified of length BITS.
#The program also computes fingerprints for interal nodes of a tree based on XORing the fingerprints of leafs.
#The program calculates distances of trees based on the Robinson Foulds distance, 
import random
from Caesal import *
import copy
import numpy as np
BITS = 64

def xor(string1,string2):
#xors the strings and outputs the exored string!
    xorstring = ""
    if len(string1) != len(string2):
        return "strings not same length!"
        
    for value in range(len(string1)):
        if list(string1)[value] == list(string2)[value]:
            xorstring += "0"
        else:
            xorstring += "1"
    return xorstring
    
def fingerprint(Tree,mydict):
#creates a dicitonary that associates every species in a newick tree with a 
#randomly created sequence of length BITS.
  if type(Tree)== int:
    randomint = random.randint(0,99999999999999999999999999999999999999999)
    mydict[Tree] = randomint
    return mydict
  else:
    fingerprint(Tree[0],mydict)
    fingerprint(Tree[1],mydict)
    return mydict

def exoredges(Tree,mydict,edgedict):
#populates an edge dictionary, where the binary sequences of the edges are stored.
    if type(Tree)== int:
      return mydict[Tree]
    else:
      FP1 = exoredges(Tree[0],mydict,edgedict)
      FP2 = exoredges(Tree[1],mydict,edgedict)
      exor = FP1 ^ FP2
      edgedict[exor] = "Blah"
      return exor

def distance(Tree1,Tree2):
#takes in two trees, and computes the distance between the two using exor.
  distance = 0
  newdict1 = {}
  newdict2 = {}
  fingerprinttree1 = fingerprint(Tree1,{})
  fingerprinttree2 = fingerprinttree1
  exoredges(Tree1,fingerprinttree1,newdict1)
  exoredges(Tree2,fingerprinttree2,newdict2)
  internalvalues1 = list(newdict1.keys())
  internalvalues2 = list(newdict2.keys())
  for index in range(len(internalvalues1)):
    if internalvalues1[index] not in newdict2:
        distance += 1
    else:
        distance = distance
  return distance



def initialize(treeList):
#returns a dictionary that is filled with clusters of trees (values) and distances between them (keys)
  distdict = {}
  for cluster in treeList:
    for cluster1 in treeList:
      if type(cluster) == tuple and type(cluster1) == tuple:
        if cluster != cluster1:
          mydist = distance(cluster,cluster1)
          distdict[mydist] = [cluster,cluster1]
      elif type(cluster) == list and type(cluster1) == tuple:
          if cluster != cluster1:
            sumtot = 0
            for tree in cluster:
              pairdist = distance(tree,cluster1)
              sumtot += pairdist
            truedistance = sumtot / (len(cluster)*len(cluster1))
            distdict[truedistance] = [cluster,cluster1]
      elif type(cluster) == tuple and type(cluster1) == list:
        if cluster != cluster1:
          sumtot = 0
          for tree in cluster1:
            pairdist = distance(tree,cluster1)
            sumtot += pairdist
          truedistance = sumtot / (len(cluster)*len(cluster1))
          distdict[truedistance] = [cluster,cluster1]
      else:
        if cluster != cluster1:
          for tree in cluster:
            for tree1 in cluster1:
              sumtot = 0
              if tree != tree1:
                pairdist = distance(tree,tree1)
                sumtot += pairdist
          truedistance = sumtot / (len(cluster)*len(cluster1))
          distdict[truedistance] = [cluster,cluster1]
  return distdict


def decrimentbyone(treeList,distdict):
#decriments the size of treeList by one by determining the minimum distance 
#and removing those individual trees and adding a 2-tree cluster.
  minclusterdistance = min(list(distdict.keys()))
  cluster1 = distdict[minclusterdistance][0]
  cluster2 = distdict[minclusterdistance][1]
  index1 = treeList.index(cluster1)
  index2 = treeList.index(cluster2)
  if index2 > index1:
    del treeList[index2]
    del treeList[index1]
  else:
    del treeList[index1]
    del treeList[index2]
  treeList.append([cluster1,cluster2])
  return treeList

def cluster(treeList, minClusters, maxClusters):
#for every clustersize within the range of minClusters and maxClusters, 
#a list is generated which provides the size of the cluster, length of each cluster, 
#diameter of each cluster, and average diameter of each cluster.
  dist = 0
  clusterdata = []
  for clustersize in range(minClusters,maxClusters):
    tempList = copy.copy(treeList)
    while len(tempList) > clustersize:
      distdict = initialize(tempList)
      tempList = decrimentbyone(tempList,distdict)
    data = [clustersize]
    for cluster in tempList:
      data.append(analyze(cluster))
    clusterdata.append(data)
  return clusterdata

def analyze(cluster):
#for every cluster in the treeList, the diameter (max distance of all pairwise comparisons) 
#and average diameter is calculated.
  distancelist = []
  if type(cluster) == tuple:
    diameter = 0
    average = 0
    length = 1
    return [1,0,0]
  else:
    for tree1 in cluster:
      for tree2 in cluster:
        if tree1 != tree2:
          mydistance = distance(tree1,tree2)
          distancelist.append(mydistance)
    diameter = max(distancelist)
    average = sum(distancelist[0:len(distancelist)]) / len(distancelist)
    length = len(cluster)
  return [length,diameter,average]

    

def bipartitions(tree):
#determines the bipartitions for a given tree
  treebipartitions = []
  edgedict = {}
  speciesdict = fingerprint(tree,{})
  speciesvals = list(speciesdict.values())
  xor = exoredges(tree,speciesdict,edgedict)
  edgeseqs = list(edgedict.keys())
  for sequence in speciesvals:
      for sequence2 in speciesvals:
         tempspeciesvals = copy.copy(speciesvals)
         XOR = sequence ^ sequence2
         if XOR in edgedict:
            index1 = speciesvals.index(sequence)
            index2 = speciesvals.index(sequence2)
            species1 = speciesvals[index1]
            species2 = speciesvals[index2]
            left = [species1,species2]
            if index1 > index2:
                del tempspeciesvals[index1]
                del tempspeciesvals[index2]
            else:
                del tempspeciesvals[index2]
                del tempspeciesvals[index1]
            left.append(tuple(tempspeciesvals))
            treebipartitions.append(left)


#==============================================================================
# 
# def changetobinary(tree):
# #takes in a tree and converts the created  bipartitions to binary associations.
#     newpartition = []
#     treebipartitions = bipartitions(tree)
#     for partition in treebipartitions:
#         binary = []
#         for speciesset in partition:
#             if type(speciesset) == int:
#                 binary.append(1)
#             else:
#                 for species in speciesset:
#                     binary.append(0)
#         newpartition.append(binary)
#     return newpartition
#==============================================================================
            
def consensuspartitionhelper(cluster,threshold):
# takes as input a cluster (a list of phylogenetic trees) and a threshold 
#and returns all of the consensus trees in a list that meet the threshold value.
  totalpartitions = []
  
  for tree in cluster:  
      mybipartitions = bipartitions(tree)
      
      for partition in mybipartitions:
        totalpartitions.append(partition)
  
  frequencydict = frequency(totalpartitions)
  print(len(totalpartitions),len(list(frequencydict.keys())))
  frequencies = list(frequencydict.values())
  partitionoptions = list(frequencydict.keys())
  consensuspartitionlist = []
  for myfrequency in frequencies:
      if myfrequency >= threshold:
          index = frequencies.index(myfrequency)
          consensuspartition = partitionoptions[index]
          consensuspartitionlist.append(consensuspartition)
  

#def specificity(cluster, threshold):

def frequency(totalpartitions):
#calculates the frequency of each binary partition
    frequencydict = {}
    for element in totalpartitions:
        if tuple(element) in list(frequencydict.keys()):
            frequencydict[tuple(element)] += 1
        else:
            frequencydict[tuple(element)] = 1
    for frequency in list(frequencydict.values()):
        frequency = frequency / len(totalpartitions)
    return frequencydict



def test1():
#testing the distance function on two trees in Caesal.
  tree1 = treeList[0]
  tree2 = treeList[42]
  return distance(tree1,tree2)

def test2():
  return initialize(treeList,{})

def test3():
#testing the cluster algorithm on a smaller dataset.
  return cluster(smallList,2,10)

def test4():
  mydict = initialize(smallList)
  return decrimentbyone(smallList,mydict)

def test5():
  return initialize(treeList,{})

def test6():
    return consensuspartitionhelper(smallerlist,0.1)









  
        
    
    
    
    
  
                
            