# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 18:54:16 2016

@author: mattthewong
"""
#Matt Wong
#this program has a few tests for NNIheuristic, which offer the expected score
#from a provided catMapping and sampleSize
import parsimony

def test1(dictionary):

    dictionary = {"Matt": "AAA" , "Mike": "AAT", "Jasmin": "AAG"}
    expectedcost = 2
    expectedtree = ('Anc', ('Matt', (), ()), ('Anc', ('Jasmin', (), ()), ('Mike', (), ())))
    if parsimony.NNIheuristic(dictionary,0) == (expectedtree,expectedcost):
        return True
    else:
        return False

def test2(dictionary):
    dictionary = {"BRCA1": "ATGT" , "P53": "AATC", "RAS": "CGTA", "MIIA" : "TTCC"}
    expectedcost = 5
    expectedtree = (('Anc', ('RAS', (), ()), ('Anc', ('BRCA1', (), ()), ('Anc', ('P53', (), ()), ('MIIA', (), ())))))
    if parsimony.NNIheuristic(dictionary,0) == (expectedtree,expectedcost):
        return True
    else:
        return False

