#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 08:55:11 2020

@author: Abdel
"""
import numpy as np


#Returns the dynamic programming table of the alignment of S1 into S2. At position [i,j], the table
#indicates the best alignment possible as well as its metrics, for S1[:i] and S2[:j].
def align(alignement_type,S1,S2):
    
    m = len(S1)+1
    n = len(S2)+1
    
    def init_res():

        res = np.zeros((m,n), dtype=object)


        #Initialization
        if alignement_type == "global":
            res[0] = [ ("L",i) for i in range(n) ]
            res[:,0] = [("_",0)] + [ ("U",i) for i in range(m-1) ]


        elif alignement_type == "local":
            res[0] = [ ("_",0) for i in range(n) ]
            res[:,0] =[ ("_",0) for i in range(m) ]


        elif alignement_type == "motif_search":#aka approximative occurence of s1 in s2
            res[0] = [ ("_",0) for i in range(n) ]
            res[:,0] = [("_",0)] + [ ("U",i) for i in range(m-1) ]
        
        return res


    res = init_res()
    
    for i in range(1,m):
        for j in range(1,n):
            if alignement_type == "global" or alignement_type == "motif_search":
                #table de distance d'edition
                tmp = [res[i-1,j-1][1] + bool(S1[i-1]!=S2[j-1]), res[i-1,j][1] + 1, res[i,j-1][1] + 1 ]
                min_idx = np.argmin(tmp)
                
                #D: best path from diagonal, U: best path from above, L: best path from left 
                if min_idx == 0:
                    fleche = "D"
                elif min_idx == 1:
                    fleche = "U"
                else :
                    fleche = "L"
                
                
                res[i,j] =  fleche, tmp[min_idx]
                
            
            elif alignement_type == "local":
                #similarity table
                
                tmp2 = 2 if  S1[i-1]==S2[j-1] else -1
                
                tmp =  [ 0,res[i-1,j-1][1] +tmp2, res[i-1,j][1] - 2, res[i,j-1][1] - 2 ]
                max_idx = np.argmax(tmp)

                
                if max_idx == 0:
                    fleche = "_"
                elif max_idx == 1:
                    fleche = "D"
                elif max_idx == 2 :
                    fleche = "U"
                else:
                    fleche = "L"
                
                res[i,j] =  fleche, tmp[max_idx]
                
    
    return res



#Classes and functions to compute the automata found in https://link.springer.com/chapter/10.1007/3-540-60268-2_315

class PretreatmentAutomataNode(object):
    
    def __init__(self,colonne):
        self.colonne = colonne
        self.pointors = dict()
    
    
    def __eq__(self,other):
        return self.colonne == other.colonne
    
    def __str__(self):
        return str(self.colonne)
    
    def __getitem__(self,i):
        return self.pointors[i]
        

class PretreatementAutomata(object):
    
    def __init__(self, word, alphabet):
        
        initial_node = PretreatmentAutomataNode([i for i in range(len(word)+1)])
        
        nodes = [initial_node]
        done = [] # Nodes which pointors have been created
        
        while len(done) != len(nodes): # automata dont change in size:

            for node in nodes:
                #for each node in the automata, creates the pointors and add newly created nodes in automata
                if node not in done:
                        for letter in alphabet:
                            #create new node
                            new_node_colonne = [0]
                            for i in range(1,len(word)+1):
                                #dynamic programming to get the error vector components
                                new_node_colonne.append(min(  node.colonne[i-1] + bool(word[i-1] != letter),
                                                              node.colonne[i] + 1, 
                                                              new_node_colonne[i-1]+1))
                            
                            new_node =  PretreatmentAutomataNode(new_node_colonne)
                        
                            if new_node in nodes:
                                new_node = nodes[nodes.index(new_node)]
                            else:
                                nodes.append(new_node)
                            
                            node.pointors[letter] = new_node
                       
                        done.append(node)
        
        
        self.nodes = nodes
        
        
    
    def plot(self):
        import networkx as nx
        from networkx.drawing.nx_agraph import to_agraph 
        
        G=nx.DiGraph()
        
        #add nodes
        for idx,node in enumerate(self.nodes):  
            G.add_node(str(node), color='green' if idx==0 else 'black' )
        
        #add arcs
        for node in self.nodes:
            for letter,adjacent in node.pointors.items():
                G.add_edge(str(node), str(adjacent), label=letter )
                
        
        A = to_agraph(G)
        A.layout('dot')         
        A.draw('automata.png')                                                       
        
        
    



    
   
P = "abaa"
alphabet = ["a","b"]
automata = PretreatementAutomata(P,alphabet)
automata.plot()
        
        
        
        
        
                                
                           
                                
        
        
        
            
        
        
        
        
        
        



    
    

    
