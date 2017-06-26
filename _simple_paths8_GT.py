#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#    Copyright (C) 2012 by
#    Sergio Nery Simoes <sergionery@gmail.com>
#	 <Graph-tool implementation> Pedro Queirós <pdqueiros@gmail.com>
#    All rights reserved.
#    BSD license (see second license in LICENSE.txt).

import graph_tool as _gt
import graph_tool.topology as _gt_tp

__author__ = """\n""".join(['Sérgio Nery Simões <sergionery@gmail.com>',
                             'Aric Hagberg <aric.hagberg@gmail.com>',
							 'Pedro Queirós <pdqueiros@gmail>'])
__all__ = ['all_circular_paths_edges', 'all_circular_paths_edges']


def all_simple_paths_edges(G, source, target, cutoff=None, data=False):
    all_paths=[]
    if data == True:
        edge_data = lambda u,v,n,E,I: (u,v,dict([E[n][I[n]]]))
    else:
        edge_data = lambda u,v,n,E,I: (u,v)
    for path in _gt_tp.all_paths(G,source,target,cutoff=cutoff):
        if path.tolist() not in all_paths:all_paths.append(path.tolist())    
    for path in all_paths:
        edges = list(zip(path[:-1],path[1:]))
        E = []  # list: items of each edge
        N = []  # list: number of items of each edge
        for u,v in edges:
            c=0
            edge_items=[]
            for parallel_edges in G.edge(u,v,all_edges=True):
                frag = G.edge_properties["frag_Property"][parallel_edges]
                iterProperty = G.edge_properties["iter_Property"][parallel_edges]
                weight = G.edge_properties["weight"][parallel_edges]                
                edge_items +=  [(c,{'frag':frag, 'weight': weight, 'i': iterProperty})]
                c+=1
            E += [edge_items]
            N += [len(edge_items)]
        I = [0 for n in N]
        idx = [i for i in reversed(list(range(len(I))))]
        while True:
            path_edges = []
            for n,(u,v) in enumerate(edges):
                path_edges += [edge_data(u,v,n,E,I)]
            yield path_edges
            for i in idx:
                I[i] = (I[i] + 1) % N[i]
                if I[i] != 0:
                    break
            if i == 0 and I[0] == 0:
                break
    
def all_circular_paths_edges(G):
    for pathArray in sorted(_gt_tp.all_circuits(G), key=len, reverse=True):
        path = list(pathArray)
        edges = list(zip(path, path[1:]+[path[0]]))
        N = []
        for u,v in edges:
            n=0
            for parallel_edges in G.edge(u,v,all_edges=True):   #this finds parallel edges
                n+=1
            N += [n]   
        I = [0 for n in N]
        idx = [i for i in reversed(list(range(len(I))))]
        while True:
            path_edges = []
            for i,(u,v) in enumerate(edges):
                frag = G.edge_properties["frag_Property"][G.edge(u,v)]
                iterProperty = G.edge_properties["iter_Property"][G.edge(u,v)]
                weight = G.edge_properties["weight"][G.edge(u,v)]
                path_edges += [(u,v,{'frag':frag, 'weight': weight, 'i': iterProperty})]
            yield path_edges
            for i in idx:
                I[i] = (I[i] + 1) % N[i]
                if I[i] != 0:
                    break
            if i == 0 and I[0] == 0:
                break
        

