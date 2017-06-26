#!/usr/bin/env python3
'''This module provides functions for assembly of sequences by homologous recombination and other
related techniques. Given a list of sequences (Dseqrecords), all sequences will be analyzed for
overlapping regions of DNA (common substrings).

The assembly algorithm is based on graph theory where each overlapping region forms a node and
sequences separating the overlapping regions form edges.

'''

import logging as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

import itertools as _itertools
import operator as _operator
import random as _random
import os as _os

from copy import copy as _copy
from collections import defaultdict as _defaultdict

from Bio.SeqFeature import FeatureLocation as _FeatureLocation
from Bio.SeqFeature import SeqFeature as _SeqFeature

#to reload packages
import pydna
from pydna.dseq  import Dseq as _Dseq
from pydna.dseqrecord import Dseqrecord as _Dseqrecord

from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
from pydna.common_sub_strings import terminal_overlap as _terminal_overlap
from pydna.contig import Contig as _Contig
from ordered_set   import OrderedSet as _OrderedSet
from pydna.utils   import memorize   as _memorize
import time

#####################################
#<Graph-tool implementation> Pedro Queirós <pdqueiros@gmail.com>
#gc= graph creation package, Graph-tool as gt whenever available or NetworkX as nx as a fallback

try:
    import graph_tool as _gt
    from pydna._simple_paths8_GT import all_simple_paths_edges   as _all_simple_paths_edges
    from pydna._simple_paths8_GT import all_circular_paths_edges as _all_circular_paths_edges
    gc="graph-tool"
except ImportError:
    import networkx as _nx
    from pydna._simple_paths8_NX import all_simple_paths_edges   as _all_simple_paths_edges
    from pydna._simple_paths8_NX import all_circular_paths_edges as _all_circular_paths_edges
    gc="networkx"


class _Fragment(_Dseqrecord):
    '''This class holds information about a DNA fragment in an assembly.
    This class is instantiated by the :class:`Assembly` class and is not
    meant to be instantiated directly.

    '''
    
    def __init__(self, record, *args,
                               start1    = 0,
                               end1      = 0,
                               start2    = 0,
                               end2      = 0,
                               alignment = 0,
                               i         = 0, 
                               **kwargs):

        super().__init__(record, *args, **kwargs)

        self.start1             = start1
        self.end1               = end1
        self.left_overlap_size  = end1-start1
        self.start2             = start2
        self.end2               = end2
        self.right_overlap_size = end2-start2
        self.alignment          = alignment
        self.i                  = i

    def __str__(self):
        return ("Fragment alignment {}\n").format(self.alignment)+super().__str__()

class Memoize(type):
    @_memorize("Assembly")
    def __call__(cls, *args, **kwargs):
        return super().__call__(*args, **kwargs)

class Assembly(object, metaclass = Memoize):
    '''Assembly of a list of linear DNA fragments into linear or circular constructs.
    The Assembly is meant to replace the Assembly method as it is easier to use.
    Accepts a list of Dseqrecords (source fragments) to initiate an Assembly object.
    Several methods are available for analysis of overlapping sequences, graph construction
    and assembly.

    Parameters
    ----------

    dsrecs : list
        a list of Dseqrecord objects.

    Examples
    --------

    >>> from pydna.assembly import Assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
    >>> b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
    >>> c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
    >>> x = Assembly((a,b,c), limit=14)
    >>> x
    Assembly:
    Sequences........................: [33] [34] [35]
    Sequences with shared homologies.: [33] [34] [35]
    Homology limit (bp)..............: 14
    Number of overlaps...............: 3
    Nodes in graph(incl. 5' & 3')....: 5
    Only terminal overlaps...........: No
    Circular products................: [59]
    Linear products..................: [74] [73] [73] [54] [54] [53] [15] [14] [14]
    >>> x.circular_products
    [Contig(o59)]
    >>> x.circular_products[0].seq.watson
    'CCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGtacgatgctatactg'

    '''

    def __init__(self, dsrecs, limit = 25, only_terminal_overlaps=False, max_nodes=None):

        self.dsrecs = dsrecs
        ''' Sequences fed to this class is stored in this property'''
        self.only_terminal_overlaps = only_terminal_overlaps
        ''' Consider only terminal overlaps?'''
        self.limit = limit
        ''' The shortest common sub strings to be considered '''
        self.max_nodes = max_nodes or len(self.dsrecs)
        ''' The max number of nodes allowed. This can be reset to some other value'''
        self.only_terminal_overlaps = only_terminal_overlaps
        self.vertices={}
        
        for dr in self.dsrecs:
            if dr.name in ("",".", "<unknown name>", None):
                dr.name = "frag{}".format(len(dr))

        if self.only_terminal_overlaps:
            algorithm = _terminal_overlap
        else:
            algorithm = _common_sub_strings

        # analyze_overlaps
        cols = {}
        for dsrec in self.dsrecs:
            dsrec.features = [f for f in dsrec.features if f.type!="overlap"]
            dsrec.seq = _Dseq(dsrec.seq.todata)
        rcs = {dsrec:dsrec.rc() for dsrec in self.dsrecs}
        matches=[]
        dsset=_OrderedSet()

        for a, b in _itertools.combinations(self.dsrecs, 2):
            match = algorithm( str(a.seq).upper(),
                               str(b.seq).upper(),
                               self.limit)
            if match:
                matches.append((a, b, match))
                dsset.add(a)
                dsset.add(b)
            match = algorithm( str(a.seq).upper(),
                               str(rcs[b].seq).upper(),
                               self.limit)
            if match:
                matches.append((a, rcs[b], match))
                dsset.add(a)
                dsset.add(rcs[b])
                matches.append((rcs[a], b, [(len(a)-sa-le,len(b)-sb-le,le) for sa,sb,le in match]))
                dsset.add(b)
                dsset.add(rcs[a])

        self.no_of_olaps=0

        for a, b, match in matches:
            for start_in_a, start_in_b, length in match:
                self.no_of_olaps+=1
                chksum = a[start_in_a:start_in_a+length].seguid()
                #assert chksum == b[start_in_b:start_in_b+length].seguid()

                try:
                    fcol, revcol = cols[chksum]
                except KeyError:
                    fcol = '#%02X%02X%02X' % (_random.randint(175,255),_random.randint(175,255),_random.randint(175,255))
                    rcol = '#%02X%02X%02X' % (_random.randint(175,255),_random.randint(175,255),_random.randint(175,255))
                    cols[chksum] = fcol,rcol

                qual      = {"note"             : ["olp_{}".format(chksum)],
                             "chksum"           : [chksum],
                             "ApEinfo_fwdcolor" : [fcol],
                             "ApEinfo_revcolor" : [rcol]}

                if not chksum in [f.qualifiers["chksum"][0] for f in a.features if f.type == "overlap"]:
                    a.features.append( _SeqFeature( _FeatureLocation(start_in_a,
                                                                   start_in_a + length),
                                                                   type = "overlap",
                                                                   qualifiers = qual))
                if not chksum in [f.qualifiers["chksum"][0] for f in b.features if f.type == "overlap"]:
                    b.features.append( _SeqFeature( _FeatureLocation(start_in_b,
                                                                   start_in_b + length),
                                                                   type = "overlap",
                                                                   qualifiers = qual))
        for ds in dsset:
            ds.features = sorted([f for f in ds.features], key = _operator.attrgetter("location.start"))

        self.analyzed_dsrecs = list(dsset)
        
################################################This is where we either use nx or gt
#####################################GRAPH-TOOL#####################################
        if gc=="graph-tool":
            timegt=time.time()
            # Create graph
            self.G=_gt.Graph(directed=True)
            #Edge properties
            self.G.edge_properties["frag_Property"]= self.G.new_edge_property("python::object")
            self.G.edge_properties["iter_Property"]= self.G.new_edge_property("float")
            self.G.edge_properties["weight"]= self.G.new_edge_property("float")
            
            #Vertex properties
            self.G.vertex_properties["vertex"]=self.G.new_vertex_property("string")
            #5' will be v5, 3' will be v3
            v5=self.G.add_vertex()
            self.G.vertex_properties["vertex"][v5]="5"
            self.vertices['5']=v5
            v3=self.G.add_vertex()
            self.G.vertex_properties["vertex"][v3]="3"
            self.vertices['3']=v3
        
            for i, dsrec in enumerate(self.analyzed_dsrecs):
    
                overlaps = sorted( list({f.qualifiers['chksum'][0]:f for f in dsrec.features
                                    if f.type=='overlap'}.values()),
                                   key = _operator.attrgetter('location.start'))
    
                if overlaps:
                    overlaps = ([_SeqFeature(_FeatureLocation(0, 0),
                                 type = 'overlap',
                                 qualifiers = {'chksum':['5']})]+
                                 overlaps+
                                [_SeqFeature(_FeatureLocation(len(dsrec),len(dsrec)),
                                            type = 'overlap',
                                            qualifiers = {'chksum':['3']})])
    
                    for olp1, olp2 in _itertools.combinations(overlaps, 2):
    
                        n1 = olp1.qualifiers['chksum'][0]
                        if n1 not in self.vertices:
                            ni1=self.G.add_vertex()
                            self.vertices[n1]=ni1                    
                            self.G.vertex_properties["vertex"][ni1]=n1
                        n2 = olp2.qualifiers['chksum'][0]
                        if n2 not in self.vertices:
                            ni2=self.G.add_vertex()
                            self.vertices[n2]=ni2
                            self.G.vertex_properties["vertex"][ni2]=n2
                        if n1 == '5' and n2=='3':
                            continue
    
                        s1,e1,s2,e2 = (olp1.location.start.position,
                                       olp1.location.end.position,
                                       olp2.location.start.position,
                                       olp2.location.end.position,)
    
                        source_fragment = _Fragment(dsrec, start1=s1,end1=e1,start2=s2,end2=e2,i=i) #_Fragment(dsrec,s1,e1,s2,e2,i)
                        #Adding the edge
                        edge=self.G.add_edge( self.vertices[n1], self.vertices[n2])
                        #Adding edge properties
                        self.G.edge_properties["frag_Property"][edge]=source_fragment
                        self.G.edge_properties["iter_Property"][edge]=i
                        self.G.edge_properties["weight"][edge]=s1-e1
            print("time graph gt-", time.time()-timegt)    

#####################################NETWORKX#####################################			
		elif gc=="networkx":
            timenx=time.time()
            self.G=_nx.MultiDiGraph(multiedges=True, name ="original graph" , selfloops=False)
            self.G.add_node( '5' )
            self.G.add_node( '3' )
    
            for i, dsrec in enumerate(self.analyzed_dsrecs):
    
                overlaps = sorted( list({f.qualifiers['chksum'][0]:f for f in dsrec.features
                                    if f.type=='overlap'}.values()),
                                   key = _operator.attrgetter('location.start'))
    
                if overlaps:
                    overlaps = ([_SeqFeature(_FeatureLocation(0, 0),
                                 type = 'overlap',
                                 qualifiers = {'chksum':['5']})]+
                                 overlaps+
                                [_SeqFeature(_FeatureLocation(len(dsrec),len(dsrec)),
                                            type = 'overlap',
                                            qualifiers = {'chksum':['3']})])
    
                    for olp1, olp2 in _itertools.combinations(overlaps, 2):
    
                        n1 = olp1.qualifiers['chksum'][0]
                        n2 = olp2.qualifiers['chksum'][0]
    
                        if n1 == '5' and n2=='3':
                            continue
    
                        s1,e1,s2,e2 = (olp1.location.start.position,
                                       olp1.location.end.position,
                                       olp2.location.start.position,
                                       olp2.location.end.position,)
    
                        source_fragment = _Fragment(dsrec, start1=s1,end1=e1,start2=s2,end2=e2,i=i) #_Fragment(dsrec,s1,e1,s2,e2,i)
    
                        self.G.add_edge( n1, n2,
                                         frag=source_fragment,
                                         weight = s1-e1,
                                         i = i)
            print("time graph nx-", time.time()-timenx) 
        linear_products=_defaultdict(list)
        node5= 0 if gc=="graph-tool" else '5'
        node3= 1 if gc=="graph-tool" else '3'
        lineartime=time.time()
        for path in _all_simple_paths_edges(self.G, node5, node3, data=True, cutoff=self.max_nodes):
            pred_frag = _copy(list(path[0][2].values()).pop()['frag'])
            source_fragments = [pred_frag, ]

            if pred_frag.start2<pred_frag.end1:
                result=pred_frag[pred_frag.start2+(pred_frag.end1-pred_frag.start2):pred_frag.end2]
            else:
                result=pred_frag[pred_frag.end1:pred_frag.end2]

            for first_node, second_node, edgedict in path[1:]:

                edgedict = list(edgedict.values()).pop()

                f  = _copy(edgedict['frag'])

                f.alignment =  pred_frag.alignment + pred_frag.start2- f.start1
                source_fragments.append(f)

                if f.start2>f.end1:
                    result+=f[f.end1:f.end2]
                else:
                    result+=f[f.start2+(f.end1-f.start2):f.end2]

                pred_frag = f
            

            add=True
            for lp in linear_products[len(result)]:
                if (str(result.seq).lower() == str(lp.seq).lower()
                    or
                    str(result.seq).lower() == str(lp.seq.reverse_complement()).lower()):
                    add=False
            for dsrec in self.dsrecs:
                if (str(result.seq).lower() == str(dsrec.seq).lower()
                    or
                    str(result.seq).lower() == str(dsrec.seq.reverse_complement()).lower()):
                    add=False
            if add:
                linear_products[len(result)].append(_Contig( result, source_fragments = source_fragments))
        print("linear path",time.time()-lineartime)
        self.linear_products = list(_itertools.chain.from_iterable(linear_products[size] for size in sorted(linear_products, reverse=True)))

###############################################GC checkpoint for circular assembly
        # circular assembly
        self.cG = self.G.copy()
        if gc=="graph-tool": self.cG.remove_vertex([0,1])
        elif gc=="networkx": self.cG.remove_nodes_from(('5','3'))
        
        circular_products={}
                  
        circulartime=time.time()
        for pth in _all_circular_paths_edges(self.cG):

            ns = min( enumerate(pth), key = lambda x:x[1][2]['i'] )[0]

            path = pth[ns:]+pth[:ns]

            pred_frag = _copy(path[0][2]['frag'])

            source_fragments = [pred_frag, ]
            
            if pred_frag.start2<pred_frag.end1:
                result=pred_frag[pred_frag.start2+(pred_frag.end1-pred_frag.start2):pred_frag.end2]
            else:
                result=pred_frag[pred_frag.end1:pred_frag.end2]
                
            result.seq = _Dseq(str(result.seq))

            for first_node, second_node, edgedict in path[1:]:

                f  = _copy(edgedict['frag'])

                f.alignment =  pred_frag.alignment + pred_frag.start2- f.start1
                source_fragments.append(f)

                if f.start2>f.end1:
                    nxt = f[f.end1:f.end2]
                else:
                    nxt =f[f.start2+(f.end1-f.start2):f.end2]
                nxt.seq = _Dseq(str(nxt.seq))
                result+=nxt

                pred_frag = f

            r = _Dseqrecord(result, circular=True)
            circular_products[r.cseguid()] = _Contig(r, source_fragments = source_fragments )
        print("circular time", time.time()-circulartime)
        
        import functools
        def comp(item1, item2):
            item1_len = len(item1.source_fragments)
            item2_len = len(item2.source_fragments)
            if item1_len < item2_len:
                return -1
            if item1_len > item2_len:
                return 1
            return 0
        
        self.circular_products = sorted(circular_products.values(), key=functools.cmp_to_key(comp))
        self.circular_products.sort(key=len, reverse=True)

    def __repr__(self):
        return   ( "Assembly:\n"
                   "Sequences........................: {sequences}\n"
                   "Sequences with shared homologies.: {analyzed_dsrecs}\n"
                   "Homology limit (bp)..............: {limit}\n"
                   "Number of overlaps...............: {no_of_olaps}\n"
                   "Nodes in graph(incl. 5' & 3')....: {nodes}\n"
                   "Only terminal overlaps...........: {pr}\n"
                   "Circular products................: {cp}\n"
                   "Linear products..................: {lp}"    ).format(sequences       = " ".join("[{}]".format(len(x)) for x in self.dsrecs),
                                                                         analyzed_dsrecs = " ".join("[{}]".format(len(x)) for x in self.analyzed_dsrecs),
                                                                         limit           = self.limit,
                                                                         no_of_olaps     = self.no_of_olaps,
                                                                         nodes           = len(self.G.get_vertices()) if gc=="graph-tool" else len(self.G.nodes()) ,
                                                                         pr              = {True:"Yes",False:"No"}[self.only_terminal_overlaps],
                                                                         cp              = " ".join("[{}]".format(len(x)) for x in self.circular_products),
                                                                         lp              = " ".join("[{}]".format(len(x)) for x in self.linear_products))

																		 


#These 3 functions are useful to profile the assembly module
def example_seqs(n_seqs):
	res=[]
	a = _Dseqrecord("ttctagaactagtggatcccccgggctgcagatgagtgaaggccccgtcaaattcgaaaaaaataccgtcatatctgtctttggtgcgtcaggtgatctggcaaagaagaagacttttcccgccttatttgggcttttcagagaaggttaccttgatccatctaccaagatcttcggttatgcccggtccaaattgtccatggaggaggacctgaagtcccgtgtcctaccccacttgaaaaaacctcacggtgaagccgatgactctaaggtcgaacagttcttcaagatggtcagctacatttcgggaaattacgacacagatgaaggcttcgacgaattaagaacgcagatcgagaaattcgagaaaagtgccaacgtcgatgtcccacaccgtctcttctatctggccttgccgccaagcgtttttttgacggtggccaagcagatcaagagtcgtgtgtacgcagagaatggcatcacccgtgtaatcgtagagaaacctttcggccacgacctggcctctgccagggagctgcaaaaaaacctggggcccctctttaaagaagaagagttgtacagaattgaccattacttgggtaaagagttggtcaagaatcttttagtcttgaggttcggtaaccagtttttgaatgcctcgtggaatagagacaacattcaaagcgttcagat")
	b = _Dseqrecord("tatcgataagcttgatatcgaattcctgcagctaattatccttcgtatcttctggcttagtcacgggccaagcgtaagggtgcttttcgggcataacatacttgtgtttttgcatatattccttcaatccctttggacctcttgatccgtaggggtaaatttccggtgttggaccgtccggacgctctatgtgcttcagtaatggggtgaatatgccccaactgatatccaattcgtcatctctgacaaagttggaatggtcacccagtagggcgtctcttatcaacacctcgtaagcctctggaatccaaaagtcttggtacctgcttgcgtaagttagattcagatctgtgacttgggtagcatttgacagaccaggggtcttagcattaaactttaggtacacagcggcatcgggctgcactctgatgaccagttcgttatttggaatgtctttgaagacacccgatgcgaccgctttgtactgcagtctgatctccaccttggactcattcaaagccttaccggcacgcatcatgatggggacgccctcccaacgctcgttttcgatgttgaaagtcattgctgcaaaagtgacacatttagagtccttgtctacagtgtcatcatccacgtaggcgggcttagacccgtcctcagatttaccgtactggcccaagaggacgtcgtccgtgtcgatgggggccacggcctttagaaccttaaccttttcgtcacgaatagattccgggtcaaaagacaccggtctttccatagtcaagagagtcatgatttgtaacagatggttctgcatcacgtctctgattatgcctatagagtcgaaatagccgccacggccttcggtgccgaacctctctttaaacgaaatctgaacgctttgaatgttgtctctattccacgaggcattcaaaaa")
	c = _Dseqrecord("gaattcgatatcaagcttatcgataccgtcgacctcgagtcatgtaattagttatgtcacgcttacattcacgccctccccccacatccgctctaaccgaaaaggaaggagttagacaacctgaagtctaggtccctatttatttttttatagttatgttagtattaagaacgttatttatatttcaaatttttcttttttttctgtacagacgcgtgtacgcatgtaacattatactgaaaaccttgcttgagaaggttttgggacgctcgaaggctttaatttgcggccggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgcgacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgtttacaatttcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatcgacggtcgaggagaacttctagtatatccacatacctaatattattgccttattaaaaatggaatcccaacaattacatcaaaatccacattctcttcaaaatcaattgtcctgtacttccttgttcatgtgtgttcaaaaacgttatatttataggataattatactctatttctcaacaagtaattggttgtttggccgagcggtctaaggcgcctgattcaagaaatatcttgaccgcagttaactgtgggaatactcaggtatcgtaagatgcaagagttcgaatctcttagcaaccattatttttttcctcaacataacgagaacacacaggggcgctatcgcacagaatcaaattcgatgactggaaattttttgttaatttcagaggtcgcctgacgcatatacctttttcaactgaaaaattgggagaaaaaggaaaggtgagaggccggaaccggcttttcatatagaatagagaagcgttcatgactaaatgcttgcatcacaatacttgaagttgacaatattatttaaggacctattgttttttccaataggtggttagcaatcgtcttactttctaacttttcttaccttttacatttcagcaatatatatatatatttcaaggatataccattctaatgtctgcccctatgtctgcccctaagaagatcgtcgttttgccaggtgaccacgttggtcaagaaatcacagccgaagccattaaggttcttaaagctatttctgatgttcgttccaatgtcaagttcgatttcgaaaatcatttaattggtggtgctgctatcgatgctacaggtgtcccacttccagatgaggcgctggaagcctccaagaaggttgatgccgttttgttaggtgctgtggctggtcctaaatggggtaccggtagtgttagacctgaacaaggtttactaaaaatccgtaaagaacttcaattgtacgccaacttaagaccatgtaactttgcatccgactctcttttagacttatctccaatcaagccacaatttgctaaaggtactgacttcgttgttgtcagagaattagtgggaggtatttactttggtaagagaaaggaagacgatggtgatggtgtcgcttgggatagtgaacaatacaccgttccagaagtgcaaagaatcacaagaatggccgctttcatggccctacaacatgagccaccattgcctatttggtccttggataaagctaatcttttggcctcttcaagattatggagaaaaactgtggaggaaaccatcaagaacgaattccctacattgaaggttcaacatcaattgattgattctgccgccatgatcctagttaagaacccaacccacctaaatggtattataatcaccagcaacatgtttggtgatatcatctccgatgaagcctccgttatcccaggttccttgggtttgttgccatctgcgtccttggcctctttgccagacaagaacaccgcatttggtttgtacgaaccatgccacggttctgctccagatttgccaaagaataaggttgaccctatcgccactatcttgtctgctgcaatgatgttgaaattgtcattgaacttgcctgaagaaggtaaggccattgaagatgcagttaaaaaggttttggatgcaggtatcagaactggtgatttaggtggttccaacagtaccaccgaagtcggtgatgctgtcgccgaagaagttaagaaaatccttgcttaaaaagattctctttttttatgatatttgtacataaactttataaatgaaattcataatagaaacgacacgaaattacaaaatggaatatgttcatagggtagacgaaactatatacgcaatctacatacatttatcaagaaggagaaaaaggaggatagtaaaggaatacaggtaagcaaattgatactaatggctcaacgtgataaggaaaaagaattgcactttaacattaatattgacaaggaggagggcaccacacaaaaagttaggtgtaacagaaaatcatgaaactacgattcctaatttgatattggaggattttctctaaaaaaaaaaaaatacaacaaataaaaaacactcaatgacctgaccatttgatggagtttaagtcaataccttcttgaagcatttcccataatggtgaaagttccctcaagaattttactctgtcagaaacggccttacgacgtagtcgatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagtatgatccaatatcaaaggaaatgatagcattgaaggatgagactaatccaattgaggagtggcagcatatagaacagctaaagggtagtgctgaaggaagcatacgataccccgcatggaatgggataatatcacaggaggtactagactacctttcatcctacataaatagacgcatataagtacgcatttaagcataaacacgcactatgccgttcttctcatgtatatatatatacaggcaacacgcagatataggtgcgacgtgaacagtgagctgtatgtgcgcagctcgcgttgcattttcggaagcgctcgttttcggaaacgctttgaagttcctattccgaagttcctattctctagaaagtataggaacttcagagcgcttttgaaaaccaaaagcgctctgaagacgcactttcaaaaaaccaaaaacgcaccggactgtaacgagctactaaaatattgcgaataccgcttccacaaacattgctcaaaagtatctctttgctatatatctctgtgctatatccctatataacctacccatccacctttcgctccttgaacttgcatctaaactcgacctctacattttttatgtttatctctagtattactctttagacaaaaaaattgtagtaagaactattcatagagtgaatcgaaaacaatacgaaaatgtaaacatttcctatacgtagtatatagagacaaaatagaagaaaccgttcataattttctgaccaatgaagaatcatcaacgctatcactttctgttcacaaagtatgcgcaatccacatcggtatagaatataatcggggatgcctttatcttgaaaaaatgcacccgcagcttcgctagtaatcagtaaacgcgggaagtggagtcaggctttttttatggaagagaaaatagacaccaaagtagccttcttctaaccttaacggacctacagtgcaaaaagttatcaagagactgcattatagagcgcacaaaggagaaaaaaagtaatctaagatgctttgttagaaaaatagcgctctcgggatgcatttttgtagaacaaaaaagaagtatagattctttgttggtaaaatagcgctctcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttttacaaaaatgaagcacagattcttcgttggtaaaatagcgctttcgcgttgcatttctgttctgtaaaaatgcagctcagattctttgtttgaaaaattagcgctctcgcgttgcatttttgttctacaaaatgaagcacagatgcttcgttcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttacctcactcattaggcaccccaggctttacactttatgcttccggctcctatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcgcgcaattaaccctcactaaagggaacaaaagctggagctcagtttatcattatcaatactcgccatttcaaagaatacgtaaataattaatagtagtgattttcctaactttatttagtcaaaaaattagccttttaattctgctgtaacccgtacatgcccaaaatagggggcgggttacacagaatatataacatcgtaggtgtctgggtgaacagtttattcctggcatccactaaatataatggagcccgctttttaagctggcatccagaaaaaaaaagaatcccagcaccaaaatattgttttcttcaccaaccatcagttcataggtccattctcttagcgcaactacagagaacaggggcacaaacaggcaaaaaacgggcacaacctcaatggagtgatgcaacctgcctggagtaaatgatgacacaaggcaattgacccacgcatgtatctatctcattttcttacaccttctattaccttctgctctctctgatttggaaaaagctgaaaaaaaaggttgaaaccagttccctgaaattattcccctacttgactaataagtatataaagacggtaggtattgattgtaattctgtaaatctatttcttaaacttcttaaattctacttttatagttagtcttttttttagttttaaaacaccagaacttagtttcgacggattctagaactagtggatcccccgggctgcag")
	d = _Dseqrecord("tcctgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgacgtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgg")
	e = _Dseqrecord("tcctgacgggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatgtcgaggaacgccaggttgcccactttctcactagtgacctgcagccgg")
	f = _Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
	g = _Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
	h = _Dseqrecord("tattctggctgtatcagtgaaggccccgtcaaattcgaaaaaaataccgtcatatctgtctttggtgcgtcaggtgatctggcaaagaagaagacttttcccgccttatttgggcttttGGGGGtacgatgctatactg")
	k = _Dseqrecord("tattctggctgtatcagtgaaggccccgtcaaattgacgggtaattttgatttgcatgccgtccgggtgagtcatagcgtctggttgttttgccagattcagcagagtctgtgcaatgcggccgctgacgtatctcagttcggtgtaggtccgaaaaaaataccgtcatatctgtctttgg")
	j = _Dseqrecord("tgtgctgtgctctaTTTTTtattctgccgtcgacctcgagtcatgtaattagttatgtcacgcttacattcacgccctccccccacatccgctctaaccgaaaaggagctgtatc")
	available_seqs=[a,b,c,d,e,f,g,h,k,j]
	for i in range(n_seqs):
		res.append(available_seqs[i])
	return tuple(res)
																		 
def testseqs(tuple_of_seqs,output_file,limit=14, average="F",average_iter=10):
	x=Assembly(tuple_of_seqs,limit)
	import sys
	import time
    sys.stdout=open(output_file,"w")
	if average='T':
		for i in range(average_iter):
			startTime=time.time()
			x=Assembly(tuple_of_seqs,limit)
			print("total time-", time.time()-startTime)
	else:
		startTime=time.time()
		x=Assembly(tuple_of_seqs,limit)
		print("total time-", time.time()-startTime)
	print("Analysis complete")

def testcProfile(tuple_of_seqs,output_file,limit=14, top_nLines=100):
	import cProfile
	import pstats
    import sys
    cProfile.run("Assembly((a,b,c,d,e,f,g,h,k,j), limit=14)","restats")
	sys.stdout=open("profilingGT_10.txt","w")
    p = pstats.Stats('restats')
    p.strip_dirs().sort_stats("tottime").print_stats(top_nLines)
	
if __name__=="__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True)
    _os.environ["pydna_cache"]=cache  
	#testseqs(example_seqs(3),'average_of_3',14,T)
	#testcProfile(example_seqs(3),'cProfile_of_3',14,50)




