# Copyright (c) 2012, CyberPoint International, LLC
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the CyberPoint International, LLC nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL CYBERPOINT INTERNATIONAL, LLC BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''This module provides tools for creating and using an individual
factorized representation of a node. See description of factorized
representations in :doc:`tablecpdfactorization`.

'''

def prod(l):
    """ Calculate the product of the iterable l """
    r = 1
    for x in l:
        r *= x
    return r

class TableCPDFactor:
    """Factorized representation of a CPD table.

    This class represents a factorized representation of a conditional
    probability distribution table.

    """
    def __init__(self, vertex, bn):
        '''Construct a factorized CPD table from a vertex in a discrete
        Bayesian network.

        This class is constructed with a
        :doc:`DiscreteBayesianNetwork <discretebayesiannetwork>`
        instance and a *vertex* name as arguments. First it stores
        these inputs in *inputvertex* and *inputbn*. Then, it creates
        a factorized representation of *vertex*, storing the values in
        *vals*, the names of the variables involved in *scope* the
        cardinality of each of these variables in *card* and the
        stride of each of these variables in *stride*.

        '''
        self.inputvertex = vertex
        '''The name of the vertex.'''
        self.inputbn = bn
        '''The :doc:`DiscreteBayesianNetwork <discretebayesiannetwork>` instance that the vertex lives in.'''
        
        self.vals = []
        '''A flat array of all the values from the CPD.'''
        self.stride = {}
        '''A dict of {vertex: value} pairs for each vertex in *self.scope*, where vertex is the name of the vertex and value is the self.stride of that vertex in the *self.vals* array.'''
        self.card = []
        '''A list of the self.cardinalities of each vertex in self.scope, where cardinality is the number of values that the vertex may take. The cardinalities are indexed according to the vertex's index in *scope*.'''
        self.scope = []
        '''An array of vertices that affect the self.vals found in *vals*. Normally, this is the node itself and its parents.'''
        root = bn.Vdata[vertex]["cprob"]

        parents = bn.Vdata[vertex]["parents"]
        # add values
        def explore(_dict, key, depth, totaldepth):
            if depth == totaldepth:
                for x in _dict[key]:
                    self.vals.append(x)
                return
            else:
                for val in bn.Vdata[parents[depth]]["vals"]:
                    ckey = key + (val,)
                    explore(_dict, ckey, depth+1, totaldepth)
                    
        if not parents:
            self.vals = bn.Vdata[vertex]["cprob"]
            assert len(self.vals) == len(bn.Vdata[vertex]["vals"])
            assert abs(sum(self.vals) - 1) < 1e-8
        else: 
            td = len(parents)
            explore(root, (), 0, td)
        
        # add self.cardinalities
        self.card.append(len(bn.Vdata[vertex]["vals"]))
        if (bn.Vdata[vertex]["parents"] != None):
            for parent in reversed(bn.Vdata[vertex]["parents"]):
                self.card.append(len(bn.Vdata[parent]["vals"]))
            
        # add self.scope
        self.scope.append(vertex)
        if (bn.Vdata[vertex]["parents"] != None):
            for parent in reversed(bn.Vdata[vertex]["parents"]):
                self.scope.append(parent)
        
        
        # add self.self.stride
        t_stride = 1
        self.stride = dict()
        for x in range(len(self.scope)):
            self.stride[self.scope[x]] = (t_stride)
            t_stride *= len(bn.Vdata[self.scope[x]]["vals"])
        
        

    def multiplyfactor(self, other):  # cf. PGM 359 
        '''Multiply this factor by another Factor

        Multiplying factors means taking the union of the scopes, and
        for each combination of variables in the scope, multiplying
        together the probabilities from each factor that that
        combination will be found.
        
        Arguments:
            1. *other* -- An instance of :doc:`TableCPDFactor <tablecpdfactor>` class representing the factor to multiply by.
                 
        Attributes modified: 
            *vals*, *scope*, *stride*, *t_card* -- Modified to reflect the data of the new product factor.
                                                         
        For more information cf. Koller et al. 359.

        '''

        # merge t_scopes
        scope = self.scope
        card = self.card
        for t_scope, t_card in zip(other.scope, other.card):
            try:
                scope.index(t_scope)
            except: 
                scope.append(t_scope)
                card.append(t_card)
    
        # algorithm (see book)
        assignment = {}
        vals = []
        j = 0
        k = 0
        for _ in range(prod(card)):
            vals.append(
                self.vals[j] * other.vals[k])
            
            for t_card, t_scope in zip(card, scope):
                assignment[t_scope] = assignment.get(t_scope, 0) + 1
                if (assignment[t_scope] == t_card):
                    assignment[t_scope] = 0
                    if t_scope in self.stride:
                        j = j - (t_card - 1) * self.stride[t_scope]
                    if t_scope in other.stride:
                        k = k - (t_card - 1) * other.stride[t_scope]
                else:
                    if t_scope in self.stride:
                        j = j + self.stride[t_scope]
                    if t_scope in other.stride:
                        k = k + other.stride[t_scope]
                    break
            
        # add strides
        t_stride = 1 
        stride = {}
        for t_card, t_scope in zip(card, scope):
            stride[t_scope] = (t_stride)
            t_stride *= t_card
    
        self.vals = vals
        self.scope = scope
        self.card = card
        self.stride = stride

    def reducefactor(self, vertex, value=None):
        '''Sum out the variable specified by *vertex* from the factor.

        Summing out means summing all sets of entries together where
        *vertex* is the only variable changing in the set. Then
        *vertex* is removed from the scope of the factor.
        
        Arguments:
            1. *vertex* -- The name of the variable to be summed out.
        
        Attributes modified: 
            *vals*, *scope*, *stride*, *card* -- Modified to reflect the data of the summed-out product factor.
        
        For more information see Koller et al. 297.

        '''
        vscope = self.scope.index(vertex)
        vstride = self.stride[vertex]
        vcard = self.card[vscope]

        result = [0 for i in range(len(self.vals)//self.card[vscope])]
        
        # machinery that calculates values in summed out factor
        k = 0
        lcardproduct = prod(self.card[:vscope])
        for i, entry in enumerate(result):
            if value is None:
                for h in range(vcard):
                    result[i] += self.vals[k + vstride * h]
            else:
                index = self.inputbn.Vdata[vertex]['vals'].index(value)
                result[i] += self.vals[k + vstride * index]
                
            k += 1
            if (k % lcardproduct == 0):
                k += (lcardproduct * (vcard - 1))
        self.vals = result
        
        # modify scope, card, and stride in new factor
        self.scope.remove(vertex)
        del(self.card[vscope])
        for i in range(vscope, len(self.stride)-1):
            self.stride[self.scope[i]] //= vcard
        del(self.stride[vertex])
        
    sumout = reducefactor

    def copy(self):
        '''Return a copy of the factor.'''
        copy = type(self)(self.inputvertex, self.inputbn)
        copy.vals = self.vals[:]
        copy.stride = self.stride.copy()
        copy.scope = self.scope[:]
        copy.card = self.card[:]
        return copy
