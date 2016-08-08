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
'''
This module provides tools for creating and using factorized representations of Bayesian networks. Factorized representations of Bayesian networks are discrete CPDs whose values have been flattened into a single array, while the cardinalities and strides of each variable represented are kept track of separately. With the proper setup, these flattened structures can be more easily multiplied together, reduced, and operated on. For more information on factors cf. Koller et al. Ch. 4.

'''

from .oldtablecpdfactorization import TableCPDFactorization as old

class TableCPDFactorization (old):
    '''Factorized discrete CPD Bayesian Network.

    This class represents a factorized Bayesian network with discrete
    CPD tables. 
    '''
    def sumproductve(self, vertices):
        '''
        Eliminate each vertex in *vertices* from *factorlist* using *sumproducteliminatevar*.
        
        Arguments:
            1. *vertices* -- A list of UUIDs of vertices to be eliminated.
        
        Returns:
            1. *factor* -- the resulting single TableCPDFactor

        '''
    
        # eliminate one by one
        for vertex in vertices:
            self.sumproducteliminatevar(vertex)
        
        # multiply together if many factors remain 
        for i in range(1, len(self.factorlist)):
            self.factorlist[0].multiplyfactor(self.factorlist[i])
        
        self.factorlist = [self.factorlist[0]]
        return self.factorlist[0]
  
    def condprobve(self, query, evidence={}):
        '''Calculate the conditional probabilities for *query*.
        
        Eliminate all variables in *factorlist* except for the ones
        queried. Adjust all distributions for the evidence
        given. Return the probability distribution over a set of
        variables given by the keys of *query* given *evidence*.
        
        Arguments:
            1. *query* -- A dict containing (key: value) pairs reflecting (variable: value) that represents what outcome to calculate the probability of. 
            2. *evidence* -- A dict containing (key: value) pairs reflecting (variable: value) that represents what is known about the system.
                    
        Attributes modified:
            1. *factorlist* -- Modified to be one factor representing the probability distribution of the query variables given the evidence.
                           
        The function returns *factorlist* after it has been modified as above.
        
        Usage example: this code would return the distribution over a queried node, given evidence::

            import json

            from libpgm.graphskeleton import GraphSkeleton
            from libpgm.nodedata import NodeData
            from libpgm.discretebayesiannetwork import DiscreteBayesianNetwork
            from libpgm.tablecpdfactorization import TableCPDFactorization

            # load nodedata and graphskeleton
            nd = NodeData()
            skel = GraphSkeleton()
            nd.load("../tests/unittestdict.txt")
            skel.load("../tests/unittestdict.txt")

            # toporder graph skeleton
            skel.toporder()

            # load evidence
            evidence = dict(Letter='weak')
            query = dict(Grade='A')

            # load bayesian network
            bn = DiscreteBayesianNetwork(skel, nd)

            # load factorization
            fn = TableCPDFactorization(bn)

            # calculate probability distribution
            result = fn.condprobve(query, evidence)

            # output
            print json.dumps(result.vals, indent=2)
            print json.dumps(result.scope, indent=2)
            print json.dumps(result.card, indent=2)
            print json.dumps(result.stride, indent=2)

        '''

        eliminate = [vertex for vertex in self.bn.V
                     if vertex not in query
                     if vertex not in evidence]

        factorlist = self.factorlist[:]
        
        # modify factors to account for the evidence
        for vertex, value in evidence.items():
            factorlist = [
                factor.reducefactor(vertex, value)
                if (factor.scope.count(vertex) > 0)
                else factor
                for factor in factorlist]
            # Eliminate skope-free vertices
            factorlist = [factor for factor in factorlist
                          if factor.scope]
                    
        # eliminate all necessary variables in the new factor set to produce result
        factor = self.sumproductve(eliminate)
        
        # normalize result
        norm = sum(factor.vals)
        for x, val in enumerate(factor.vals):
            factor.vals[x] /= norm
            
        # return table
        return factor

