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
    def sumproductve(self,
                     vertices: "A sequence of UUIDs of vertices to be eliminated."
    ) -> "the resulting single TableCPDFactor":
        '''Eliminate each vertex in *vertices* from *factorlist*

        Using *sumproducteliminatevar*, remove all vertices in the
        sequence from self.factorlist

        '''
    
        # eliminate one by one
        for vertex in vertices:
            self.sumproducteliminatevar(vertex)
        
        # multiply together if many factors remain 
        result = self.factorlist[0].copy()
        for i in range(1, len(self.factorlist)):
            result.multiplyfactor(self.factorlist[i])
        
        return result

    def condition(self,
                  evidence: "A dictionary of {node name: value} mappings",
                  in_place: "If True, write the changes back to self.factorlist." = False,
                  reset_before: "If True, reset all evidence before adding this, otherwise keep old evidence." = False
    ) -> "the resulting factorlist":
        """Create a factorlist reflecting evidence.
        
        Adjust all distributions in self.factorlist for the evidence
        given, and return the resulting new factorlist.

        """

        if reset_before:
            self.reset()

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

        if in_place:
            self.factorlist = factorlist

        return factorlist
  
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

        self.condition(evidence, in_place=True)
                   
        eliminate = [vertex for vertex in self.bn.V
                     if vertex not in query
                     if vertex not in evidence]

        # eliminate all necessary variables in the new factor set to produce result
        factor = self.sumproductve(eliminate)
        
        # normalize result
        norm = sum(factor.vals)
        for x, val in enumerate(factor.vals):
            factor.vals[x] /= norm
            
        # return table
        return factor
    
    def specificquery(self, query, evidence=None):
        '''
        Eliminate all variables except for the ones specified by *query*. Adjust all distributions to reflect *evidence*. Return the entry that matches the exact probability of a specific event, as specified by *query*.
        
        Arguments:
            1. *query* -- A dict containing (key: value) pairs reflecting (variable: value) that represents what outcome to calculate the probability of. The value must be a list of values (for ordinary queries do a list of length one).
            2. *evidence* -- A dict containing (key: value) pairs reflecting (variable: value) evidence that is known about the system.
                    
        Attributes modified:
            1. *factorlist* -- Modified as in *condprobve*.
                           
        The function then chooses the entries of *factorlist* that match the queried event or events. It then operates on them to return the probability that the event (or events) specified will occur, represented as a float between 0 and 1.

        Note that in this function, queries of the type P((x=A or x=B) and (y=C or y=D)) are permitted. They are executed by formatting the *query* dictionary like so::

            {
                "x": ["A", "B"],
                "y": ["C", "D"]
            }
        
        Usage example: this code would answer the specific query that vertex ``Grade`` gets outcome ``A`` given that ``Letter`` has outcome ``weak``, in :doc:`this Bayesian network <unittestdict>`::

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
            query = dict(Grade=['A'])

            # load bayesian network
            bn = DiscreteBayesianNetwork(skel, nd)

            # load factorization
            fn = TableCPDFactorization(bn)

            # calculate probability distribution
            result = fn.specificquery(query, evidence)

            # output
            print result

        '''
        condprob = self.condprobve(query, evidence)

        # now self.factorlist contains the joint distribution across the
        # variables designated in query. next, we have to locate the entries
        # where the variables have values matching the query (e.g., where "Grade"
        # is "A" and "Intelligence" is "High"). because must loop once for each 
        # variable, and we don't know how many variables there are, we use 
        # recursion to iterate through the variables
        visited = {}
        rindices = {}
        findices = []
        
        # find corresponding numbers to possible values, store in rindices
        for var, alternatives in query.items():
            rindices[var] = []
            visited[var] = False
            for alternative_value in alternatives:
                # Look the alternative_value up in the Bayesian network and get its index
                rindices[var].append(self.bn.Vdata[var]["vals"].index(alternative_value))
        
        # define function to help iterate recursively through all combinations of variables
        def findentry(var, index):
            visited[var] = True 
        
            for x in range(len(rindices[var])):
                newindex = index + rindices[var][x] * condprob.stride[var]
                if (list(visited.values()).count(False) > 0):
                    i = list(visited.values()).index(False)
                    nextvar = list(visited.keys())[i]
                    findentry(nextvar, newindex)
                else:
                    # we've accounted for all variable assignments and found an entry
                    findices.append(newindex)
            visited[var] = False
            return
        
        # calculate all relevant entries
        findentry(list(visited.keys())[0], 0)
            
        # sum entries
        fanswer = 0
        for findex in findices:
            fanswer += condprob.vals[findex]
            
        # return result
        return fanswer

