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
A module for creating and managing node data. Node data in this library can have many types, dependent on whether the conditional probability distributions are discrete, linear Gaussian, or hybrid, and on whether the Bayesian network is static or dynamic. For example input files, see :doc:`unittestdict`, :doc:`unittesthdict`, :doc:`unittestlgdict`, and :doc:`unittestdyndict`.

'''

import re
import json
from . import graphskeleton

class NodeData:
    '''This class represents the node data for each node in a graph.

    If the Bayesian network is static, it contains the attribute
    *Vdata*.

    If the Bayesian network is dynamic, it contains two attributes,
    *initial_Vdata* and *twotbn_Vdata*.

    If the Bayesian network has hybrid CPDs, it contains the
    additional attribute *nodes*.

    '''
    @classmethod
    def load(k, path):
        '''
        Load node data from an input file located at *path*. Input file must be a plaintext .txt file with a JSON-style representation of a dict. The dict must have the top-level key ``Vdata`` or two top-level keys, ``initial_Vdata`` and ``twotbn_Vdata``. For example:: 

            {
                "Vdata": {
                    "<vertex 1>": <dict containing vertex 1 data>,
                    ...
                    "<vertex n>": <dict containing vertex n data>
                }
            }

        or::

            {
                "initial_Vdata": {
                    "<vertex 1>": <dict containing vertex 1 data>,
                    ...
                    "<vertex n>": <dict containing vertex n data>
                }
                "twotbn_Vdata": {
                    "<vertex 1>": <dict containing vertex 1 data>,
                    ...
                    "<vertex n>": <dict containing vertex n data>
                }
            }
            
        The function takes the following arguments:
            1. *path* -- The path to the text file that contains input data (e.g., "mydictionary.txt")
        
        In the static case, it modifies *Vdata* to hold the dictionary found at path. In the dynamic case, it modifies the *initial_Vdata* and *twotbn_Vdata* attributes to hold the dictionaries found at path.
            
        '''
        alldata = json.load(open(path))

        # try to load both for normal and dynamic cases
        self = k()
        try: 
            self.Vdata = self.parse_cprob(alldata["Vdata"])
        except KeyError:
            self.initial_Vdata = self.parse_cprob(
                alldata["initial_Vdata"])
            self.twotbn_Vdata = self.parse_cprob(
                alldata["twotbn_Vdata"])
        return self

    def parse_cprob(self, vdata):
        """Parse a cprob dictionary.

        Convert a cprob dictionary provided in json into a dictionary
        indexed by parent-value tuples.

        """

        inlist = re.compile(r"^\w*\['(.*)'\]\w*$")
        comma = re.compile(r"', *'")
        for node, props in vdata.items():
            new_cprob = {}
            if "cprob" in props:
                # It is a discrete node
                cprob = props["cprob"]
                try:
                    for key, ps in cprob.items():
                        content = inlist.match(key).group(1)
                        key = tuple(comma.split(content))
                        new_cprob[key] = ps
                    props["cprob"] = new_cprob
                except AttributeError:
                    pass
        return vdata

    @property
    def V(self):
        try:
            return self._graphskeleton.V
        except AttributeError:
            try:
                V = list(self.Vdata)
            except (TypeError, AttributeError):
                V = list(self.initialVdata)
            self._graphskeleton = graphskeleton.GraphSkeleton(V, self.E)
            self._graphskeleton.toporder()
            return self._graphskeleton.V
            

    @property
    def E(self):
        try:
            return self._E
        except AttributeError:
            self._E = []
            for (node, properties) in self.Vdata.items():
                parents = properties["parents"]
                for parent in parents:
                    self._E.append((parent, node))
            return self._E
            
class StaticNodeData(NodeData):
    def __init__(self, Vdata={}):
        self.Vdata = Vdata
        '''A dictionary of node data.'''

class DynamicNodeData(NodeData):
    def __init__(self, initial_Vdata, dwotbn_Vdata):
        self.initial_Vdata = {}
        '''In dynamic graphs, a dictionary containing node data for the initial time interval.'''
        self.twotbn_Vdata = {}
        '''In dynamic graphs, a dictionary containing node data for every time step interval after the first one.'''

class HybridNodeData(NodeData):
    def __init__(self, nodes={}):
        self.nodes = nodes
        '''In hybrid graphs, a dictionary of {key:value} pairs linking the name of each node (the key) to a clas instance (the value) which represents the node, its data, and its sampling function.'''
        
    def entriestoinstances(self):
        '''
        For each node, convert dictionary entry to class instance.

        This method is used only when dealing with Hybrid Bayesian networks as found in the :doc:`hybayesiannetwork` module. 

        The type of the node must be located in the 'type' attribute of the node's dictionary entry. To see an example of such a dictionary, see :doc:`unittesthdict`. This type is used to instantiate a corresponding class from libpgm/CPDtypes/, and store the node's dictionary info in that class. Thus we lose none of the dictionary data, yet we gain the ability to use the instantiated class's built-in function to choose its own outcome based on the outcomes of its parents. 

        In order for this method to be called, the self.Vdata attribute must have dictionary entries of the following form::

            <vertex name>: {
                'type': <type of node -- must correspond to module in /CPDtypes>,
                'parents': <array of parents of node>,
                'children': <array of children of node>,
                <whatever further entries are required by the type*> 
            }

        For instance, type "discrete" requires a "cprob" entry, while type "lg"
        requires "mean_base", "mean_scal", and "variance" entries.

        The function draws on the data in the *Vdata* attribute, and instantiates the attribute *nodes*, which is a dictionary of {name: instance} pairs where 'name' is the name of the node and 'instance' is a class instance containing the node data and the proper sampling function.

        '''
        # declare result dict
        rarray = dict()

        # transform into class instances
        for entry in self.Vdata.keys():
            
            # import module containing class
            path = str(self.Vdata[entry]["type"])
            exec("from libpgm.CPDtypes import " + path)
            
            # instantiate class
            exec("tmpnode = " + path + "." + str.capitalize(path) + "(self.Vdata[entry])")

            # append to array
            exec("rarray['" + str(entry) + "'] = tmpnode")

        self.nodes = rarray
