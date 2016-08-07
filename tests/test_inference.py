"""Construct a small Bayes network and perform some inferences on it"""

import unittest

import libpgm.discretebayesiannetwork as discrete
import libpgm.nodedata as nodedata
import libpgm.tablecpdfactorization as factorization

class TestInference(unittest.TestCase):
    def setUp(self):
        nd = nodedata.NodeData()
        nd.Vdata = {
            "C": {
                "parents": [],
                "vals": ["C1", "C2", "C3"],
                "cprob": [0.1, 0.6, 0.3]},
            "Q1": {
                "parents": ["C"],
                "vals": ["A1.1", "A1.2", "A1.3"],
                "cprob": {
                    ("C1",): [0.5, 0.0, 0.5],
                    ("C2",): [0.0, 0.5, 0.5],
                    ("C3",): [0.5, 0.5, 0.0]}},
            "Q2": {
                "parents": ["C"],
                "vals": ["A2.1", "A2.2"],
                "cprob": {
                    ("C1",): [0.9, 0.1],
                    ("C2",): [0.1, 0.9],
                    ("C3",): [0.5, 0.5]}},
            "Q3": {
                "parents": ["C"],
                "vals": ["A", "B", "Both"],
                'cprob': {
                    ("C1",): [0.5,0.3,0.2],
                    ("C2",): [0.9,0.1,0.0],
                    ("C3",): [0.6,0.3,0.1]}},
            "Q4": {
                "parents": ["Q3"],
                "vals": ["Yes", "No", "N/A"],
                'cprob': {
                    ("A",):    [0.5,0.5,0.0],
                    ("Both",): [0.5,0.5,0.0],
                    ("B",):    [0.0,0.0,1.0]}}}

        self.tan = discrete.DiscreteBayesianNetwork(nodedata=nd)

    def test_condprobve(self):
        fact = factorization.TableCPDFactorization(self.tan)
        result = fact.condprobve(
            {"C"},
            {"Q4": "Yes"})
        ps = result.vals
        self.assertAlmostEqual(ps[0], 0.0853658536585365)
        self.assertAlmostEqual(ps[1], 0.6585365853658537)
        self.assertAlmostEqual(ps[2], 0.2560975609756097)

