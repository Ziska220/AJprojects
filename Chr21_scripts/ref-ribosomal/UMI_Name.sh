"""
import os
import sys
import doctest
from re import findall
from pysam import Samfile
from toolshed import nopen
from itertools import islice, izip
from collections import Counter, defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from ._version import __version__

os.environ['TERM'] = 'linux'

IUPAC = {"A":"A","T":"T","C":"C","G":"G","R":"GA","Y":"TC",
         "M":"AC","K":"GT","S":"GC","W":"AT","H":"ACT",
         "B":"GTC","V":"GCA","D":"GAT","N":"GATC"}


class Fastq(object):
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual)

    def __repr__(self):
        return "Fastq({name})".format(name=self.name)

    def __str__(self):
        return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                seq=self.seq, qual=self.qual)


def umi_from_name(name):
    """
    extract the UMI sequence from the read name.
    >>> umi_from_name("cluster_1017333:UMI_GCCGCA")
    'GCCGCA'
    """
    return findall(r'UMI_([\w]*)', name)[0].strip()
