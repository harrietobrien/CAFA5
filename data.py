import os
import json
import time
from PIL import Image
from typing import Dict
from collections import Counter
import random
# import cv2
import obonet
import networkx
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from Bio import SeqIO
from pyvis.network import Network

CAFA5 = "/Users/harrietobrien/Desktop/CAFA5"
TRAIN = "{home}/Train".format(home=CAFA5)


class CFG:
    gobasic = "{path}/go-basic.obo".format(path=TRAIN)
    sequences = "{path}/train_sequences.fasta".format(path=TRAIN)
    terms = "{path}/train_terms.tsv".format(path=TRAIN)
    taxonomy = "{path}/train_taxonomy.tsv".format(path=TRAIN)
    ia = "{path}/IA.txt".format(path=CAFA5)


class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


# obonet package to load this file
graph = obonet.read_obo(CFG.gobasic)
print(f"Number of nodes: {len(graph)}")
print(f"Number of edges: {graph.number_of_edges()}")
term = "GO:0034655"
# nucleobase-containing compound catabolic process node properties
print(graph.nodes[term])
print(type(graph))  # --> <class 'networkx.classes.multidigraph.MultiDiGraph'>
print(type(graph.nodes[term]))  # --> <class 'dict'>
# plt.plot_dag(graph, term, radius=1)
# plot_dag(graph, term, radius=1000)

# analyze  protein sequences from sequences file with Biopython package
print("Sequence example:\n\n", next(iter(SeqIO.parse(CFG.sequences, "fasta"))))
#  count the number of sequences
sequences = SeqIO.parse(CFG.sequences, "fasta")
num_sequences = sum(1 for seq in sequences)
print("Number of sequences:", num_sequences)

# plot the length distribution of the protein sequences
sequences = SeqIO.parse(CFG.sequences, "fasta")
# get the length of each sequence
lengths = [len(seq) for seq in sequences]
fig = px.histogram(x=lengths, nbins=1000, color_discrete_sequence=['goldenrod'])
fig.update_layout(
    title={
        'text': "Distribution of protein sequence lengths",
        'y': 0.95,
        'x': 0.5,
        'xanchor': 'center',
        'yanchor': 'top'
    },
    xaxis_title="Sequence length", yaxis_title="Count"
)
fig.show()
np.percentile(lengths, 99)

# calculate the amino acid composition of each protein sequence
records = SeqIO.parse(CFG.sequences, "fasta")
# create a list of all amino acids in the sequences
aa_list = [aa for record in records for aa in record.seq]
# count the frequency of each amino acid
aa_count = Counter(aa_list)
fig = px.bar(
    x=list(aa_count.values()), y=list(aa_count.keys()),
    color_discrete_sequence=['darkslateblue'],
    orientation='h', height=700
)
fig.update_layout(
    title={
        'text': "Amino Acid Composition",
        'y': 0.95,
        'x': 0.5,
        'xanchor': 'center',
        'yanchor': 'top'
    },
    xaxis_title="Frequency", yaxis_title="Amino Acid"
)
fig.show()
# Load the train terms dataframe
train_terms = pd.read_csv(CFG.terms, sep="\t")
train_terms.head()
train_terms.describe()
# TODO: Join train_terms_df.EntryID on train_taxonomy_df.EntryID
# Information accretion
limit = 10

with open(CFG.ia) as f:
    ia_weights = [x.replace("\n", "").split("\t") for x in f.readlines()]

ia_weights[:limit]