##The purpose of this algorithm is to sort through the frames and find the best 
##combination of frames that cover all the atoms in the reference frames 
## and find the representative frames for each conformational state going forward


#### First let us import the necessary libraries and modules####
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx 

### Now we need to import the dataset of RMSD comparisons 
### between all frames

RMSD_matrix = pd.read_csv("RMSD_matrix.csv", index_col=0)

#### Rename the columns such that starting at column 1 
### we have number 1, then number 2, then number 3, etc. and column 0 is not labeled as a number 
# but rather as the index of the dataframe.
RMSD_matrix.columns = [f'Frame_{i}' for i in range(1, RMSD_matrix.shape[1]+1)]

frames = RMSD_matrix 
#### Once we have done this let us find 
#### the best combination of frames that cover all the atoms in 
# the reference frames and find the representative frames for each 
# conformational state going forward

def representative_frames(frames):
    ##### What we want to do is parse through the frame comparisons in each column
    ### and filter based on our threshold of 4.5 Angstroms. 
    # This will allow us to generate the first category/class with its representative
    # The goal is to then iterate and go to the next frame which
    # WAS NOT selected in the first class, use that as a new representative 
    # and generate the next class based on the same threshold of 4.5 Angstroms.
    # let u start with column number 1 (in this iteration)
    categories = []
    skipped_columns = []
    col = frames.shape[1]
    for i in range(col):
        eligible_values = frames[(col[i] <= 4.5)]

        eligible_values.columns = ["Frame", "RMSD"]
        categories[f'Class_{i}'] = eligible_values

        for j in range(col):
            if frames.shape[0] == col[j]:
                skipped_columns.append(col[j])
                continue
            else:
                representative_frames(frames[frames.index != col[j]])

    for name, categories, in categories.items():
        print(name)
        print(categories)
        print() 
    #### We also want to generate a graph 
    ### of the representative frames and their 
    ### relative distance from each other based on 
    ### the RMSD values. This is achieved using 
    ### networkx package which was installed earlier
    G = nx.Graph()
    G.add_nodes_from(categories)
    
    return categories, G 