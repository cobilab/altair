import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import sys
import argparse
import pandas as pd
import numpy as np
from numba import jit, prange
from skbio.tree import nj
from skbio import DistanceMatrix
from ete3 import Tree as EteTree, TreeStyle, NodeStyle, AttrFace, faces

@jit(nopython=True, parallel=True)
def create_ncd_matrix_fast(ncd_values, n):
    matrix = np.zeros((n, n), dtype=np.float32)
    for i in prange(n):
        for j in prange(i + 1, n):
            matrix[i, j] = matrix[j, i] = ncd_values[i]
    return matrix

def read_data(file_path, num_closest_nodes=None):
    data = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1, 2], dtype={0: str, 1: np.float32, 2: str})
    data.columns = ['Sequence_Number', 'NCD', 'Location_Date']
    # Replace spaces with underscores in labels to comply with Newick format
    data['Label'] = data['Sequence_Number'].str.replace(' ', '') + '|' + data['Location_Date'].str.replace(' ', '')
    if num_closest_nodes is not None and num_closest_nodes < len(data):
        return find_closest_nodes(data, num_closest_nodes)
    return data

def add_origin_entry(data):
    # Create a fake origin entry
    origin_entry = pd.DataFrame([['Origin', 0, '']], columns=['Sequence_Number', 'NCD', 'Location_Date'])
    origin_entry['Label'] = 'Origin||'
    return pd.concat([origin_entry, data], ignore_index=True)


def find_closest_nodes(data, x):
    return data.sort_values(by='NCD').head(x)

def construct_tree(matrix, data):
    labels = data['Label'].tolist()
    dist_matrix = DistanceMatrix(matrix, labels)
    tree = nj(dist_matrix)
    return tree

def plot_tree(tree_newick, data, output_file):
    ete_tree = EteTree(tree_newick, format=1)
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = True
    ts.show_branch_support = False
    ts.scale = 120
    ts.mode = "c"
    ts.layout_fn = lambda node: custom_layout(node,data)

    # Render and save to file directly without trying to display
    ete_tree.render(output_file, tree_style=ts, w=183, units="mm")


def custom_layout(node, data):
    if node.is_leaf():
        
        node_name = node.name
        parts = node_name.split('|')

        location = parts[1] if len(parts) > 1 else "Unknown Location"
        date = parts[2] if len(parts) > 2 else "Unknown Date"
        sequence_label = f"Seq: {parts[0]}, From: {location}, {date}"
        
        # Include NCD value in label
        ncd_value = data.loc[data['Label'] == node_name, 'NCD'].iloc[0]
        ncd_value_print = format(data.loc[data['Label'] == node_name, 'NCD'].iloc[0], '.5f')
        sequence_label += f", NCD: {ncd_value_print}"
        # Change text color for closest branch leaf nodes
        text_color = "black"  # Default color
        if ncd_value == data['NCD'].min() and ncd_value != 0:
            text_color = "green"

        sequence_face = faces.TextFace(sequence_label, fgcolor=text_color)
        node.add_face(sequence_face, column=0, position="branch-right")




def main(file_path, num_closest_nodes=None):
    try:
        data = read_data(file_path, num_closest_nodes)
        n = len(data)
        ncd_values = data['NCD'].values
        ncd_matrix = create_ncd_matrix_fast(ncd_values, n)
        tree = construct_tree(ncd_matrix, data)
        tree_newick = tree.__str__()
        output_file = os.path.splitext(file_path)[0] + '-tree.pdf'
        plot_tree(tree_newick, data, output_file)
        print(f"Tree saved to {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a phylogenetic tree based on NCD values.")
    parser.add_argument("file_path", type=str, help="Path to the input data file")
    parser.add_argument("-N", "--nodes", type=int, help="Number of closest nodes to use", default=None)
    args = parser.parse_args()
    main(args.file_path, args.nodes)
