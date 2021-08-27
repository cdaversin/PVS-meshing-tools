import meshio
import numpy as np
import argparse
import gmsh
from os import path
import sys

# Configuration from arguments parsing ###########################################
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output_file', required=True, default="", help="Name of the output file(s) [without extension]") # XDMF
args = parser.parse_args()
output_filename = args.output_file
##################################################################################

# Load centerline coords
centerline_points = np.load(output_filename + "_centerline_points.npy")
# Load cells array
cells_array= np.load(output_filename + "_centerline_cells_array.npy")

def create_line_mesh(points, cells_array, filename) :
    centerline_cells = {"line": cells_array}
    meshio.write_points_cells(
        filename,
        points,
        centerline_cells,
        # Optional
        # point_data=point_data,
        # cell_data=cell_data,
    )

create_line_mesh(centerline_points, cells_array, output_filename + "_centerline_mesh.msh")

### ------------ Remove duplicates if any (GMSH API) ----------- ###
gmsh.initialize()
gmsh.open(output_filename + "_centerline_mesh.msh")
# Node tags
initial_node_tags = np.array(gmsh.model.mesh.getNodes()[0])
# Node coords (re-arranged as numpy array [ [x1,x2,x3]_1, ..., [x2,x3,x3]_N])
initial_node_coords = np.array(gmsh.model.mesh.getNodes()[1])
initial_node_coords = np.split(initial_node_coords, len(initial_node_coords)/3.0)
# Removing duplicate nodes
gmsh.model.mesh.removeDuplicateNodes()
# Updated node tags
node_tags = np.array(gmsh.model.mesh.getNodes()[0])
# Updated node coords (re-arranged as numpy array [ [x1,x2,x3]_1, ..., [x2,x3,x3]_N])
node_coords = np.array(gmsh.model.mesh.getNodes()[1])
node_coords = np.split(node_coords, len(node_coords)/3.0)
# Getting array of removed nodes to update centerline data
_, common_idx_init, _ = np.intersect1d(initial_node_tags, node_tags, return_indices=True)
removed_node_tags = np.array([]) # Tags of remote nodes
removed_idx_init = np.array([]) # Indices of these removed nodes in the initial array
if len(node_tags) != len(initial_node_tags): # Some nodes have been removed
    removed_node_tags = np.delete(initial_node_tags, common_idx_init)
    _, removed_idx_init, _ = np.intersect1d(initial_node_tags, removed_node_tags, return_indices=True)

# Writing updated mesh
gmsh.write(output_filename + "_centerline_mesh-noduplicates.msh")
gmsh.finalize()
# Converting to xdmf
meshio._cli.convert([output_filename + "_centerline_mesh-noduplicates.msh", output_filename + "_centerline_mesh.xdmf", "--input-format", "gmsh", "--output-format", "xdmf"])
### ------------------------------------------------------------- ###

### ------- Update centerline data (mesh without duplicates) ---- ###
if len(removed_node_tags) != 0:
    radius_data = np.load(output_filename + "_centerline_max_radius.npy")
    radius_data = np.delete(radius_data, removed_idx_init)
    np.save(path.join(output_filename + "_centerline_max_radius.npy"), radius_data)

    torsion_data = np.load(output_filename + "_centerline_torsion.npy")
    torsion_data = np.delete(torsion_data, removed_idx_init)
    np.save(path.join(output_filename +  "_centerline_torsion.npy"), torsion_data)

    curvature_data = np.load(output_filename + "_centerline_curvature.npy")
    curvature_data = np.delete(curvature_data, removed_idx_init)
    np.save(path.join(output_filename + "_centerline_curvature.npy"), curvature_data)

    normal_data = np.load(output_filename + "_centerline_frenet_normal.npy")
    normal_data = np.delete(normal_data, removed_idx_init, axis=0)
    np.save(path.join(output_filename + "_centerline_frenet_normal.npy"), normal_data)

    binormal_data = np.load(output_filename + "_centerline_frenet_binormal.npy")
    binormal_data = np.delete(binormal_data, removed_idx_init, axis=0)
    np.save(path.join(output_filename + "_centerline_frenet_binormal.npy"), binormal_data)

    tangent_data = np.load(output_filename + "_centerline_frenet_tangent.npy")
    tangent_data = np.delete(tangent_data, removed_idx_init, axis=0)
    np.save(path.join(output_filename + "_centerline_frenet_tangent.npy"), tangent_data)

# Update inlet/outlets indices
inlets_arr = []
outlets_arr = []
with open(path.join(output_filename + "_centerline_inlets.txt"),'r') as inlet_file:
    inlets = np.loadtxt(inlet_file, dtype='int')
    inlets = np.atleast_1d(inlets)
with open(path.join(output_filename + "_centerline_inlets.txt"),'ab') as inlet_file:
    inlet_file.truncate(0)
    for inlet in inlets:
        inlet_index = np.where(node_coords == initial_node_coords[inlet])[0][0]
        if inlet_index not in inlets_arr:
            np.savetxt(inlet_file, [inlet_index], fmt='%i')
            inlets_arr.append(inlet_index)

with open(path.join(output_filename + "_centerline_outlets.txt"),'r') as outlet_file:
    outlets = np.loadtxt(outlet_file, dtype='int')
    outlets = np.atleast_1d(outlets)
with open(path.join(output_filename + "_centerline_outlets.txt"),'ab') as outlet_file:
    outlet_file.truncate(0)
    for outlet in outlets:
        outlet_index = np.where(node_coords == initial_node_coords[outlet])[0][0]
        if outlet_index not in outlets_arr:
            np.savetxt(outlet_file, [np.where(node_coords == initial_node_coords[outlet])[0][0]], fmt='%i')
            outlets_arr.append(outlet_index)

### ------------------------------------------------------------- ###
