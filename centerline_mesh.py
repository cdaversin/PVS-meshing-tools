import meshio
import numpy as np
import argparse
import gmsh
from os import path

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

# Use GMSH to remove duplicates if any
gmsh.initialize()
gmsh.open(output_filename + "_centerline_mesh.msh")

# Before removing duplicates
# Return nodeTags, coord, parametricCoord
initial_node_tags = np.array(gmsh.model.mesh.getNodes()[0])
initial_node_coords = np.array(gmsh.model.mesh.getNodes()[1])
gmsh.model.mesh.removeDuplicateNodes()
node_tags = np.array(gmsh.model.mesh.getNodes()[0])
node_coords = np.array(gmsh.model.mesh.getNodes()[1])
common_node_tags = np.intersect1d(initial_node_tags, node_tags)
removed_node_tags = np.delete(initial_node_tags, common_node_tags)
#print("removed = ", removed_node_tags)

gmsh.write(output_filename + "_centerline_mesh-noduplicates.msh")
gmsh.finalize()

# Convert to xdmf
meshio._cli.convert([output_filename + "_centerline_mesh-noduplicates.msh", output_filename + "_centerline_mesh.xdmf", "--input-format", "gmsh", "--output-format", "xdmf"])

# Update centerline data after duplicated nods removal
radius_data = np.load(output_filename + "_centerline_max_radius.npy")
radius_data = np.delete(radius_data, removed_node_tags)
np.save(path.join(output_filename + "_centerline_max_radius.npy"), radius_data)

torsion_data = np.load(output_filename + "_centerline_torsion.npy")
torsion_data = np.delete(torsion_data, removed_node_tags)
np.save(path.join(output_filename +  "_centerline_torsion.npy"), torsion_data)

curvature_data = np.load(output_filename + "_centerline_curvature.npy")
curvature_data = np.delete(curvature_data, removed_node_tags)
np.save(path.join(output_filename + "_centerline_curvature.npy"), curvature_data)

normal_data = np.load(output_filename + "_centerline_frenet_normal.npy")
normal_data = np.delete(normal_data, removed_node_tags, axis=0)
np.save(path.join(output_filename + "_centerline_frenet_normal.npy"), normal_data)

binormal_data = np.load(output_filename + "_centerline_frenet_binormal.npy")
binormal_data = np.delete(binormal_data, removed_node_tags, axis=0)
np.save(path.join(output_filename + "_centerline_frenet_binormal.npy"), binormal_data)

tangent_data = np.load(output_filename + "_centerline_frenet_tangent.npy")
tangent_data = np.delete(tangent_data, removed_node_tags, axis=0)
np.save(path.join(output_filename + "_centerline_frenet_tangent.npy"), tangent_data)

# Update inlet/outlets indices
with open(path.join(output_filename + "_centerline_inlets.txt"),'r') as inlet_file:
    inlets = np.loadtxt(inlet_file, dtype='int')
    inlets = np.atleast_1d(inlets)
with open(path.join(output_filename + "_centerline_inlets.txt"),'ab') as inlet_file:
    inlet_file.truncate(0)
    for inlet in inlets:
        np.savetxt(inlet_file, np.where(node_tags == inlet), fmt='%i')
        #np.savetxt(inlet_file, np.where(node_coords == initial_node_coords[inlet]), fmt='%i')

with open(path.join(output_filename + "_centerline_outlets.txt"),'r') as outlet_file:
    outlets = np.loadtxt(outlet_file, dtype='int')
    outlets = np.atleast_1d(outlets)
with open(path.join(output_filename + "_centerline_outlets.txt"),'ab') as outlet_file:
    outlet_file.truncate(0)
    for outlet in outlets:
        np.savetxt(outlet_file, np.where(node_tags == outlet), fmt='%i')
        #np.savetxt(outlet_file, np.where(node_coords == initial_node_coords[outlet]), fmt='%i')
