import meshio
import numpy as np
import argparse
import gmsh

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
gmsh.model.mesh.removeDuplicateNodes()
#gmsh.model.mesh.generate(1) # Remeshing renumbers the nodes
gmsh.write(output_filename + "_centerline_mesh-noduplicates.msh")
gmsh.finalize()

# Convert to xdmf
meshio._cli.convert([output_filename + "_centerline_mesh-noduplicates.msh", output_filename + "_centerline_mesh.xdmf", "--input-format", "gmsh", "--output-format", "xdmf"])
