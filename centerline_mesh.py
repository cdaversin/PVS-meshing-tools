import meshio
import numpy as np
import argparse

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

create_line_mesh(centerline_points, cells_array, output_filename + "_centerline_mesh.xdmf")
