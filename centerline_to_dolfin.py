from dolfin import *
import numpy as np
import argparse

# Configuration from arguments parsing ###########################################
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file_path', required=True, default="", help="Path to .xml.gz file (without extension))")
args = parser.parse_args()
prefix = args.input_file_path
##################################################################################

# Read mesh
mesh_file = prefix + "_centerline_mesh.xdmf"
mesh = Mesh()
with XDMFFile(MPI.comm_world, mesh_file) as xdmf:
    xdmf.read(mesh)

# Mark inlets and outlets
inlet_tag = 10
outlet_tag = 20
mf = MeshFunction('size_t',mesh, mesh.topology().dim()-1)
# Tag inlets
with open(prefix + "_centerline_inlets.txt", "r") as inlet_file:
    inlets = np.loadtxt(inlet_file, dtype='int')
for i, inlet in enumerate(inlets):
    mf[inlet] = int(inlet_tag) + i
# Tag outlets
with open(prefix + "_centerline_outlets.txt", "r") as outlet_file:
    outlets = np.loadtxt(outlet_file, dtype='int')
for i, outlet in enumerate(outlets):
    mf[outlet] = int(outlet_tag) + i

# Export markers (check)
mfile = XDMFFile(MPI.comm_world, "./C0075/markers.xdmf")
mfile.write(mf)
mfile.close()

# Load data
# Scalar
radius_data = np.load(prefix + "_centerline_max_radius.npy")
torsion_data = np.load(prefix + "_centerline_torsion.npy")
curvature_data = np.load(prefix + "_centerline_curvature.npy")
# Vectorial data
normal_data = np.load(prefix + "_centerline_frenet_normal.npy")
normal_data = np.concatenate(normal_data, axis = None)
binormal_data = np.load(prefix + "_centerline_frenet_binormal.npy")
binormal_data = np.concatenate(binormal_data, axis = None)
tangent_data = np.load(prefix + "_centerline_frenet_tangent.npy")
tangent_data = np.concatenate(tangent_data, axis = None)

# Function spaces
V = FunctionSpace(mesh, "CG", 1)
W = VectorFunctionSpace(mesh, "CG", 1)

# Radius
u = Function(V)
u.vector().set_local(radius_data)
# Torsion
t = Function(V)
t.vector().set_local(torsion_data)
# Curvature
c = Function(V)
c.vector().set_local(curvature_data)
# Normal (Frenet)
w = Function(W)
w.vector().set_local(normal_data)
# Binormal
w2 = Function(W)
w2.vector().set_local(binormal_data)
# Tangent
w3 = Function(W)
w3.vector().set_local(tangent_data)


# Export
field_file = XDMFFile(MPI.comm_world, prefix + "_centerline_radius.xdmf")
field_file.write(u,0)
field_file.close()
field_file = XDMFFile(MPI.comm_world, prefix + "_centerline_curvature.xdmf")
field_file.write(c,0)
field_file.close()
field_file = XDMFFile(MPI.comm_world, prefix + "_centerline_torsion.xdmf")
field_file.write(t,0)
field_file.close()

field_file = XDMFFile(MPI.comm_world, prefix + "_centerline_normal.xdmf")
field_file.write(w,0)
field_file.close()
field_file = XDMFFile(MPI.comm_world, prefix + "_centerline_binormal.xdmf")
field_file.write(w2,0)
field_file.close()
field_file = XDMFFile(MPI.comm_world, prefix + "_centerline_tangent.xdmf")
field_file.write(w3,0)
field_file.close()
