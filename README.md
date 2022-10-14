
# Mesher for arteries and perivascular spaces (PVS) from VTP data

PVS-meshing-tools is a python-based package for meshing image-based blood vessels (VTP format) and possibly the perivascular space
surrounding them. PVS-meshing-tools is based on :
* [VMTK](http://www.vmtk.org/), a software dedicated to 3D reconstruction, geometric analysis, mesh generation and surface data analysis for image-based modeling of blood vessels
* [FEniCS](https://fenicsproject.org/download/archive/), an open-source computing platform for solving partial differential equations (PDEs).

## Installation
We provide a [Docker container](https://github.com/cdaversin/PVS-meshing-tools/pkgs/container/pvs-meshing-tools) including the aformentioned dependencies. VMTK uses an interactive GUI, which requires the use of specific display options with Docker.
On Linux, the appropriate command to run the Docker is:
```
docker run -ti --network=host -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix --privileged ghcr.io/cdaversin/pvs-meshing-tools:latest
```
The PVS-meshing-tools repository is included into the container.

## Meshing 

Once in the Docker container, the main bash script `generate_mesh.sh` can be configured with
the correct input data :
- Input/output files
- Mesh resolution
- Additional extensions, allowing to get "smoother" inlets and outlets surfaces if needed
- Perivascular space generation options (thickness factor, mesh resolution)
- Centerline mesh (optional)

### Input file

The input file is given as a VTP file (see example [here](https://github.com/cdaversin/PVS-meshing-tools/blob/master/C0075/C0075_clip2.vtp)), obtained from the (former) [AneuriskWeb](http://ecm2.mathcs.emory.edu/aneuriskweb/index) repository.
Note that [Paraview](https://www.paraview.org/) offers the possibility to clip geometries and generate the resulting VTP file : 
* Open <input_file>.vtp in Paraview
* Use "Clip" feature to select the part of interest
* Apply "Extract Surface" filter on that "Clip"
* Save Data - VTP format
    
The variable `VTP_FILE` takes the path of your input geometry (with VTP format) :
```
VTP_FILE=<path_to_input_file>/<input_file>.vtp
```

### Extensions

VMTK offers the possibility to add extensions to the inlets and outlets, to make those surfaces smoother when needed.
To do so, you can set the `EXT` variable to `true` and set the name of the output VTP file (including the extensions) as `VTP_FILE_EXT`
```
EXT=<true or false>
VTP_FILE_EXT=<path_to_input_file>/<input_file_ext>.vtp
```

### Output

You need to specify a path and a prefix (without file extension) for the generated output files (meshes and post-processing data), for example:
```
`OFILE=$PWD/myresults`
```
The mesh file will be generated from this prefix, in XDMF format, together with tag data containing volume and surface tags.
These can be visualized with [Paraview](https://www.paraview.org/).
* Volume tags (`<prefix>_mc.xdmf`)
   + Vessel tag : 0
   + PVS tag : 1
   
* Surface tags (`<prefix>_mf.xdmf`)
   + Vessel inlets and outlets : 10, 11, 12, ... 
   + PVS inlets and outlets : 20, 21, 22, ...
   + Vessel wall : 30
   + PVS wall : 40
   
### Mesh resolution

The parameters `EDGE_LENGTH` and `EDGE_FACTOR` can be modified to setup the mesh resolution, and the optimal parameter varies
from geometry to geometry.
```
EDGE_LENGTH=0.5
EDGE_FACTOR=0.5
```

### Perivascular mesh

`PVS-meshing-tools` generates the mesh of both the vessel and the perivascular space (PVS).
You can specify the desired thickness of the PVS mesh, as a ratio of the vessel diameter. Note that the `PVS_THICKNESS` takes a list of
thicknesses, and will generate as many meshes as PVS thicknesses given in the list. If `NON_CONST_THICKNESS` is set to true, the thickness of the PVS will be proportinal to the vessel diameter (at every cross section), with the given proportional factor. If `NON_CONST_THICKNESS` is set to false, the PVS thickness will be constant 

The number of sublayers (impacting the number of elements in the PVS thickness), can also be customized with `NB_SUB_LAYERS` Example :
```
PVS_THICKNESS=(0.95)
NB_SUB_LAYERS=3
NON_CONST_THICKNESS=true
```

It is possible to export only the PVS mesh as output (i.e. remove the vessel mesh), with the option
```
PVS_ONLY=<true or false>
```

### Centerline mesh (optional, default:false)

The variable `CENTERLINE_MESH` configure the building of the 1D centerline mesh of the input vessel.
```
CENTERLINE_MESH=<true or false>
```
If set to true, a 1D mesh of the centerline will be generated, together with related data :
* centerline_markers (`<prefix>_centerline_markers.xdmf`):
   + inlet tags : `10, 11, 12, ...`
   + outlet tags : `20, 21, 22, ...`
* Other relevant data (in `<prefix>_XDMF` directory)
   + Maximum inscribed sphere radius
   + Normal / binormal, and tangent (Frenet)
   + Curvature
   + Torsiom

Once everything configured in the main script `generate_mesh.sh`, just use the command :
```
bash generate_meshes.sh
```
