#!/bin/bash

VTP_FILE=$(pwd)/C0075/C0075_clip2.vtp
VTP_FILE_EXT=$(pwd)/C0075/C0075_clip2_ext.vtp
OFILE=$(pwd)/C0075/C0075_clip2_mesh1

# VTP_FILE=$(pwd)/C0092/C0092_clip1_0.vtp
# VTP_FILE_EXT=$(pwd)/C0092/C0092_clip1_ext.vtp
# OFILE=$(pwd)/C0092/C0092_clip1_mesh1

# Extensions to get "smoother" inlet/outlets surfaces
#EXT=true
EXT=false

# Extraction of the PVS mesh
PVS_ONLY=true
#PVS_ONLY=false

# Centerline mesh
CENTERLINE_MESH=true
REFINE=1

# Mesh size
NON_CONST_THICKNESS=true
# Coarse
# EDGE_LENGTH=0.6
# EDGE_FACTOR=0.6
# Fine
EDGE_LENGTH=0.4
EDGE_FACTOR=0.4

PVS_THICKNESS=(0.95)
# Fine : 3 / coarse : 2
#NB_SUB_LAYERS=2
NB_SUB_LAYERS=3

if [ "$EXT" != false ]; then
    vmtksurfacereader -ifile $VTP_FILE \
                      --pipe vmtkcenterlines -seedselector openprofiles \
                      --pipe vmtkflowextensions -adaptivelength 1 -extensionratio 1.5 -normalestimationratio 1 -interactive 1 \
                      --pipe vmtksurfacewriter -ofile $VTP_FILE_EXT || exit 1
    VTP_FILE=$VTP_FILE_EXT;
fi

for i in "${PVS_THICKNESS[@]}"
do
   : 
   echo " --> Generating meshes from $VTP_FILE with pvs thickness = $i <--"
   if [ "$NON_CONST_THICKNESS" == true ];then
       OFILE="${OFILE}_${i}_ratio"
       if [ "$CENTERLINE_MESH" == true ];then
           python2 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --non_const_thickness --centerline_mesh --target_edge_factor $EDGE_FACTOR --nb_sub_layers $NB_SUB_LAYERS || exit 1
           #python3 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --non_const_thickness --centerline_mesh --target_edge_factor $EDGE_FACTOR --nb_sub_layers $NB_SUB_LAYERS || exit 1
       else
           python2 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --non_const_thickness --target_edge_factor $EDGE_FACTOR --nb_sub_layers $NB_SUB_LAYERS || exit 1
           #python3 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --non_const_thickness --target_edge_factor $EDGE_FACTOR --nb_sub_layers $NB_SUB_LAYERS || exit 1
       fi
   else
       OFILE="${OFILE}_${i}_cst"
       if [ "$CENTERLINE_MESH" == true ];then
           python2 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --centerline_mesh --refine_centerline_mesh $REFINE --target_edge_length $EDGE_LENGTH --nb_sub_layers $NB_SUB_LAYERS || exit 1
           #python3 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --centerline_mesh --target_edge_length $EDGE_LENGTH --nb_sub_layers $NB_SUB_LAYERS || exit 1
       else
           python2 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --target_edge_length $EDGE_LENGTH --nb_sub_layers $NB_SUB_LAYERS || exit 1
           #python3 $PWD/generate_xml.py -i $VTP_FILE -o $OFILE -t $i --target_edge_length $EDGE_LENGTH --nb_sub_layers $NB_SUB_LAYERS || exit 1
       fi
   fi

   if [ "$CENTERLINE_MESH" == true ];then
       echo " --> Generating centerline mesh (1D) <--"
       python3 $PWD/centerline_mesh.py -o ${OFILE} || exit 1
       python3 $PWD/centerline_to_dolfin.py -i $OFILE || exit 1
   fi

   echo " --> Converting into XDMF format <--"
   if [ "$PVS_ONLY" == true ];then
       python3 $PWD/xml_to_xdmf.py -i $OFILE.xml.gz -o ${OFILE}_PVS --pvs_only || exit 1
   else
       python3 $PWD/xml_to_xdmf.py -i $OFILE.xml.gz -o $OFILE || exit 1
   fi
done
