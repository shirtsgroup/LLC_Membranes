#!/usr/bin/env bash

# Do all relevant analysis on finished simulations

trajectory="wiggle.trr"
gro="wiggle.gro"
tpr="wiggle.tpr"
pore_components="C C1 C2 C3 C4 C5"  # components used to estimate pore size
xrd_frames=50

while getopts "b:x:y:z:r:t:x" opt; do
    case $opt in
    b) BUILD_MON=$OPTARG;;
    x) x=$OPTARG;;
    y) y=$OPTARG;;
    z) z=$OPTARG;;
    r) ring_restraints=$OPTARG;;
    t) tail_restraints=$OPTARG;;
    x) xrd_frames=$OPTARG;;
    esac
done

echo "All of the analysis!" > analysis.log
echo "" >> analysis.log
# modify the trajectory so all monomers stay whole
echo 0 | gmx trjconv -f ${trajectory} -o traj_whole.xtc -s ${tpr} -pbc whole
echo "Trajectory Converted" >> analysis.log
echo "" >> analysis.log
ConvertTraj.py -l ${xrd_frames} -t ${trajectory} -g ${gro} --avg_dims
echo "Last ${xrd_frames} frames of trajectory converted for X-ray Diffraction calculation using ConvertTraj.py" >> analysis.log
echo "Average dimension written to dims.txt" >> analysis.log
echo "" >> analysis.log
# Possibly need to extract the last frame and use it as the .gro but I think the molecules are made whole in the output
# .gro file so it should be fine. If there is a problem in the future it might be worth trying this
echo "Structure_char.py output" >> analysis.log
Structure_char.py -t ${trajectory} -g ${gro} --noshow --auto_exclude >> analysis.log
echo "" >> analysis.log
echo "poresize.py output" >> analysis.log
poresize.py -t traj_whole.xtc -g ${gro} -c ${pore_components} >> analysis.log
echo "" >> analysis.log
