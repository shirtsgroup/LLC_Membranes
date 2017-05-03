#!/bin/bash

# Do all relevant analysis on finished simulations

trajectory="wiggle.trr"
gro="wiggle.gro"
tpr="wiggle.tpr"
pore_components="C C1 C2 C3 C4 C5"  # components used to estimate pore size
xrd_frames=50
build_mon='NAcarb11V'

while getopts "t:g:r:p:x:b:" opt; do
    case $opt in
    t) trajectory=$OPTARG;;
    g) gro=$OPTARG;;
    r) tpr=$OPTARG;;
    p) pore_components=$OPTARG;;
    x) xrd_frames=$OPTARG;;
    b) build_mon=$OPTARG;;
    esac
done

## modify the trajectory so all monomers stay whole
#echo 0 | gmx trjconv -f ${trajectory} -o traj_whole.xtc -s ${tpr} -pbc whole
#echo "Trajectory Converted"
#
## Convert the trajectory so that it is usable with XrayDiffraction.exe
#ConvertTraj.py -l ${xrd_frames} -t ${trajectory} -g ${gro} --avg_dims
#echo "Last ${xrd_frames} frames of trajectory converted for X-ray Diffraction calculation using ConvertTraj.py"
#echo "Average dimension written to dims.txt" >> analysis.log

# Structure_char.py
echo 'Calculating pore to pore distances'
p2p_output="$(Structure_char.py -t ${trajectory} -g ${gro} --noshow --auto_exclude --save)"  # capture all stdout in one variable
# extract all information we care about from p2p_ouput. http://stackoverflow.com/questions/19399238/under-bash-how-to-extract-two-numbers-by-parsing-a-string#comment28753888_19399302
autocorrelation="$(grep -Po '(?<=Maximum Autocorrelation Time: )[0-9]+\.[0-9]+' <<< ${p2p_output})"
p2p_avg="$(grep -Po '(?<=Pore to Pore distance: )[0-9]+\.[0-9]+' <<< ${p2p_output})"
p2p_std="$(grep -Po '(?<=Standard Deviation of Pore to Pore distances: )[0-9]+\.[0-9]+' <<< ${p2p_output})"
equil_frame_p2p="$(grep -Po '(?<=frame )\d+' <<< ${p2p_output})"
equil_percent="$(grep -Po '(?<=\()[0-9]+\.[0-9]+' <<< ${p2p_output})"

echo 'Calculating thickness'
# Thickness.py
thick="$(Thickness.py -g ${gro} -t traj_whole.xtc --noshow --save --trajectory)"
avg_thick="$(grep -Po '(?<=Average membrane thickness: )[0-9]+\.[0-9]+' <<< ${thick})"
std_thick="$(grep -Po '(?<=- )[0-9]+\.[0-9]+' <<< ${thick})"
equil_frame_thick="$(grep -Po '(?<=Equilibration detected after )\d+' <<< ${thick})"

echo 'Calculating overlap and pi-stacking distance'
eclipse="$(eclipse.py -g ${gro} -t traj_whole.xtc -b ${build_mon} --noshow --save)"
pistack="$(grep -Po '(?<=Average stacking distance: )[0-9]+\.[0-9]+' <<< ${eclipse})"
lowerlimit_stack="$(grep -Po '(?<=CI \[)[0-9]+\.[0-9]+' <<< ${eclipse})"
limit_stack="$(grep -Po '(?<=, )[0-9]+\.[0-9]+' <<< ${eclipse})"
avg_overlap="$(grep -Po '(?<=Average overlap: )[0-9]+\.[0-9]+' <<< ${eclipse})"
std_overlap="$(grep -Po '(?<=- )[0-9]+\.[0-9]+' <<< ${eclipse})"
pistack_equil="$(grep -Po '(?<=Pi-stacking distance equilibrated after )\d+' <<< ${eclipse})"
overlap_equil="$(grep -Po '(?<=Overlap equilibrated after )\d+' <<< ${eclipse})"

echo 'Calculating pore size and order parameter'
poresize="$(poresize.py -t traj_whole.xtc -g ${gro} -c ${pore_components} --save --noshow)"
poresize_equil="$(grep -Po '(?<=Pore size equilibrated after )\d+' <<< ${poresize})"
order_equil="$(grep -Po '(?<=Order parameter equilibrated after )\d+' <<< ${poresize})"
avg_poresize="$(grep -Po '(?<=Average Pore Size: )[0-9]+\.[0-9]+' <<< ${poresize})"
std_poresize="$(grep -Po '(?<=- )[0-9]+\.[0-9]+' <<< ${poresize})"
avg_order="$(grep -Po '(?<=Average Order Parameter: )[0-9]+\.[0-9]+' <<< ${poresize})"

vmd wiggle.gro -e input.txt  # take a few pictures of the membrane

# Convert tga files to 
for i in top_full top side; do  # convert the pictures to .png files using imagemagick
	convert ${i}.tga ${i}.png;
done

rm *.tga  # get rid of the .tga files

# write a .tex file organizing everything
echo '\documentclass{article}' > analysis.tex
echo '\usepackage{graphicx}' >> analysis.tex
echo '\usepackage{geometry}' >> analysis.tex
echo '\usepackage{subcaption}' >> analysis.tex
echo '\title{Analysis of LLC membrane simulation}' >> analysis.tex
echo '\author{Benjamin J. Coscia}' >> analysis.tex
echo '\date{\today}' >> analysis.tex
echo '\geometry{legalpaper, margin=0.5in}' >> analysis.tex
echo '\begin{document}' >> analysis.tex
echo '\maketitle' >> analysis.tex
# VMD screenshots
echo '\begin{figure}[h]' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\begin{subfigure}{0.3\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{top_full.png}' >> analysis.tex
echo '\end{subfigure}'>> analysis.tex
echo '\centering' >> analysis.tex
echo '\begin{subfigure}{0.3\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{top.png}' >> analysis.tex
echo '\end{subfigure}'>> analysis.tex
echo '\centering' >> analysis.tex
echo '\begin{subfigure}{0.3\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{side.png}' >> analysis.tex
echo '\end{subfigure}'>> analysis.tex
echo '\end{figure}' >> analysis.tex
# Pore-to-pore and Thickness
echo '\begin{figure}[h]' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\includegraphics[width=\textwidth]{p2p.png}' >> analysis.tex
echo "\caption*{Average pore-to-Pore distance : ${p2p_avg} $\pm$ ${p2p_std} nm. Equilibration detected after frame
        ${equil_frame} , ${equil_percent} \% into simulation}\label{fig:p2p}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{thickness.png}' >> analysis.tex
echo "\caption*{Average thickness : ${avg_thick} $\pm$ ${std_thick} nm. Equilibrationd detected after ${equil_frame_thick} ns}\label{fig:thickness}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\end{figure}' >> analysis.tex
# Overlap and pi-stacking
echo '\begin{figure}[h]' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\includegraphics[width=\textwidth]{overlap.png}' >> analysis.tex
echo "\caption*{Average overlap : ${avg_overlap} $\pm$ ${std_overlap} nm. Equilibration detected after
        ${overlap_equil} ns}\label{fig:overlap}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{pistack.png}' >> analysis.tex
echo "\caption*{Average stacking distance: ${pistack} 95 \% CI [${lowerlimit_stack}, ${upperlimit_stack}].
        Equilibration detected after ${pistack_equil} ns}\label{fig:pistack}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\end{figure}' >> analysis.tex
# Pore size and order parameter
echo '\begin{figure}[h]' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\centering' >> analysis.tex
echo '\includegraphics[width=\textwidth]{poresize.png}' >> analysis.tex
echo "\caption*{Average pore size : ${avg_poresize} $\pm$ ${std_poresize} nm. Equilibration detected after
        ${poresize_equil} ns}\label{fig:poresize}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{order.png}' >> analysis.tex
echo "\caption*{Average order parameter: ${avg_order}. Equilibration detected after ${order_equil} ns}\label{fig:order}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\end{figure}' >> analysis.tex

echo '\end{document}' >> analysis.tex

pdflatex analysis.tex # compile the .tex document

gnome-open analysis.pdf  # open the compiled document for viewing
