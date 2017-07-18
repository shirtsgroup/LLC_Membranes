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

# modify the trajectory so all monomers stay whole
echo 0 | gmx trjconv -f ${trajectory} -o traj_whole.xtc -s ${tpr} -pbc whole
echo "Trajectory Converted"

# Simulate X-ray Diffraction
#main_gromacs.py -top ${gro} -traj traj_whole.xtc --lcscale 1.25 --cscale 0.05

# Structure_char.py
echo 'Calculating pore to pore distances'
p2p_output="$(Structure_char.py -t traj_whole.xtc -g ${gro} --noshow --auto_exclude --save --plot_avg)"  # capture all stdout in one variable
# extract all information we care about from p2p_ouput. http://stackoverflow.com/questions/19399238/under-bash-how-to-extract-two-numbers-by-parsing-a-string#comment28753888_19399302
autocorrelation="$(grep -Po '(?<=Maximum Autocorrelation Time: )[0-9]+\.[0-9]+' <<< ${p2p_output})"
p2p_avg="$(grep -Po '(?<=Pore to Pore distance: )[0-9]+\.[0-9]+' <<< ${p2p_output})"
p2p_std="$(grep -Po '(?<=Standard Deviation of Pore to Pore distances: )[0-9]+\.[0-9]+' <<< ${p2p_output})"
equil_p2p="$(grep -Po '(?<=Equilibration detected after )\d+' <<< ${p2p_output})"

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
upperlimit_stack="$(grep -Po '(?<=, )[0-9]+\.[0-9]+' <<< ${eclipse})"
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

# Create an input file that will cause vmd to produce the desired output
echo "color Display Background white" >> input.txt  # change background color to white
echo "axes location off" >> input.txt  # turn off axes
echo "display height 5" >> input.txt  # set size of window
echo "render TachyonInternal top_full.tga" >> input.txt  # render the scene
echo "mol selection name ${pore_components}" >> input.txt  # make an atom selection with only benzene rings
echo "mol modrep 0 0" >> input.txt  # display the new representation
echo "mol modstyle 0 0 CPK 0.5 0.5 10 10" >> input.txt  # make the drawing mode CPK with bond width = 0.5 and sphere size=0.5
echo "display height 4" >> input.txt  # change the display window size again
echo "render TachyonInternal top.tga" >> input.txt  # render the scene
echo "rotate x by 90" >> input.txt  # rotate the system to look at a side view
echo "render TachyonInternal side.tga" >> input.txt  # render the scene
echo "exit" >> input.txt  # exit vmdi

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
echo "\caption*{Average pore-to-pore distance : ${p2p_avg} $\pm$ ${p2p_std} nm. Equilibration detected after
        ${equil_p2p} ns}\label{fig:p2p}" >> analysis.tex
echo '\end{subfigure}' >> analysis.tex
echo '\begin{subfigure}{0.45\textwidth}' >> analysis.tex
echo '\includegraphics[width=\textwidth]{thickness.png}' >> analysis.tex
echo "\caption*{Average thickness : ${avg_thick} $\pm$ ${std_thick} nm. Equilibration detected after ${equil_frame_thick} ns}\label{fig:thickness}" >> analysis.tex
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
# Summarize results in a table
echo '\begin{center}' >> analysis.tex
echo '\begin{tabular}{|c|c|c|}' >> analysis.tex
echo '\hline' >> analysis.tex
echo '\bf{Measurement} & \bf{Value} & \bf{Equil Time} \\' >> analysis.tex
echo '\hline' >> analysis.tex
echo "Average Pore-to-Pore distance & ${p2p_avg} $\pm$ ${p2p_std} nm & ${equil_p2p} ns \\\\" >> analysis.tex  # need 4 backslashes in double quotes since each printed backslash needs to be escaped with a backslash
echo '\hline' >> analysis.tex
echo "Average Thickness & ${avg_thick} $\pm$ ${std_thick} nm & ${equil_frame_thick} ns \\\\" >> analysis.tex
echo '\hline' >> analysis.tex
echo "Average Overlap & ${avg_overlap} $\pm$ ${std_overlap} nm & ${overlap_equil} ns \\\\" >> analysis.tex
echo '\hline' >> analysis.tex
echo "Average stacking distance & ${pistack} [${lowerlimit_stack}, ${upperlimit_stack}] & ${pistack_equil} ns \\\\" >> analysis.tex
echo '\hline' >> analysis.tex
echo "Average pore size & ${avg_poresize} $\pm$ ${std_poresize} & ${poresize_equil} ns \\\\" >> analysis.tex
echo '\hline' >> analysis.tex
echo "Average order parameter & ${avg_order} & ${order_equil} ns \\\\" >> analysis.tex
echo '\hline' >> analysis.tex
echo '\end{tabular}' >> analysis.tex
echo '\end{center}' >> analysis.tex

echo '\end{document}' >> analysis.tex

pdflatex analysis.tex # compile the .tex document

gnome-open analysis.pdf  # open the compiled document for viewing
