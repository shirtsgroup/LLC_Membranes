#!/bin/bash

scontrol show job $1 > job
line_number=$(awk '/WorkDir=/{ print NR; exit}' job)
str=$(sed "${line_number}q;d" job)
rm job
dir=${str#*=}
cd $dir
