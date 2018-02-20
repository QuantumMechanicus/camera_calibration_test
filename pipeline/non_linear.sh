#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/A
out_dir=/home/danielbord/CLionProjects/camera_calibration/pipeline/bin/automatic_solver_results_A/
list_l=${out_dir}list.left
list_r=${out_dir}list.right
find $dir -name "points_*_left" | sed 's/^\(.*\)$/\1/g' | sort > $list_l
find $dir -name "points_*_right" | sed 's/^\(.*\)$/\1/g' | sort > $list_r
fundamental=`find $out_dir -name "*.f" | sort`
lambdas=`find $out_dir -name "*.l" | sort`
s_left=$(paste -s $list_l)
s_right=$(paste -s $list_r)
s_left="--lp "${s_left}
s_right="--rp "${s_right}
bin/non_linear_optimizer --i 1 $s_left $s_right --ff $fundamental --lf $lambdas --nl 2 --h 4912 --w 7360 --q 0.1
