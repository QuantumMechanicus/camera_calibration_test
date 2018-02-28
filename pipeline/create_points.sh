#!/usr/bin/env bash



camera_name=A
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/${camera_name}
out_dir=/home/danielbord/CLionProjects/camera_calibration/pipeline/bin/${camera_name}_data/
keypoints_left=${out_dir}path.keypoints.left
keypoints_right=${out_dir}path.keypoints.right
intr_left=${out_dir}path.left.intrinsics
intr_right=${out_dir}path.right.intrinsics
extr_left=${out_dir}path.left.extrinsics
extr_right=${out_dir}path.right.extrinsics
fund=${out_dir}path.fundamental
motion=${out_dir}path.motion
info_left=${out_dir}path.left.info
info_right=${out_dir}path.right.info
paste_all=${out_dir}paste.all
temp=${out_dir}temp.file

mkdir -p $out_dir
find $dir -name "points_*_left" | sed 's/^\(.*\)$/Left keypoints: \1#/g' | sort > $keypoints_left
find $dir -name "points_*_*_right" | sed 's/^\(.*\)$/Right keypoints: \1#/g' | sort > $keypoints_right
find $dir -name "points_*_*_left" | grep -o '[0-9]\+_[0-9]\+' | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Left intrinsics:${out_dir}\1_intrinsics# |g" | sort > $intr_left
find $dir -name "points_*_*_right" | grep -o '[0-9]\+_[0-9]\+'| sed "s|^\([0-9]\+\)_\([0-9]\+\)|Right intrinsics:${out_dir}\2_intrinsics# |g" | sort > $intr_right
find $dir -name "points_*_*_left" | grep -o '[0-9]\+_[0-9]\+'  | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Left extrinsics:${out_dir}\1_extrinsics# |g" | sort > $extr_left
find $dir -name "points_*_*_right" | grep -o '[0-9]\+_[0-9]\+' | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Right extrinsics:${out_dir}\2_extrinsics# |g" | sort > $extr_right
find $dir -name "points_*_*_left" | grep -o '[0-9]\+_[0-9]\+' | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Fundamental matrix:${out_dir}\1_\2_fundamental_matrix# |g" | sort > $fund
find $dir -name "points_*_*_right" | grep -o '[0-9]\+_[0-9]\+' | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Relative motion:${out_dir}\1_\2_relative_motion# |g" | sort > $motion
find $dir -name "points_*_*_left" | grep -o '[0-9]\+_[0-9]\+' | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Left camera info:${out_dir}\1_info# | g" | sort > $info_left
find $dir -name "points_*_*_right" | grep -o '[0-9]\+_[0-9]\+' | sed "s|^\([0-9]\+\)_\([0-9]\+\)|Right camera info:${out_dir}\2_info |g"| sort > $info_right


paste $keypoints_left $keypoints_right $intr_left $intr_right $extr_left $extr_right $fund $motion $info_left $info_right > $paste_all
find $dir -name "points_*_left" | grep -o '[0-9]\+' |sort -u > $temp

while read p1<&3 && read p2<&4; do
  data_f=${out_dir}$p2.data
  echo $p1 | sed 's/#\s/\n/g' > $data_f
done 3<$paste_all 4<$temp

rm $keypoints_left
rm $keypoints_right
rm $intr_left
rm $intr_right
rm $extr_left
rm $extr_right
rm $fund
rm $motion
rm $info_left
rm $info_right
rm $paste_all

sed -i '$ d' $temp
sed -ri "s|$|.data|" $temp
sed -ri "s|^|--d $out_dir|" $temp

paste $temp| sed 's#^\(.*\)$#./bin/groebner_automatic_solver \1 --i 1000000#g' | parallel

rm $temp

