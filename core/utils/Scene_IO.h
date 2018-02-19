//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_IO_H
#define CAMERA_CALIBRATION_SCENE_IO_H

#include <utility>

#include "../scene/Two_View.h"

namespace scene_serialization {

    static const int CompileTimeKnown = -1;
    template<typename IntrinsicsArchiver>
    class SimpleSceneArchiver {
        std::string left_keypoints_file_name_;
        std::string right_keypoints_file_name_;
        std::string left_intrinsics_parameters_file_name_;
        std::string right_intrinsics_parameters_file_name_;
        std::string left_extrinsics_parameters_file_name_;
        std::string right_extrinsics_parameters_file_name_;
        std::string fundamental_matrix_file_name_;

    public:
        explicit SimpleSceneArchiver(std::string fundamental_matrix_file_name = "",
                                     std::string left_intrinsics_parameters_file_name = "",
                                     std::string right_intrinsics_parameters_file_name = "",
                                     std::string left_keypoints_file_name = "",
                                     std::string right_keypoints_file_name = "",
                                     std::string left_extrinsics_parameters_file_name = "",
                                     std::string right_extrinsics_parameters_file_name = ""
        ) : left_keypoints_file_name_(std::move(
                left_keypoints_file_name)), right_keypoints_file_name_(std::move(right_keypoints_file_name)),
            left_intrinsics_parameters_file_name_(
                    std::move(
                            left_intrinsics_parameters_file_name)),
            right_intrinsics_parameters_file_name_(
                    std::move(
                            right_intrinsics_parameters_file_name)),
            left_extrinsics_parameters_file_name_(
                    std::move(
                            left_extrinsics_parameters_file_name)),
            right_extrinsics_parameters_file_name_(
                    std::move(
                            right_extrinsics_parameters_file_name)),
            fundamental_matrix_file_name_(std::move(
                    fundamental_matrix_file_name)) {}


        void serialize(const scene::TwoView<typename IntrinsicsArchiver::Model> &stereo_pair) const {


            utils::saveMatrix(fundamental_matrix_file_name_, stereo_pair.getFundamentalMatrix(), true);
            utils::saveMatrix(left_keypoints_file_name_, stereo_pair.getLeftKeypoints().transpose().eval(), true);
            utils::saveMatrix(right_keypoints_file_name_, stereo_pair.getRightKeypoints().transpose().eval(), true);

            Sophus::SE3d left_motion(stereo_pair.getLeftRotation(), stereo_pair.getLeftTranslation());
            Sophus::SE3d right_motion(stereo_pair.getRightRotation(), stereo_pair.getRightTranslation());
            utils::saveMatrix(left_extrinsics_parameters_file_name_, left_motion.matrix());
            utils::saveMatrix(right_extrinsics_parameters_file_name_, right_motion.matrix());
            IntrinsicsArchiver intrinsics_archiver;

            intrinsics_archiver = IntrinsicsArchiver(left_intrinsics_parameters_file_name_);

            intrinsics_archiver.serialize(stereo_pair.getLeftIntrinsics());
            if (left_intrinsics_parameters_file_name_ != right_intrinsics_parameters_file_name_) {
                intrinsics_archiver = IntrinsicsArchiver(right_intrinsics_parameters_file_name_);
                intrinsics_archiver.serialize(stereo_pair.getRightIntrinsics());
            }
        }

        void deserialize(scene::TwoView<typename IntrinsicsArchiver::Model> &stereo_pair) const {
            Eigen::Matrix3d fundamental_matrix;
            Eigen::Matrix<double, 4, 4> left_motion_matrix, right_motion_matrix;
            scene::ImagePoints left_points, right_points;
            Sophus::SE3d left_motion, right_motion;


            utils::loadMatrix(fundamental_matrix_file_name_, fundamental_matrix);
            utils::loadMatrix(left_keypoints_file_name_, left_points, true);
            utils::loadMatrix(right_keypoints_file_name_, right_points, true);

            utils::loadMatrix(left_extrinsics_parameters_file_name_, left_motion_matrix);
            utils::loadMatrix(right_extrinsics_parameters_file_name_, right_motion_matrix);

            if (!left_motion_matrix.isZero())
                left_motion = Sophus::SE3d(Sophus::SE3d::fitToSE3(left_motion_matrix));
            if (!right_motion_matrix.isZero())
                right_motion = Sophus::SE3d(Sophus::SE3d::fitToSE3(right_motion_matrix));

            typename IntrinsicsArchiver::Model left_intrinsics, right_intrinsics;

            IntrinsicsArchiver intrinsics_archiver;

            intrinsics_archiver = IntrinsicsArchiver(left_intrinsics_parameters_file_name_);

            intrinsics_archiver.deserialize(left_intrinsics);
            if (left_intrinsics_parameters_file_name_ != right_intrinsics_parameters_file_name_) {
                intrinsics_archiver = IntrinsicsArchiver(right_intrinsics_parameters_file_name_);
                intrinsics_archiver.deserialize(right_intrinsics);

            } else {
                right_intrinsics = left_intrinsics;
            }
            scene::Camera<typename IntrinsicsArchiver::Model> left_camera, right_camera;
            if (left_intrinsics == right_intrinsics) {
                left_camera = scene::Camera<typename IntrinsicsArchiver::Model>(
                        std::make_shared<typename IntrinsicsArchiver::Model>(left_intrinsics),
                        left_motion.rotationMatrix(),
                        left_motion.translation());
                right_camera = scene::Camera<typename IntrinsicsArchiver::Model>(left_camera.getIntrinsicsPointer(),
                                                                                          right_motion.rotationMatrix(),
                                                                                          right_motion.translation());


            } else {
                left_camera = scene::Camera<typename IntrinsicsArchiver::Model>(
                        std::make_shared<typename IntrinsicsArchiver::Model>(left_intrinsics),
                        left_motion.rotationMatrix(),
                        left_motion.translation());
                right_camera = scene::Camera<typename IntrinsicsArchiver::Model>(
                        std::make_shared<typename IntrinsicsArchiver::Model>(right_intrinsics),
                        right_motion.rotationMatrix(),
                        right_motion.translation());

            }
            stereo_pair = scene::TwoView<typename IntrinsicsArchiver::Model>(left_camera, right_camera, left_points,
                                                                             right_points,
                                                                             fundamental_matrix);
        }
    };

    template<int N = 1>
    class SimpleDivisionModelArchiver {
        std::string intrinsics_parameters_file_name_;


    public:
        typedef intrinsics::DivisionModelIntrinsic<N> Model;

        SimpleDivisionModelArchiver() {
            intrinsics_parameters_file_name_ = "";
        };

        explicit SimpleDivisionModelArchiver(std::string intrinsics_parameters_file_name)
                : intrinsics_parameters_file_name_(std::move(intrinsics_parameters_file_name)) {}



        void serialize(const intrinsics::DivisionModelIntrinsic<N> &intrinsics) const {

            Eigen::VectorXd vector_of_parameters(intrinsics.getNumberOfCofficients() + 5);
            vector_of_parameters.head(intrinsics.getNumberOfCofficients()) = intrinsics.getDistortionCoefficients();
            vector_of_parameters.tail(5)
                    << intrinsics.getFocalLength(), intrinsics.getPrincipalPointX(), intrinsics.getPrincipalPointY(), intrinsics.getHeight(), intrinsics.getWidth();
            utils::saveMatrix(intrinsics_parameters_file_name_, vector_of_parameters, true);

        }

        void deserialize(intrinsics::DivisionModelIntrinsic<N> &intrinsics) {
            Eigen::VectorXd vector_of_parameters;
            utils::loadMatrix(intrinsics_parameters_file_name_, vector_of_parameters);
            Eigen::VectorXd distortion_coefficients = vector_of_parameters.head(vector_of_parameters.size()-5);
            double ppx, ppy, w, h, f;
            w = vector_of_parameters[vector_of_parameters.size()-1];
            h = vector_of_parameters[vector_of_parameters.size()-2];
            ppy = vector_of_parameters[vector_of_parameters.size()-3];
            ppx = vector_of_parameters[vector_of_parameters.size()-4];
            f = vector_of_parameters[vector_of_parameters.size()-5];
            intrinsics = intrinsics::DivisionModelIntrinsic<N>(distortion_coefficients,
                                                                   static_cast<unsigned int>(w),
                                                                   static_cast<unsigned int>(h), f, ppx, ppy);
        }
    };
}


#endif //CAMERA_CALIBRATION_SCENE_IO_H
