//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_IO_H
#define CAMERA_CALIBRATION_SCENE_IO_H

#include <utility>

#include "../scene/Two_View.h"

namespace scene_serialization {

    template<typename IntrinsicsArchiver>
    class SimpleSceneArchiver {
        std::string left_keypoints_file_name_;
        std::string right_keypoints_file_name_;
        std::string left_intrinsics_parameters_file_name_;
        std::string right_intrinsics_parameters_file_name_;
        std::string left_extrinsics_parameters_file_name_;
        std::string right_extrinsics_parameters_file_name_;
        std::string fundamental_matrix_file_name_;
        int number_of_intrinsic_parameters_;

    public:
        const int CompileTimeKnown = -1;

        explicit SimpleSceneArchiver(std::string fundamental_matrix_file_name = "",
                                     std::string left_intrinsics_parameters_file_name = "",
                                     std::string right_intrinsics_parameters_file_name = "",
                                     std::string left_extrinsics_parameters_file_name = "",
                                     std::string right_extrinsics_parameters_file_name = "",
                                     std::string left_keypoints_file_name = "",
                                     std::string right_keypoints_file_name = ""
        ) : number_of_intrinsic_parameters_(CompileTimeKnown), left_keypoints_file_name_(std::move(
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

        explicit SimpleSceneArchiver(int number_of_intrinsic_parameters,
                                     std::string fundamental_matrix_file_name = "",
                                     std::string left_intrinsics_parameters_file_name = "",
                                     std::string right_intrinsics_parameters_file_name = "",
                                     std::string left_extrinsics_parameters_file_name = "",
                                     std::string right_extrinsics_parameters_file_name = "",
                                     std::string left_keypoints_file_name = "",
                                     std::string right_keypoints_file_name = ""
        ) : number_of_intrinsic_parameters_(number_of_intrinsic_parameters), left_keypoints_file_name_(std::move(
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
            if (number_of_intrinsic_parameters_ != CompileTimeKnown)
                intrinsics_archiver = IntrinsicsArchiver (number_of_intrinsic_parameters_,
                                                          left_intrinsics_parameters_file_name_);
            else
                intrinsics_archiver = IntrinsicsArchiver(left_intrinsics_parameters_file_name_);

            intrinsics_archiver.serialize(stereo_pair.getLeftIntrinsics());
            if (left_intrinsics_parameters_file_name_ != right_intrinsics_parameters_file_name_) {
                intrinsics_archiver = IntrinsicsArchiver(right_intrinsics_parameters_file_name_);
                intrinsics_archiver.serialize(stereo_pair.getRightIntrinsics());
            }
        }

        void deserialize(scene::TwoView<typename IntrinsicsArchiver::Model> &stereo_pair) {
            Eigen::Matrix3d fundamental_matrix;
            Eigen::Matrix<double, 4, 3> left_motion_matrix, right_motion_matrix;
            scene::ImagePoints left_points, right_points;
            Sophus::SE3d left_motion, right_motion;


            utils::loadMatrix(fundamental_matrix_file_name_, fundamental_matrix);
            utils::loadMatrix(left_keypoints_file_name_, left_points, true);
            utils::loadMatrix(right_keypoints_file_name_, right_points, true);

            utils::loadMatrix(left_extrinsics_parameters_file_name_, left_motion_matrix);
            utils::loadMatrix(right_extrinsics_parameters_file_name_, right_motion_matrix);

            //TODO check it's ok
            left_motion = Sophus::SE3d(left_motion_matrix);
            right_motion = Sophus::SE3d(right_motion);

            typename IntrinsicsArchiver::Model left_intrinsics, right_intrinsics;

            IntrinsicsArchiver intrinsics_archiver;
            if (number_of_intrinsic_parameters_ != CompileTimeKnown)
                intrinsics_archiver = IntrinsicsArchiver(number_of_intrinsic_parameters_,
                                                         left_intrinsics_parameters_file_name_);
            else
                intrinsics_archiver = IntrinsicsArchiver(left_intrinsics_parameters_file_name_);

            intrinsics_archiver.deserialize(left_intrinsics);
            if (left_intrinsics_parameters_file_name_ != right_intrinsics_parameters_file_name_) {
                intrinsics_archiver = IntrinsicsArchiver(right_intrinsics_parameters_file_name_);
                intrinsics_archiver.deserialize(right_intrinsics);

            } else {
                right_intrinsics = left_intrinsics;
            }
            internal_scene::Camera<typename IntrinsicsArchiver::Model> left_camera, right_camera;
            if (left_intrinsics == right_intrinsics) {
                left_camera = internal_scene::Camera<typename IntrinsicsArchiver::Model>(
                        std::make_shared<typename IntrinsicsArchiver::Model>(left_intrinsics),
                        left_motion.rotationMatrix(),
                        left_motion.translation());
                right_camera = internal_scene::Camera<typename IntrinsicsArchiver::Model>(left_camera.getIntrinsics(),
                                                                                          right_motion.rotationMatrix(),
                                                                                          right_motion.translation());


            } else {
                left_camera = internal_scene::Camera<typename IntrinsicsArchiver::Model>(
                        std::make_shared<typename IntrinsicsArchiver::Model>(left_intrinsics),
                        left_motion.rotationMatrix(),
                        left_motion.translation());
                right_camera = internal_scene::Camera<typename IntrinsicsArchiver::Model>(
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
        unsigned int number_of_parameters;

    public:
        typedef intrinsics::DivisionModelIntrinsic<N> Model;

        SimpleDivisionModelArchiver() {
            intrinsics_parameters_file_name_ = "";
            number_of_parameters = 0;
        };

        explicit SimpleDivisionModelArchiver(std::string intrinsics_parameters_file_name)
                : intrinsics_parameters_file_name_(std::move(intrinsics_parameters_file_name)) {
            assert(N != Eigen::Dynamic && "Use should use dynamic constructor");
            number_of_parameters = static_cast<unsigned int>(N);
        }

        explicit SimpleDivisionModelArchiver(unsigned int n, std::string intrinsics_parameters_file_name)
                : intrinsics_parameters_file_name_(std::move(intrinsics_parameters_file_name)) {
            number_of_parameters = n;
        }

        void serialize(const intrinsics::DivisionModelIntrinsic<N> &intrinsics) const {
            Eigen::VectorXd vector_of_parameters(number_of_parameters + 5);
            vector_of_parameters.head(number_of_parameters) = intrinsics.getDistortionCoefficients();
            vector_of_parameters.tail(5)
                    << intrinsics.getFocalLength(), intrinsics.getPrincipalPointX(), intrinsics.getPrincipalPointY(), intrinsics.getHeight(), intrinsics.getWidth();
            utils::saveMatrix(intrinsics_parameters_file_name_, vector_of_parameters, true);

        }

        void deserialize(intrinsics::DivisionModelIntrinsic<N> &intrinsics) {
            Eigen::VectorXd vector_of_parameters(number_of_parameters + 5);
            utils::loadMatrix(intrinsics_parameters_file_name_, vector_of_parameters);
            Eigen::VectorXd distortion_coefficients = vector_of_parameters.head(number_of_parameters);
            double ppx, ppy, w, h, f;
            w = vector_of_parameters[number_of_parameters + 4];
            h = vector_of_parameters[number_of_parameters + 3];
            ppy = vector_of_parameters[number_of_parameters + 2];
            ppx = vector_of_parameters[number_of_parameters + 1];
            f = vector_of_parameters[number_of_parameters];
            if (N != Eigen::Dynamic)
                intrinsics = intrinsics::DivisionModelIntrinsic<N>(distortion_coefficients,
                                                                   static_cast<unsigned int>(w),
                                                                   static_cast<unsigned int>(h), f, ppx, ppy);
            else
                intrinsics = intrinsics::DivisionModelIntrinsic<N>(number_of_parameters, distortion_coefficients,
                                                                   static_cast<unsigned int>(w),
                                                                   static_cast<unsigned int>(h), f, ppx, ppy);

        }
    };
}


#endif //CAMERA_CALIBRATION_SCENE_IO_H
