//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_IO_H
#define CAMERA_CALIBRATION_SCENE_IO_H

#include <utility>

#include "../scene/Two_View.h"

namespace scene_serialization {

    template<typename IntrinsicsArchiver, typename IntrinsicsModel>
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
                                     std::string left_extrinsics_parameters_file_name = "",
                                     std::string right_extrinsics_parameters_file_name = "",
                                     std::string left_keypoints_file_name = "",
                                     std::string right_keypoints_file_name = ""
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

        void serialize(const scene::TwoView<IntrinsicsModel> &stereo_pair) const {


            utils::saveMatrix(fundamental_matrix_file_name_, stereo_pair.getFundamentalMatrix(), true);
            utils::saveMatrix(left_keypoints_file_name_, stereo_pair.getLeftKeypoints().transpose().eval(), true);
            utils::saveMatrix(right_keypoints_file_name_, stereo_pair.getRightKeypoints().transpose().eval(), true);

            Sophus::SE3d left_motion(stereo_pair.getLeftRotation(), stereo_pair.getLeftTranslation());
            Sophus::SE3d right_motion(stereo_pair.getRightRotation(), stereo_pair.getRightTranslation());
            utils::saveMatrix(left_extrinsics_parameters_file_name_, left_motion.matrix());
            utils::saveMatrix(right_extrinsics_parameters_file_name_, right_motion.matrix());

            IntrinsicsArchiver intrinsics_archiver(left_intrinsics_parameters_file_name_);
            intrinsics_archiver.serialize(stereo_pair.getLeftIntrinsics());
            if (left_intrinsics_parameters_file_name_ != right_intrinsics_parameters_file_name_) {
                intrinsics_archiver = IntrinsicsArchiver(right_intrinsics_parameters_file_name_);
                intrinsics_archiver.serialize(stereo_pair.getRightIntrinsics());
            }
        }

        void deserialize(scene::TwoView<IntrinsicsModel> &stereo_pair) {
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

            IntrinsicsModel left_intrinsics, right_intrinsics;


            IntrinsicsArchiver intrinsics_archiver(left_intrinsics_parameters_file_name_);
            intrinsics_archiver.deserialize(left_intrinsics);
            if (left_intrinsics_parameters_file_name_ != right_intrinsics_parameters_file_name_) {
                intrinsics_archiver = IntrinsicsArchiver(right_intrinsics_parameters_file_name_);
                intrinsics_archiver.deserialize(right_intrinsics);

            } else {
                right_intrinsics = left_intrinsics;
            }
            internal_scene::Camera<IntrinsicsModel> left_camera, right_camera;
            if (left_intrinsics == right_intrinsics) {
                left_camera = internal_scene::Camera<IntrinsicsModel>(
                        std::make_shared<IntrinsicsModel>(left_intrinsics), left_motion.rotationMatrix(),
                        left_motion.translation());
                right_camera = internal_scene::Camera<IntrinsicsModel>(left_camera.getIntrinsics(),
                                                                       right_motion.rotationMatrix(),
                                                                       right_motion.translation());


            } else {
                left_camera = internal_scene::Camera<IntrinsicsModel>(
                        std::make_shared<IntrinsicsModel>(left_intrinsics), left_motion.rotationMatrix(),
                        left_motion.translation());
                right_camera = internal_scene::Camera<IntrinsicsModel>(
                        std::make_shared<IntrinsicsModel>(right_intrinsics), right_motion.rotationMatrix(),
                        right_motion.translation());

            }
            stereo_pair = scene::TwoView<IntrinsicsModel>(left_camera, right_camera, left_points, right_points,
                                                          fundamental_matrix);
        }
    };

    template<int N = 1>
    class SimpleDivisionModelArchiver {
        std::string intrinsics_parameters_file_name_;


    public:
        explicit SimpleDivisionModelArchiver(std::string intrinsics_parameters_file_name = "")
                : intrinsics_parameters_file_name_(std::move(intrinsics_parameters_file_name)) {}

        void serialize(const intrinsics::DivisionModelIntrinsic<N> &intrinsics) const {
            Eigen::VectorXd vector_of_parameters(N + 5);
            vector_of_parameters.head(N) = intrinsics.getDistortionCoefficients();
            vector_of_parameters.tail(5)
                    << intrinsics.getFocalLength(), intrinsics.getPrincipalPointX(), intrinsics.getPrincipalPointY(), intrinsics.getHeight(), intrinsics.getWidth();
            utils::saveMatrix(intrinsics_parameters_file_name_, vector_of_parameters, true);

        }

        void deserialize(intrinsics::DivisionModelIntrinsic<N> &intrinsics) {
            Eigen::VectorXd vector_of_parameters(N + 5);
            utils::loadMatrix(intrinsics_parameters_file_name_, vector_of_parameters);
            Eigen::VectorXd distortion_coefficients = vector_of_parameters.head(N);
            double ppx, ppy, w, h, f;
            w = vector_of_parameters[N + 4];
            h = vector_of_parameters[N + 3];
            ppy = vector_of_parameters[N + 2];
            ppx = vector_of_parameters[N + 1];
            f = vector_of_parameters[N];
            intrinsics = intrinsics::DivisionModelIntrinsic<N>(distortion_coefficients, static_cast<unsigned int>(w),
                                                               static_cast<unsigned int>(h), f, ppx, ppy);

        }
    };
}


#endif //CAMERA_CALIBRATION_SCENE_IO_H