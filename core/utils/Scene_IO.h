//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_IO_H
#define CAMERA_CALIBRATION_SCENE_IO_H

#include <utility>
#include <regex>
#include <fstream>
#include <random>
#include <boost/algorithm/string.hpp>
#include "../scene/Two_View.h"
#include "Utilities.h"

namespace scene_serialization {

    const static int CompileTimeKnown = -1;

    struct SceneFiles {
        std::string left_keypoints_file_name_;
        std::string right_keypoints_file_name_;
        std::string relative_motion_file_name_;
        std::string left_intrinsics_parameters_file_name_;
        std::string right_intrinsics_parameters_file_name_;
        std::string left_extrinsics_parameters_file_name_;
        std::string right_extrinsics_parameters_file_name_;
        std::string fundamental_matrix_file_name_;
        std::string left_camera_info_file_name_;
        std::string right_camera_info_file_name_;

        const static std::regex rgx_left_keypoints;

        const static std::regex rgx_right_keypoints;

        const static std::regex rgx_left_info;

        const static std::regex rgx_right_info;

        const static std::regex rgx_left_intr;

        const static std::regex rgx_right_intr;

        const static std::regex rgx_left_extr;

        const static std::regex rgx_right_extr;

        const static std::regex rgx_fundamental;

        const static std::regex rgx_motion;

        SceneFiles();

        explicit SceneFiles(std::string fundamental_matrix_file_name = "",
                            std::string left_intrinsics_parameters_file_name = "",
                            std::string right_intrinsics_parameters_file_name = "",
                            std::string left_keypoints_file_name = "",
                            std::string right_keypoints_file_name = "",
                            std::string left_camera_info_file_name = "",
                            std::string right_camera_info_file_name = "",
                            std::string relative_motion_file_name = "",
                            std::string left_extrinsics_parameters_file_name = "",
                            std::string right_extrinsics_parameters_file_name = "");

        bool parse(const std::string &info_file);

    private:
        bool handleRegex(const std::string &lines, const std::regex &rgx, std::string &result_string);
    };

    template<typename CameraArchiver>
    class SimpleSceneArchiver {
        SceneFiles files_;
        std::vector<bool> overwrites_;

    public:

        explicit SimpleSceneArchiver(
                std::vector<bool> overwrites = {true, true, true, false, false, false, false, true, true, true},
                std::string fundamental_matrix_file_name = "",
                std::string left_intrinsics_parameters_file_name = "",
                std::string right_intrinsics_parameters_file_name = "",
                std::string left_keypoints_file_name = "",
                std::string right_keypoints_file_name = "",
                std::string left_camera_info_file_name = "",
                std::string right_camera_info_file_name = "",
                std::string relative_motion_file_name = "",
                std::string left_extrinsics_parameters_file_name = "",
                std::string right_extrinsics_parameters_file_name = "") :
                files_(std::move(fundamental_matrix_file_name),
                       std::move(left_intrinsics_parameters_file_name),
                       std::move(right_intrinsics_parameters_file_name),
                       std::move(left_keypoints_file_name),
                       std::move(right_keypoints_file_name),
                       std::move(left_camera_info_file_name),
                       std::move(right_camera_info_file_name),
                       std::move(relative_motion_file_name),
                       std::move(left_extrinsics_parameters_file_name),
                       std::move(right_extrinsics_parameters_file_name)), overwrites_(std::move(overwrites)) {}

        bool parse(const std::string &info_file) {
            return files_.parse(info_file);
        }

        void serialize(const scene::TwoView<typename CameraArchiver::Model_t> &stereo_pair) const {


            utils::saveMatrix(files_.fundamental_matrix_file_name_, stereo_pair.getFundamentalMatrix(), overwrites_[0]);
            utils::saveMatrix(files_.left_keypoints_file_name_, stereo_pair.getLeftKeypoints().transpose().eval(),
                              overwrites_[3]);
            utils::saveMatrix(files_.right_keypoints_file_name_, stereo_pair.getRightKeypoints().transpose().eval(),
                              overwrites_[4]);

            const Sophus::SO3d &relative_rotation = stereo_pair.getRelativeRotation();
            const Eigen::Vector3d &relative_translation = stereo_pair.getRelativeTranslation();
            Sophus::SE3d relative_motion(relative_rotation, relative_translation);
            utils::saveMatrix(files_.relative_motion_file_name_, relative_motion.matrix(), overwrites_[7]);

            CameraArchiver archiver;

            archiver = CameraArchiver(files_.left_camera_info_file_name_, files_.left_intrinsics_parameters_file_name_,
                                      files_.left_extrinsics_parameters_file_name_, overwrites_[1], overwrites_[8],
                                      overwrites_[5]);
            archiver.serialize(stereo_pair.getStartVertex());
            archiver = CameraArchiver(files_.right_camera_info_file_name_,
                                      files_.right_intrinsics_parameters_file_name_,
                                      files_.right_extrinsics_parameters_file_name_, overwrites_[2], overwrites_[9],
                                      overwrites_[6]);
            archiver.serialize(stereo_pair.getFinishVertex());
        }

        void deserialize(scene::TwoView<typename CameraArchiver::Model_t> &stereo_pair,
                         std::shared_ptr<typename scene::TwoView<typename CameraArchiver::Model_t>::VertexMap_t>
                         ptr_to_list_of_vertices) const {
            std::cout << "A" << std::endl;
            Eigen::Matrix3d fundamental_matrix;
            scene::ImagePoints left_points, right_points;
            Eigen::Matrix4d relative_motion_matrix;
            Sophus::SE3d relative_motion;
            utils::loadMatrix(files_.fundamental_matrix_file_name_, fundamental_matrix);
            utils::loadMatrix(files_.left_keypoints_file_name_, left_points, true);
            utils::loadMatrix(files_.right_keypoints_file_name_, right_points, true);
            utils::loadMatrix(files_.relative_motion_file_name_, relative_motion_matrix);

            if (!relative_motion_matrix.isZero())
                relative_motion = Sophus::SE3d::fitToSE3(relative_motion_matrix);

            typename CameraArchiver::Camera_t left_camera, right_camera;

            CameraArchiver archiver;

            archiver = CameraArchiver(files_.left_camera_info_file_name_, files_.left_intrinsics_parameters_file_name_,
                                      files_.left_extrinsics_parameters_file_name_);

            archiver.deserialize(left_camera);
            archiver = CameraArchiver(files_.right_camera_info_file_name_,
                                      files_.right_intrinsics_parameters_file_name_,
                                      files_.right_extrinsics_parameters_file_name_);
            archiver.deserialize(right_camera);
            if (files_.left_intrinsics_parameters_file_name_ == files_.right_intrinsics_parameters_file_name_) {
                right_camera = typename CameraArchiver::Camera_t(right_camera.getLabel(),
                                                                 left_camera.getIntrinsicsPointer(),
                                                                 right_camera.getRotation(),
                                                                 right_camera.getTranslation());


            }
            //TODO check if key left or right camera exsits
            (*ptr_to_list_of_vertices)[left_camera.getLabel()] = left_camera;
            (*ptr_to_list_of_vertices)[right_camera.getLabel()] = right_camera;

            stereo_pair = scene::TwoView<typename CameraArchiver::Model_t, typename CameraArchiver::Label_t>(
                    ptr_to_list_of_vertices,
                    left_camera.getLabel(),
                    right_camera.getLabel(), left_points,
                    right_points,
                    fundamental_matrix, relative_motion.rotationMatrix(), relative_motion.translation());
        }
    };

    template<int N = 1, typename TInfo = std::string>
    class SimpleDivisionModelArchiver {
        std::string intrinsics_parameters_file_name_;
        std::string absolute_motion_file_name_;
        std::string camera_info_file_name_;
        bool overwrite_intrinsics_;
        bool overwrite_extrinsics_;
        bool overwrite_info_;

    public:
        using Camera_t = scene::Camera<intrinsics::DivisionModel<N>, TInfo>;
        using Model_t =  intrinsics::DivisionModel<N>;
        using Label_t = TInfo;


        explicit SimpleDivisionModelArchiver(std::string camera_info_file_name = "",
                                             std::string intrinsics_parameters_file_name = "",
                                             std::string absolute_motion_file_name = "",
                                             bool overwrite_intrinsics = true,
                                             bool overwrite_extrinsics = true, bool overwrite_info = false)
                : overwrite_extrinsics_(overwrite_extrinsics),
                  overwrite_intrinsics_(overwrite_intrinsics),
                  overwrite_info_(overwrite_info),
                  camera_info_file_name_(std::move(camera_info_file_name)),
                  intrinsics_parameters_file_name_(std::move(intrinsics_parameters_file_name)),
                  absolute_motion_file_name_(std::move(absolute_motion_file_name)) {}


        void serialize(const scene::Camera<intrinsics::DivisionModel<N>, TInfo> &camera) const {
            auto intrinsics_ptr = camera.getIntrinsicsPointer();
            Eigen::VectorXd vector_of_parameters(intrinsics_ptr->getNumberOfCoefficients() + 5);
            vector_of_parameters.head(
                    intrinsics_ptr->getNumberOfCoefficients()) = intrinsics_ptr->getDistortionCoefficients();
            vector_of_parameters.tail(5)
                    << intrinsics_ptr->getFocalLength(), intrinsics_ptr->getPrincipalPointX(), intrinsics_ptr->getPrincipalPointY(), camera.getHeight(), camera.getWidth();
            utils::saveMatrix(intrinsics_parameters_file_name_, vector_of_parameters, true);
            std::fstream info(camera_info_file_name_, std::fstream::out);
            info << camera.getLabel();
            utils::saveMatrix(absolute_motion_file_name_, camera.getMotion().matrix(), true);

        }

        void deserialize(scene::Camera<intrinsics::DivisionModel<N>, TInfo> &camera) {
            intrinsics::DivisionModel<N> intrinsics;
            Eigen::VectorXd vector_of_parameters;
            utils::loadMatrix(intrinsics_parameters_file_name_, vector_of_parameters);
            if (vector_of_parameters.size() >= 5) {
                Eigen::VectorXd distortion_coefficients_dynamic = vector_of_parameters.head(
                        vector_of_parameters.size() - 5);
                double ppx, ppy, w, h, f;
                w = vector_of_parameters[vector_of_parameters.size() - 1];
                h = vector_of_parameters[vector_of_parameters.size() - 2];
                ppy = vector_of_parameters[vector_of_parameters.size() - 3];
                ppx = vector_of_parameters[vector_of_parameters.size() - 4];
                f = vector_of_parameters[vector_of_parameters.size() - 5];
                if (N != Eigen::Dynamic) {
                    intrinsics = intrinsics::DivisionModel<N>(distortion_coefficients_dynamic.head(N),
                                                                       static_cast<unsigned int>(w),
                                                                       static_cast<unsigned int>(h), f, ppx, ppy);
                } else
                    intrinsics = intrinsics::DivisionModel<N>(distortion_coefficients_dynamic,
                                                                       static_cast<unsigned int>(w),
                                                                       static_cast<unsigned int>(h), f, ppx, ppy);

            }
            std::fstream info(camera_info_file_name_, std::fstream::in);
            TInfo label;
            info >> label;
            if (label == "") {
                std::random_device rd;
                std::mt19937 gen(rd());
                label = std::to_string(gen());
            }
            Sophus::SE3d motion;
            Eigen::Matrix<double, 4, 4> motion_matrix;
            motion_matrix.setZero();
            utils::loadMatrix(absolute_motion_file_name_, motion_matrix);
            if (!motion_matrix.isZero())
                motion = Sophus::SE3d::fitToSE3(motion_matrix);
            camera = scene::Camera<intrinsics::DivisionModel<N>, TInfo>(label, intrinsics,
                                                                                 motion.rotationMatrix(),
                                                                                 motion.translation());

        }
    };

    using StandartSceneArchiver = SimpleSceneArchiver<SimpleDivisionModelArchiver<>>;
    using DynamicSceneArchiver = SimpleSceneArchiver<SimpleDivisionModelArchiver<Eigen::Dynamic>>;

}
#endif //CAMERA_CALIBRATION_SCENE_IO_H
