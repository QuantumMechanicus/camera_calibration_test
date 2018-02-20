//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_IO_H
#define CAMERA_CALIBRATION_SCENE_IO_H

#include <utility>
#include <regex>
#include <fstream>
#include "../scene/Two_View.h"

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

        /*const static std::regex rgx_left_keypoints;

        const static  std::regex rgx_right_keypoints;

        const static std::regex rgx_left_info;

        const static  std::regex rgx_right_info;

        const static  std::regex rgx_left_intr;

        const static  std::regex rgx_right_intr;

        const static std::regex rgx_left_extr;

        const static std::regex rgx_right_extr;

        const static  std::regex rgx_fundamental;

        const static  std::regex rgx_motion;*/

        SceneFiles() = default;

        explicit SceneFiles(std::string fundamental_matrix_file_name = "",
                            std::string left_intrinsics_parameters_file_name = "",
                            std::string right_intrinsics_parameters_file_name = "",
                            std::string left_keypoints_file_name = "",
                            std::string right_keypoints_file_name = "",
                            std::string left_camera_info_file_name = "",
                            std::string right_camera_info_file_name = "",
                            std::string relative_motion_file_name = "",
                            std::string left_extrinsics_parameters_file_name = "",
                            std::string right_extrinsics_parameters_file_name = ""
        ) : left_camera_info_file_name_(std::move(left_camera_info_file_name)),
            right_camera_info_file_name_(std::move(right_camera_info_file_name)),
            relative_motion_file_name_(std::move(relative_motion_file_name)),
            left_keypoints_file_name_(std::move(left_keypoints_file_name)),
            right_keypoints_file_name_(std::move(right_keypoints_file_name)),
            left_intrinsics_parameters_file_name_(std::move(left_intrinsics_parameters_file_name)),
            right_intrinsics_parameters_file_name_(std::move(right_intrinsics_parameters_file_name)),
            left_extrinsics_parameters_file_name_(std::move(left_extrinsics_parameters_file_name)),
            right_extrinsics_parameters_file_name_(std::move(right_extrinsics_parameters_file_name)),
            fundamental_matrix_file_name_(std::move(fundamental_matrix_file_name)) {}

        bool parse(const std::string &info_file) {
            std::ifstream info(info_file);
            if (info.good())
            {   std::string lines;
                std::string line;
                while (getline(info, line))
                {
                    lines = lines + lines + "\n";
                }



                handleRegex(lines, rgx_left_keypoints, left_keypoints_file_name_);
                handleRegex(lines, rgx_right_keypoints, right_keypoints_file_name_);

                handleRegex(lines, rgx_left_info, left_camera_info_file_name_);
                handleRegex(lines, rgx_right_info, right_camera_info_file_name_);

                handleRegex(lines, rgx_left_intr, left_intrinsics_parameters_file_name_);
                handleRegex(lines, rgx_right_intr, right_intrinsics_parameters_file_name_);

                handleRegex(lines, rgx_left_extr, left_extrinsics_parameters_file_name_);
                handleRegex(lines, rgx_right_extr, right_extrinsics_parameters_file_name_);

                handleRegex(lines, rgx_fundamental, fundamental_matrix_file_name_);
                handleRegex(lines, rgx_motion, relative_motion_file_name_);

                return true;
            }
            else
                return false;
        }

    private:
        bool handleRegex(const std::string &lines, const std::regex &rgx, std::string &result_string)
        {
            std::smatch result;
            bool found = std::regex_search(lines, result, rgx);
            if (found) {
                result_string = result[0].str();
                result_string = result_string.substr(result_string.find(':') + 1);
            }
            return found;
        }
    };

    template<typename CameraArchiver>
    class SimpleSceneArchiver {
        SceneFiles files_;


    public:
        SimpleSceneArchiver() = default;

        explicit SimpleSceneArchiver(std::string fundamental_matrix_file_name = "",
                                     std::string left_intrinsics_parameters_file_name = "",
                                     std::string right_intrinsics_parameters_file_name = "",
                                     std::string left_keypoints_file_name = "",
                                     std::string right_keypoints_file_name = "",
                                     std::string left_camera_info_file_name = "",
                                     std::string right_camera_info_file_name = "",
                                     std::string relative_motion_file_name = "",
                                     std::string left_extrinsics_parameters_file_name = "",
                                     std::string right_extrinsics_parameters_file_name = ""
        ) : files_(std::move(left_camera_info_file_name),
                   std::move(right_camera_info_file_name),
                   std::move(relative_motion_file_name),
                   std::move(right_keypoints_file_name),
                   std::move(left_intrinsics_parameters_file_name),
                   std::move(right_intrinsics_parameters_file_name),
                   std::move(left_extrinsics_parameters_file_name),
                   std::move(right_extrinsics_parameters_file_name),
                   std::move(fundamental_matrix_file_name)) {}

        bool parse(const std::string &info_file) {
            return files_.parse(info_file);
        }

        void serialize(const scene::TwoView<typename CameraArchiver::TModel> &stereo_pair) const {


            utils::saveMatrix(files_.fundamental_matrix_file_name_, stereo_pair.getFundamentalMatrix(), true);
            utils::saveMatrix(files_.left_keypoints_file_name_, stereo_pair.getLeftKeypoints().transpose().eval(),
                              true);
            utils::saveMatrix(files_.right_keypoints_file_name_, stereo_pair.getRightKeypoints().transpose().eval(),
                              true);

            const Sophus::SO3d &relative_rotation = stereo_pair.getRelativeRotation();
            const Eigen::Vector3d &relative_translation = stereo_pair.getRelativeTranslation();
            Sophus::SE3d relative_motion(relative_rotation, relative_translation);
            utils::saveMatrix(files_.relative_motion_file_name_, relative_motion.matrix(), true);

            CameraArchiver archiver;

            archiver = CameraArchiver(files_.left_camera_info_file_name_, files_.left_intrinsics_parameters_file_name_,
                                      files_.left_extrinsics_parameters_file_name_);
            archiver.serialize(stereo_pair.getStartVertex());
            archiver = CameraArchiver(files_.right_camera_info_file_name_,
                                      files_.right_intrinsics_parameters_file_name_,
                                      files_.right_extrinsics_parameters_file_name_);
            archiver.serialize(stereo_pair.getFinishVertex());
        }

        void deserialize(scene::TwoView<typename CameraArchiver::TCamera> &stereo_pair,
                         std::shared_ptr<std::map<typename CameraArchiver::TLabel,
                                 typename CameraArchiver::TCamera>> ptr_to_list_of_vertices) const {

            Eigen::Matrix3d fundamental_matrix;
            scene::ImagePoints left_points, right_points;


            utils::loadMatrix(files_.fundamental_matrix_file_name_, fundamental_matrix);
            utils::loadMatrix(files_.left_keypoints_file_name_, left_points, true);
            utils::loadMatrix(files_.right_keypoints_file_name_, right_points, true);


            typename CameraArchiver::TCamera left_camera, right_camera;

            CameraArchiver archiver;

            archiver = CameraArchiver(files_.left_camera_info_file_name_, files_.left_intrinsics_parameters_file_name_,
                                      files_.left_extrinsics_parameters_file_name_);

            archiver.deserialize(left_camera);
            archiver = CameraArchiver(files_.right_camera_info_file_name_,
                                      files_.right_intrinsics_parameters_file_name_,
                                      files_.right_extrinsics_parameters_file_name_);
            archiver.deserialize(right_camera);
            if (files_.left_intrinsics_parameters_file_name_ == files_.right_intrinsics_parameters_file_name_) {
                right_camera = CameraArchiver::TCamera(right_camera.getLabel(),
                                                       left_camera.getIntrinsicsPointer(),
                                                       right_camera.getRotation(), right_camera.getTranslation());


            }
            //TODO check if key left or right camera exsits
            (*ptr_to_list_of_vertices)[left_camera.getLabel()] = left_camera;
            (*ptr_to_list_of_vertices)[right_camera.getLabel()] = right_camera;

            stereo_pair = scene::TwoView<typename CameraArchiver::TCamera>(ptr_to_list_of_vertices,
                                                                           left_camera.getLabel(),
                                                                           right_camera.getLabel(), left_points,
                                                                           right_points,
                                                                           fundamental_matrix);
        }
    };

    template<typename TInfo = std::string, int N = 1>
    class SimpleDivisionModelArchiver {
        std::string intrinsics_parameters_file_name_;
        std::string absolute_motion_file_name_;
        std::string camera_info_file_name_;
    public:
        typedef scene::Camera<intrinsics::DivisionModelIntrinsic<N>, TInfo> TCamera;
        typedef intrinsics::DivisionModelIntrinsic<N> TModel;
        typedef TInfo TLabel;


        SimpleDivisionModelArchiver() {
            intrinsics_parameters_file_name_ = "";
            absolute_motion_file_name_ = "";
            camera_info_file_name_ = "";
        };

        explicit SimpleDivisionModelArchiver(std::string camera_info_file_name,
                                             std::string intrinsics_parameters_file_name,
                                             std::string absolute_motion_file_name = "")
                : camera_info_file_name_(std::move(camera_info_file_name)),
                  intrinsics_parameters_file_name_(std::move(intrinsics_parameters_file_name)),
                  absolute_motion_file_name_(std::move(absolute_motion_file_name)) {}


        void serialize(const scene::Camera<intrinsics::DivisionModelIntrinsic<N>, TInfo> &camera) const {
            auto intrinsics_ptr = camera.getIntrinsicsPointer();
            Eigen::VectorXd vector_of_parameters(intrinsics_ptr->getNumberOfCofficients() + 5);
            vector_of_parameters.head(
                    intrinsics_ptr->getNumberOfCofficients()) = intrinsics_ptr->getDistortionCoefficients();
            vector_of_parameters.tail(5)
                    << intrinsics_ptr->getFocalLength(), intrinsics_ptr->getPrincipalPointX(), intrinsics_ptr->getPrincipalPointY(), camera.getHeight(), camera.getWidth();
            utils::saveMatrix(intrinsics_parameters_file_name_, vector_of_parameters, true);
            std::fstream info(camera_info_file_name_, std::fstream::out);
            info << camera.getLabel();

        }

        void deserialize(scene::Camera<intrinsics::DivisionModelIntrinsic<N>, TInfo> &camera) {
            Eigen::VectorXd vector_of_parameters;
            utils::loadMatrix(intrinsics_parameters_file_name_, vector_of_parameters);
            Eigen::VectorXd distortion_coefficients = vector_of_parameters.head(vector_of_parameters.size() - 5);
            double ppx, ppy, w, h, f;
            w = vector_of_parameters[vector_of_parameters.size() - 1];
            h = vector_of_parameters[vector_of_parameters.size() - 2];
            ppy = vector_of_parameters[vector_of_parameters.size() - 3];
            ppx = vector_of_parameters[vector_of_parameters.size() - 4];
            f = vector_of_parameters[vector_of_parameters.size() - 5];

            auto intrinsics = intrinsics::DivisionModelIntrinsic<N>(distortion_coefficients,
                                                                    static_cast<unsigned int>(w),
                                                                    static_cast<unsigned int>(h), f, ppx, ppy);
            std::fstream info(camera_info_file_name_, std::fstream::in);
            TInfo label;
            info >> label;
            Sophus::SE3d motion;
            Eigen::Matrix<double, 4, 4> motion_matrix;
            motion_matrix.setZero();
            utils::loadMatrix(absolute_motion_file_name_, motion_matrix);
            if (!motion_matrix.isZero())
                motion = Sophus::SE3d::fitToSE3(motion_matrix);
            camera = scene::Camera<intrinsics::DivisionModelIntrinsic<N>, TInfo>(label, intrinsics,
                                                                                 motion.rotationMatrix(),
                                                                                 motion.translation());

        }
    };
}


#endif //CAMERA_CALIBRATION_SCENE_IO_H
