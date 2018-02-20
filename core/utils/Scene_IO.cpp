//
// Created by danielbord on 1/26/18.
//
#include "Scene_IO.h"

namespace scene_serialization {
    const std::regex SceneFiles::rgx_left_keypoints = std::regex("(Left keypoints:|left keypoints:)(.*)");

    const std::regex SceneFiles::rgx_right_keypoints = std::regex("(Right keypoints:|right keypoints:)(.*)");

    const std::regex SceneFiles::rgx_left_info = std::regex("(Left camera info:|left camera info:)(.*)");

    const std::regex SceneFiles::rgx_right_info = std::regex("(Right camera info:|right camera info:)(.*)");

    const std::regex SceneFiles::rgx_left_intr = std::regex("(Left intrinsics:|left intrinsics:)(.*)");

    const std::regex SceneFiles::rgx_right_intr = std::regex("(Right intrinsics:|right intrinsics:)(.*)");

    const std::regex SceneFiles::rgx_left_extr = std::regex("(Left extrinsics:|left extrinsics:)(.*)");

    const std::regex SceneFiles::rgx_right_extr = std::regex("(Right extrinsics:|Right extrinsics:)(.*)");

    const std::regex SceneFiles::rgx_fundamental = std::regex("(Fundamental matrix:|fundamental matrix:)(.*)");

    const std::regex SceneFiles::rgx_motion = std::regex("(Relative motion:|relative motion:)(.*)");

    SceneFiles::SceneFiles() = default;

    SceneFiles::SceneFiles(std::string fundamental_matrix_file_name,
                           std::string left_intrinsics_parameters_file_name,
                           std::string right_intrinsics_parameters_file_name,
                           std::string left_keypoints_file_name,
                           std::string right_keypoints_file_name,
                           std::string left_camera_info_file_name,
                           std::string right_camera_info_file_name,
                           std::string relative_motion_file_name,
                           std::string left_extrinsics_parameters_file_name,
                           std::string right_extrinsics_parameters_file_name)
            : left_camera_info_file_name_(std::move(left_camera_info_file_name)),
              right_camera_info_file_name_(std::move(right_camera_info_file_name)),
              relative_motion_file_name_(std::move(relative_motion_file_name)),
              left_keypoints_file_name_(std::move(left_keypoints_file_name)),
              right_keypoints_file_name_(std::move(right_keypoints_file_name)),
              left_intrinsics_parameters_file_name_(std::move(left_intrinsics_parameters_file_name)),
              right_intrinsics_parameters_file_name_(std::move(right_intrinsics_parameters_file_name)),
              left_extrinsics_parameters_file_name_(std::move(left_extrinsics_parameters_file_name)),
              right_extrinsics_parameters_file_name_(std::move(right_extrinsics_parameters_file_name)),
              fundamental_matrix_file_name_(std::move(fundamental_matrix_file_name)) {}

    bool SceneFiles::parse(const std::string &info_file) {
        std::ifstream info(info_file);
        if (info.good()) {
            std::string lines;
            std::string line;
            while (getline(info, line)) {
                lines.append(line + "\n");
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
        } else
            return false;
    }

    bool SceneFiles::handleRegex(const std::string &lines, const std::regex &rgx, std::string &result_string) {
        std::smatch result;
        bool found = std::regex_search(lines, result, rgx);
        if (found) {
            result_string = result[0].str();
            result_string = result_string.substr(result_string.find(':') + 1);
            boost::trim(result_string);
        }
        return found;
    }


}
