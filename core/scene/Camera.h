//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_CAMERA_H
#define CAMERA_CALIBRATION_CAMERA_H

#include <memory>
#include <utility>
#include "Sophus/sophus/so3.hpp"
#include "Sophus/sophus/se3.hpp"
#include "Camera_Intrinsics.h"
#include "../interfaces/ICamera.h"
#include "../interfaces/INode.h"

namespace scene {
    /**
     * @brief Class describing position (L = RW + t, where L denotes local coordinates system and W --- world coordinates system) and camera inner parameters
     * @tparam IntrinsicsModel Parametrization of camera intrinsics model (pinhole, fisheye, etc.)
     */
    template<typename IntrinsicsModel>
    class Camera : public ICamera<IntrinsicsModel>, public graph::INode<> {
        std::shared_ptr<IntrinsicsModel> intrinsics_;
        Sophus::SO3d world_rotation_;
        Eigen::Vector3d world_translation_;

    public:
        /**
         * @brief Constructor
         */
        Camera() = default;

        Camera(Camera &&rhs) noexcept = default;

        Camera(const Camera &rhs) = default;

        Camera &operator=(const Camera &rhs) = default;

        Camera &operator=(Camera &&rhs) noexcept = default;

        /**
         * @brief Constructor
         * @param intrinsics Pointer to intrinsic parameters
         * @param rotation  Rotation element R of transform from world coordinates to local camera coordinates
         * @param translation Translation element t of transform from world coordinates to local camera coordinates
         */
        Camera(const std::shared_ptr<IntrinsicsModel> intrinsics, Sophus::SO3d rotation,
               Eigen::Vector3d translation) : intrinsics_(std::move(intrinsics)), world_rotation_(std::move(rotation)),
                                              world_translation_(std::move(translation)) {}


        /**
         * @brief Constructor
         * @param intrinsics Intrinsic parameters
         * @param rotation  Rotation element R of transform from world coordinates to local camera coordinates
         * @param translation Translation element t of transform from world coordinates to local camera coordinates
         */
        Camera(const IntrinsicsModel intrinsics, Sophus::SO3d rotation,
               Eigen::Vector3d translation) : intrinsics_(std::make_shared<IntrinsicsModel>(intrinsics)), world_rotation_(std::move(rotation)),
                                              world_translation_(std::move(translation)) {}


        /**
         * @brief Another version of constructor
         */
        explicit Camera(const std::shared_ptr<IntrinsicsModel> intrinsics) : intrinsics_(std::move(intrinsics)),
                                                                       world_rotation_() {
            world_translation_.setZero();
        }

        /**
         * @brief Another version of constructor
         */
        explicit Camera(const std::shared_ptr<const IntrinsicsModel> intrinsics) : intrinsics_(std::move(intrinsics)),
                                                                             world_rotation_() {
            world_translation_.setZero();
        }

        /**
         * @brief Method to estimate intrinsics parameters of camera
         * @tparam IntrinsicsEstimator Class with 'estimate' method
         * @param estimator Instance of intrinsic estimator
         */


        void estimateIntrinsics(estimators::AbstractEstimator<IntrinsicsModel> &estimator) override
        {
            intrinsics_->estimateParameters(estimator);
        }

        void estimatedExtrinsicsRotation(estimators::AbstractEstimator<Sophus::SO3d> &estimator) override {
            world_rotation_ = estimator.getEstimation();
        }

        void estimatedExtrinsicsTranslation(estimators::AbstractEstimator<Eigen::Vector3d> &estimator) override {
            world_translation_ = estimator.getEstimation();
        }

        /*void estimatedExtrinsicsMotion(estimators::AbstractEstimator<Sophus::SE3d> &estimator) override {

        }*/

        /**
         * @brief Getter for intrinsic parameters
         * @return Pointer to intrinsics
         */
        const std::shared_ptr<IntrinsicsModel> getIntrinsicsPointer() const {
            return intrinsics_;
        }

        /**
         * @brief Getter for intrinsic parameters
         * @return Intrinsic object
         */
        const IntrinsicsModel &getIntrinsics() const {
            return *intrinsics_;
        }

        /**
         * @brief Getter for R
         * @return Sophus representation of SO3 (3D-rotation) group
         */
        const Sophus::SO3d &getRotation() const {
            return world_rotation_;
        }

        /**
         * @brief Getter for t
         * @return 3D vector
         */
        const Eigen::Vector3d &getTranslation() const {
            return world_translation_;
        }

        /**
         * @brief getter for height of image produced by camera
         * @return Height
         */
        unsigned int getHeight() const {
            return (*intrinsics_).getHeight();
        }

        /**
         * @brief getter for width of image produced by camera
         * @return Width
         */
        unsigned int getWidth() const {
            return (*intrinsics_).getWidth();
        }

    };

}
#endif //CAMERA_CALIBRATION_CAMERA_H
