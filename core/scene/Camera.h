//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_CAMERA_H
#define CAMERA_CALIBRATION_CAMERA_H

#include <memory>
#include <utility>
#include "Sophus/sophus/so3.hpp"
#include "Sophus/sophus/se3.hpp"
#include "Intrinsics.h"
#include "../interfaces/ICamera.h"
#include "../interfaces/INode.h"

namespace scene {
    /**
     * @brief Class describing position (L = RW + t, where L denotes local coordinates system and W --- world coordinates system) and camera inner parameters
     * @tparam TIntrinsicsModel Parametrization of camera intrinsics model (pinhole, fisheye, etc.)
     */
    template<typename TIntrinsicsModel, typename TLabel = std::string>
    class Camera
            : public ICamera<Camera<TIntrinsicsModel, TLabel>>,
              public graph::INode<Camera<TIntrinsicsModel, TLabel>, TLabel> {

        friend class ICamera<Camera<TIntrinsicsModel, TLabel>>;
        friend graph::INode<Camera<TIntrinsicsModel, TLabel>, TLabel>;

        std::shared_ptr<TIntrinsicsModel> intrinsics_;
        Sophus::SO3d world_rotation_;
        Eigen::Vector3d world_translation_;
        TLabel label_;

    protected:
        TLabel getLabelImpl() const {
            return label_;
        }

        template<typename TEstimator>
        void estimateImpl(TEstimator &estimator) {
            intrinsics_->estimateParameter(estimator);
        }

        void estimateImpl(estimators::AbstractEstimator<Sophus::SO3d> &estimator) {
            world_rotation_ = estimator.getEstimation();
        }

        void estimateImpl(std::shared_ptr<TIntrinsicsModel> simple_estimator) {
            intrinsics_ = simple_estimator;
        }

        //TODO more simple estimators == setters


    public:

        using Model_t = TIntrinsicsModel;



        Camera() : world_rotation_{}, world_translation_{}, label_{} {
            intrinsics_ = std::make_shared<TIntrinsicsModel>(TIntrinsicsModel());
        }

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
        Camera(TLabel label, std::shared_ptr<TIntrinsicsModel> intrinsics, Sophus::SO3d rotation,
               Eigen::Vector3d translation) : label_(std::move(label)), intrinsics_(std::move(intrinsics)),
                                              world_rotation_(std::move(rotation)),
                                              world_translation_(std::move(translation)) {}


        /**
         * @brief Constructor
         * @param intrinsics Intrinsic parameters
         * @param rotation  Rotation element R of transform from world coordinates to local camera coordinates
         * @param translation Translation element t of transform from world coordinates to local camera coordinates
         */
        Camera(TLabel label, TIntrinsicsModel intrinsics, Sophus::SO3d rotation,
               Eigen::Vector3d translation) : label_(std::move(label)),
                                              intrinsics_(std::make_shared<TIntrinsicsModel>(intrinsics)),
                                              world_rotation_(std::move(rotation)),
                                              world_translation_(std::move(translation)) {}


        /**
         * @brief Another version of constructor
         */
        Camera(TLabel label, std::shared_ptr<TIntrinsicsModel> intrinsics) : label_(std::move(label)),
                                                                             intrinsics_(std::move(intrinsics)),
                                                                             world_rotation_() {
            world_translation_.setZero();
        }

        /**
         * @brief Another version of constructor
         */
        Camera(TLabel label, std::shared_ptr<const TIntrinsicsModel> intrinsics) : label_(std::move(label)),
                                                                                   intrinsics_(std::move(intrinsics)),
                                                                                   world_rotation_() {
            world_translation_.setZero();
        }


        /**
         * @brief Getter for intrinsic parameters
         * @return Pointer to intrinsics
         */
        const std::shared_ptr<TIntrinsicsModel> getIntrinsicsPointer() const {
            return intrinsics_;
        }

        /**
         * @brief Getter for intrinsic parameters
         * @return Intrinsic object
         */
        const TIntrinsicsModel &getIntrinsics() const {
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
         * @brief Getter for camera local coordinates
         * @return element of SE3d
         */
        Sophus::SE3d getMotion() const {

            return Sophus::SE3d(world_rotation_, world_translation_);
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

    using DynamicDivisionModelCamera = Camera<intrinsics::DynamicDivisionModel >;
    using StandartDivisionModelCamera = Camera<intrinsics::StandardDivisionModel>;
    using PinholeCamera = Camera<intrinsics::PinholeModel >;

}
#endif //CAMERA_CALIBRATION_CAMERA_H
