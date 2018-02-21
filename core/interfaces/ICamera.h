//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ICAMERA_H
#define CAMERA_CALIBRATION_ICAMERA_H

#include "Sophus/sophus/so3.hpp"
#include "Sophus/sophus/se3.hpp"
#include "Abstract_Estimator.h"

namespace scene {

    template<typename T>
    using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;

    template<typename T>
    using TImagePoint = Eigen::Matrix<T, 2, 1>;

    template<typename T>
    using TImagePoints = Eigen::Matrix<T, 2, Eigen::Dynamic>;

    template<typename T>
    using THomogenousImagePoint = Eigen::Matrix<T, 3, 1>;

    template<typename T>
    using THomogenousImagePoints = Eigen::Matrix<T, 3, Eigen::Dynamic>;

    template<typename T>
    using TFundamentalMatrix = Eigen::Matrix<T, 3, 3>;

    typedef Eigen::Matrix<double, 2, 1> ImagePoint;

    typedef Eigen::Matrix<double, 3, 1> HomogenousImagePoint;

    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> ImagePoints;

    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> HomogenousImagePoints;

    typedef Eigen::Matrix3d FundamentalMatrix;

    template<typename Intrinsics>
    struct ICamera {

        virtual ~ICamera() = default;

        virtual void estimateIntrinsics(estimators::AbstractEstimator<Intrinsics> &estimator) = 0;

        virtual void estimateExtrinsicsRotation(estimators::AbstractEstimator<Sophus::SO3d> &estimator) = 0;

        virtual void estimateExtrinsicsTranslation(estimators::AbstractEstimator<Eigen::Vector3d> &estimator) = 0;

        virtual void estimateIntrinsics(const Intrinsics &simple_estimation) = 0;

        virtual void estimateExtrinsicsRotation(const Sophus::SO3d &simple_estimation) = 0;

        virtual void estimateExtrinsicsTranslation(const Eigen::Vector3d &simple_estimation) = 0;

        //virtual void estimatedExtrinsicsMotion(estimators::AbstractEstimator<Sophus::SE3d> &estimator) = 0;

    };
}


#endif //CAMERA_CALIBRATION_ICAMERA_H
