//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ICAMERA_H
#define CAMERA_CALIBRATION_ICAMERA_H

#include "Sophus/sophus/so3.hpp"
#include "Sophus/sophus/se3.hpp"
#include "Abstract_Estimator.h"

template<typename Intrinsics>
struct ICamera {

    virtual ~ICamera() = default;

    virtual void estimateIntrinsics(estimators::AbstractEstimator<Intrinsics> &estimator) = 0;

    virtual void estimatedExtrinsicsRotation(estimators::AbstractEstimator<Sophus::SO3d> &estimator) = 0;

    virtual void estimatedExtrinsicsTranslation(estimators::AbstractEstimator<Eigen::Vector3d> &estimator) = 0;

    //virtual void estimatedExtrinsicsMotion(estimators::AbstractEstimator<Sophus::SE3d> &estimator) = 0;

};

#endif //CAMERA_CALIBRATION_ICAMERA_H
