//
// Created by danielbord on 2/28/18.
//

#ifndef CAMERA_CALIBRATION_FOCAL_ESTIMATOR_H
#define CAMERA_CALIBRATION_FOCAL_ESTIMATOR_H

#include "Core.h"

namespace estimators {
    class SimpleFocalEstimator : public estimators::AbstractEstimator<intrinsics::FocalLength> {

        unsigned int w_, h_;
        bool is_estimated_;
        intrinsics::FocalLength f_;
        scene::FundamentalMatrices fundamental_matrices_;

        Eigen::Vector2d makeLinearEquation1(const Eigen::Matrix3d &F);

        Eigen::Vector2d makeLinearEquation2(const Eigen::Matrix3d &F);

        Eigen::Vector3d makeQuadricEquation(const Eigen::Matrix3d &F);

        Eigen::Matrix<double, 8, 1> makeSquaredQuadricEquation(const Eigen::Matrix3d &F);

    protected:
        void estimateImpl() override;

        void getEstimationImpl(intrinsics::FocalLength &result) override;

    public:
        bool isEstimated() const override;

        explicit SimpleFocalEstimator(scene::FundamentalMatrices fundamental_matrices_, unsigned int w, unsigned int h);

    };
}
#endif //CAMERA_CALIBRATION_FOCAL_ESTIMATOR_H
