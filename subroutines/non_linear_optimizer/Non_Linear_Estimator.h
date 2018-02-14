//
// Created by danielbord on 2/14/18.
//

#ifndef CAMERA_CALIBRATION_NON_LINEAR_ESTIMATOR_H
#define CAMERA_CALIBRATION_NON_LINEAR_ESTIMATOR_H

#include "Core.h"

namespace non_linear_optimization {
    class NonLinearEstimator
            : public estimators::internal::DivisionModelIntrinsicsEstimator<Eigen::Dynamic> {
    public:
        void estimate() override;

        bool isEstimated() const override;
    };
}

#endif //CAMERA_CALIBRATION_NON_LINEAR_ESTIMATOR_H
