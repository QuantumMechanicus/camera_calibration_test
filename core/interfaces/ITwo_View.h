//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ITWO_VIEW_H
#define CAMERA_CALIBRATION_ITWO_VIEW_H

#include "Eigen/Dense"
#include "ICamera.h"
#include "Abstract_Estimator.h"

namespace scene {
    template<typename TDerived>
    struct ITwoView {
    public:
        //TODO forward
        template <typename TEstimator>
        void estimateFundamentalMatrix(TEstimator &estimator)
        {
            static_cast<TDerived*>(this)->estimateFundamentalMatrixImpl(estimator);
        }

        template <typename TEstimator>
        void estimateLeftIntrinsics(TEstimator &estimator) {
            static_cast<TDerived*>(this)->estimateLeftIntrinsicsImpl(estimator);
        }

        template <typename TEstimator>
        void estimateRightIntrinsics(TEstimator &estimator) {
            static_cast<TDerived*>(this)->estimateRightIntrinsicsImpl(estimator);
        }

    };
}


#endif //CAMERA_CALIBRATION_ITWO_VIEW_H
