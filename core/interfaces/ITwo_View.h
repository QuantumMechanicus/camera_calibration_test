//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ITWO_VIEW_H
#define CAMERA_CALIBRATION_ITWO_VIEW_H

#include "Eigen/Dense"
#include "ICamera.h"
#include "Abstract_Estimator.h"

namespace scene {
    /**
     * @brief Base class for stereo pair
     * @tparam TDerived --- CRTP
     */
    template<typename TDerived>
    struct ITwoView {
    public:
        //TODO forward
        /**
         * @brief Estimate fundamental matrix via estimator
         */
        template <typename TEstimator>
        void estimateFundamentalMatrix(TEstimator &estimator)
        {
            static_cast<TDerived*>(this)->estimateFundamentalMatrixImpl(estimator);
        }

        /**
         * @brief Estimate any parameter of left camera (e.g. intrinsics, extrinsic), see implementation for details
         */
        template <typename TEstimator>
        void estimateLeftCamera(TEstimator &estimator) {
            static_cast<TDerived*>(this)->estimateLeftCameraImpl(estimator);
        }

        /**
         * @brief Estimate any parameter of left camera (e.g. intrinsics, extrinsic), see implementation for details
         */
        template <typename TEstimator>
        void estimateRightCamera(TEstimator &estimator) {
            static_cast<TDerived*>(this)->estimateRightCameraImpl(estimator);
        }

    };

}


#endif //CAMERA_CALIBRATION_ITWO_VIEW_H
