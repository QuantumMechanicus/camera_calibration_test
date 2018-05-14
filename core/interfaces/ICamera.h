//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ICAMERA_H
#define CAMERA_CALIBRATION_ICAMERA_H

#include "sophus/so3.hpp"
#include "sophus/se3.hpp"
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

    typedef Eigen::Matrix<double, 3, 1> WorldPoint;

    typedef Eigen::Matrix<double, 4, 1> HomogenousWorldPoint;

    typedef Eigen::Matrix<double, 3, 1> HomogenousImagePoint;

    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> ImagePoints;

    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> HomogenousImagePoints;

    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> WorldPoints;

    typedef Eigen::Matrix<double, 4, Eigen::Dynamic> HomogenousWorldPoints;

    typedef Eigen::Matrix3d FundamentalMatrix;

    typedef StdVector<FundamentalMatrix> FundamentalMatrices;

    /**
     * @brief Base class for camera
     * @tparam TDerived --- CRTP
     */
    template<typename TDerived>
    struct ICamera {

        /**
         * @brief Interface to estimates camera parameters depends on type of TEstimator (intrinsics or its part, world coordinates, etc), for implemented types (see Camera.h)
         */
        template<typename TEstimator>
        void estimate(TEstimator &estimator) {
            static_cast<TDerived *>(this)->estimateImpl(estimator);
        }
    };
}


#endif //CAMERA_CALIBRATION_ICAMERA_H
