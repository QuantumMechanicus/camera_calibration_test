//
// Created by danielbord on 2/21/18.
//

#ifndef CAMERA_CALIBRATION_ISCENE_H
#define CAMERA_CALIBRATION_ISCENE_H


#include "ITwo_View.h"
#include <memory>

namespace scene {
    template<typename TCamera, typename TTwoView>
    struct IScene {

        virtual void estimateExtrinsicsTranslation(const typename TCamera::Label &label,
                                                   estimators::AbstractEstimator<Eigen::Vector3d> &estimator) = 0;


        virtual void estimateExtrinsicsRotation(const typename TCamera::Label &label,
                                                estimators::AbstractEstimator<Sophus::SO3d> &estimator) = 0;


        virtual void estimateIntrinsics(const typename TCamera::Label &label,
                                        estimators::AbstractEstimator<typename TCamera::Model> &estimator) = 0;

        virtual void
        estimateFundamentalMatrix(size_t k, estimators::AbstractEstimator<FundamentalMatrix> &estimator) = 0;

        virtual void
        estimateFundamentalMatrices(estimators::AbstractEstimator<scene::StdVector<FundamentalMatrix>> &estimator) = 0;
    };
}
#endif //CAMERA_CALIBRATION_ISCENE_H
