//
// Created by danielbord on 2/21/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_H
#define CAMERA_CALIBRATION_SCENE_H

#include "Two_View.h"
#include "../interfaces/IScene.h"

namespace scene {

    template<typename TCamera, typename TTwoView = TwoView<typename TCamera::Model>>
    class Scene : public IScene<TCamera, TTwoView> {
        std::shared_ptr<typename TTwoView::VertexMap> ptr_to_map_;
        std::vector<TTwoView> list_of_stereo_pairs_;

    public:

        void estimateExtrinsicsTranslation(const typename TCamera::Label &label,
                                           estimators::AbstractEstimator<Eigen::Vector3d> &estimator) override {
            (*ptr_to_map_)[label].estimateExtrinsicsTranslation(estimator);
        }


        void estimateExtrinsicsRotation(const typename TCamera::Label &label,
                                        estimators::AbstractEstimator<Sophus::SO3d> &estimator) override {
            (*ptr_to_map_)[label].estimateExtrinsicsRotation(estimator);
        }


        void estimateIntrinsics(const typename TCamera::Label &label,
                                estimators::AbstractEstimator<typename TCamera::Model> &estimator) override {
            (*ptr_to_map_)[label].estimateIntrinsics(estimator);
        }


        void estimateFundamentalMatrix(size_t k, estimators::AbstractEstimator<FundamentalMatrix> &estimator) override {
            list_of_stereo_pairs_[k].estimateFundamentalMatrix(estimator);
        }


        void
        estimateFundamentalMatrices(
                estimators::AbstractEstimator<scene::StdVector<FundamentalMatrix>> &estimator) override {
            auto result = estimator.getEstimation();
            assert(result.size() >= list_of_stereo_pairs_.size() &&
                   "Number of estimators should me no less than number of stereo pairs");
            for (size_t k = 0; k < list_of_stereo_pairs_.size(); ++k)
                list_of_stereo_pairs_[k].estimateFundamentalMatrix(result[k]);


        }
    };

}

#endif //CAMERA_CALIBRATION_SCENE_H
