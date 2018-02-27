//
// Created by danielbord on 2/21/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_H
#define CAMERA_CALIBRATION_SCENE_H

#include "Two_View.h"
#include "../interfaces/IScene.h"

namespace scene {

    template<typename TCamera, typename TTwoView = TwoView<typename TCamera::Model>>
    class Scene : public IScene<Scene<TCamera, TTwoView>> {
        friend class IScene<Scene<TCamera, TTwoView>>;

        std::shared_ptr<typename TTwoView::VertexMap> ptr_to_map_;
        std::vector<TTwoView> list_of_stereo_pairs_;

    protected:
        void estimateCameraImpl(const typename TCamera::Label &label,
                                estimators::AbstractEstimator<Eigen::Vector3d> &estimator) {
            (*ptr_to_map_)[label].estimate(estimator);
        }


        void estimateCameraImpl(const typename TCamera::Label &label,
                                estimators::AbstractEstimator<Sophus::SO3d> &estimator) {
            (*ptr_to_map_)[label].estimate(estimator);
        }


        void estimateCameraImpl(const typename TCamera::Label &label,
                                estimators::AbstractEstimator<typename TCamera::Model> &estimator) {
            (*ptr_to_map_)[label].estimate(estimator);
        }


        void estimateStereoPairImpl(size_t k, estimators::AbstractEstimator<FundamentalMatrix> &estimator) {
            list_of_stereo_pairs_[k].estimateFundamentalMatrix(estimator);
        }


        void
        estimateStereoPairsImpl(
                estimators::AbstractEstimator<scene::StdVector<FundamentalMatrix>> &estimator) {
            auto result = estimator.getEstimation();
            assert(result.size() >= list_of_stereo_pairs_.size() &&
                   "Number of estimators should me no less than number of stereo pairs");
            for (size_t k = 0; k < list_of_stereo_pairs_.size(); ++k)
                list_of_stereo_pairs_[k].estimateFundamentalMatrix(result[k]);


        }

    public:

        using Camera_t = TCamera;

    };

}

#endif //CAMERA_CALIBRATION_SCENE_H
