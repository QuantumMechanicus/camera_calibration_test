//
// Created by danielbord on 2/21/18.
//

#ifndef CAMERA_CALIBRATION_SCENE_H
#define CAMERA_CALIBRATION_SCENE_H

#include "Two_View.h"
#include "../interfaces/IScene.h"

namespace scene {

    template<typename TCamera, typename TTwoView = TwoView<typename TCamera::Model_t>>
    class Scene : public IScene<Scene<TCamera, TTwoView>> {
        friend class IScene<Scene<TCamera, TTwoView>>;

        std::shared_ptr<typename TTwoView::VertexMap_t> ptr_to_map_;
        std::vector<TTwoView> list_of_stereo_pairs_;
        //TODO approve TwoViews ptr_to_map and scene ptr_to_map

    protected:
        void estimateCameraImpl(const typename TCamera::Label_t &label,
                                estimators::AbstractEstimator<Eigen::Vector3d> &estimator) {
            (*ptr_to_map_)[label].estimate(estimator);
        }


        void estimateCameraImpl(const typename TCamera::Label_t &label,
                                estimators::AbstractEstimator<Sophus::SO3d> &estimator) {
            (*ptr_to_map_)[label].estimate(estimator);
        }


        void estimateCameraImpl(const typename TCamera::Label_t &label,
                                estimators::AbstractEstimator<typename TCamera::Model_t> &estimator) {
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

        Scene() = default;

        Scene(std::shared_ptr<typename TTwoView::VertexMap_t> ptr_to_map, std::vector<TTwoView> list_of_stereo_pairs)
                : ptr_to_map_(std::move(ptr_to_map)), list_of_stereo_pairs_(std::move(list_of_stereo_pairs)) {}

        template<typename SceneArchiver>
        void saveScene(const std::vector<SceneArchiver> &serializators) const {
            assert(serializators.size() >= list_of_stereo_pairs_.size() &&
                   "Number of serializators should me no less than number of stereo pairs");
            for (size_t k = 0; k < list_of_stereo_pairs_.size(); ++k)
                serializators[k].serialize(list_of_stereo_pairs_[k]);
        }

        template<typename SceneArchiver>
        void loadScene(const std::vector<SceneArchiver> &serializators,
                       std::shared_ptr<typename TTwoView::VertexMap_t> ptr_to_map) {
            ptr_to_map_ = ptr_to_map_;
            list_of_stereo_pairs_.resize(serializators.size());
            for (size_t k = 0; k < list_of_stereo_pairs_.size(); ++k)
                serializators[k].deserialize(list_of_stereo_pairs_[k], ptr_to_map);
        }

        template<typename SceneArchiver>
        void loadScene(const std::vector<SceneArchiver> &serializators) {
            list_of_stereo_pairs_.resize(serializators.size());
            for (size_t k = 0; k < list_of_stereo_pairs_.size(); ++k)
                serializators[k].deserialize(list_of_stereo_pairs_[k], ptr_to_map_);
        }


    };

    using StandartDivisionModelScene = Scene<StandartDivisionModelCamera>;
    using DynamicDivisionModelScene = Scene<DynamicDivisionModelCamera>;
}

#endif //CAMERA_CALIBRATION_SCENE_H
