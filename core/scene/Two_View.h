//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_TWO_VIEW_H
#define CAMERA_CALIBRATION_TWO_VIEW_H

#include <utility>
#include "Camera.h"
#include "../interfaces/ITwo_View.h"
#include "../interfaces/IEdge.h"

namespace scene {


    template<typename TIntrinsicsModel, typename TLabel = std::string>
    class TwoView
            : public graph::AbstractEdge<scene::Camera<TIntrinsicsModel, TLabel>>, public ITwoView<
                    TwoView<TIntrinsicsModel, TLabel>> {

        friend class ITwoView<TwoView<TIntrinsicsModel, TLabel>>;

        Eigen::Vector3d relativeTranslation_{};
        Sophus::SO3d relativeRotation_{};
        ImagePoints left_keypoints_{};
        ImagePoints right_keypoints_{};
        FundamentalMatrix bifocal_tensor_{};
        //TODO add recompute f_matrix
        long number_of_points_{};

    protected:

        template <typename TEstimator>
        void estimateLeftCameraImpl(TEstimator &estimator)
        {
            this->ptr_to_list_of_vertices_->at(this->start_vertex_label_).estimate(estimator);
        }

        template <typename TEstimator>
        void estimateRightCameraImpl(TEstimator &estimator) {
            this->ptr_to_list_of_vertices_->at(this->end_vertex_label_).estimate(estimator);
        }


        void estimateFundamentalMatrixImpl(estimators::AbstractEstimator<FundamentalMatrix> &estimator) {
            bifocal_tensor_ = estimator.getEstimation();
        }

        void estimateFundamentalMatrixImpl(const Eigen::Matrix3d &simple_estimation) {
            bifocal_tensor_ = simple_estimation;
        }


    public:


        using VertexMap_t = typename graph::AbstractEdge<scene::Camera<TIntrinsicsModel>>::VertexMap_t;

        TwoView() = default;

        TwoView(std::shared_ptr<VertexMap_t> cameras, TLabel left_camera_label,
                TLabel right_camera_label,
                ImagePoints left_keypoints,
                ImagePoints right_keypoints,
                FundamentalMatrix bifocal_tensor = FundamentalMatrix::Zero(),
                Sophus::SO3d relativeRotation = Sophus::SO3d(),
                Eigen::Vector3d relativeTranslation = Eigen::Vector3d::Zero())
                : graph::AbstractEdge<scene::Camera<TIntrinsicsModel>>(
                std::move(left_camera_label),
                std::move(right_camera_label), cameras),
                  left_keypoints_(std::move(left_keypoints)),
                  right_keypoints_(std::move(right_keypoints)),
                  bifocal_tensor_(std::move(bifocal_tensor)),
                  relativeRotation_(std::move(relativeRotation)),
                  relativeTranslation_(std::move(relativeTranslation)) {
            number_of_points_ = TwoView::left_keypoints_.cols();

        }


        bool normalizeLeftKeypoints() {

            auto &left_camera_ = this->getStartVertex();
            double w = left_camera_.getWidth();
            double h = left_camera_.getHeight();
            if (w > 0 && h > 0) {
                double r = std::sqrt(w * w + h * h) / 2.0;
                left_keypoints_.row(0) =
                        (left_keypoints_.row(0) - (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;
                left_keypoints_.row(1) =
                        (left_keypoints_.row(1) - (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;
                return true;
            }
            return false;

        }

        bool normalizeRightKeypoints() {

            auto &right_camera_ = this->getFinishVertex();
            double w = right_camera_.getWidth();
            double h = right_camera_.getHeight();
            if (w > 0 && h > 0) {
                double r = std::sqrt(w * w + h * h) / 2.0;
                right_keypoints_.row(0) =
                        (right_keypoints_.row(0) - (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) /
                        r;
                right_keypoints_.row(1) =
                        (right_keypoints_.row(1) - (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) /
                        r;
                return true;
            }
            return false;

        }

        bool denormalizeLeftKeypoints() {

            auto &left_camera_ = this->getStartVertex();
            double w = left_camera_.getWidth();
            double h = left_camera_.getHeight();
            if (w > 0 && h > 0) {
                double r = std::sqrt(w * w + h * h) / 2.0;
                left_keypoints_.row(0) =
                        r * left_keypoints_.row(0) + (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
                left_keypoints_.row(1) =
                        r * left_keypoints_.row(1) + (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
                return true;
            }
            return false;
        }

        bool denormalizeRightKeypoints() {

            auto &right_camera_ = this->getFinishVertex();
            double w = right_camera_.getWidth();
            double h = right_camera_.getHeight();
            if (w > 0 && h > 0) {
                double r = std::sqrt(w * w + h * h) / 2.0;
                right_keypoints_.row(0) =
                        r * right_keypoints_.row(0) + (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
                right_keypoints_.row(1) =
                        r * right_keypoints_.row(1) + (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
                return true;
            }
            return false;
        }


        template<typename SceneArchiver>
        void saveScene(const SceneArchiver &serializator) const {
            serializator.serialize(*this);
        }

        template<typename SceneArchiver>
        void loadScene(const SceneArchiver &serializator,
        std::shared_ptr<VertexMap_t> ptr_to_list_of_vertices) {
            serializator.deserialize(*this, ptr_to_list_of_vertices);
        }

        template<typename SceneArchiver>
        void loadScene(const SceneArchiver &serializator) {
            serializator.deserialize(*this, this->ptr_to_list_of_vertices_);
        }

        const ImagePoints &getLeftKeypoints() const {
            return left_keypoints_;
        }

        const ImagePoints &getRightKeypoints() const {
            return right_keypoints_;
        }

        const FundamentalMatrix &getFundamentalMatrix() const {
            return bifocal_tensor_;
        }

        const std::shared_ptr<TIntrinsicsModel> getLeftIntrinsicsPointer() const {
            return this->getStartVertex().getIntrinsicsPointer();
        }

        const std::shared_ptr<TIntrinsicsModel> getRightIntrinsicsPointer() const {
            return this->getFinishVertex().getIntrinsicsPointer();
        }

        const TIntrinsicsModel &getLeftIntrinsics() const {
            return this->getStartVertex().getIntrinsics();
        }

        const TIntrinsicsModel &getRightIntrinsics() const {
            return this->getFinishVertex().getIntrinsics();
        }

        const Sophus::SO3d &getRelativeRotation() const {
            return relativeRotation_;
        }


        const Eigen::Vector3d &getRelativeTranslation() const {
            return relativeTranslation_;
        }

        const Sophus::SO3d &getLeftAbsoluteRotation() const {
            return this->getStartVertex().getRotation();
        }

        const Sophus::SO3d &getRightAbsoluteRotation() const {
            return this->getFinishVertex().getRotation();
        }

        const Eigen::Vector3d &getLeftAbsoluteTranslation() const {
            return this->getStartVertex().getTranslation();
        }

        const Eigen::Vector3d &getRightAbsoluteTranslation() const {
            return this->getFinishVertex().getTranslation();
        }

    };

    using StandartDivisionModelStereoPair = TwoView<intrinsics::StandardDivisionModel>;
    using DynamicDivisionModelStereoPair = TwoView<intrinsics::DynamicDivisionModel>;
}
#endif //CAMERA_CALIBRATION_TWO_VIEW_H
