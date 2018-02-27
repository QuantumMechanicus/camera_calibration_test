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

    template<typename TLabel, typename IntrinsicsModel>
    using MapLabelToCamera = std::map<TLabel, scene::Camera<IntrinsicsModel>>;

    template<typename IntrinsicsModel, typename TLabel = std::string>
    class TwoView
            : public graph::AbstractEdge<scene::Camera<IntrinsicsModel, TLabel>>, public ITwoView<
                    TwoView<IntrinsicsModel, TLabel>> {

        friend class ITwoView<TwoView<IntrinsicsModel, TLabel>>;

        Eigen::Vector3d relativeTranslation_{};
        Sophus::SO3d relativeRotation_{};
        ImagePoints left_keypoints_{};
        ImagePoints right_keypoints_{};
        FundamentalMatrix bifocal_tensor_{};

        long number_of_points_{};

    protected:

        void estimateLeftIntrinsicsImpl(estimators::AbstractEstimator<IntrinsicsModel> &estimator) {
            //if (doesExist())
            //st_vertex_.lock()->estimateIntrinsics(estimator);
            this->ptr_to_list_of_vertices_->at(this->start_vertex_label_).estimate(estimator);
        }


        void estimateRightIntrinsicsImpl(estimators::AbstractEstimator<IntrinsicsModel> &estimator) {
            //if (doesExist())
            //end_vertex_.lock()->estimateIntrinsics(estimator);
            this->ptr_to_list_of_vertices_->at(this->end_vertex_label_).estimate(estimator);
        }


        void estimateFundamentalMatrixImpl(estimators::AbstractEstimator<FundamentalMatrix> &estimator) {
            bifocal_tensor_ = estimator.getEstimation();
        }

        void estimateFundamentalMatrixImpl(const Eigen::Matrix3d &simple_estimation) {
            bifocal_tensor_ = simple_estimation;
        }


    public:


        using VertexMap = typename graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::VertexMap_t;

        TwoView() = default;

        TwoView(std::shared_ptr<VertexMap> cameras, TLabel left_camera_label,
                TLabel right_camera_label,
                ImagePoints left_keypoints,
                ImagePoints right_keypoints,
                FundamentalMatrix bifocal_tensor = FundamentalMatrix::Zero(),
                Sophus::SO3d relativeRotation = Sophus::SO3d(),
                Eigen::Vector3d relativeTranslation = Eigen::Vector3d::Zero())
                : graph::AbstractEdge<scene::Camera<IntrinsicsModel>>(
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
            //if (doesExist()) {
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
            //if (doesExist()) {
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
        void loadScene(const SceneArchiver &serializator) {
            serializator.deserialize(*this);
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

        const std::shared_ptr<IntrinsicsModel> getLeftIntrinsicsPointer() const {
            /*if (doesExist())
                return st_vertex_.lock()->getIntrinsicsPointer();*/
            return this->getStartVertex().getIntrinsicsPointer();
        }

        const std::shared_ptr<IntrinsicsModel> getRightIntrinsicsPointer() const {
            /*if (doesExist())
                return end_vertex_.lock()->getIntrinsicsPointer();*/
            return this->getFinishVertex().getIntrinsicsPointer();
        }

        const IntrinsicsModel &getLeftIntrinsics() const {
            /*if (doesExist())
                return st_vertex_.lock()->getIntrinsics();*/
            return this->getStartVertex().getIntrinsics();
        }

        const IntrinsicsModel &getRightIntrinsics() const {
            /*if (doesExist())
                return end_vertex_.lock()->getIntrinsics();*/
            return this->getFinishVertex().getIntrinsics();
        }

        const Sophus::SO3d &getRelativeRotation() const {
            /*if (doesExist())
                return st_vertex_.lock()->getRotation();*/
            return relativeRotation_;
        }


        const Eigen::Vector3d &getRelativeTranslation() const {
            /*if (doesExist())
                return st_vertex_.lock()->getRotation();*/
            return relativeTranslation_;
        }

        const Sophus::SO3d &getLeftAbsoluteRotation() const {
            /*if (doesExist())
                return st_vertex_.lock()->getRotation();*/
            return this->getStartVertex().getRotation();
        }

        const Sophus::SO3d &getRightAbsoluteRotation() const {
            /*if (doesExist())
                return end_vertex_.lock()->getRotation();*/
            return this->getFinishVertex().getRotation();
        }

        const Eigen::Vector3d &getLeftAbsoluteTranslation() const {
            /*if (doesExist())
                return st_vertex_.lock()->getTranslation();*/
            return this->getStartVertex().getTranslation();
        }

        const Eigen::Vector3d &getRightAbsoluteTranslation() const {
            /*if (doesExist())
                return end_vertex_.lock()->getTranslation();*/
            return this->getFinishVertex().getTranslation();
        }


    };

}
#endif //CAMERA_CALIBRATION_TWO_VIEW_H
