//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_TWO_VIEW_H
#define CAMERA_CALIBRATION_TWO_VIEW_H

#include <utility>
#include "Camera.h"
#include "../interfaces/ITwoView.h"
#include "../interfaces/AbstractEdge.h"

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

    template<typename TLabel, typename IntrinsicsModel>
    using MapLabelToCamera = std::map<TLabel, scene::Camera<IntrinsicsModel>>;


    typedef Eigen::Matrix<double, 2, 1> ImagePoint;

    typedef Eigen::Matrix<double, 3, 1> HomogenousImagePoint;

    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> ImagePoints;

    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> HomogenousImagePoints;

    typedef Eigen::Matrix3d FundamentalMatrix;

    template<typename IntrinsicsModel, typename TLabel = int>
    class TwoView
            : public graph::AbstractEdge<scene::Camera<IntrinsicsModel>, TLabel>, public ITwoView<IntrinsicsModel> {

        Eigen::Vector3d relativeTranslation_{};
        Sophus::SO3d relativeRotation_{};
        ImagePoints left_keypoints_{};
        ImagePoints right_keypoints_{};
        FundamentalMatrix bifocal_tensor_{};

        long number_of_points_{};

    public:

        //using graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::doesExist;
        using graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::start_vertex_label_;
        using graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::end_vertex_label_;
        using graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::ptr_to_list_of_vertices_;
        using graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::getStartVertex;
        using graph::AbstractEdge<scene::Camera<IntrinsicsModel>>::getFinishVertex;


        TwoView() = default;

        TwoView(std::shared_ptr<MapLabelToCamera<TLabel, IntrinsicsModel>> cameras, TLabel left_camera_label,
                TLabel right_camera_label,
                ImagePoints left_keypoints,
                ImagePoints right_keypoints,
                FundamentalMatrix bifocal_tensor) : graph::AbstractEdge<scene::Camera<IntrinsicsModel>>(
                std::move(left_camera_label),
                std::move(right_camera_label), cameras),
                                                    left_keypoints_(std::move(left_keypoints)),
                                                    right_keypoints_(std::move(right_keypoints)),
                                                    bifocal_tensor_(std::move(bifocal_tensor)) {
            number_of_points_ = TwoView::left_keypoints_.cols();
            relativeRotation_ = Sophus::SO3d();
            relativeTranslation_.setZero();
        }

        TwoView(std::shared_ptr<MapLabelToCamera<TLabel, IntrinsicsModel>> cameras,
                scene::Camera<IntrinsicsModel> left_camera,
                scene::Camera<IntrinsicsModel> right_camera,
                ImagePoints left_keypoints,
                ImagePoints right_keypoints) : graph::AbstractEdge<scene::Camera<IntrinsicsModel>>(
                left_camera.getLabel(), right_camera.getLabel(), cameras),
                                               left_keypoints_(std::move(left_keypoints)),
                                               right_keypoints_(std::move(right_keypoints)) {
            bifocal_tensor_.setZero();
            number_of_points_ = TwoView::left_keypoints_.cols();
        }

        void normalizeLeftKeypoints() {
            //if (doesExist()) {
            auto &left_camera_ = getStartVertex();//.lock();
            double w = left_camera_.getWidth();
            double h = left_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            left_keypoints_.row(0) =
                    (left_keypoints_.row(0) - (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;
            left_keypoints_.row(1) =
                    (left_keypoints_.row(1) - (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;


            //}
        }

        void normalizeRightKeypoints() {
            //if (doesExist()) {
            auto &right_camera_ = getFinishVertex();//.lock();
            double w = right_camera_.getWidth();
            double h = right_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            right_keypoints_.row(0) =
                    (right_keypoints_.row(0) - (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) /
                    r;
            right_keypoints_.row(1) =
                    (right_keypoints_.row(1) - (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) /
                    r;
            // }

        }

        void denormalizeLeftKeypoints() {
            //if (doesExist()) {
            auto &left_camera_ = getStartVertex();//.lock();
            double w = left_camera_.getWidth();
            double h = left_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            left_keypoints_.row(0) =
                    r * left_keypoints_.row(0) + (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
            left_keypoints_.row(1) =
                    r * left_keypoints_.row(1) + (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();

            //}
        }

        void denormalizeRightKeypoints() {
            //if (doesExist()) {
            auto &right_camera_ = getFinishVertex();//.lock();
            double w = right_camera_.getWidth();
            double h = right_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            right_keypoints_.row(0) =
                    r * right_keypoints_.row(0) + (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
            right_keypoints_.row(1) =
                    r * right_keypoints_.row(1) + (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
            //}
        }


        void estimateLeftIntrinsics(estimators::AbstractEstimator<IntrinsicsModel> &estimator) {
            //if (doesExist())
            //st_vertex_.lock()->estimateIntrinsics(estimator);
            ptr_to_list_of_vertices_->at(start_vertex_label_).estimateIntrinsics(estimator);
        }


        void estimateRightIntrinsics(estimators::AbstractEstimator<IntrinsicsModel> &estimator) {
            //if (doesExist())
            //end_vertex_.lock()->estimateIntrinsics(estimator);
            ptr_to_list_of_vertices_->at(end_vertex_label_).estimateIntrinsics(estimator);
        }


        void estimateFundamentalMatrix(estimators::AbstractEstimator<FundamentalMatrix> &estimator) override {
            bifocal_tensor_ = estimator.getEstimation();
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
            return ptr_to_list_of_vertices_->at(start_vertex_label_).getIntrinsicsPointer();
        }

        const std::shared_ptr<IntrinsicsModel> getRightIntrinsicsPointer() const {
            /*if (doesExist())
                return end_vertex_.lock()->getIntrinsicsPointer();*/
            return ptr_to_list_of_vertices_->at(end_vertex_label_).getIntrinsicsPointer();
        }

        const IntrinsicsModel &getLeftIntrinsics() const {
            /*if (doesExist())
                return st_vertex_.lock()->getIntrinsics();*/
            return ptr_to_list_of_vertices_->at(start_vertex_label_).getIntrinsics();
        }

        const IntrinsicsModel &getRightIntrinsics() const {
            /*if (doesExist())
                return end_vertex_.lock()->getIntrinsics();*/
            return ptr_to_list_of_vertices_->at(end_vertex_label_).getIntrinsics();
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
            return ptr_to_list_of_vertices_->at(start_vertex_label_).getRotation();
        }

        const Sophus::SO3d &getRightAbsoluteRotation() const {
            /*if (doesExist())
                return end_vertex_.lock()->getRotation();*/
            return ptr_to_list_of_vertices_->at(end_vertex_label_).getRotation();
        }

        const Eigen::Vector3d &getLeftAbsoluteTranslation() const {
            /*if (doesExist())
                return st_vertex_.lock()->getTranslation();*/
            return ptr_to_list_of_vertices_->at(start_vertex_label_).getTranslation();
        }

        const Eigen::Vector3d &getRightAbsoluteTranslation() const {
            /*if (doesExist())
                return end_vertex_.lock()->getTranslation();*/
            return ptr_to_list_of_vertices_->at(end_vertex_label_).getTranslation();
        }


    };

}
#endif //CAMERA_CALIBRATION_TWO_VIEW_H
