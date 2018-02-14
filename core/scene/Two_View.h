//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_TWO_VIEW_H
#define CAMERA_CALIBRATION_TWO_VIEW_H

#include <utility>
#include "Camera.h"

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

    typedef Eigen::Matrix<double, 3, 1> HomogenousImagePoint;

    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> ImagePoints;

    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> HomogenousImagePoints;

    typedef Eigen::Matrix3d FundamentalMatrix;

    template<typename IntrinsicsModel>
    class TwoView {
        internal_scene::Camera<IntrinsicsModel> left_camera_;
        internal_scene::Camera<IntrinsicsModel> right_camera_;

        ImagePoints left_keypoints_{};
        ImagePoints right_keypoints_{};
        FundamentalMatrix bifocal_tensor_{};

        long number_of_points_{};

    public:

        TwoView() = default;

        TwoView(internal_scene::Camera<IntrinsicsModel> left_camera,
                internal_scene::Camera<IntrinsicsModel> right_camera,
                ImagePoints left_keypoints,
                ImagePoints right_keypoints,
                FundamentalMatrix bifocal_tensor) : left_camera_(std::move(left_camera)),
                                                    right_camera_(std::move(right_camera)),
                                                    left_keypoints_(std::move(left_keypoints)),
                                                    right_keypoints_(std::move(right_keypoints)),
                                                    bifocal_tensor_(std::move(bifocal_tensor)) {
            number_of_points_ = TwoView::left_keypoints_.cols();
        }

        TwoView(internal_scene::Camera<IntrinsicsModel> left_camera,
                internal_scene::Camera<IntrinsicsModel> right_camera,
                ImagePoints left_keypoints,
                ImagePoints right_keypoints) : left_camera_(std::move(left_camera)),
                                               right_camera_(std::move(right_camera)),
                                               left_keypoints_(std::move(left_keypoints)),
                                               right_keypoints_(std::move(right_keypoints)) {
            bifocal_tensor_.setZero();
            number_of_points_ = TwoView::left_keypoints_.cols();
        }

        void normalizeLeftKeypoints() {
            double w = left_camera_.getWidth();
            double h = left_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            left_keypoints_.row(0) =
                    (left_keypoints_.row(0) - (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;
            left_keypoints_.row(1) =
                    (left_keypoints_.row(1) - (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;


        }

        void normalizeRightKeypoints() {
            double w = right_camera_.getWidth();
            double h = right_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            right_keypoints_.row(0) =
                    (right_keypoints_.row(0) - (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;
            right_keypoints_.row(1) =
                    (right_keypoints_.row(1) - (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose()) / r;

        }

        void denormalizeLeftKeypoints() {
            double w = left_camera_.getWidth();
            double h = left_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            left_keypoints_.row(0) =
                    r * left_keypoints_.row(0) + (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
            left_keypoints_.row(1) =
                    r * left_keypoints_.row(1) + (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();


        }

        void denormalizeRightKeypoints() {
            double w = right_camera_.getWidth();
            double h = right_camera_.getHeight();
            double r = std::sqrt(w * w + h * h) / 2.0;
            right_keypoints_.row(0) =
                    r * right_keypoints_.row(0) + (w / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();
            right_keypoints_.row(1) =
                    r * right_keypoints_.row(1) + (h / 2.0) * Eigen::VectorXd::Ones(number_of_points_).transpose();

        }


        template<typename IntrinsicsEstimator>
        void estimateLeftIntrinsics(IntrinsicsEstimator &estimator) {
            left_camera_.estimateIntrinsics(estimator);
        }

        template<typename IntrinsicsEstimator>
        void estimateRightIntrinsics(IntrinsicsEstimator &estimator) {
            right_camera_.estimateIntrinsics(estimator);
        }

        template<typename FundamentalMatrixEstimator>
        void estimateFundamentalMatrix(FundamentalMatrixEstimator &estimator) {
            if (!estimator.isEstimated())
                estimator.estimate();
            bifocal_tensor_ = estimator.getFundamentalMatrix();
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
            return left_camera_.getIntrinsicsPointer();
        }

        const std::shared_ptr<IntrinsicsModel> getRightIntrinsicsPointer() const {
            return right_camera_.getIntrinsicsPointer();
        }

        const IntrinsicsModel &getLeftIntrinsics() const {
            return left_camera_.getIntrinsics();
        }

        const IntrinsicsModel &getRightIntrinsics() const {
            return right_camera_.getIntrinsics();
        }

        const Sophus::SO3d &getLeftRotation() const {

            return left_camera_.getRotation();
        }

        const Sophus::SO3d &getRightRotation() const {

            return right_camera_.getRotation();
        }

        const Eigen::Vector3d &getLeftTranslation() const {

            return left_camera_.getTranslation();
        }

        const Eigen::Vector3d &getRightTranslation() const {

            return right_camera_.getTranslation();
        }


    };

}
#endif //CAMERA_CALIBRATION_TWO_VIEW_H
