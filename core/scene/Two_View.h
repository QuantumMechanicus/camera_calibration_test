//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_TWO_VIEW_H
#define CAMERA_CALIBRATION_TWO_VIEW_H

#include "../interfaces/ITwo_View.h"
#include "../interfaces/IEdge.h"
#include "../utils/Utilities.h"

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
        Eigen::Matrix<double, 2, Eigen::Dynamic> i1d{}, i2d{};
        //TODO add recompute f_matrix
        long number_of_points_{};




    protected:

        template<typename TEstimator>
        void estimateLeftCameraImpl(TEstimator &estimator) {
            this->ptr_to_list_of_vertices_->at(this->start_vertex_label_).estimate(estimator);
        }

        template<typename TEstimator>
        void estimateRightCameraImpl(TEstimator &estimator) {
            this->ptr_to_list_of_vertices_->at(this->end_vertex_label_).estimate(estimator);
        }


        void estimateFundamentalMatrixImpl(estimators::AbstractEstimator<FundamentalMatrix> &estimator) {
            bifocal_tensor_ = estimator.getEstimation();
        }

        void estimateFundamentalMatrixImpl(const Eigen::Matrix3d &simple_estimation) {
            bifocal_tensor_ = simple_estimation;
        }

        struct costEssentialFunctor {
            Eigen::Matrix3d fm;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            costEssentialFunctor(const Eigen::Matrix3d &f) : fm(f) {}

            template<typename T>
            bool operator()(const T *f_ptr, const T *cx_ptr, const T *cy_ptr, T *residuals) const {
                Eigen::Matrix<T, 3, 3> K;
                K.setIdentity();
                K(0, 0) = K(1, 1) = *f_ptr;
                K(0, 2) = *cx_ptr;
                K(1, 2) = *cy_ptr;
                auto essential_matrix = (K.transpose() * fm.template cast<T>() * K).eval();

                Eigen::Map<Eigen::Matrix<T, 3, 3> > e(residuals);

                e = ((essential_matrix * essential_matrix.transpose() * essential_matrix
                      - 0.5 * (essential_matrix * essential_matrix.transpose()).trace() * essential_matrix) /
                     ceres::pow(essential_matrix.norm(), 3));
                return *f_ptr > T(0.0);

            }
        };

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

        const FundamentalMatrix getEssentialMatrix() {
            auto left_K = getLeftIntrinsicsPointer()->getCalibrationMatrix();
            auto right_K = getRightIntrinsicsPointer()->getCalibrationMatrix();


            /*double f, cx, cy;
            f = left_K(0, 0);
            cx = left_K(0, 2);
            cy = left_K(1, 2);
            ceres::Problem problem;
            problem.AddParameterBlock(&f, 1);
            problem.AddParameterBlock(&cx, 1);
            problem.AddParameterBlock(&cy, 1);

            problem.AddResidualBlock(new ceres::AutoDiffCostFunction<costEssentialFunctor, 9, 1, 1, 1>(
                    new costEssentialFunctor(bifocal_tensor_)),
                                     nullptr, &f, &cx, &cy);
            std::cout << "Intitial K: \n" << left_K << std::endl;
            ceres::Solver::Options options;
            options.max_num_iterations = 500;
            options.linear_solver_type = ceres::DENSE_QR;
            options.num_threads = 8;
            options.function_tolerance = 1e-16;
            options.parameter_tolerance = 1e-16;
            options.minimizer_progress_to_stdout = true;
            options.preconditioner_type = ceres::IDENTITY;
            options.jacobi_scaling = false;

            // Solve
            ceres::Solver::Summary summary;
            ceres::Solve(options, &problem, &summary);
            std::cout << summary.BriefReport() << std::endl;

            std::cout << "Final: " << f << " " << cx << " " << cy << std::endl;*/

            /*Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(bifocal_tensor_,
                                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;

            bifocal_tensor_= fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                     fmatrix_svd.matrixV().transpose();
            bifocal_tensor_ /= bifocal_tensor_(2, 2);
            std::cout << "K1: \n" << left_K << std::endl;
            std::cout << "K2: \n" << right_K << std::endl;*/
            return right_K.transpose() * bifocal_tensor_ * left_K;
        }

        void recoverRelativeMotion() {


            Eigen::Matrix3d matrix_D;

            matrix_D.setZero();

            Eigen::Matrix3d essential_matrix(getEssentialMatrix());

            //essential_matrix.normalize();

            matrix_D(0, 1) = 1;
            matrix_D(1, 0) = -1;
            matrix_D(2, 2) = 1;
            std::cout << (essential_matrix * essential_matrix.transpose() * essential_matrix
                          - 0.5 * (essential_matrix * essential_matrix.transpose()).trace() * essential_matrix).norm() /
                         std::pow(essential_matrix.norm(), 3) << std::endl;

            Eigen::JacobiSVD<Eigen::Matrix3d> svd(essential_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d matrixU = svd.matrixU();
            Eigen::Matrix3d matrixV = svd.matrixV();
            Eigen::Matrix3d singVal = svd.singularValues().asDiagonal();
            std::cout << "Sing values: " << singVal.transpose() << std::endl;
            singVal(0, 0) = singVal(1, 1) = 1;
            singVal(2, 2) = 0;
            Eigen::Matrix3d original_E = essential_matrix;
            essential_matrix = matrixU * singVal * matrixV.transpose();
            std::cout << (essential_matrix * essential_matrix.transpose() * essential_matrix
                          - 0.5 * (essential_matrix * essential_matrix.transpose()).trace() * essential_matrix).norm() /
                         std::pow(essential_matrix.norm(), 3) << std::endl;
            //std::cout << svd.singularValues().transpose() << std::endl;
            // std::cout << essential_matrix - (matrixU*svd.singularValues().asDiagonal()*matrixV.transpose()).normalized() << std::endl;
            if (matrixU.determinant() < 0)
                matrixU = -matrixU;
            if (matrixV.determinant() < 0)
                matrixV = -matrixV;

            Eigen::Matrix3d rotation_matrix = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d translation_matrix = Eigen::Matrix3d::Zero();
            Eigen::Vector3d translation_vector = matrixU.col(2);
            //translation_vector.normalize();
            translation_matrix = utils::screw_hat(translation_vector);

            Eigen::Matrix3d current_rotation = Eigen::Matrix3d::Identity();
            Eigen::Matrix3d current_translation = translation_matrix;

            int max_counter = 0;
            for (std::size_t k = 0; k < 4; ++k) {
                int counter = 0;
                switch (k) {
                    case 0:
                        current_rotation = matrixU * matrix_D * matrixV.transpose();
                        break;
                    case 1:
                        current_rotation = matrixU * matrix_D * matrixV.transpose();
                        current_translation.transposeInPlace();
                        break;
                    case 2:
                        current_rotation = matrixU * matrix_D.transpose() * matrixV.transpose();
                        current_translation.transposeInPlace();
                        break;
                    case 3:
                        current_rotation = matrixU * matrix_D.transpose() * matrixV.transpose();
                        current_translation.transposeInPlace();
                        break;
                    default:
                        break;
                }
                std::cout << std::min((essential_matrix.normalized() - (current_translation *
                                                                        current_rotation).normalized()).norm(),
                                      (essential_matrix.normalized() + (current_translation *
                                                                        current_rotation).normalized()).norm())
                          << " check" << std::endl;
                std::cout << std::min((original_E.normalized() - (current_translation *
                                                                  current_rotation).normalized()).norm(),
                                      (original_E.normalized() + (current_translation *
                                                                  current_rotation).normalized()).norm())
                          << " check_OE" << std::endl;

                std::vector<size_t> inliers_ind;
                Sophus::SE3d leftToRight(current_rotation, utils::inverted_screw_hat(current_translation));

                double interval = utils::distortion_problem::findInliers(left_keypoints_,
                                                                         right_keypoints_,
                                                                         getLeftIntrinsicsPointer()->getDistortionCoefficients(),
                                                                         leftToRight,
                                                                         getLeftIntrinsicsPointer()->getCalibrationMatrix(),
                                                                         0.1,
                                                                         inliers_ind, Eigen::Vector2d(getLeftIntrinsicsPointer()->getWidth(), getLeftIntrinsicsPointer()->getHeight()).norm() / 2.0);



                counter = inliers_ind.size();
                LOG(INFO) << "Number of points in front of camera " << counter << " Max: " << max_counter;

                if (counter > max_counter) {
                    rotation_matrix = current_rotation;
                    translation_matrix = current_translation;
                    max_counter = counter;
                }
            }
            relativeRotation_ = Sophus::SO3d(rotation_matrix);
            relativeTranslation_ = utils::inverted_screw_hat(translation_matrix).normalized();

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

        long getNumberOfPoints() const {
            return number_of_points_;
        }

    };

    using StandartDivisionModelStereoPair = TwoView<intrinsics::StandardDivisionModel>;
    using DynamicDivisionModelStereoPair = TwoView<intrinsics::DynamicDivisionModel>;
}
#endif //CAMERA_CALIBRATION_TWO_VIEW_H
