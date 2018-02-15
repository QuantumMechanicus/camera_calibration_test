//
// Created by danielbord on 2/14/18.
//

#include "Non_Linear_Estimator.h"

#include <utility>
#include "../distortion_groebner_estimator/Groebner_Estimator.h"

namespace non_linear_optimization {
    void NonLinearEstimator::estimate() {
        if (!estimated_) {
            estimated_ = true;
            int residuals = 0;
            long number_of_distortion_coefficients = lambdas_.rows();
            std::vector<bool> skip(number_of_pairs, false);
            for (size_t iters = 0; iters < options_.number_of_non_linear_iters_; ++iters) {
                ceres::Problem problem;
                double *lambda_ptr = lambdas_.data();
                problem.AddParameterBlock(lambda_ptr, static_cast<int>(number_of_distortion_coefficients));

                for (size_t kth_pair = 0; kth_pair < number_of_pairs; ++kth_pair) {
                    if (skip[kth_pair])
                        continue;
                    double *f_ptr = two_view_fundamental_estimators_[kth_pair].getRawFundamentalMatrix();
                    problem.AddParameterBlock(f_ptr, 8);

                    std::vector<int> inliers_ind;
                    double interval = utils::distortion_problem::findInliers(left_pictures_keypoints_[kth_pair],
                                                                             right_pictures_keypoints_[kth_pair],
                                                                             lambdas_,
                                                                             two_view_fundamental_estimators_[kth_pair].getFundamentalMatrix(),
                                                                             options_.quantile_to_minimize_,
                                                                             inliers_ind, options_.image_radius_);
                    std::cout << "Interval: " << interval << " " << options_.image_radius_<< std::endl;
                    if (interval > options_.max_interval_) {
                        skip[kth_pair] = true;
                        std::cout << "Skip " << kth_pair + 1 << "-nth stereo pair with high confidence interval --- "
                                  << interval << std::endl;
                        continue;
                    }
                    Eigen::Matrix<double, 2, Eigen::Dynamic> i1d, i2d;
                    i1d.resize(Eigen::NoChange, inliers_ind.size());
                    i2d.resize(Eigen::NoChange, inliers_ind.size());
                    for (size_t kth_inlier = 0; kth_inlier < inliers_ind.size(); ++kth_inlier) {
                        i1d.col(kth_inlier) = left_pictures_keypoints_[kth_pair].col(inliers_ind[kth_inlier]);
                        i2d.col(kth_inlier) = right_pictures_keypoints_[kth_pair].col(inliers_ind[kth_inlier]);
                    }

                    for (size_t k = 0; k < i1d.cols(); ++k) {
                        Eigen::Vector3d left, right;
                        left.template block<2, 1>(0, 0) = i1d.col(k);
                        right.template block<2, 1>(0, 0) = i2d.col(k);
                        left[2] = right[2] = 1.0;

                        auto fun = new ceres::DynamicAutoDiffCostFunction<ErrorFunctor, 10>(
                                new ErrorFunctor(left, right, number_of_distortion_coefficients, options_.image_radius_));
                        fun->AddParameterBlock(static_cast<int>(number_of_distortion_coefficients));
                        fun->AddParameterBlock(8);
                        fun->SetNumResiduals(2);
                        problem.AddResidualBlock(fun, /*new ceres::HuberLoss(15) */nullptr, lambda_ptr, f_ptr);

                        ++residuals;
                    }

                }
                std::cout << lambdas_.transpose() << " --- coefficients before estimation" << std::endl;
                ceres::Solver::Options options;
                //options.max_trust_region_radius = 0.01;
                options.max_num_iterations = 500;
                options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
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
                std::cout << lambdas_.transpose() << " --- estimated coefficients" << std::endl;

            }


        }
    }

    bool NonLinearEstimator::isEstimated() const {
        return estimated_;
    }

    NonLinearEstimator::NonLinearEstimator(
            scene::StdVector<scene::ImagePoints> left_pictures_keypoints,
            scene::StdVector<scene::ImagePoints> right_pictures_keypoints,
            const scene::StdVector<scene::FundamentalMatrix> &fundamental_matrices,
            const intrinsics::DivisionModelIntrinsic<-1> &intrinsic_, NonLinearEstimatorOptions options)
            : estimators::internal::DivisionModelIntrinsicsEstimator<Eigen::Dynamic>(
            intrinsic_.getDistortionCoefficients(),
            intrinsic_.getFocalLength(),
            intrinsic_.getPrincipalPointX(),
            intrinsic_.getPrincipalPointY()),
              left_pictures_keypoints_(std::move(left_pictures_keypoints)),
              right_pictures_keypoints_(std::move(right_pictures_keypoints)),
              number_of_pairs(fundamental_matrices.size()),
              estimated_(false), options_(options) {
        two_view_fundamental_estimators_.resize(fundamental_matrices.size());
        for (size_t k = 0; k < fundamental_matrices.size(); ++k) {
            two_view_fundamental_estimators_[k] = SimpleFundamentalMatrixEstimator(fundamental_matrices[k]);
        }
    }

    void SimpleFundamentalMatrixEstimator::estimate() {
        Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(fundamental_matrix_, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
        singular_values[2] = 0.0;
        fundamental_matrix_ = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                              fmatrix_svd.matrixV().transpose();
    }

    bool SimpleFundamentalMatrixEstimator::isEstimated() const {
        return estimated_;
    }

    SimpleFundamentalMatrixEstimator::SimpleFundamentalMatrixEstimator() {
        estimated_ = false;
    }

    SimpleFundamentalMatrixEstimator::SimpleFundamentalMatrixEstimator(
            const scene::FundamentalMatrix &estimated_fundamental_matrix) {
        fundamental_matrix_ = estimated_fundamental_matrix;
        estimated_ = true;

    }

    double *SimpleFundamentalMatrixEstimator::getRawFundamentalMatrix() {
        return fundamental_matrix_.data();
    }

    NonLinearEstimatorOptions::NonLinearEstimatorOptions(int number_of_non_linear_iters, double quantile_to_minimize,
                                                         double max_interval,
                                                         double image_radius) :
            number_of_non_linear_iters_(number_of_non_linear_iters),
            max_interval_(max_interval),
            quantile_to_minimize_(quantile_to_minimize),
            image_radius_(image_radius) {}

    ErrorFunctor::ErrorFunctor(const Eigen::Vector3d &left_point, const Eigen::Vector3d &right_point,
                               int number_of_distortion_coefficients, double image_radius)
            : left_point_(left_point),
              right_point_(right_point), number_of_distortion_coefficients_(number_of_distortion_coefficients), image_radius_(image_radius) {}
}