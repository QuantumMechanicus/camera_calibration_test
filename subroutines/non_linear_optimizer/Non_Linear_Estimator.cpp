//
// Created by danielbord on 2/14/18.
//

#include "Non_Linear_Estimator.h"


namespace non_linear_optimization {
    NonLinearEstimatorOptions::NonLinearEstimatorOptions(int number_of_non_linear_iters, double quantile_to_minimize,
                                                         double max_interval,
                                                         double image_radius) :
            number_of_non_linear_iters_(number_of_non_linear_iters),
            max_interval_(max_interval),
            quantile_to_minimize_(quantile_to_minimize),
            image_radius_(image_radius) {}

    void NonLinearEstimator::estimateImpl() {


        std::cout << "Iters: " << options_.number_of_non_linear_iters_ << std::endl;

        is_estimated_ = true;
        int residuals = 0;
        long number_of_distortion_coefficients = lambdas_.cols();
        std::vector<bool> skip(number_of_pairs_, false);
        for (size_t iters = 0; iters < options_.number_of_non_linear_iters_; ++iters) {
            ceres::Problem problem;
            double *lambda_ptr = lambdas_.data();
            problem.AddParameterBlock(lambda_ptr, static_cast<int>(number_of_distortion_coefficients));

            for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
                if (skip[kth_pair])
                    continue;
                auto &kth_fundamental_matrix = fundamental_matrices_[kth_pair];
                Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(kth_fundamental_matrix,
                                                              Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
                singular_values[2] = 0.0;
                kth_fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                         fmatrix_svd.matrixV().transpose();

                double *f_ptr = kth_fundamental_matrix.data();
                problem.AddParameterBlock(f_ptr, 8);

                std::vector<size_t> inliers_ind;
                double interval = utils::distortion_problem::findInliers(left_pictures_keypoints_[kth_pair],
                                                                         right_pictures_keypoints_[kth_pair],
                                                                         lambdas_,
                                                                         kth_fundamental_matrix,
                                                                         options_.quantile_to_minimize_,
                                                                         inliers_ind, options_.image_radius_);
                std::cout << "Interval: " << interval << " " << options_.image_radius_ << std::endl;
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

                    auto fun = new ceres::DynamicAutoDiffCostFunction<DivisionDistortionAndFundamentalMatrixOptimizerFunctor<>>(
                            new DivisionDistortionAndFundamentalMatrixOptimizerFunctor<>(left, right,
                                                                                         static_cast<int>(number_of_distortion_coefficients),
                                                                                         options_.image_radius_));
                    fun->AddParameterBlock(static_cast<int>(number_of_distortion_coefficients));
                    fun->AddParameterBlock(8);
                    fun->SetNumResiduals(2);
                    problem.AddResidualBlock(fun, nullptr, lambda_ptr, f_ptr);

                    ++residuals;
                }

            }
            std::cout << lambdas_ << " --- coefficients before estimation" << std::endl;
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
            std::cout << lambdas_ << " --- estimated coefficients" << std::endl;

        }

    }

    void NonLinearEstimator::getEstimationImpl(Eigen::RowVectorXd &result) {
        result = lambdas_;
    }

    void NonLinearEstimator::getEstimationImpl(scene::FundamentalMatrices &result) {
        for (auto &kth_fundamental_matrix : fundamental_matrices_) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(kth_fundamental_matrix,
                                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;
            kth_fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                     fmatrix_svd.matrixV().transpose();
            kth_fundamental_matrix /= kth_fundamental_matrix(2,2);

        }
        result = fundamental_matrices_;
    }

    NonLinearEstimator::NonLinearEstimator(scene::StdVector<scene::ImagePoints> left_pictures_keypoints,
                                           scene::StdVector<scene::ImagePoints> right_pictures_keypoints,
                                           scene::StdVector<scene::FundamentalMatrix> fundamental_matrices,
                                           Eigen::RowVectorXd distortion_coefficients,
                                           NonLinearEstimatorOptions options) :
            left_pictures_keypoints_(std::move(left_pictures_keypoints)),
            right_pictures_keypoints_(std::move(right_pictures_keypoints)),
            lambdas_(std::move(distortion_coefficients)),
            fundamental_matrices_(std::move(fundamental_matrices)),
            is_estimated_(false),
            options_(options) {
        number_of_pairs_ = left_pictures_keypoints_.size();
    }

    bool NonLinearEstimator::isEstimated() const {
        return is_estimated_;
    }
}