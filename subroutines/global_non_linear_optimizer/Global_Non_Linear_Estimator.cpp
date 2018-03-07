//
// Created by danielbord on 3/1/18.
//

#include "Global_Non_Linear_Estimator.h"

namespace non_linear_optimization {

    GlobalNonLinearEstimatorOptions::GlobalNonLinearEstimatorOptions(double quantile_to_minimize,
                                                                     double max_interval,
                                                                     double image_radius) :
            max_interval_(max_interval),
            quantile_to_minimize_(quantile_to_minimize),
            image_radius_(image_radius) {}

    void GlobalNonLinearEstimator::estimateImpl() {
        is_estimated_ = true;
        int residuals = 0;
        long number_of_distortion_coefficients = lambdas_.cols();

        ceres::Problem problem;
        double *lambda_ptr = lambdas_.data();
        double *focal_ptr = &focal_length_;
        double pp_ptr[2];
        pp_ptr[0] = ppx_;
        pp_ptr[1] = ppy_;

        problem.AddParameterBlock(lambda_ptr, static_cast<int>(number_of_distortion_coefficients));
        //problem.SetParameterBlockConstant(lambda_ptr);
        problem.AddParameterBlock(focal_ptr, 1);
        problem.AddParameterBlock(pp_ptr, 2);
        //problem.SetParameterBlockConstant(focal_ptr);
        //problem.SetParameterBlockConstant(pp_ptr);
        Eigen::Matrix3d calibration_matrix;
        calibration_matrix.setIdentity();
        calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_ptr;
        calibration_matrix(0, 2) = pp_ptr[0];
        calibration_matrix(1, 2) = pp_ptr[1];
        std::cout << *focal_ptr << " " << focal_length_ << " Check" << std::endl;

        for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
            auto &kth_translation = translations_[kth_pair];
            auto &kth_rotation = rotations_[kth_pair];

            //TODO to function essintial to fund
            fundamental_matrices_[kth_pair] =
                    calibration_matrix.inverse().transpose() * utils::screw_hat(kth_translation) *
                    kth_rotation.matrix() * calibration_matrix.inverse();
            std::cout << fundamental_matrices_[kth_pair] << std::endl;

            double *tr_ptr = kth_translation.data();
            double *rot_ptr = kth_rotation.data();
            problem.AddParameterBlock(tr_ptr, 3, new local_parametrization::LocalParameterizationSphere(1));
            problem.AddParameterBlock(rot_ptr, 4, new local_parametrization::LocalParameterizationSO3());


            std::vector<size_t> inliers_ind;
            double interval = utils::distortion_problem::findInliers(left_pictures_keypoints_[kth_pair],
                                                                     right_pictures_keypoints_[kth_pair],
                                                                     lambdas_,
                                                                     fundamental_matrices_[kth_pair],
                                                                     options_.quantile_to_minimize_,
                                                                     inliers_ind, options_.image_radius_);

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

                auto fun = new ceres::DynamicAutoDiffCostFunction<GlobalOptimizerFunctor<>>(
                        new GlobalOptimizerFunctor<>(left, right, static_cast<int>(number_of_distortion_coefficients),
                                                     options_.image_radius_));
                fun->AddParameterBlock(static_cast<int>(number_of_distortion_coefficients));
                fun->AddParameterBlock(1);
                fun->AddParameterBlock(2);
                fun->AddParameterBlock(3);
                fun->AddParameterBlock(4);
                fun->SetNumResiduals(2);
                problem.AddResidualBlock(fun, nullptr, lambda_ptr, focal_ptr, pp_ptr, tr_ptr, rot_ptr);
                double ress[2];
                double *ddd[] = {lambda_ptr, focal_ptr, pp_ptr, tr_ptr, rot_ptr};
                std::cout << fun->Evaluate(ddd, ress, NULL) << " " << residuals << std::endl;
                std::cout << ress[0] << " " << ress[1] << std::endl;
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
        std::cout << focal_length_ << " --- focal length" << std::endl;
        std::cout << ppx_ << " " << ppy_ << std::endl;
        std::cout << 2*atan(1/focal_length_)*180/M_PI << std::endl;
        std::cout << std::sqrt(summary.final_cost*2/residuals) << std::endl;
        for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
            calibration_matrix.setIdentity();
            calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_ptr;
            calibration_matrix(0, 2) = pp_ptr[0];
            calibration_matrix(1, 2) = pp_ptr[1];
            fundamental_matrices_[kth_pair] =
                    calibration_matrix.inverse().transpose() * utils::screw_hat(translations_[kth_pair]) *
                    rotations_[kth_pair].matrix() * calibration_matrix.inverse();

        }

    }


    void GlobalNonLinearEstimator::getEstimationImpl(Eigen::RowVectorXd &result) {
        result = lambdas_;
    }

    void GlobalNonLinearEstimator::getEstimationImpl(scene::FundamentalMatrices &result) {
        for (auto &kth_fundamental_matrix : fundamental_matrices_) {
            Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(kth_fundamental_matrix,
                                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
            singular_values[2] = 0.0;
            kth_fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                     fmatrix_svd.matrixV().transpose();
            kth_fundamental_matrix /= kth_fundamental_matrix(2, 2);

        }
        result = fundamental_matrices_;
    }


    bool GlobalNonLinearEstimator::isEstimated() const {
        return is_estimated_;
    }

    GlobalNonLinearEstimator::GlobalNonLinearEstimator(scene::StdVector<scene::ImagePoints> left_pictures_keypoints,
                                                       scene::StdVector<scene::ImagePoints> right_pictures_keypoints,
                                                       Eigen::RowVectorXd lambdas,
                                                       scene::StdVector<Sophus::SO3d> rotations,
                                                       scene::StdVector<Eigen::Vector3d> translations,
                                                       double focal_length, double ppx, double ppy,
                                                       GlobalNonLinearEstimatorOptions options)
            : rotations_(std::move(rotations)), translations_(std::move(translations)),
              focal_length_(focal_length), ppx_(ppx), ppy_(ppy),
              lambdas_(std::move(lambdas)), left_pictures_keypoints_(std::move(left_pictures_keypoints)),
              right_pictures_keypoints_(std::move(right_pictures_keypoints)), options_(options),
              is_estimated_(false) {
        number_of_pairs_ = rotations_.size();
        fundamental_matrices_.resize(number_of_pairs_);

    }


}