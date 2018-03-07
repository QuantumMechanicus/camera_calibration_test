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
        long number_of_distortion_coefficients = lambdas_.rows();

        ceres::Problem problem;
        double *lambda_ptr = lambdas_.data();
        double *focal_ptr = &focal_length_;
        double pp_ptr[2];
        pp_ptr[0] = ppx_;
        pp_ptr[1] = ppy_;

        problem.AddParameterBlock(lambda_ptr, static_cast<int>(number_of_distortion_coefficients));
        problem.AddParameterBlock(focal_ptr, 1);
        problem.AddParameterBlock(pp_ptr, 2);
        for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
            auto &kth_fundamental_matrix = fundamental_matrices_[kth_pair];
            auto &kth_translation = translations_[kth_pair];
            auto &kth_rotation = rotation_matrices_[kth_pair];
            double *tr_ptr = kth_translation.data();
            double *rot_ptr = Sophus::SO3d(kth_rotation).data();
            problem.AddParameterBlock(tr_ptr, 3, new local_parametrization::LocalParameterizationSphere(1));
            problem.AddParameterBlock(rot_ptr, 4, new local_parametrization::LocalParameterizationSO3());

            std::vector<size_t> inliers_ind;
            double interval = utils::distortion_problem::findInliers(left_pictures_keypoints_[kth_pair],
                                                                     right_pictures_keypoints_[kth_pair],
                                                                     lambdas_,
                                                                     kth_fundamental_matrix,
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

                ++residuals;
            }

        }
        std::cout << *focal_ptr << " " << focal_length_ << " Check" << std::endl;
        for (size_t kth_pair = 0; kth_pair < number_of_pairs_; ++kth_pair) {
            Eigen::Matrix3d calibration_matrix;
            calibration_matrix.setIdentity();
            calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_ptr;
            calibration_matrix(0, 2) = pp_ptr[0];
            calibration_matrix(1, 2) = pp_ptr[1];
            fundamental_matrices_[kth_pair] =
                    utils::screw_hat(translations_[kth_pair]) * rotation_matrices_[kth_pair];

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
                                                       scene::StdVector<Eigen::Matrix3d> rotation_matrices,
                                                       scene::StdVector<Eigen::Vector3d> translations,
                                                       double focal_length, double ppx, double ppy,
                                                       GlobalNonLinearEstimatorOptions options)
            : rotation_matrices_(std::move(rotation_matrices)), translations_(std::move(translations)),
              focal_length_(focal_length), ppx_(ppx), ppy_(ppy),
              lambdas_(std::move(lambdas)), left_pictures_keypoints_(std::move(left_pictures_keypoints)),
              right_pictures_keypoints_(std::move(right_pictures_keypoints)), options_(std::move(options)),
              is_estimated_(false) {
        number_of_pairs_ = rotation_matrices_.size();


    }


}