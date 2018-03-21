//
// Created by danielbord on 2/14/18.
//

#ifndef CAMERA_CALIBRATION_NON_LINEAR_ESTIMATOR_H
#define CAMERA_CALIBRATION_NON_LINEAR_ESTIMATOR_H

#include "Core.h"


namespace non_linear_optimization {

    template<template<typename> class TCostFunctor = utils::distortion_problem::EpipolarCurveDistanceError>
    class DivisionDistortionAndFundamentalMatrixOptimizerFunctor {
        Eigen::Vector3d left_point_, right_point_;
        int number_of_distortion_coefficients_;
        double image_radius_;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        DivisionDistortionAndFundamentalMatrixOptimizerFunctor(
                Eigen::Vector3d left_point, Eigen::Vector3d right_point,
                int number_of_distortion_coefficients, double image_radius = 1)
                : left_point_(std::move(left_point)),
                  right_point_(std::move(right_point)),
                  number_of_distortion_coefficients_(number_of_distortion_coefficients),
                  image_radius_(image_radius) {}


        template<typename T>
        bool operator()(T const *const *parameters, T *residuals) const {
            const T *lambda_ptr = parameters[0];
            const T *f_ptr = parameters[1];


            using Vector3T = Eigen::Matrix<T, 3, 1>;
            using Matrix3T = Eigen::Matrix<T, 3, 3>;
            using VectorNL = Eigen::Matrix<T, Eigen::Dynamic, 1>;
            scene::TImagePoint<T> left_point_T = left_point_.topRows(2).cast<T>();
            scene::TImagePoint<T> right_point_T = right_point_.topRows(2).cast<T>();
            scene::TImagePoint<T> center;


            VectorNL lambdas = Eigen::Map<const VectorNL>(lambda_ptr, number_of_distortion_coefficients_);


            bool is_invertible = utils::distortion_problem::checkUndistortionInvertibility<T>(lambdas);

            if (!is_invertible)
                return false;

            Matrix3T F;
            for (size_t k = 0; k < 3; ++k)
                for (size_t j = 0; j < 3; ++j) {
                    if (k < 2 || j < 2)
                        F(k, j) = f_ptr[3 * j + k];
                }

            F(2, 2) = T(1);
            Eigen::JacobiSVD<Matrix3T> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Vector3T singular_values = fmatrix_svd.singularValues();
            singular_values[2] = T(0.0);
            F = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                fmatrix_svd.matrixV().transpose();
            F = F / F(2, 2);


            std::pair<T, T> res;
            bool is_correct = TCostFunctor<T>()(left_point_T,
                                                right_point_T,
                                                F,
                                                lambdas,
                                                res.first,
                                                res.second);


            residuals[0] = image_radius_ * res.first;
            residuals[1] = image_radius_ * res.second;
            return is_correct;
        }
    };

    struct NonLinearEstimatorOptions {
        int number_of_non_linear_iters_;
        double quantile_to_minimize_;
        double image_radius_;
        double max_interval_;

        explicit NonLinearEstimatorOptions(int number_of_non_linear_iters = 1, double quantile_to_minimize = 0.1,
                                           double max_interval = 10, double image_radius = 1);
    };

    class NonLinearEstimator
            : public estimators::AbstractEstimator<Eigen::VectorXd>,
              public estimators::AbstractEstimator<scene::FundamentalMatrices> {

        scene::FundamentalMatrices fundamental_matrices_;
        Eigen::VectorXd lambdas_;
        scene::StdVector<scene::ImagePoints> left_pictures_keypoints_, right_pictures_keypoints_;
        size_t number_of_pairs_;
        bool is_estimated_;
        NonLinearEstimatorOptions options_;

    protected:

        void estimateImpl() override;

        void getEstimationImpl(Eigen::VectorXd &result) override;

        void getEstimationImpl(scene::FundamentalMatrices &result) override;

    public:
        NonLinearEstimator(const scene::StdVector<scene::ImagePoints> &left_pictures_keypoints,
                           const scene::StdVector<scene::ImagePoints> &right_pictures_keypoints,
                           const scene::StdVector<scene::FundamentalMatrix> &fundamental_matrices,
                           const Eigen::VectorXd &distortion_coefficients,
                           NonLinearEstimatorOptions options = NonLinearEstimatorOptions());


        bool isEstimated() const override;
    };

}

#endif //CAMERA_CALIBRATION_NON_LINEAR_ESTIMATOR_H
