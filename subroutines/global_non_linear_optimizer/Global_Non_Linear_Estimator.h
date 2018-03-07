//
// Created by danielbord on 3/1/18.
//

#ifndef CAMERA_CALIBRATION_GLOBAL_NON_LINEAR_ESTIMATOR_H
#define CAMERA_CALIBRATION_GLOBAL_NON_LINEAR_ESTIMATOR_H

#include "Core.h"

namespace non_linear_optimization {
    template<template<typename> class TCostFunctor = utils::distortion_problem::EpipolarCurveDistanceError>
    class GlobalOptimizerFunctor {
        Eigen::Vector3d left_point_, right_point_;
        int number_of_distortion_coefficients_;
        double image_radius_;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        GlobalOptimizerFunctor(
                Eigen::Vector3d left_point, Eigen::Vector3d right_point,
                int number_of_distortion_coefficients, double image_radius = 1)
                : left_point_(std::move(left_point)),
                  right_point_(std::move(right_point)),
                  number_of_distortion_coefficients_(number_of_distortion_coefficients),
                  image_radius_(image_radius) {}


        template<typename T>
        bool operator()(T const *const *parameters, T *residuals) const {
            const T *lambda_ptr = parameters[0];
            const T *translation_ptr = parameters[3];
            const T *rotation_ptr = parameters[4];
            const T *focal_length_ptr = parameters[1];
            const T *principal_point_ptr = parameters[2];
            //TODO points ot parameters?

            using Vector3T = Eigen::Matrix<T, 3, 1>;
            using Matrix3T = Eigen::Matrix<T, 3, 3>;
            using VectorNL = Eigen::Matrix<T, Eigen::Dynamic, 1>;


            VectorNL lambdas = Eigen::Map<const VectorNL>(lambda_ptr, number_of_distortion_coefficients_);
            bool is_invertible = utils::distortion_problem::checkUndistortionInvertibility<T>(lambdas);
            if (!is_invertible)
                return false;

            scene::TImagePoint<T> left_point_T = left_point_.topRows(2).cast<T>();
            scene::TImagePoint<T> right_point_T = right_point_.topRows(2).cast<T>();

            const Eigen::Map<const Matrix3T> rotation(rotation_ptr);
            const Eigen::Map<const Eigen::Matrix<T, 3, 1>> translation(translation_ptr);
            Matrix3T translation_matrix = utils::screw_hat<T>(translation);
            Matrix3T calibration_matrix;
            calibration_matrix.setIdentity();
            calibration_matrix(0, 0) = calibration_matrix(1, 1) = *focal_length_ptr;
            calibration_matrix(0, 2) = principal_point_ptr[0];
            calibration_matrix(1, 2) = principal_point_ptr[1];

            Matrix3T fundamental_matrix =
                    calibration_matrix.transpose().inverse() * translation_matrix * rotation *
                    calibration_matrix.inverse();


            Eigen::JacobiSVD<Matrix3T> fmatrix_svd(fundamental_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Vector3T singular_values = fmatrix_svd.singularValues();
            singular_values[2] = T(0.0);
            fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                 fmatrix_svd.matrixV().transpose();
            fundamental_matrix = fundamental_matrix / fundamental_matrix(2, 2);


            std::pair<T, T> res;
            bool is_correct = TCostFunctor<T>()(left_point_T,
                                                right_point_T,
                                                fundamental_matrix,
                                                lambdas,
                                                res.first,
                                                res.second);


            residuals[0] = image_radius_ * res.first;
            residuals[1] = image_radius_ * res.second;
            return is_correct;
        }
    };

    struct GlobalNonLinearEstimatorOptions {
        double quantile_to_minimize_;
        double image_radius_;
        double max_interval_;

        explicit GlobalNonLinearEstimatorOptions(double quantile_to_minimize = 0.1,
                                           double max_interval = 10, double image_radius = 1);
    };

    class GlobalNonLinearEstimator
            : public estimators::AbstractEstimator<Eigen::RowVectorXd>,
              public estimators::AbstractEstimator<scene::FundamentalMatrices> {

        scene::StdVector<Eigen::Matrix3d> rotation_matrices_;
        scene::StdVector<Eigen::Vector3d> translations_;
        scene::FundamentalMatrices fundamental_matrices_;
        double focal_length_;
        double ppx_, ppy_;
        Eigen::RowVectorXd lambdas_;
        scene::StdVector<scene::ImagePoints> left_pictures_keypoints_, right_pictures_keypoints_;
        size_t number_of_pairs_;
        bool is_estimated_;
        GlobalNonLinearEstimatorOptions options_;


    protected:

        void estimateImpl() override;

        void getEstimationImpl(Eigen::RowVectorXd &result) override;

        void getEstimationImpl(scene::FundamentalMatrices &result) override;

    public:
        GlobalNonLinearEstimator(scene::StdVector<scene::ImagePoints> left_pictures_keypoints,
                                 scene::StdVector<scene::ImagePoints> right_pictures_keypoints,
                                 Eigen::RowVectorXd lambdas,
                                 scene::StdVector<Eigen::Matrix3d> rotation_matrices,
                                 scene::StdVector<Eigen::Vector3d> translations,
                                 double focal_length, double ppx, double ppy, GlobalNonLinearEstimatorOptions options = GlobalNonLinearEstimatorOptions());


        bool isEstimated() const override;
    };

}
#endif //CAMERA_CALIBRATION_GLOBAL_NON_LINEAR_ESTIMATOR_H
