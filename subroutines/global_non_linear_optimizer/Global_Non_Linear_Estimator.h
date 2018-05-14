//
// Created by danielbord on 3/1/18.
//

#ifndef CAMERA_CALIBRATION_GLOBAL_NON_LINEAR_ESTIMATOR_H
#define CAMERA_CALIBRATION_GLOBAL_NON_LINEAR_ESTIMATOR_H

#include "Core.h"

namespace non_linear_optimization {
    template<template<typename> class TCostFunctor = utils::distortion_problem::EpipolarCurveDistanceError>
    class GlobalOptimizerFunctor {
        Eigen::Vector2d left_point_, right_point_;
        int number_of_distortion_coefficients_;
        double image_radius_;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        GlobalOptimizerFunctor(
                const Eigen::Vector2d &left_point, const Eigen::Vector2d &right_point,
                int number_of_distortion_coefficients, double image_radius = 1)
                : left_point_(left_point),
                  right_point_(right_point),
                  number_of_distortion_coefficients_(number_of_distortion_coefficients),
                  image_radius_(image_radius) {}


        template<typename T>
        bool operator()(T const *const *parameters, T *residuals) const {
            const T *lambda_ptr = parameters[0];
            const T *translation_ptr = parameters[3];
            const T *rotation_ptr = parameters[4];
            const T *focal_length_ptr = parameters[1];
            const T *principal_point_ptr = parameters[2];
            const T *world_point_ptr = parameters[5];
            //TODO points ot parameters?
            using Vector2T = Eigen::Matrix<T, 2, 1>;
            using Vector3T = Eigen::Matrix<T, 3, 1>;
            using Matrix3T = Eigen::Matrix<T, 3, 3>;
            using VectorNT = Eigen::Matrix<T, Eigen::Dynamic, 1>;


            VectorNT lambdas = Eigen::Map<const VectorNT>(lambda_ptr, number_of_distortion_coefficients_);
            Vector3T wp = Eigen::Map<const Vector3T>(world_point_ptr);
            Matrix3T calibration = Matrix3T::Identity();
            Sophus::SO3<T> rotation = Eigen::Map<const Sophus::SO3<T>>(rotation_ptr);
            Vector3T translation = Eigen::Map<const Vector3T>(translation_ptr);
            Eigen::Map<Eigen::Matrix<T, 4, 1>> residual(residuals);

            if (wp(2) < T(0.0)){
                residual = T(std::numeric_limits<double>::max())*Eigen::Matrix<T, 4, 1>::Ones();

                LOG(INFO) << "____It happend____\n";
                LOG(INFO) << wp;
                LOG(INFO) << "_________________"<<std::endl;
                return false;
            }

            calibration(0, 0) = calibration(1, 1) = *focal_length_ptr;
            calibration(0, 2) = principal_point_ptr[0];
            calibration(1, 2) = principal_point_ptr[1];
            Vector2T lip = (calibration * wp).hnormalized();
            Vector2T rip = (calibration * (rotation * wp + translation)).hnormalized();

            Vector2T ldp = utils::distortion_problem::distortion(lip, lambdas);
            Vector2T rdp = utils::distortion_problem::distortion(rip, lambdas);

            residual.template head<2>() = (ldp - left_point_.cast<T>()) * image_radius_;
            residual.template tail<2>() = (rdp - right_point_.cast<T>()) * image_radius_;

            return true;
        }
    };

    struct GlobalNonLinearEstimatorOptions {
        double quantile_to_minimize_;
        double image_radius_;
        double max_interval_;
        double w_, h_;
        explicit GlobalNonLinearEstimatorOptions(double quantile_to_minimize = 0.1,
                                                 double max_interval = 10, double image_radius = 1, double w = 0, double h = 0);
    };

    class GlobalNonLinearEstimator
            : public estimators::AbstractEstimator<Eigen::VectorXd>,
              public estimators::AbstractEstimator<scene::FundamentalMatrices>,
    public estimators::AbstractEstimator<intrinsics::FocalLength >{

        scene::StdVector<Sophus::SO3d> rotations_;
        scene::StdVector<Eigen::Vector3d> translations_;
        scene::FundamentalMatrices fundamental_matrices_;
        double focal_length_;
        double ppx_, ppy_;
        Eigen::VectorXd lambdas_;
        scene::StdVector<scene::ImagePoints> left_pictures_keypoints_, right_pictures_keypoints_;
        size_t number_of_pairs_;
        bool is_estimated_;
        GlobalNonLinearEstimatorOptions options_;


    protected:

        void estimateImpl() override;

        void getEstimationImpl(Eigen::VectorXd &result) override;

        void getEstimationImpl(scene::FundamentalMatrices &result) override;

        void getEstimationImpl(intrinsics::FocalLength &result) override;

    public:
        GlobalNonLinearEstimator(const scene::StdVector<scene::ImagePoints> &left_pictures_keypoints,
                                                           const scene::StdVector<scene::ImagePoints> &right_pictures_keypoints,
                                                           const Eigen::VectorXd &lambdas,
                                                           const scene::StdVector<Sophus::SO3d> &rotations,
                                                           const scene::StdVector<Eigen::Vector3d> &translations,
                                                           double focal_length, double ppx, double ppy,
                                                           GlobalNonLinearEstimatorOptions options = GlobalNonLinearEstimatorOptions());


        bool isEstimated() const override;
    };

}
#endif //CAMERA_CALIBRATION_GLOBAL_NON_LINEAR_ESTIMATOR_H
