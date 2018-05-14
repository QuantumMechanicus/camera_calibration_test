//
// Created by danielbord on 3/1/18.
//

#ifndef CAMERA_CALIBRATION_LOCAL_PARAMETRIZATION_SO3_H
#define CAMERA_CALIBRATION_LOCAL_PARAMETRIZATION_SO3_H

#include "ceres/local_parameterization.h"
#include "sophus/so3.hpp"

namespace local_parametrization {

    class LocalParameterizationSO3 : public ceres::LocalParameterization {
    public:
        ~LocalParameterizationSO3() final {}

        // SO3 plus operation for Ceres
        //
        //  T * exp(x)
        //
        bool Plus(double const *T_raw, double const *delta_raw,
                  double *T_plus_delta_raw) const final {
            Eigen::Map<Sophus::SO3d const> const T(T_raw);
            Eigen::Map<Eigen::Vector3d const> const delta(delta_raw);
            Eigen::Map<Sophus::SO3d> T_plus_delta(T_plus_delta_raw);
            T_plus_delta = T * Sophus::SO3d::exp(delta);
            return true;
        }

        // Jacobian of SO3 plus operation for Ceres
        //
        // dx T * exp(x)  with  x=0
        //
        bool ComputeJacobian(double const *T_raw,
                             double *jacobian_raw) const final {
            Eigen::Map<Sophus::SO3d const> T(T_raw);
            Eigen::Map<Eigen::Matrix<double, Sophus::SO3d::num_parameters, Sophus::SO3d::DoF,  Eigen::RowMajor> > jacobian(jacobian_raw);
            jacobian = T.Dx_this_mul_exp_x_at_0();
            return true;

        }


        int GlobalSize() const final { return Sophus::SO3d::num_parameters; }

        int LocalSize() const final { return Sophus::SO3d::DoF; }
    };
}
#endif //CAMERA_CALIBRATION_LOCAL_PARAMETRIZATION_SO3_H
