//
// Created by danielbord on 3/1/18.
//

#ifndef CAMERA_CALIBRATION_LOCAL_PARAMETRIZATION_SPHERE_H
#define CAMERA_CALIBRATION_LOCAL_PARAMETRIZATION_SPHERE_H
#include <ceres/local_parameterization.h>
#include "Eigen/Dense"
namespace local_parametrization {
    class LocalParameterizationSphere : public ceres::LocalParameterization {
    private:
        const double radius_;

        static void calculateBasis(const Eigen::Vector3d &p, Eigen::Vector3d &du, Eigen::Vector3d &dv) {
            Eigen::Vector3d x, y, z, c;
            x << 1, 0, 0;
            y << 0, 1, 0;
            z << 0, 0, 1;
            double maxNonCollinearity = x.cross(p).norm();
            c = x;
            if (y.cross(p).norm() > maxNonCollinearity) {
                c = y;
                maxNonCollinearity = y.cross(p).norm();
            }
            if (z.cross(p).norm() > maxNonCollinearity) {
                c = z;
            }

            du = (p - c).cross(p);
            dv = du.cross(p);
            du.normalize();
            dv.normalize();
        }

    public:
        explicit LocalParameterizationSphere(double norm)
                : radius_(norm) {}

        ~LocalParameterizationSphere() final {}


        bool Plus(const double *x_raw, const double *delta_raw, double *x_plus_delta_raw) const final {

            Eigen::Map<Eigen::Vector3d const> x(x_raw);
            Eigen::Map<Eigen::Vector3d> x_plus_delta(x_plus_delta_raw);

            Eigen::Vector3d du, dv;
            calculateBasis(x, du, dv);

            x_plus_delta = x + delta_raw[0] * du + delta_raw[1] * dv;

            x_plus_delta *= radius_ / x_plus_delta.norm();

            return true;
        }

        bool ComputeJacobian(const double *x_raw, double *jacobian_raw) const final {
            Eigen::Map<Eigen::Vector3d const> x(x_raw);
            Eigen::Map<Eigen::Matrix<double, 3, 2, Eigen::RowMajor> > jacobian(jacobian_raw);
            Eigen::Vector3d du, dv;
            calculateBasis(x, du, dv);

            jacobian(0, 0) = du[0];
            jacobian(1, 0) = du[1];
            jacobian(2, 0) = du[2];

            jacobian(0, 1) = dv[0];
            jacobian(1, 1) = dv[1];
            jacobian(2, 1) = dv[2];

            return true;
        }

        int GlobalSize() const final {
            return 3;
        }

        int LocalSize() const final {
            return 2;
        }

    };
}

#endif //CAMERA_CALIBRATION_LOCAL_PARAMETRIZATION_SPHERE_H
