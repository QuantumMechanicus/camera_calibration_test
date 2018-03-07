//
// Created by danielbord on 2/28/18.
//
#include "Focal_Estimator.h"

#include <utility>

bool estimators::SimpleFocalEstimator::isEstimated() const {
    return is_estimated_;
}

void estimators::SimpleFocalEstimator::estimateImpl() {
    is_estimated_ = true;

    Eigen::Matrix<double, 8, 1> sum_of_derivatives_of_squared_equations;
    Eigen::Vector3d sum_of_quadric_equations;
    Eigen::Vector2d sum_of_linear_equations;

    sum_of_derivatives_of_squared_equations.setZero();
    sum_of_quadric_equations.setZero();
    sum_of_linear_equations.setZero();
    for (const auto &fundamental_matrix : fundamental_matrices_) {

        sum_of_derivatives_of_squared_equations = sum_of_derivatives_of_squared_equations + makeSquaredQuadricEquation(
                fundamental_matrix);
        sum_of_quadric_equations = sum_of_quadric_equations + makeQuadricEquation(fundamental_matrix);
        sum_of_linear_equations =
                sum_of_linear_equations + makeLinearEquation1(fundamental_matrix) + makeLinearEquation2(
                        fundamental_matrix);
    }

    std::cout << "Poly sum_of_derivatives_of_squared_equations: " << sum_of_derivatives_of_squared_equations.transpose()
              << "\n";

    sum_of_derivatives_of_squared_equations =
            sum_of_derivatives_of_squared_equations / sum_of_derivatives_of_squared_equations(7);

    std::cout << "Poly sum_of_derivatives_of_squared_equations divided: "
              << sum_of_derivatives_of_squared_equations.transpose() << "\n";

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> companion(7, 7);
    companion.setZero();
    companion.col(6) = -1 * sum_of_derivatives_of_squared_equations.topRows(7);
    companion.template block<6, 6>(1, 0).setIdentity();

    double f = std::numeric_limits<double>::max();
    auto r = companion.eigenvalues();
    double h2 = std::numeric_limits<double>::max();
    for (int kk = 0; kk < r.rows(); ++kk) {
        double real = r[kk].real();
        double imag = r[kk].imag();

        std::cout << r[kk] << std::endl;
        if (/*std::abs(imag) < 1e-9 &&*/ real > 1e-5) {
            double cur_f2 = real * real;
            double cur_h2 = (sum_of_quadric_equations(2) * cur_f2 * cur_f2 + sum_of_quadric_equations(1) * cur_f2 +
                             sum_of_quadric_equations(0));
            cur_h2 = cur_h2 * cur_h2;

            cur_h2 += (sum_of_linear_equations(1) * cur_f2 + sum_of_linear_equations(0)) *
                      (sum_of_linear_equations(1) * cur_f2 + sum_of_linear_equations(0));
            if (cur_h2 < h2) {
                h2 = cur_h2;
                f = real;
            }
            std::cout << "real root: " << real << " " << cur_h2 << std::endl;
        }

    }
    f_ = f;
}

void estimators::SimpleFocalEstimator::getEstimationImpl(intrinsics::FocalLength &result) {
    result = f_;
}

Eigen::Vector2d estimators::SimpleFocalEstimator::makeLinearEquation1(const Eigen::Matrix3d &F) {
    Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    double a, b, u13, u23, v13, v23;
    a = fmatrix_svd.singularValues()[0];
    b = fmatrix_svd.singularValues()[1];

    u13 = fmatrix_svd.matrixU()(2, 0);
    u23 = fmatrix_svd.matrixU()(2, 1);
    v13 = fmatrix_svd.matrixV()(2, 0);
    v23 = fmatrix_svd.matrixV()(2, 1);

    Eigen::Vector2d res;

    res[1] = a * u13 * u23 * (1 - v13 * v13) + b * v13 * v23 * (1 - u23 * u23);
    res[0] = u23 * v13 * (a * u13 * v13 + b * u23 * v23);
    res = res / res[1];
    return res;

}

Eigen::Vector2d estimators::SimpleFocalEstimator::makeLinearEquation2(const Eigen::Matrix3d &F) {
    Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    double a, b, u13, u23, v13, v23;
    a = fmatrix_svd.singularValues()[0];
    b = fmatrix_svd.singularValues()[1];

    u13 = fmatrix_svd.matrixU()(2, 0);
    u23 = fmatrix_svd.matrixU()(2, 1);
    v13 = fmatrix_svd.matrixV()(2, 0);
    v23 = fmatrix_svd.matrixV()(2, 1);

    Eigen::Vector2d res;

    res[1] = a * v13 * v23 * (1 - u13 * u13) + b * u13 * u23 * (1 - v23 * v23);
    res[0] = v23 * u13 * (a * u13 * v13 + b * u23 * v23);
    res = res / res[1];


    return res;

}

Eigen::Vector3d estimators::SimpleFocalEstimator::makeQuadricEquation(const Eigen::Matrix3d &F) {
    Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    double a, b, u13, u23, v13, v23;
    a = fmatrix_svd.singularValues()[0];
    b = fmatrix_svd.singularValues()[1];

    u13 = fmatrix_svd.matrixU()(2, 0);
    u23 = fmatrix_svd.matrixU()(2, 1);
    v13 = fmatrix_svd.matrixV()(2, 0);
    v23 = fmatrix_svd.matrixV()(2, 1);

    Eigen::Vector3d res;
    res[2] = a * a * (1 - u13 * u13) * (1 - v13 * v13)
             - b * b * (1 - u23 * u23) * (1 - v23 * v23);
    res[1] = a * a * (u13 * u13 + v13 * v13 - 2 * u13 * u13 * v13 * v13)
             - b * b * (u23 * u23 + v23 * v23 - 2 * u23 * u23 * v23 * v23);
    res[0] = a * a * u13 * u13 * v13 * v13 - b * b * u23 * u23 * v23 * v23;
    res = res / res[2];

    return res;

}

Eigen::Matrix<double, 8, 1> estimators::SimpleFocalEstimator::makeSquaredQuadricEquation(const Eigen::Matrix3d &F) {
    Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Vector3d coeff = makeQuadricEquation(F);
    double a, b, c;
    a = coeff(2);
    b = coeff(1);
    c = coeff(0);
    Eigen::Matrix<double, 8, 1> res;
    res.setZero();
    res(7) = 8 * a * a;
    res(5) = 12 * a * b;
    res(3) = 4 * b * b + 8 * a * c;
    res(1) = 4 * b * c;
    res = res / res(7);
    return res;

}

estimators::SimpleFocalEstimator::SimpleFocalEstimator(scene::FundamentalMatrices fundamental_matrices_, unsigned int w,
                                                       unsigned int h)
        : fundamental_matrices_(std::move(fundamental_matrices_)), f_{}, is_estimated_(false), w_(w), h_(h) {}

