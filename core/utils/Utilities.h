//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_UTILITIES_H
#define CAMERA_CALIBRATION_UTILITIES_H

#include <fstream>
#include <boost/math/special_functions/erf.hpp>
#include <boost/filesystem.hpp>
#include "ceres/ceres.h"
#include "../scene/Camera.h"
#include "Local_Parametrization_SO3.h"
#include "Local_Parametrization_Sphere.h"

#include <glog/logging.h>

namespace utils {
    static const double EPS = 1e-10;

    template<typename Scalar = double, int RowsAtCompileTime = Eigen::Dynamic, int ColsAtCompileTime = Eigen::Dynamic>
    inline bool loadMatrix(const std::string &filename, Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> &m,
                           bool transposed = false) {
        if (filename.empty())
            return false;
        std::ifstream input(filename.c_str());
        if (input.fail()) {
            std::cerr << "Cannot find file '" << filename << "'." << std::endl;
            m.setZero();
            return false;
        }
        std::string line;
        Scalar d;

        std::vector<Scalar> v;
        std::size_t n_rows = 0;
        while (getline(input, line)) {
            ++n_rows;
            std::stringstream input_line(line);
            while (!input_line.eof()) {
                input_line >> d;
                v.push_back(d);
            }
        }
        input.close();

        std::size_t n_cols = v.size() / n_rows;
        if (transposed)
            std::swap(n_cols, n_rows);

        if (RowsAtCompileTime == Eigen::Dynamic)
            m.resize(n_rows, Eigen::NoChange);
        if (ColsAtCompileTime == Eigen::Dynamic)
            m.resize(Eigen::NoChange, n_cols);

        CHECK_EQ(m.cols(), n_cols) << "Invalid column count";
        CHECK_EQ(m.rows(), n_rows) << "Invalid row count";

        if (transposed)
            std::swap(n_cols, n_rows);
        for (int i = 0; i < n_rows; i++)
            for (int j = 0; j < n_cols; j++)
                if (transposed)
                    m(j, i) = v[i * n_cols + j];
                else
                    m(i, j) = v[i * n_cols + j];


        return true;
    }

    template<typename Scalar = double, int RowsAtCompileTime = Eigen::Dynamic, int ColsAtCompileTime = Eigen::Dynamic>
    inline bool saveMatrix(const std::string &directory, std::string filename,
                           Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> matrix, bool overwrite = false) {

        if (directory.empty() or filename.empty())
            return false;


        if (!boost::filesystem::exists(directory)) {
            // Directory doesn't exist. Try to create it.
            if (!boost::filesystem::create_directories(directory)) {
                std::cerr << "Couldn't make directory '" << directory << "'. Not saving data." << std::endl;
                return false;
            }
        }

        filename = directory + "/" + filename;
        saveMatrix(filename, matrix, overwrite);

        return true;
    }

    template<typename Scalar = double, int RowsAtCompileTime = Eigen::Dynamic, int ColsAtCompileTime = Eigen::Dynamic>
    inline bool
    saveMatrix(const std::string &filename, Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> matrix,
               bool overwrite = false) {
        if (filename.empty())
            return false;
        if (boost::filesystem::exists(filename)) {
            if (!overwrite) {
                // File exists, but overwriting is not allowed. Abort.
                std::cerr << "File '" << filename << "' already exists. Not saving data." << std::endl;
                return false;
            }
        }


        std::ofstream file;
        file.open(filename.c_str());
        if (!file.is_open()) {
            std::cerr << "Couldn't open file '" << filename << "' for writing." << std::endl;
            return false;
        }

        file << std::fixed;
        file << matrix;
        file.close();

        return true;

    }

    template<typename T>
    auto solvePoly(Eigen::Matrix<T, Eigen::Dynamic, 1> &coefficients) {
        auto deg = coefficients.size() - 1;
        coefficients /= coefficients[deg];
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> companion(deg, deg);
        companion.setZero();
        companion.col(deg - 1) = -T(1) * coefficients.topRows(deg);
        companion.block(1, 0, deg - 1, deg - 1).setIdentity();

        return companion.eigenvalues();
    }

    template<typename Scalar>
    Eigen::Matrix<Scalar, 3, 3> screw_hat(const Eigen::Matrix<Scalar, 3, 1> &t) {
        Eigen::Matrix<Scalar, 3, 3> t_hat = Eigen::Matrix<Scalar, 3, 3>::Zero(3, 3);
        t_hat << Scalar(0), -t(2), t(1),
                t(2), Scalar(0), -t(0),
                -t(1), t(0), Scalar(0);
        return t_hat;
    }

    template<typename Scalar>
    Eigen::Matrix<Scalar, 3, 1> inverted_screw_hat(const Eigen::Matrix<Scalar, 3, 3> &t_hat) {
        Eigen::Matrix<Scalar, 3, 1> t;
        t << t_hat(2, 1), t_hat(0, 2), t_hat(1, 0);
        return t;
    }





    namespace distortion_problem {

        double estimateQuantile(std::vector<double> errors,
                                double expected_percent_of_inliers);

        double estimateConfidenceInterval(double quantile, double expected_percent_of_inliers);

        double findInliers(const scene::ImagePoints &u1d,
                           const scene::ImagePoints &u2d,
                           const Eigen::Matrix<double, Eigen::Dynamic, 1> &distortion_coefficients,
                           const scene::FundamentalMatrix &fundamental_matrix, double expected_percent_of_inliers,
                           std::vector<size_t> &inliers_indices, double image_r = 1.0);

        double findInliers(const scene::ImagePoints &u1d,
                           const scene::ImagePoints &u2d,
                           const Eigen::Matrix<double, Eigen::Dynamic, 1> &distortion_coefficients,
                           const Sophus::SE3d &leftToRight, const Eigen::Matrix3d &calibration,
                           double expected_percent_of_inliers,
                           std::vector<size_t> &inliers_indices, double image_r = 1.0);


        template<typename T>
        T undistortionDenominator(const T &r_distorted,
                                  const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients) {
            T denominator(1.0);
            T r_distorted2 = r_distorted * r_distorted;
            T r_distorted2_pow = r_distorted2;
            for (int i = 0; i < distortion_coefficients.rows(); ++i) {
                denominator += distortion_coefficients[i] * r_distorted2_pow;
                r_distorted2_pow *= r_distorted2;
            }
            return denominator;
        }

        template<typename T>
        T undistortionDenominatorDerivative(const T &r_distorted,
                                            const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients) {
            T denominator(0.0);
            T c(2.0);
            T r_distorted2 = r_distorted * r_distorted;
            T r_distorted_pow = r_distorted;
            for (int i = 0; i < distortion_coefficients.rows(); ++i) {
                denominator += c * distortion_coefficients[i] * r_distorted_pow;
                c = c + T(2.0);
                r_distorted_pow *= r_distorted2;
            }

            return denominator;
        }

        template<typename T>
        bool checkUndistortionInvertibility(const Eigen::Matrix<T, Eigen::Dynamic, 1> &lambdas) {
            T k1 = lambdas(0);
            T k2 = T(0);
            if (lambdas.rows() > 1) {
                k2 = lambdas(1);
            }
            return (k1 > T(-2) and ((k1 > T(2) and (T(-1) - k1 < k2 and k2 < -(k1 * k1 / T(12)))) or
                                    (k1 <= T(2) and (T(-1) - k1 < k2 and k2 < (T(1) - k1) / T(3.0)))));

        }

        template<typename T>
        scene::TImagePoint<T> distortion(const scene::TImagePoint<T> &u, T rd,
                                         const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients) {
            T denominator = undistortionDenominator<T>(rd, distortion_coefficients);
            return u * denominator;
        }


        template<typename T>
        scene::TImagePoint<T> undistortion(const scene::TImagePoint<T> &ud,
                                           const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients) {
            T rd = ud.norm();
            T denominator = undistortionDenominator<T>(rd, distortion_coefficients);
            return ud / denominator;
        }

        template<typename T>
        T findDistortedRadius(const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients, T r) {
            Eigen::Matrix<T, Eigen::Dynamic, 1> checked_distortion_coefficients;
            int count_until_non_zero = static_cast<int>(distortion_coefficients.size() - 1);
            while (count_until_non_zero > 0 && ceres::abs(distortion_coefficients[count_until_non_zero]) < T(1e-9)) {
                --count_until_non_zero;
            }
            //std::cout << distortion_coefficients << std::endl;
            //std::cout << count_until_non_zero << std::endl;
            checked_distortion_coefficients = distortion_coefficients.head(count_until_non_zero + 1);
            auto deg = 2 * checked_distortion_coefficients.size();
            Eigen::Matrix<T, Eigen::Dynamic, 1> coeff = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(deg + 1);

            for (auto k = static_cast<size_t>(checked_distortion_coefficients.size()); k > 0; --k)
                coeff(2 * k, 0) = r * checked_distortion_coefficients[k - 1];

            coeff(1, 0) = T(-1);
            coeff(0, 0) = r;
            coeff /= coeff[deg];

            T rd = T(std::numeric_limits<double>::max());


            auto eigenvalues = solvePoly<T>(coeff);
            for (size_t j = 0; j < eigenvalues.rows(); ++j) {
                T real = eigenvalues[j].real();
                T imag = eigenvalues[j].imag();
                if (ceres::abs(imag) < 1e-9 && real > T(1e-9) && rd > real) {
                    rd = real;
                }

            }
            return rd;
        }

        template<typename T>
        scene::TImagePoint<T> distortion(const scene::TImagePoint<T> &u,
                                         const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients) {
            return distortion<T>(u, findDistortedRadius(distortion_coefficients, u.norm()), distortion_coefficients);
        }

        template<typename T>
        T solveQuadric(T distortion_coefficient, T r) {
            if (ceres::abs(distortion_coefficient) < T(EPS) or r < T(1e-9)) {
                return r;

            } else {

                auto tmp = T(4.0) * distortion_coefficient * r * r;

                if (T(1) - tmp < T(0)) {
                    return T(std::numeric_limits<double>::max());
                } else {
                    return (T(1) - ceres::sqrt(T(1.0) - tmp)) /
                           (T(2.0) * distortion_coefficient * r);
                }
            }
        }


        template<template<typename> class ErrorFunctor, typename T>
        void computeErrors(const scene::TImagePoints<T> &u1d,
                           const scene::TImagePoints<T> &u2d,
                           const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients,
                           const scene::TFundamentalMatrix<T> &fundamental_matrix, std::vector<T> &left_residuals,
                           std::vector<T> &right_residuals, T image_r = T(1.0)) {
            assert(u1d.cols() == u2d.cols() && u1d.cols() > 0 &&
                   "Numbers of left and right keypoints should be equal and more than 0");
            auto number_of_points = static_cast<size_t>(u1d.cols());
            left_residuals.resize(number_of_points, T(std::numeric_limits<double>::max()));
            right_residuals.resize(number_of_points, T(std::numeric_limits<double>::max()));
            for (size_t k = 0; k < number_of_points; ++k) {

                ErrorFunctor<T> cost;
                T left_residual, right_residual;
                bool is_correct = cost(u1d.col(k), u2d.col(k), fundamental_matrix,
                                       distortion_coefficients, left_residual, right_residual);
                left_residuals[k] = image_r * left_residual;
                right_residuals[k] = image_r * right_residual;
            }
        }

        template<typename T>
        struct EpipolarCurveDistanceError {

            bool operator()(const scene::TImagePoint<T> &u1d, const scene::TImagePoint<T> &u2d,
                            const scene::TFundamentalMatrix<T> &fundamental_matrix,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1> &distortion_coefficients, T &left_residual,
                            T &right_residual) {

                bool is_correct = true;
                scene::TImagePoint<T> u1, u2;
                u1 = undistortion(u1d, distortion_coefficients);
                u2 = undistortion(u2d, distortion_coefficients);

                scene::THomogenousImagePoint<T> homogeneous_u1, homogeneous_u2;
                Eigen::Matrix<T, 2, 1> line_point1, line_point2;
                homogeneous_u1 = u1.homogeneous();
                homogeneous_u2 = u2.homogeneous();

                Eigen::Matrix<T, 3, 1> l2 = fundamental_matrix * homogeneous_u1;
                Eigen::Matrix<T, 3, 1> l1 = fundamental_matrix.transpose() * homogeneous_u2;

                T n1 = l1.template block<2, 1>(0, 0).norm();
                T n2 = l2.template block<2, 1>(0, 0).norm();
                T err = l1.dot(homogeneous_u1);


                left_residual = err / n1;
                right_residual = err / n2;
                line_point1 = u1 - left_residual * l1.template block<2, 1>(0, 0) / n1;
                line_point2 = u2 - right_residual * l2.template block<2, 1>(0, 0) / n2;


                T root_r1d_estimation, root_r2d_estimation;


                T r1u = line_point1.norm();
                T r2u = line_point2.norm();
                //TODO try to change this to Eigen polynomial solver (need fix for std::Complex with Ceres::Jet)
                //TODO fix when root_r1d_estimation is bad
                root_r1d_estimation = findDistortedRadius(distortion_coefficients, r1u);
                root_r2d_estimation = findDistortedRadius(distortion_coefficients, r2u);

                if (root_r1d_estimation == T(std::numeric_limits<double>::max()) or
                    root_r2d_estimation == T(std::numeric_limits<double>::max())) {
                    left_residual = T(std::numeric_limits<double>::max());
                    right_residual = T(std::numeric_limits<double>::max());
                    is_correct = false;
                    return is_correct;
                } else {

                    scene::TImagePoint<T> curve_point1 = distortion(line_point1, root_r1d_estimation,
                                                                    distortion_coefficients);
                    scene::TImagePoint<T> curve_point2 = distortion(line_point2, root_r2d_estimation,
                                                                    distortion_coefficients);

                    left_residual = (u1d - curve_point1).norm();
                    right_residual = (u2d - curve_point2).norm();

                    if (ceres::IsNaN(left_residual) or ceres::IsNaN(right_residual)) {
                        left_residual = T(std::numeric_limits<double>::max());
                        right_residual = T(std::numeric_limits<double>::max());
                        is_correct = false;
                    }
                    return is_correct;
                }
                /*
                if (distortion_coefficients.rows() > 1) {

                    root_r1d_estimation = findDistortedRadius(distortion_coefficients, r1u);
                    root_r2d_estimation = findDistortedRadius(distortion_coefficients, r2u);
                    if (root_r1d_estimation < T(0))
                        root_r1d_estimation += T(2) * epsilon1;
                    if (root_r2d_estimation < T(0))
                        root_r2d_estimation += T(2) * epsilon2;

                    int newtone_number_of_iterations = 0;
                    while (newtone_number_of_iterations < 200) {
                        newtone_number_of_iterations++;

                        root_r1d_estimation = root_r1d_estimation -
                                              (r1u *
                                               undistortionDenominator(root_r1d_estimation, distortion_coefficients) -
                                               root_r1d_estimation) /
                                              (r1u * undistortionDenominatorDerivative(root_r1d_estimation,
                                                                                       distortion_coefficients) - T(1));

                        root_r2d_estimation = root_r2d_estimation -
                                              (r2u *
                                               undistortionDenominator(root_r2d_estimation, distortion_coefficients) -
                                               root_r2d_estimation) /
                                              (r2u * undistortionDenominatorDerivative(root_r2d_estimation,
                                                                                       distortion_coefficients) - T(1));

                    }
                } else {
                    root_r1d_estimation = solveQuadric(distortion_coefficients(0), r1u);
                    root_r2d_estimation = solveQuadric(distortion_coefficients(0), r2u);
                    std::cout << root_r1d_estimation  << std::endl;
                    std::cout << findDistortedRadius(distortion_coefficients, r1u) << "!!\n" << std::endl;
                    std::cout << root_r2d_estimation  << std::endl;
                    std::cout << findDistortedRadius(distortion_coefficients, r2u) << "!!!\n" << std::endl;
                }

                if (root_r1d_estimation == T(std::numeric_limits<double>::max()) or
                    root_r2d_estimation == T(std::numeric_limits<double>::max())) {
                    left_residual = T(std::numeric_limits<double>::max());
                    right_residual = T(std::numeric_limits<double>::max());
                    is_correct = false;
                    return is_correct;
                } else {

                    scene::TImagePoint<T> curve_point1 = distortion(line_point1, root_r1d_estimation,
                                                                    distortion_coefficients);
                    scene::TImagePoint<T> curve_point2 = distortion(line_point2, root_r2d_estimation,
                                                                    distortion_coefficients);

                    left_residual = (u1d - curve_point1).norm();
                    right_residual = (u2d - curve_point2).norm();

                    if (ceres::IsNaN(left_residual) or ceres::IsNaN(right_residual)) {
                        left_residual = T(std::numeric_limits<double>::max());
                        right_residual = T(std::numeric_limits<double>::max());
                        is_correct = false;
                    }
                    return is_correct;
                }*/
            }
        };
    }
    struct CostFunction {
        Sophus::SE3d leftToRight;
        Eigen::VectorXd distortion_c;
        Eigen::Matrix3d calibration;
        Eigen::Vector2d lft;
        Eigen::Vector2d rht;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        CostFunction(const Sophus::SE3d &l2r, const Eigen::VectorXd &dst, const Eigen::Matrix3d &clb,

                     const Eigen::Vector2d &l, const Eigen::Vector2d &r) : leftToRight(l2r), distortion_c(dst),
                                                                           calibration(clb), lft(l), rht(r) {}

        template<typename T>
        bool operator()(const T *point_ptr, T *residuals) const {
            Eigen::Map<const Eigen::Matrix<T, 3, 1>> point(point_ptr);
            Eigen::Map<Eigen::Matrix<T, 4, 1>> rs(residuals);

            rs.template head<2>() = distortion_problem::distortion<T>( (calibration.template cast<T>() * point).hnormalized(),
                                                                      distortion_c.template cast<T>()) -
                                    lft.template cast<T>();
            rs.template tail<2>() = distortion_problem::distortion<T>(
                    (calibration.template cast<T>() * (leftToRight.template cast<T>() * point)).hnormalized(),
                    distortion_c.template cast<T>()) - rht.template cast<T>();

            return true;
        }
    };

    inline void triangulate(const Eigen::Matrix3d &bifocal_tensor,
                            const Eigen::Matrix3d &rotation_matrix,
                            const Eigen::Matrix3d &translation_matrix,
                            const scene::HomogenousImagePoint &left_keypoint,
                            const scene::HomogenousImagePoint &right_keypoint,
                            scene::HomogenousWorldPoint &left_backprojected,
                            scene::HomogenousWorldPoint &right_backprojected,
                            const Eigen::VectorXd &distortion_coefficients = Eigen::VectorXd(),
                            const Eigen::Matrix3d &calibration = Eigen::Matrix3d::Zero()) {

#if 1
        Eigen::Vector3d t = utils::inverted_screw_hat(translation_matrix);
        Sophus::SO3d so3(rotation_matrix);
        Sophus::SE3d leftToRight = Sophus::SE3d(so3, t);


        Eigen::Vector3d dir_left = (leftToRight.so3() * left_keypoint).normalized(),
                dir_right = right_keypoint.normalized();
        Eigen::Matrix<double, 3, 2> A;
        A.col(0) = dir_left;
        A.col(1) = -dir_right;
        Eigen::Vector3d b = -leftToRight.translation();
        Eigen::Vector2d alphas = A.fullPivHouseholderQr().solve(b);
        Eigen::Vector3d pt = (alphas[0] * dir_left + t + alphas[1] * dir_right) / 2.0;
        right_backprojected = pt.homogeneous();
        left_backprojected = (leftToRight.inverse() * pt).homogeneous();

        ceres::Problem problem;
        problem.AddParameterBlock(left_backprojected.data(), 3);
        problem.AddResidualBlock(new ceres::AutoDiffCostFunction<CostFunction, 4, 3>(
                new CostFunction(leftToRight, distortion_coefficients, calibration,

                                 left_keypoint.hnormalized(), right_keypoint.hnormalized())), nullptr,
                                 left_backprojected.data());
        ceres::Solver::Options options;
        //options.max_trust_region_radius = 0.01;
        options.max_num_iterations = 500;
        options.linear_solver_type = ceres::DENSE_QR;
        options.num_threads = 8;
        options.function_tolerance = 1e-16;
        options.parameter_tolerance = 1e-16;
        options.minimizer_progress_to_stdout = true;
        options.preconditioner_type = ceres::IDENTITY;
        options.jacobi_scaling = false;

        // Solve
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        right_backprojected.template head<3>() = leftToRight * left_backprojected.template head<3>();


#else
        Eigen::Matrix<double, 3, 4> projection_matrix;
        Eigen::Vector3d translation_vector = utils::inverted_screw_hat(translation_matrix);
        projection_matrix << rotation_matrix, translation_vector;
        Eigen::Matrix3d matrixH = Eigen::Matrix3d::Zero();
        matrixH(0, 0) = 1;
        matrixH(1, 1) = 1;

        Eigen::Vector3d a = bifocal_tensor.transpose() * right_keypoint;
        Eigen::Vector3d h1 = matrixH * a;
        Eigen::Vector3d h2 = matrixH * bifocal_tensor * left_keypoint;

        Eigen::Vector3d b = left_keypoint.cross(h1);
        Eigen::Vector3d c = right_keypoint.cross(h2);

        a.normalize();
        b.normalize();
        Eigen::Vector3d d = a.cross(b);
        d.normalize();

        Eigen::Vector4d pC = (projection_matrix.transpose() * c);

        left_backprojected[0] = pC[3] * d[0];
        left_backprojected[1] = pC[3] * d[1];
        left_backprojected[2] = pC[3] * d[2];
        left_backprojected[3] = -d.cwiseProduct(pC.template block<3, 1>(0, 0)).sum();
        Eigen::Vector3d pQ = projection_matrix * left_backprojected / left_backprojected[3];
        right_backprojected[0] = pQ[0];
        right_backprojected[1] = pQ[1];
        right_backprojected[2] = pQ[2];
        right_backprojected[3] = 1;
#endif
    }

    inline void triangulate(const Eigen::Matrix3d &bifocal_tensor,
                            const Eigen::Matrix3d &rotation_matrix,
                            const Eigen::Matrix3d &translation_matrix,
                            const scene::ImagePoint &left_keypoint,
                            const scene::ImagePoint &right_keypoint,
                            scene::WorldPoint &left_backprojected,
                            scene::WorldPoint &right_backprojected,
                            const Eigen::VectorXd &distortion_coefficients = Eigen::VectorXd(),
                            const Eigen::Matrix3d &calibration = Eigen::Matrix3d::Zero()) {
        auto left_backprojected_h = left_backprojected.homogeneous().eval();
        auto right_backprojected_h = right_backprojected.homogeneous().eval();
        std::cout << "FM: " << right_keypoint.transpose().homogeneous() * bifocal_tensor * left_keypoint.homogeneous()
                  << std::endl;
        triangulate(bifocal_tensor, rotation_matrix, translation_matrix, left_keypoint.homogeneous(),
                    right_keypoint.homogeneous(),
                    left_backprojected_h, right_backprojected_h, distortion_coefficients, calibration);
        left_backprojected = left_backprojected_h.hnormalized();
        right_backprojected = right_backprojected_h.hnormalized();
        std::cout << left_backprojected(2) << " " << right_backprojected(2) << std::endl;
        std::cout << "E " << right_backprojected.normalized().transpose() * translation_matrix * rotation_matrix *
                             left_backprojected.normalized() << std::endl;


    }

    inline bool
    chiralityTest(const Eigen::Matrix3d &fundamental_matrix,
                  const Eigen::Matrix3d &rotation_matrix,
                  const Eigen::Matrix3d &translation_matrix,
                  const scene::ImagePoint &left_keypoint, const scene::ImagePoint &right_keypoint,
                  const Eigen::VectorXd &distortion_coefficients = Eigen::VectorXd(),
                  const Eigen::Matrix3d &calibration = Eigen::Matrix3d::Zero()) {
        scene::HomogenousWorldPoint left_backprojected, right_backprojected;

        triangulate(fundamental_matrix, rotation_matrix, translation_matrix, left_keypoint.homogeneous(),
                    right_keypoint.homogeneous(),
                    left_backprojected, right_backprojected, distortion_coefficients, calibration);
        bool c1 = left_backprojected[2] * left_backprojected[3] > 0;
        bool c2 = right_backprojected[2] * right_backprojected[3] > 0;
        return (c1 && c2);
    }
}
#endif //CAMERA_CALIBRATION_UTILITIES_H
