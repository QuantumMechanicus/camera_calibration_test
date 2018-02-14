//
// Created by danielbord on 1/26/18.
//

#ifndef CAMERA_CALIBRATION_UTILITIES_H
#define CAMERA_CALIBRATION_UTILITIES_H

#include <fstream>
#include <boost/math/special_functions/erf.hpp>
#include <boost/filesystem.hpp>
#include "ceres/ceres.h"
#include "../scene/Two_View.h"

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

        if (m.cols() != n_cols or m.rows() != n_rows) {
            std::cerr << "Invalid matrix size\n";
            m.setZero();
            return false;
        }
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

    namespace distortion_problem {

        double estimateQuantile(std::vector<double> &errors,
                                double expected_percent_of_inliers);

        double estimateConfidenceInterval(double quantile, double expected_percent_of_inliers);


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
                           std::vector<T> &right_residuals) {
            assert(u1d.cols() == u2d.cols() && "Numbers of left and right keypoints should bee equal");
            auto number_of_points = static_cast<size_t>(u1d.cols());
            left_residuals.resize(number_of_points, T(std::numeric_limits<double>::max()));
            right_residuals.resize(number_of_points, T(std::numeric_limits<double>::max()));
            for (size_t k = 0; k < number_of_points; ++k) {

                ErrorFunctor<T> cost;
                T left_residual, right_residual;
                bool is_correct = cost(u1d.col(k), u2d.col(k), fundamental_matrix,
                                       distortion_coefficients, left_residual, right_residual);
                left_residuals[k] = left_residual;
                right_residuals[k] = right_residual;
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

                Eigen::Matrix<T, 3, 1> l1 = fundamental_matrix * homogeneous_u1;
                Eigen::Matrix<T, 3, 1> l2 = fundamental_matrix.transpose() * homogeneous_u2;

                T n1 = l1.template block<2, 1>(0, 0).norm();
                T n2 = l2.template block<2, 1>(0, 0).norm();
                T err = l1.dot(homogeneous_u2);


                left_residual = err / n1;
                right_residual = err / n2;
                line_point1 = u1 + left_residual * l1.template block<2, 1>(0, 0) / n1;
                line_point2 = u2 + right_residual * l1.template block<2, 1>(0, 0) / n1;


                T root_r1d_estimation, root_r2d_estimation;

                T epsilon1, epsilon2;
                epsilon1 = std::min(ceres::abs(left_residual), T(0.5));
                epsilon2 = std::min(ceres::abs(right_residual), T(0.5));


                T r1u = line_point1.norm();
                T r2u = line_point2.norm();

                if (distortion_coefficients.rows() > 1) {

                    root_r1d_estimation = u1d.norm() - epsilon1;
                    root_r2d_estimation = u2d.norm() - epsilon2;
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

                    if (std::isnan(left_residual) or std::isnan(right_residual)) {
                        left_residual = T(std::numeric_limits<double>::max());
                        right_residual = T(std::numeric_limits<double>::max());
                        is_correct = false;
                    }
                    return is_correct;
                }
            }
        };
    }
}
#endif //CAMERA_CALIBRATION_UTILITIES_H
