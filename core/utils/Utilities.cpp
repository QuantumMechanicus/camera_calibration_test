//
// Created by danielbord on 1/26/18.
//
#include "Utilities.h"

namespace utils {
    namespace distortion_problem {

        double estimateQuantile(std::vector<double> errors,
                                double expected_percent_of_inliers) {

            std::nth_element(errors.begin(), errors.begin() + int(errors.size() * expected_percent_of_inliers),
                             errors.end());
            double quantile = errors[int(errors.size() * expected_percent_of_inliers) + 1];
            return quantile;
        }

        double estimateConfidenceInterval(double quantile, double expected_percent_of_inliers) {
            return quantile * boost::math::erfc_inv((0.95 + 1.0)) /
                   boost::math::erfc_inv((expected_percent_of_inliers + 1.0));

        }

        double findInliers(const scene::ImagePoints &u1d,
                           const scene::ImagePoints &u2d,
                           const Eigen::Matrix<double, Eigen::Dynamic, 1> &distortion_coefficients,
                           const Sophus::SE3d &leftToRight, const Eigen::Matrix3d &calibration,
                           double expected_percent_of_inliers,
                           std::vector<size_t> &inliers_indices, double image_r) {
            Eigen::Matrix3d translation_matrix = screw_hat(
                    leftToRight.translation()), rotation_matrix = leftToRight.so3().matrix();
            Eigen::Matrix3d fundamental_matrix =
                    calibration.inverse().transpose() * translation_matrix * rotation_matrix * calibration.inverse();

            std::vector<double> left_residuals, right_residuals, errors(static_cast<unsigned long>(u1d.cols()));
            utils::distortion_problem::computeErrors<utils::distortion_problem::EpipolarCurveDistanceError, double>(u1d,
                                                                                                                    u2d,
                                                                                                                    distortion_coefficients,
                                                                                                                    fundamental_matrix,
                                                                                                                    left_residuals,
                                                                                                                    right_residuals,
                                                                                                                    image_r);
            for (std::size_t k = 0; k < u1d.cols(); ++k) {
                errors[k] = std::abs(left_residuals[k]) + std::abs(right_residuals[k]);
            }
            double quantile = estimateQuantile(errors, expected_percent_of_inliers);
            double interval = estimateConfidenceInterval(quantile, expected_percent_of_inliers);
            inliers_indices.clear();

            for (size_t k = 0; k < u1d.cols(); ++k) {
                if (errors[k] < interval && chiralityTest(fundamental_matrix, rotation_matrix, translation_matrix,
                                                          undistortion<double>(u1d.col(k), distortion_coefficients),
                                                          undistortion<double>(u2d.col(k), distortion_coefficients))) {
                    inliers_indices.push_back(k);
                }
            }
            return interval;

        }

        double findInliers(const scene::ImagePoints &u1d,
                           const scene::ImagePoints &u2d,
                           const Eigen::Matrix<double, Eigen::Dynamic, 1> &distortion_coefficients,
                           const scene::FundamentalMatrix &fundamental_matrix, double expected_percent_of_inliers,
                           std::vector<size_t> &inliers_indices, double image_r) {
            std::vector<double> left_residuals, right_residuals, errors(static_cast<unsigned long>(u1d.cols()));
            utils::distortion_problem::computeErrors<utils::distortion_problem::EpipolarCurveDistanceError, double>(u1d,
                                                                                                                    u2d,
                                                                                                                    distortion_coefficients,
                                                                                                                    fundamental_matrix,
                                                                                                                    left_residuals,
                                                                                                                    right_residuals,
                                                                                                                    image_r);
            for (std::size_t k = 0; k < u1d.cols(); ++k) {
                errors[k] = std::abs(left_residuals[k]) + std::abs(right_residuals[k]);
            }
            double quantile = estimateQuantile(errors, expected_percent_of_inliers);
            double interval = estimateConfidenceInterval(quantile, expected_percent_of_inliers);
            inliers_indices.resize(0);

            for (size_t k = 0; k < u1d.cols(); ++k) {
                if (errors[k] < interval) {
                    inliers_indices.push_back(k);
                }
            }
            return interval;

        }


    }


}