//
// Created by danielbord on 1/26/18.
//
#include <tbb/tbb.h>
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

#if 0
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
#else
            std::vector<double> errors(u1d.cols());
#if 1
            tbb::parallel_for(tbb::blocked_range<int>(0, u1d.cols()),
                              [&](auto range) {
                                  for (int i = range.begin(); i != range.end(); ++i) {
#else
                                      for (int i = 0; i < u1d.cols(); ++i) {
#endif
                                          double err = std::numeric_limits<double>::max();

                                          chiralityTest(fundamental_matrix, rotation_matrix, translation_matrix,
                                                                      undistortion<double>(u1d.col(i), distortion_coefficients),
                                                                      undistortion<double>(u2d.col(i),
                                                                                           distortion_coefficients), distortion_coefficients, calibration, &err) ? 1 : 0;

                                          double px_err = err * image_r;
                                          errors[i] = px_err;
                                      }
#if 1
                                  }
                                  );
#endif
#endif
            double quantile = estimateQuantile(errors, expected_percent_of_inliers);
            double interval = estimateConfidenceInterval(quantile, expected_percent_of_inliers);
            inliers_indices.clear();

            std::atomic<int> invalid_threshold(0), invalid_chirality(0), invalid_reprojection_threshold(0);

            std::vector<int> validity(u1d.cols());
#if 1

            tbb::parallel_for(tbb::blocked_range<int>(0, u1d.cols()),
                              [&](auto range) {
                                  for (int i = range.begin(); i != range.end(); ++i) {
#else
            for (int i = 0; i < u1d.cols(); ++i) {
#endif
                                      if (errors[i] >= interval) {
                                        validity[i] = 0;
                                          ++invalid_threshold;
                                          continue;
                                      }
                                      double err = std::numeric_limits<double>::max();
                                      validity[i] = chiralityTest(fundamental_matrix, rotation_matrix, translation_matrix,
                                                                                          undistortion<double>(u1d.col(i), distortion_coefficients),
                                                                                          undistortion<double>(u2d.col(i),
                                                                                                               distortion_coefficients), distortion_coefficients, calibration, &err) ? 1 : 0;

                                                                 double px_err = err * image_r;
                                      if (!validity[i]) {
                                        ++invalid_chirality;
                                      }

                                      if (px_err > interval * 5.0) {
                                          validity[i] = 0;
                                          ++invalid_reprojection_threshold;
                                          LOG(WARNING) << "Expected error: " << interval << " observed error: " << px_err;
                                          continue;
                                      }
                                      LOG(INFO) << "Reprojection error: " << px_err;
                                  }
#if 1
                              }
            );
#endif

            for (size_t k = 0; k < u1d.cols(); ++k) {
                if (validity[k]) {
                    inliers_indices.push_back(k);
                }
            }
            LOG(INFO) << "Inliers: " << inliers_indices.size() << " [out of: " << u1d.cols() << "]; " << invalid_threshold << "/" << invalid_chirality << "/" << invalid_reprojection_threshold << " thresh/chir/repr";
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
    double points_on_curves(const Eigen::Matrix3d &F, const Eigen::VectorXd &distortion, const Eigen::Vector2d &distorted_left, const Eigen::Vector2d &distorted_right, Eigen::Vector2d &distorted_curve_left, Eigen::Vector2d &distorted_curve_right) {
        double residual_left, residual_right;
        distortion_problem::EpipolarCurveDistanceError2<double> cst;
        cst(distorted_left, distorted_right, F, distortion, distorted_curve_left, distorted_curve_right, residual_left, residual_right);
        return std::pow(residual_left, 2) + std::pow(residual_right, 2);
    }

    Eigen::Vector3d triangulate(const Sophus::SE3d &leftToRight,
                                const Eigen::Vector3d &dirLeft,
                                const Eigen::Vector3d &dirRight) {
    Eigen::Vector3d dir_left = (leftToRight.so3() * dirLeft).normalized(),
            dir_right = dirRight.normalized(),
            t = leftToRight.translation();
    Eigen::Matrix<double, 3, 2> A;
    A.col(0) = dir_left;
    A.col(1) = -dir_right;
    Eigen::Vector3d b = -leftToRight.translation();
    Eigen::Vector2d alphas = A.fullPivHouseholderQr().solve(b);
    return (alphas[0] * dir_left + t + alphas[1] * dir_right) / 2.0;
}

    void triangulate(const Eigen::Matrix3d &bifocal_tensor,
                     const Eigen::Matrix3d &rotation_matrix,
                     const Eigen::Matrix3d &translation_matrix,
                     const scene::HomogenousImagePoint &left_keypoint,
                     const scene::HomogenousImagePoint &right_keypoint,
                     scene::HomogenousWorldPoint &left_backprojected,
                     scene::HomogenousWorldPoint &right_backprojected,
                     const Eigen::VectorXd &distortion_coefficients,
                     const Eigen::Matrix3d &calibration,
    double* reproj_error) {

#if 1



        Eigen::FullPivHouseholderQR<Eigen::Matrix3d> calibrationQR(calibration);

        Eigen::Vector2d distorted_keypoint_left = utils::distortion_problem::distortion<double>(left_keypoint.hnormalized(), distortion_coefficients),
                distorted_keypoint_right = utils::distortion_problem::distortion<double>(right_keypoint.hnormalized(), distortion_coefficients);

        Eigen::Vector2d distorted_curve_left, distorted_curve_right;
        double expected_error = points_on_curves(bifocal_tensor, distortion_coefficients, distorted_keypoint_left, distorted_curve_right, distorted_curve_left, distorted_curve_right);

        // XXX: taking left ray as point on curve corresponding to right's point epipolar line
        Eigen::Vector3d left_ray = calibrationQR.solve(utils::distortion_problem::undistortion<double>(distorted_curve_left, distortion_coefficients).homogeneous()).normalized(),
                right_ray = calibrationQR.solve(right_keypoint).normalized();
        double left_error_initial = (utils::distortion_problem::distortion<double>((calibration * left_ray).hnormalized(), distortion_coefficients) - distorted_keypoint_left).squaredNorm();
        double right_error_initial = (utils::distortion_problem::distortion<double>((calibration * right_ray).hnormalized(), distortion_coefficients) - distorted_keypoint_right).squaredNorm();
        //CHECK_LT(left_error_initial, 2.0*expected_error) << "Left ray error mismatch";
        //CHECK_LT(right_error_initial, 1e-5) << "Right ray error mismatch";



        Eigen::Vector3d t = utils::inverted_screw_hat(translation_matrix);
        Sophus::SO3d so3(rotation_matrix);
        Sophus::SE3d leftToRight = Sophus::SE3d(so3, t);

        Eigen::Matrix3d E_expected = calibration.transpose() * bifocal_tensor * calibration;
        E_expected /= E_expected.norm();
        Eigen::Matrix3d E_observed = translation_matrix * rotation_matrix;
        E_observed /= E_observed.norm();
        double diff_E = std::min((E_expected - E_observed).norm(), (E_expected + E_observed).norm());
        //CHECK_LT(diff_E, 1e-5) << "Pose does not match fundamental matrix";
        ////CHECK_LT(right_ray.transpose() * E_observed * left_ray, 1e-2) << "Rays/translation are not coplanar";


        const Eigen::Vector3d pt_init = triangulate(leftToRight, left_ray, right_ray);
        double left_error_lsq = (utils::distortion_problem::distortion<double>((calibration * (leftToRight.inverse() * pt_init)).hnormalized(), distortion_coefficients) - distorted_keypoint_left).squaredNorm();
        double right_error_lsq = (utils::distortion_problem::distortion<double>((calibration * pt_init).hnormalized(), distortion_coefficients) - distorted_keypoint_right).squaredNorm();
        //CHECK_LT(left_error_lsq, 10.0*expected_error) << "Left point LSQ-solution error mismatch";
        //CHECK_LT(right_error_lsq, 10.0*expected_error) << "Right point LSQ-solution error mismatch";


        right_backprojected = pt_init.homogeneous();
        left_backprojected = (leftToRight.inverse() * pt_init).homogeneous();

        ceres::Problem problem;
        problem.AddParameterBlock(left_backprojected.data(), 3);
        problem.AddResidualBlock(new ceres::AutoDiffCostFunction<CostFunction, 4, 3>(
                new CostFunction(leftToRight, distortion_coefficients, calibration,
                                 distorted_keypoint_left,
                                 distorted_keypoint_right)), nullptr,
                                 left_backprojected.data());
        ceres::Solver::Options options;
        //options.max_trust_region_radius = 0.01;
        options.max_num_iterations = 500;
        options.linear_solver_type = ceres::DENSE_QR;
        options.num_threads = 8;
        options.function_tolerance = 1e-16;
        options.parameter_tolerance = 1e-16;
        options.minimizer_progress_to_stdout = false;
        options.preconditioner_type = ceres::IDENTITY;
        options.jacobi_scaling = false;


        // Solve
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        right_backprojected.template head<3>() = leftToRight * left_backprojected.template head<3>();
        if (reproj_error)
            *reproj_error = std::sqrt(summary.final_cost);

        const Eigen::Vector3d pt_nnls = leftToRight * left_backprojected.head<3>();
        double left_error_nnls = (utils::distortion_problem::distortion<double>((calibration * (leftToRight.inverse() * pt_nnls)).hnormalized(), distortion_coefficients) - distorted_keypoint_left).squaredNorm();
        double right_error_nnls = (utils::distortion_problem::distortion<double>((calibration * pt_nnls).hnormalized(), distortion_coefficients) - distorted_keypoint_right).squaredNorm();
        //CHECK_LT(left_error_nnls, 10.0*expected_error) << "Left point NNLS-solution error mismatch";
        //CHECK_LT(right_error_nnls, 10.0*expected_error) << "Right point NNLS-solution error mismatch";

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
}