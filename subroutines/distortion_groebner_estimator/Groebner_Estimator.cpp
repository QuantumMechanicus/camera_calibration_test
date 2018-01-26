//
// Created by danielbord on 1/23/18.
//
#include "Groebner_Estimator.h"

namespace estimators {
    GroebnerEstimatorOptions::GroebnerEstimatorOptions(int number_of_RANSAC_iterations, double quantile_to_minimize,
                                                       double lambda_lower_bound, double lambda_upper_bound,
                                                       double image_radius) : number_of_RANSAC_iterations(
            number_of_RANSAC_iterations), quantile_to_minimize(quantile_to_minimize),
                                                                              lambda_lower_bound(lambda_lower_bound),
                                                                              lambda_upper_bound(lambda_upper_bound),
                                                                              image_radius(image_radius) {}

    GroebnerDivisionModelEstimator::GroebnerDivisionModelEstimator(Eigen::Matrix<double, 2, -1> u1d,
                                                                   Eigen::Matrix<double, 2, -1> u2d, GroebnerEstimatorOptions opt)
            : u1d_(std::move(u1d)), u2d_(std::move(u2d)) {
        assert(u1d_.cols() == u2d_.cols() && "Numbers of left and right keypoints should bee equal");
        assert(u1d_.cols() > 7 && "Numbers of keypoints must be greater than 7");

        number_of_points_ = static_cast<size_t>(u1d_.cols());
        left_residuals_.resize(number_of_points_, 0);
        right_residuals_.resize(number_of_points_, 0);
        is_estimated_ = false;
        options_ = opt;
    }

    void GroebnerDivisionModelEstimator::setOptions(const GroebnerEstimatorOptions &options_) {
        is_estimated_ = false;
        GroebnerDivisionModelEstimator::options_ = options_;
    }

    bool GroebnerDivisionModelEstimator::isEstimated() const
    {
        return is_estimated_;
    }

    void GroebnerDivisionModelEstimator::estimate() {
        is_estimated_ = true;
        ppx_ = 0;
        ppy_ = 0;
        f_ = 0;
        fundamental_matrix_.setZero();
        lambdas_.setZero();
        double min_quantile = std::numeric_limits<double>::max();

        auto number_of_associated_points = static_cast<unsigned long>(u1d_.cols());
        std::random_device rd;
        std::mt19937 gen(rd());
        std::vector<std::size_t> numbers(number_of_associated_points);
        std::iota(numbers.begin(), numbers.end(), 0);

        AutomaticSolver groebner_automatic_solver;
        for (std::size_t k = 0; k < options_.number_of_RANSAC_iterations; ++k) {
            std::shuffle(numbers.begin(), numbers.end(), gen);
            AutomaticSolver::EightPoints subset_Q1, subset_Q2;
            subset_Q1 << u1d_.col(numbers[0]), u1d_.col(numbers[1]), u1d_.col(numbers[2]), u1d_.col(
                    numbers[3]), u1d_.col(numbers[4]), u1d_.col(numbers[5]), u1d_.col(numbers[6]), u1d_.col(
                    numbers[7]);
            subset_Q2 << u2d_.col(numbers[0]), u2d_.col(numbers[1]), u2d_.col(numbers[2]), u2d_.col(
                    numbers[3]), u2d_.col(numbers[4]), u2d_.col(numbers[5]), u2d_.col(numbers[6]), u2d_.col(
                    numbers[7]);


            AutomaticSolver::FundamentalMatricesAndDistortionCoefficients models = groebner_automatic_solver.runSolver(
                    subset_Q1, subset_Q2);

            unsigned long count = models.first.size();
            for (std::size_t i = 0; i < count; ++i) {
                if (models.second[i] > options_.lambda_upper_bound or
                    models.second[i] < options_.lambda_lower_bound or
                    models.first[i].determinant() > 1)
                    continue;

                double hyp_lambda = models.second[i];
                Eigen::Matrix3d hyp_fundamental_matrix = models.first[i];

                Eigen::JacobiSVD<Eigen::Matrix3d> fmatrix_svd(hyp_fundamental_matrix,
                                                              Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Vector3d singular_values = fmatrix_svd.singularValues();
                singular_values[2] = 0.0;

                hyp_fundamental_matrix = fmatrix_svd.matrixU() * singular_values.asDiagonal() *
                                         fmatrix_svd.matrixV().transpose();
                hyp_fundamental_matrix = hyp_fundamental_matrix / hyp_fundamental_matrix(2, 2);


                double quantile = computeErrorsAndEstimateQuantile(hyp_lambda, hyp_fundamental_matrix);

                if (quantile < min_quantile) {

                    min_quantile = quantile;
                    lambdas_(0) = hyp_lambda;
                    fundamental_matrix_ = hyp_fundamental_matrix;

                }

            }
        }
    }

    double GroebnerDivisionModelEstimator::computeErrorsAndEstimateQuantile(double distortion_coefficient,
                                                                            const scene::FundamentalMatrix &fundamental_matrix) {
        Eigen::Matrix<double, 1, 1> distortion_coefficients;
        std::vector<double> errors(number_of_points_);
        distortion_coefficients(0) = distortion_coefficient;
        utils::distortion_problem::computeErrors<utils::distortion_problem::EpipolarCurveDistanceError, double>(u1d_, u2d_, distortion_coefficients, fundamental_matrix, left_residuals_, right_residuals_);
        for (std::size_t k = 0; k < number_of_points_; ++k)
        {
            errors[k] = std::abs(left_residuals_[k]) + std::abs(right_residuals_[k]);
        }
        return utils::distortion_problem::estimateQuantile(errors, options_.quantile_to_minimize);
    }


    GroebnerDivisionModelEstimator::AutomaticSolver::FundamentalMatricesAndDistortionCoefficients
    GroebnerDivisionModelEstimator::AutomaticSolver::solver(const GPolynomial &g1, const GPolynomial &g2,
                                                            const GPolynomial &g3, const GPolynomial &g4,
                                                            const GPolynomial &g5, const GPolynomial &g6,
                                                            const GPolynomial &g7, const GPolynomial &g8) {
        Eigen::Matrix<double, 59, 1> c;
        c(0) = -g6(2);
        c(1) = g5(2) - g6(5);
        c(2) = -g6(1);
        c(3) = -g6(0);
        c(4) = g5(1) - g6(4);
        c(5) = g5(0) - g6(3);
        c(6) = g5(5) - g6(6);
        c(7) = g5(4);
        c(8) = g5(3);
        c(9) = g5(6);
        c(10) = -g8(2);
        c(11) = g7(2) - g8(5);
        c(12) = -g8(1);
        c(13) = -g8(0);
        c(14) = g7(1) - g8(4);
        c(15) = g7(0) - g8(3);
        c(16) = g7(5) - g8(6);
        c(17) = g7(4);
        c(18) = g7(3);
        c(19) = g7(6);

        c(20) = g3(0) * g6(2) - g1(2) * g8(0) - g1(0) * g8(2) + g3(2) * g6(0) + g2(1) * g8(2) + g2(2) * g8(1) -
                g4(1) * g6(2) - g4(2) * g6(1);
        c(21) = g3(0) * g6(3) - g1(3) * g8(0) - g1(0) * g8(3) + g3(3) * g6(0) + g2(0) * g8(4) + g2(1) * g8(3) +
                g2(3) * g8(1) + g2(4) * g8(0) - g4(0) * g6(4) - g4(1) * g6(3) - g4(3) * g6(1) - g4(4) * g6(0);
        c(22) = g3(0) * g6(4) - g1(1) * g8(3) - g1(3) * g8(1) - g1(4) * g8(0) - g1(0) * g8(4) + g3(1) * g6(3) +
                g3(3) * g6(1) + g3(4) * g6(0) + g2(1) * g8(4) + g2(4) * g8(1) - g4(1) * g6(4) - g4(4) * g6(1);
        c(23) = g1(0) * g4(1) + g1(1) * g4(0) - g2(0) * g3(1) - g2(1) * g3(0) - g1(0) * g8(5) - g1(2) * g8(3) -
                g1(3) * g8(2) - g1(5) * g8(0) + g3(0) * g6(5) + g3(2) * g6(3) + g3(3) * g6(2) + g3(5) * g6(0) +
                g2(1) * g8(5) + g2(2) * g8(4) + g2(4) * g8(2) + g2(5) * g8(1) - g4(1) * g6(5) - g4(2) * g6(4) -
                g4(4) * g6(2) - g4(5) * g6(1);
        c(24) = g1(0) * g4(4) + g1(1) * g4(3) + g1(3) * g4(1) + g1(4) * g4(0) - g2(0) * g3(4) - g2(1) * g3(3) -
                g2(3) * g3(1) - g2(4) * g3(0) - g1(0) * g8(6) - g1(6) * g8(0) + g3(0) * g6(6) + g3(6) * g6(0) -
                g1(3) * g8(5) - g1(5) * g8(3) + g3(3) * g6(5) + g3(5) * g6(3) + g2(1) * g8(6) + g2(6) * g8(1) -
                g4(1) * g6(6) - g4(6) * g6(1) + g2(4) * g8(5) + g2(5) * g8(4) - g4(4) * g6(5) - g4(5) * g6(4);
        c(25) = g3(0) * g6(0) - g1(0) * g8(0) + g2(0) * g8(1) + g2(1) * g8(0) - g4(0) * g6(1) - g4(1) * g6(0);
        c(26) = g3(0) * g6(1) - g1(1) * g8(0) - g1(0) * g8(1) + g3(1) * g6(0) + g2(1) * g8(1) - g4(1) * g6(1);
        c(27) = g3(1) * g6(1) - g1(1) * g8(1);
        c(28) = g3(1) * g6(2) - g1(2) * g8(1) - g1(1) * g8(2) + g3(2) * g6(1);
        c(29) = g1(1) * g4(1) - g2(1) * g3(1) - g1(1) * g8(5) - g1(2) * g8(4) - g1(4) * g8(2) - g1(5) * g8(1) +
                g3(1) * g6(5) + g3(2) * g6(4) + g3(4) * g6(2) + g3(5) * g6(1);
        c(30) = g2(0) * g8(0) - g4(0) * g6(0);
        c(31) = g2(0) * g8(2) + g2(2) * g8(0) - g4(0) * g6(2) - g4(2) * g6(0);
        c(32) = g1(0) * g4(0) - g2(0) * g3(0) + g2(0) * g8(5) + g2(2) * g8(3) + g2(3) * g8(2) + g2(5) * g8(0) -
                g4(0) * g6(5) - g4(2) * g6(3) - g4(3) * g6(2) - g4(5) * g6(0);
        c(33) = g1(2) * g4(2) - g2(2) * g3(2);
        c(34) = g1(2) * g4(5) + g1(5) * g4(2) - g2(2) * g3(5) - g2(5) * g3(2);
        c(35) = g1(2) * g4(6) + g1(6) * g4(2) - g2(2) * g3(6) - g2(6) * g3(2) + g1(5) * g4(5) - g2(5) * g3(5);
        c(36) = g3(4) * g6(4) - g1(4) * g8(4);
        c(37) = g1(4) * g4(4) - g2(4) * g3(4) - g1(4) * g8(6) - g1(6) * g8(4) + g3(4) * g6(6) + g3(6) * g6(4);
        c(38) = g2(3) * g8(3) - g4(3) * g6(3);
        c(39) = g1(3) * g4(3) - g2(3) * g3(3) + g2(3) * g8(6) + g2(6) * g8(3) - g4(3) * g6(6) - g4(6) * g6(3);
        c(40) = g3(2) * g6(2) - g1(2) * g8(2);
        c(41) = g3(1) * g6(4) - g1(4) * g8(1) - g1(1) * g8(4) + g3(4) * g6(1);
        c(42) = g1(1) * g4(2) + g1(2) * g4(1) - g2(1) * g3(2) - g2(2) * g3(1) - g1(2) * g8(5) - g1(5) * g8(2) +
                g3(2) * g6(5) + g3(5) * g6(2);
        c(43) = g2(2) * g8(2) - g4(2) * g6(2);
        c(44) = g1(1) * g4(4) + g1(4) * g4(1) - g2(1) * g3(4) - g2(4) * g3(1) - g1(1) * g8(6) - g1(6) * g8(1) +
                g3(1) * g6(6) + g3(6) * g6(1) - g1(4) * g8(5) - g1(5) * g8(4) + g3(4) * g6(5) + g3(5) * g6(4);
        c(45) = g2(0) * g8(3) + g2(3) * g8(0) - g4(0) * g6(3) - g4(3) * g6(0);
        c(46) = g1(1) * g4(5) + g1(2) * g4(4) + g1(4) * g4(2) + g1(5) * g4(1) - g2(1) * g3(5) - g2(2) * g3(4) -
                g2(4) * g3(2) - g2(5) * g3(1) - g1(2) * g8(6) - g1(6) * g8(2) + g3(2) * g6(6) + g3(6) * g6(2) -
                g1(5) * g8(5) + g3(5) * g6(5);
        c(47) = g1(0) * g4(2) + g1(2) * g4(0) - g2(0) * g3(2) - g2(2) * g3(0) + g2(2) * g8(5) + g2(5) * g8(2) -
                g4(2) * g6(5) - g4(5) * g6(2);
        c(48) = g1(0) * g4(3) + g1(3) * g4(0) - g2(0) * g3(3) - g2(3) * g3(0) + g2(0) * g8(6) + g2(6) * g8(0) -
                g4(0) * g6(6) - g4(6) * g6(0) + g2(3) * g8(5) + g2(5) * g8(3) - g4(3) * g6(5) - g4(5) * g6(3);
        c(49) = g1(0) * g4(5) + g1(2) * g4(3) + g1(3) * g4(2) + g1(5) * g4(0) - g2(0) * g3(5) - g2(2) * g3(3) -
                g2(3) * g3(2) - g2(5) * g3(0) + g2(2) * g8(6) + g2(6) * g8(2) - g4(2) * g6(6) - g4(6) * g6(2) +
                g2(5) * g8(5) - g4(5) * g6(5);
        c(50) = g3(3) * g6(3) - g1(3) * g8(3) + g2(3) * g8(4) + g2(4) * g8(3) - g4(3) * g6(4) - g4(4) * g6(3);
        c(51) = g3(3) * g6(4) - g1(4) * g8(3) - g1(3) * g8(4) + g3(4) * g6(3) + g2(4) * g8(4) - g4(4) * g6(4);
        c(52) = g1(1) * g4(6) + g1(6) * g4(1) - g2(1) * g3(6) - g2(6) * g3(1) + g1(4) * g4(5) + g1(5) * g4(4) -
                g2(4) * g3(5) - g2(5) * g3(4) - g1(5) * g8(6) - g1(6) * g8(5) + g3(5) * g6(6) + g3(6) * g6(5);
        c(53) = g1(0) * g4(6) + g1(6) * g4(0) - g2(0) * g3(6) - g2(6) * g3(0) + g1(3) * g4(5) + g1(5) * g4(3) -
                g2(3) * g3(5) - g2(5) * g3(3) + g2(5) * g8(6) + g2(6) * g8(5) - g4(5) * g6(6) - g4(6) * g6(5);
        c(54) = g1(3) * g4(4) + g1(4) * g4(3) - g2(3) * g3(4) - g2(4) * g3(3) - g1(3) * g8(6) - g1(6) * g8(3) +
                g3(3) * g6(6) + g3(6) * g6(3) + g2(4) * g8(6) + g2(6) * g8(4) - g4(4) * g6(6) - g4(6) * g6(4);
        c(55) = g1(5) * g4(6) + g1(6) * g4(5) - g2(5) * g3(6) - g2(6) * g3(5);
        c(56) = g1(4) * g4(6) + g1(6) * g4(4) - g2(4) * g3(6) - g2(6) * g3(4) - g1(6) * g8(6) + g3(6) * g6(6);
        c(57) = g1(3) * g4(6) + g1(6) * g4(3) - g2(3) * g3(6) - g2(6) * g3(3) + g2(6) * g8(6) - g4(6) * g6(6);
        c(58) = g1(6) * g4(6) - g2(6) * g3(6);


        Eigen::Matrix<double, 32, 48> M = Eigen::Matrix<double, 32, 48>::Zero();

        M(0, 4) = c[0];
        M(1, 5) = c[0];
        M(2, 6) = c[0];
        M(3, 8) = c[0];
        M(11, 13) = c[0];
        M(12, 14) = c[0];
        M(13, 15) = c[0];
        M(14, 16) = c[0];
        M(15, 17) = c[0];
        M(23, 25) = c[0];
        M(24, 26) = c[0];
        M(25, 27) = c[0];
        M(29, 37) = c[0];
        M(0, 13) = c[1];
        M(1, 14) = c[1];
        M(2, 15) = c[1];
        M(3, 17) = c[1];
        M(11, 22) = c[1];
        M(12, 23) = c[1];
        M(13, 24) = c[1];
        M(14, 25) = c[1];
        M(15, 26) = c[1];
        M(16, 27) = c[1];
        M(23, 30) = c[1];
        M(24, 31) = c[1];
        M(25, 37) = c[1];
        M(29, 43) = c[1];
        M(0, 1) = c[2];
        M(1, 2) = c[2];
        M(2, 3) = c[2];
        M(3, 6) = c[2];
        M(4, 8) = c[2];
        M(11, 10) = c[2];
        M(12, 11) = c[2];
        M(13, 12) = c[2];
        M(14, 14) = c[2];
        M(15, 15) = c[2];
        M(16, 17) = c[2];
        M(23, 23) = c[2];
        M(24, 24) = c[2];
        M(25, 26) = c[2];
        M(29, 31) = c[2];
        M(0, 0) = c[3];
        M(1, 1) = c[3];
        M(2, 2) = c[3];
        M(3, 5) = c[3];
        M(4, 7) = c[3];
        M(11, 9) = c[3];
        M(12, 10) = c[3];
        M(13, 11) = c[3];
        M(14, 13) = c[3];
        M(15, 14) = c[3];
        M(16, 16) = c[3];
        M(23, 22) = c[3];
        M(24, 23) = c[3];
        M(25, 25) = c[3];
        M(29, 30) = c[3];
        M(0, 10) = c[4];
        M(1, 11) = c[4];
        M(2, 12) = c[4];
        M(3, 15) = c[4];
        M(4, 17) = c[4];
        M(11, 19) = c[4];
        M(12, 20) = c[4];
        M(13, 21) = c[4];
        M(14, 23) = c[4];
        M(15, 24) = c[4];
        M(16, 26) = c[4];
        M(23, 29) = c[4];
        M(24, 36) = c[4];
        M(25, 31) = c[4];
        M(29, 42) = c[4];
        M(0, 9) = c[5];
        M(1, 10) = c[5];
        M(2, 11) = c[5];
        M(3, 14) = c[5];
        M(4, 16) = c[5];
        M(11, 18) = c[5];
        M(12, 19) = c[5];
        M(13, 20) = c[5];
        M(14, 22) = c[5];
        M(15, 23) = c[5];
        M(16, 25) = c[5];
        M(23, 28) = c[5];
        M(24, 29) = c[5];
        M(25, 30) = c[5];
        M(29, 41) = c[5];
        M(0, 22) = c[6];
        M(1, 23) = c[6];
        M(2, 24) = c[6];
        M(3, 26) = c[6];
        M(4, 27) = c[6];
        M(11, 28) = c[6];
        M(12, 29) = c[6];
        M(13, 36) = c[6];
        M(14, 30) = c[6];
        M(15, 31) = c[6];
        M(16, 37) = c[6];
        M(23, 41) = c[6];
        M(24, 42) = c[6];
        M(25, 43) = c[6];
        M(29, 46) = c[6];
        M(0, 19) = c[7];
        M(1, 20) = c[7];
        M(2, 21) = c[7];
        M(3, 24) = c[7];
        M(4, 26) = c[7];
        M(11, 33) = c[7];
        M(12, 34) = c[7];
        M(13, 35) = c[7];
        M(14, 29) = c[7];
        M(15, 36) = c[7];
        M(16, 31) = c[7];
        M(23, 39) = c[7];
        M(24, 40) = c[7];
        M(25, 42) = c[7];
        M(29, 45) = c[7];
        M(0, 18) = c[8];
        M(1, 19) = c[8];
        M(2, 20) = c[8];
        M(3, 23) = c[8];
        M(4, 25) = c[8];
        M(11, 32) = c[8];
        M(12, 33) = c[8];
        M(13, 34) = c[8];
        M(14, 28) = c[8];
        M(15, 29) = c[8];
        M(16, 30) = c[8];
        M(23, 38) = c[8];
        M(24, 39) = c[8];
        M(25, 41) = c[8];
        M(29, 44) = c[8];
        M(0, 28) = c[9];
        M(1, 29) = c[9];
        M(2, 36) = c[9];
        M(3, 31) = c[9];
        M(4, 37) = c[9];
        M(11, 38) = c[9];
        M(12, 39) = c[9];
        M(13, 40) = c[9];
        M(14, 41) = c[9];
        M(15, 42) = c[9];
        M(16, 43) = c[9];
        M(23, 44) = c[9];
        M(24, 45) = c[9];
        M(25, 46) = c[9];
        M(29, 47) = c[9];
        M(5, 4) = c[10];
        M(6, 5) = c[10];
        M(7, 6) = c[10];
        M(8, 7) = c[10];
        M(9, 8) = c[10];
        M(17, 13) = c[10];
        M(18, 14) = c[10];
        M(19, 15) = c[10];
        M(20, 16) = c[10];
        M(21, 17) = c[10];
        M(26, 25) = c[10];
        M(27, 26) = c[10];
        M(28, 27) = c[10];
        M(30, 37) = c[10];
        M(5, 13) = c[11];
        M(6, 14) = c[11];
        M(7, 15) = c[11];
        M(8, 16) = c[11];
        M(9, 17) = c[11];
        M(17, 22) = c[11];
        M(18, 23) = c[11];
        M(19, 24) = c[11];
        M(20, 25) = c[11];
        M(21, 26) = c[11];
        M(22, 27) = c[11];
        M(26, 30) = c[11];
        M(27, 31) = c[11];
        M(28, 37) = c[11];
        M(30, 43) = c[11];
        M(5, 1) = c[12];
        M(6, 2) = c[12];
        M(7, 3) = c[12];
        M(8, 5) = c[12];
        M(9, 6) = c[12];
        M(10, 8) = c[12];
        M(17, 10) = c[12];
        M(18, 11) = c[12];
        M(19, 12) = c[12];
        M(20, 14) = c[12];
        M(21, 15) = c[12];
        M(22, 17) = c[12];
        M(26, 23) = c[12];
        M(27, 24) = c[12];
        M(28, 26) = c[12];
        M(30, 31) = c[12];
        M(5, 0) = c[13];
        M(6, 1) = c[13];
        M(7, 2) = c[13];
        M(8, 4) = c[13];
        M(9, 5) = c[13];
        M(10, 7) = c[13];
        M(17, 9) = c[13];
        M(18, 10) = c[13];
        M(19, 11) = c[13];
        M(20, 13) = c[13];
        M(21, 14) = c[13];
        M(22, 16) = c[13];
        M(26, 22) = c[13];
        M(27, 23) = c[13];
        M(28, 25) = c[13];
        M(30, 30) = c[13];
        M(5, 10) = c[14];
        M(6, 11) = c[14];
        M(7, 12) = c[14];
        M(8, 14) = c[14];
        M(9, 15) = c[14];
        M(10, 17) = c[14];
        M(17, 19) = c[14];
        M(18, 20) = c[14];
        M(19, 21) = c[14];
        M(20, 23) = c[14];
        M(21, 24) = c[14];
        M(22, 26) = c[14];
        M(26, 29) = c[14];
        M(27, 36) = c[14];
        M(28, 31) = c[14];
        M(30, 42) = c[14];
        M(5, 9) = c[15];
        M(6, 10) = c[15];
        M(7, 11) = c[15];
        M(8, 13) = c[15];
        M(9, 14) = c[15];
        M(10, 16) = c[15];
        M(17, 18) = c[15];
        M(18, 19) = c[15];
        M(19, 20) = c[15];
        M(20, 22) = c[15];
        M(21, 23) = c[15];
        M(22, 25) = c[15];
        M(26, 28) = c[15];
        M(27, 29) = c[15];
        M(28, 30) = c[15];
        M(30, 41) = c[15];
        M(5, 22) = c[16];
        M(6, 23) = c[16];
        M(7, 24) = c[16];
        M(8, 25) = c[16];
        M(9, 26) = c[16];
        M(10, 27) = c[16];
        M(17, 28) = c[16];
        M(18, 29) = c[16];
        M(19, 36) = c[16];
        M(20, 30) = c[16];
        M(21, 31) = c[16];
        M(22, 37) = c[16];
        M(26, 41) = c[16];
        M(27, 42) = c[16];
        M(28, 43) = c[16];
        M(30, 46) = c[16];
        M(5, 19) = c[17];
        M(6, 20) = c[17];
        M(7, 21) = c[17];
        M(8, 23) = c[17];
        M(9, 24) = c[17];
        M(10, 26) = c[17];
        M(17, 33) = c[17];
        M(18, 34) = c[17];
        M(19, 35) = c[17];
        M(20, 29) = c[17];
        M(21, 36) = c[17];
        M(22, 31) = c[17];
        M(26, 39) = c[17];
        M(27, 40) = c[17];
        M(28, 42) = c[17];
        M(30, 45) = c[17];
        M(5, 18) = c[18];
        M(6, 19) = c[18];
        M(7, 20) = c[18];
        M(8, 22) = c[18];
        M(9, 23) = c[18];
        M(10, 25) = c[18];
        M(17, 32) = c[18];
        M(18, 33) = c[18];
        M(19, 34) = c[18];
        M(20, 28) = c[18];
        M(21, 29) = c[18];
        M(22, 30) = c[18];
        M(26, 38) = c[18];
        M(27, 39) = c[18];
        M(28, 41) = c[18];
        M(30, 44) = c[18];
        M(5, 28) = c[19];
        M(6, 29) = c[19];
        M(7, 36) = c[19];
        M(8, 30) = c[19];
        M(9, 31) = c[19];
        M(10, 37) = c[19];
        M(17, 38) = c[19];
        M(18, 39) = c[19];
        M(19, 40) = c[19];
        M(20, 41) = c[19];
        M(21, 42) = c[19];
        M(22, 43) = c[19];
        M(26, 44) = c[19];
        M(27, 45) = c[19];
        M(28, 46) = c[19];
        M(30, 47) = c[19];
        M(31, 14) = c[20];
        M(31, 19) = c[21];
        M(31, 20) = c[22];
        M(31, 23) = c[23];
        M(31, 29) = c[24];
        M(31, 10) = c[25];
        M(31, 11) = c[26];
        M(31, 12) = c[27];
        M(31, 15) = c[28];
        M(31, 24) = c[29];
        M(31, 9) = c[30];
        M(31, 13) = c[31];
        M(31, 22) = c[32];
        M(31, 27) = c[33];
        M(31, 37) = c[34];
        M(31, 43) = c[35];
        M(31, 35) = c[36];
        M(31, 40) = c[37];
        M(31, 32) = c[38];
        M(31, 38) = c[39];
        M(31, 17) = c[40];
        M(31, 21) = c[41];
        M(31, 26) = c[42];
        M(31, 16) = c[43];
        M(31, 36) = c[44];
        M(31, 18) = c[45];
        M(31, 31) = c[46];
        M(31, 25) = c[47];
        M(31, 28) = c[48];
        M(31, 30) = c[49];
        M(31, 33) = c[50];
        M(31, 34) = c[51];
        M(31, 42) = c[52];
        M(31, 41) = c[53];
        M(31, 39) = c[54];
        M(31, 46) = c[55];
        M(31, 45) = c[56];
        M(31, 44) = c[57];
        M(31, 47) = c[58];


        Eigen::FullPivHouseholderQR<Eigen::Matrix<long double, 32, 32>> qr(
                M.template block<32, 32>(0, 0).cast<long double>());

        Eigen::Matrix<long double, 32, 16> Mres = qr.solve(M.template block<32, 16>(0, 32).cast<long double>());

        int mrows[10] = {32, 31, 30, 29, 28, 25, 22, 21, 20, 19};
        int arows[10] = {6, 7, 9, 10, 11, 12, 13, 14, 15, 16};

        Eigen::Matrix<double, 16, 16> A;
        A.setZero();
        A(0, 1) = 1;
        A(1, 4) = 1;
        A(2, 5) = 1;
        A(3, 6) = 1;
        A(4, 10) = 1;
        A(7, 11) = 1;
        for (size_t i = 0; i < 10; ++i)
            for (size_t j = 0; j < 16; ++j) {
                A(arows[i] - 1, j) = -Mres(mrows[i] - 1, 15 - j);

            }

        Eigen::EigenSolver<Eigen::Matrix<double, 16, 16> > eigen_solver(A);


        auto V = eigen_solver.eigenvectors();
        auto sol = (V.template block<3, 16>(1, 0) *
                    Eigen::DiagonalMatrix<std::complex<double>, 16>(
                            V.template block<1, 16>(0, 0).cwiseInverse())).eval();
        sol.row(0).swap(sol.row(2));
        FundamentalMatricesAndDistortionCoefficients res;
        for (size_t k = 0; k < sol.cols(); ++k) {
            if (std::abs(sol(0, k).imag()) < utils::EPS) {
                double lambda = sol(2, k).real();

                double f31 = sol(0, k).real();
                double f32 = sol(1, k).real();


                Eigen::Matrix<double, 7, 1> mon;
                mon(0) = f31 * lambda;
                mon(1) = f32 * lambda;
                mon(2) = lambda * lambda;
                mon(3) = f31;
                mon(4) = f32;
                mon(5) = lambda;
                mon(6) = 1;
                double f11 = -mon.transpose() * g1;
                double f12 = -mon.transpose() * g2;
                double f13 = -mon.transpose() * g6;
                double f21 = -mon.transpose() * g3;
                double f22 = -mon.transpose() * g4;
                double f23 = -mon.transpose() * g8;


                Eigen::Matrix3d resF;
                resF.setZero();
                resF(0, 0) = f11;
                resF(0, 1) = f12;
                resF(0, 2) = f13;
                resF(1, 0) = f21;
                resF(1, 1) = f22;
                resF(1, 2) = f23;
                resF(2, 0) = f31;
                resF(2, 1) = f32;
                resF(2, 2) = 1;
                res.second.push_back(lambda);
                res.first.emplace_back(resF.transpose());

            }
        }

        return res;
    }

    GroebnerDivisionModelEstimator::AutomaticSolver::FundamentalMatricesAndDistortionCoefficients
    GroebnerDivisionModelEstimator::AutomaticSolver::runSolver(
            EightPoints u1d, EightPoints u2d) {
        Eigen::Matrix<double, 15, 8> C;
        Eigen::Matrix<double, 8, 15> Cfm;

        C.row(0) = u1d.row(0).cwiseProduct(u2d.row(0));
        C.row(1) = u1d.row(0).cwiseProduct(u2d.row(1));
        C.row(2) = u1d.row(1).cwiseProduct(u2d.row(0));
        C.row(3) = u1d.row(1).cwiseProduct(u2d.row(1));
        C.row(4) = u1d.row(0).cwiseProduct(u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1)));
        C.row(5) = u1d.row(0);
        C.row(6) = u1d.row(1).cwiseProduct(u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1)));
        C.row(7) = u1d.row(1);
        C.row(8) = u2d.row(0).cwiseProduct(u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1)));
        C.row(9) = u2d.row(1).cwiseProduct(u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1)));
        C.row(10) = (u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1))).cwiseProduct(
                u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1)));
        C.row(11) = u2d.row(0);
        C.row(12) = u2d.row(1);
        C.row(13) = (u2d.row(0).cwiseProduct(u2d.row(0)) + u2d.row(1).cwiseProduct(u2d.row(1))) +
                    (u1d.row(0).cwiseProduct(u1d.row(0)) + u1d.row(1).cwiseProduct(u1d.row(1)));
        C.row(14) = Eigen::Matrix<double, 1, 8>::Ones();

        Cfm = C.transpose();
        Eigen::FullPivHouseholderQR<Eigen::Matrix<double, 8, 8>> qr(Cfm.template block<8, 8>(0, 0));

        Eigen::Matrix<double, 8, 7> G;
        G = qr.solve(Cfm.template block<8, 7>(0, 8));
        GPolynomial g1 = G.row(0);
        GPolynomial g2 = G.row(1);
        GPolynomial g3 = G.row(2);
        GPolynomial g4 = G.row(3);
        GPolynomial g5 = G.row(4);
        GPolynomial g6 = G.row(5);
        GPolynomial g7 = G.row(6);
        GPolynomial g8 = G.row(7);

        return solver(g1, g2, g3, g4, g5, g6, g7, g8);
    }
}