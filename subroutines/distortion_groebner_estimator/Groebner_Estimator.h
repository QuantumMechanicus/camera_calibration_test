//
// Created by danielbord on 1/23/18.
//

#ifndef CAMERA_CALIBRATION_GROEBNER_ESTIMATOR_H
#define CAMERA_CALIBRATION_GROEBNER_ESTIMATOR_H

#include "Core.h"

namespace estimators {

    struct GroebnerEstimatorOptions {
        int number_of_RANSAC_iterations = 1000000;
        double quantile_to_minimize = 0.1;
        double lambda_lower_bound = -1;
        double lambda_upper_bound = 0.25;
        double image_radius = 1;

        explicit GroebnerEstimatorOptions(int number_of_RANSAC_iterations = 1000000, double quantile_to_minimize = 0.1,
                                 double lambda_lower_bound = -1, double lambda_upper_bound = 0.25, double image_radius = 1);
    };

    class GroebnerDivisionModelEstimator
            : public internal::DivisionModelIntrinsicsEstimator<1>, public estimators::internal::FundamentalMatrixEstimator {

        scene::ImagePoints u1d_, u2d_;
        std::size_t number_of_points_;
        GroebnerEstimatorOptions options_;
        std::vector<double> left_residuals_;
        std::vector<double> right_residuals_;
        bool is_estimated_;
        struct AutomaticSolver {

            typedef Eigen::Matrix<double, 7, 1> GPolynomial;
            typedef Eigen::Matrix<double, 2, 8> EightPoints;
            typedef std::pair<scene::StdVector<scene::FundamentalMatrix>, std::vector<double> > FundamentalMatricesAndDistortionCoefficients;

            FundamentalMatricesAndDistortionCoefficients runSolver(EightPoints u1d, EightPoints u2d);

            FundamentalMatricesAndDistortionCoefficients solver(const GPolynomial &g1,
                                                                const GPolynomial &g2,
                                                                const GPolynomial &g3,
                                                                const GPolynomial &g4,
                                                                const GPolynomial &g5,
                                                                const GPolynomial &g6,
                                                                const GPolynomial &g7,
                                                                const GPolynomial &g8);
        };

        double computeErrorsAndEstimateQuantile(double distortion_coefficient,
                                                const scene::FundamentalMatrix &fundamental_matrix);

    public:
        GroebnerDivisionModelEstimator(scene::ImagePoints u1d,
                                       scene::ImagePoints u2d, GroebnerEstimatorOptions opt = GroebnerEstimatorOptions());

        void setOptions(const GroebnerEstimatorOptions &options_);

        void estimate() override;

        bool isEstimated() const override;

    };
}

#endif //CAMERA_CALIBRATION_GROEBNER_ESTIMATOR_H
