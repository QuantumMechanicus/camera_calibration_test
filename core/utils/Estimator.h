//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_ESTIMATOR_H
#define CAMERA_CALIBRATION_ESTIMATOR_H

/*#include <Eigen/Dense>
#include <vector>

namespace estimators {

    namespace internal {
        class BaseEstimator {

        public:
            virtual void estimate() = 0;

            virtual bool isEstimated() const = 0;

            virtual ~BaseEstimator() = default;
        };

        class FundamentalMatrixEstimator : virtual public BaseEstimator {

        protected:
            Eigen::Matrix3d fundamental_matrix_;

        public:
            FundamentalMatrixEstimator() {
                fundamental_matrix_.setZero();
            }

            const Eigen::Matrix3d &getFundamentalMatrix() {
                return fundamental_matrix_;
            }


        };


        template<int N>
        class DivisionModelIntrinsicsEstimator : virtual public BaseEstimator {

        protected:
            double ppx_;
            double ppy_;
            double f_;
            Eigen::Matrix<double, N, 1> lambdas_;

        public:
            explicit DivisionModelIntrinsicsEstimator(double f = 0, double ppx = 0,
                                                      double ppy = 0) : ppx_(ppx), ppy_(ppy), f_(f) {
                lambdas_.setZero();
            }

            explicit DivisionModelIntrinsicsEstimator(unsigned int n, double f = 0, double ppx = 0,
                                                      double ppy = 0) : ppx_(ppx), ppy_(ppy), f_(f) {
                lambdas_.resize(n, Eigen::NoChange);
                lambdas_.setZero();
            }

            explicit DivisionModelIntrinsicsEstimator(Eigen::Matrix<double, N, 1> lambdas, double f = 0, double ppx = 0,
                                                      double ppy = 0) : lambdas_(lambdas), ppx_(ppx), ppy_(ppy), f_(f) {
            }


            double getPrincipalPointX() const {
                return ppx_;
            }

            double getPrincipalPointY() const {
                return ppy_;
            }

            double getFocalLength() const {
                return f_;
            }

            const Eigen::Matrix<double, N, 1> &getDistortionCoefficients() const {
                return lambdas_;
            }

        };
    }


}*/
#endif //CAMERA_CALIBRATION_ESTIMATOR_H
