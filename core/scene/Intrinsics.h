//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_CAMERA_INTRINSICS_H
#define CAMERA_CALIBRATION_CAMERA_INTRINSICS_H

#include <Eigen/Dense>
#include "../interfaces/Abstract_Estimator.h"
#include "../interfaces/IIntrinsics.h"

namespace intrinsics {


/**
 * @brief Intrinsic parameters of division model for radial distortion (see A. W. Fitzgibbon "Simultaneous linear estimation of multiple view geometry and lens distortion")
 * \f$ x_u = \frac{x_d}{1 + \lambda_1 ||x_d||^2 + ... \lambda_N ||x_d||^{2N}} \f$
 */
    using FocalLength = double;

    template<int N = 1>
    class DivisionModel : public AbstractIntrinsics<DivisionModel<N>> {
        double ppx_{};
        double ppy_{};
        double f_{};
        Eigen::Matrix<double, N, 1> lambdas_;

        friend class AbstractIntrinsics<DivisionModel<N>>;

    protected:
        /**
         * @brief See definition above
         */
        bool isEqualImpl(const AbstractIntrinsics<DivisionModel<N>> &other) const {
            const auto *other_casted = dynamic_cast<const DivisionModel<N> *>(&other);
            return other_casted != nullptr && ppx_ == other_casted->getPrincipalPointX() &&
                   ppy_ == other_casted->getPrincipalPointY() &&
                   f_ == other_casted->getFocalLength() &&
                   lambdas_ == other_casted->getDistortionCoefficients();
        }

        /**
       * @brief See base definition above
       * @param Class with 'estimate' principal point, focal length and distortion coefficients
       */

        void estimateParameterImpl(estimators::AbstractEstimator<DivisionModel<N>> &estimator) {

            *this = estimator.getEstimation();
        }

        void estimateParameterImpl(DivisionModel<N> estimator) {

            *this = std::move(estimator);
        }

        void estimateParameterImpl(estimators::AbstractEstimator<Eigen::Matrix<double, 1, N>> &estimator) {

            lambdas_ = estimator.getEstimation();
        }

        void estimateParameterImpl(estimators::AbstractEstimator<FocalLength> &estimator) {

            f_ = estimator.getEstimation();
        }

        Eigen::Vector2d undistortImpl(const  Eigen::Vector2d &p) {
            double rd = p.norm();
            double denominator(1.0);
            double r_distorted2 = rd * rd;
            double r_distorted2_pow = r_distorted2;
            for (int i = 0; i < lambdas_.rows(); ++i) {
                denominator += lambdas_[i] * r_distorted2_pow;
                r_distorted2_pow *= r_distorted2;
            }

            return p / denominator;

        }

    public:

        /**
         * @brief Constructor
         */
        DivisionModel() = default;

        /**
         * @brief Constructor
         * @param ppx X-axis coordinate of principal point
         * @param ppy Y-axis coordinate of principal point
         * @param f Focal length
         * @param lambdas Parameters of division model
         */
        explicit DivisionModel(const Eigen::Matrix<double, N, 1> &lambdas, unsigned int w = 0,
                               unsigned int h = 0,
                               double f = 0, double ppx = 0,
                               double ppy = 0)
                : AbstractIntrinsics<DivisionModel>(w, h), ppx_(ppx),
                  ppy_(ppy),
                  f_(f),
                  lambdas_(lambdas) {}

        /**
         * @brief Constructor
         * @param ppx X-axis coordinate of principal point
         * @param ppy Y-axis coordinate of principal point
         * @param f Focal length
         * @param lambdas Parameters of division model
         * @param n Number of distortion coefficients (lambdas)
         */
        explicit DivisionModel(unsigned int n, const Eigen::Matrix<double, N, 1> &lambdas, unsigned int w = 0,
                               unsigned int h = 0,
                               double f = 0, double ppx = 0,
                               double ppy = 0)
                : AbstractIntrinsics<DivisionModel>(w, h), ppx_(ppx),
                  ppy_(ppy),
                  f_(f),
                  lambdas_(lambdas) {}

        /**
         @brief Constructor for non-dynamic lambdas
         * @param ppx X-axis coordinate of principal point
         * @param ppy Y-axis coordinate of principal point
         * @param f Focal length
         */
        DivisionModel(unsigned int w, unsigned int h, double f = 0, double ppx = 0,
                      double ppy = 0) : AbstractIntrinsics<DivisionModel>(w, h), ppx_(ppx), ppy_(ppy),
                                        f_(f) {
            assert(N != Eigen::Dynamic && "You should pass number of parameters for dynamic model");
            lambdas_.setZero();
        }

        /**
         @brief Constructor for dynamic lambdas
         * @param ppx X-axis coordinate of principal point
         * @param ppy Y-axis coordinate of principal point
         * @param f Focal length
         * @param n Number of distortion coefficients (lambdas)
         */
        DivisionModel(unsigned int n, unsigned int w, unsigned int h, double f = 0, double ppx = 0,
                      double ppy = 0) : AbstractIntrinsics<DivisionModel>(w, h), ppx_(ppx), ppy_(ppy),
                                        f_(f) {
            lambdas_.resize(n, Eigen::NoChange);
            lambdas_.setZero();
        }


        /**
         * @brief Getter for principal point
         * @return X-axis coordinate of principal point
         */

        double getPrincipalPointX() const {
            return ppx_;
        }

        /**
         * @brief Getter for principal point
         * @return Y-axis coordinate of principal point
         */

        double getPrincipalPointY() const {
            return ppy_;
        }

        /**
         * @brief Getter for focal length
         * @return Focal length
         */

        double getFocalLength() const {
            return f_;
        }

        /**
         * @brief Getter for division model coefficients (see definition above)
         * @return Distortion coefficients
         */
        const Eigen::Matrix<double, N, 1> &getDistortionCoefficients() const {
            return lambdas_;
        }

        int getNumberOfCoefficients() const {
            return static_cast<int>(lambdas_.rows());
        }

        Eigen::Matrix3d getCalibrationMatrix() const {
            Eigen::Matrix3d res(Eigen::Matrix3d::Identity());
            //TODO change remove r
            auto r = std::sqrt(this->w_ * this->w_ + this->h_ * this->h_) / 2;
            res(0, 0) = res(1, 1) = f_;
            res(0, 2) = ppx_;
            res(1, 2) = ppy_;
            return res;
        }

        double getAngeleOfView() const {

            return 2 * std::atan(1 / f_) * 180.0 / M_PI;
        }

        /**
         * @param new_size
         */
        void resizeDistortionCoefficients(long new_size) {
            assert(N == Eigen::Dynamic && "You can't delete or add coefficients on static model");
            long old_size = lambdas_.rows();
            lambdas_.conservativeResize(new_size, Eigen::NoChange);
            if (new_size > old_size)
                lambdas_.bottomRows(new_size - old_size).setZero();
        }
    };

    using PinholeModel = DivisionModel<0>;

    using StandardDivisionModel = DivisionModel<1>;

    using DynamicDivisionModel = DivisionModel<Eigen::Dynamic>;

}
#endif //CAMERA_CALIBRATION_CAMERA_INTRINSICS_H
