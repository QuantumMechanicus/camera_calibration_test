//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_CAMERA_INTRINSICS_H
#define CAMERA_CALIBRATION_CAMERA_INTRINSICS_H

#include <Eigen/Dense>
#include "../utils/Estimator.h"


namespace intrinsics {
/**
 * @brief Base class to store intrinsic parameters of camera (e. g. width, height, focal length)
 */
    class IntrinsicsBase {

    protected:
        unsigned int w_;
        unsigned int h_;

    public:
        /**
         * @brief Constructor
         * @param w Width of the image
         * @param h Height of the image
         */
        explicit IntrinsicsBase(unsigned int w = 0, unsigned int h = 0);

        /**
        * @brief Destructor
        */
        virtual ~IntrinsicsBase() = default;

        /**
        * @brief Method for identifying unknown parameters of model
        * @param estimator Class with 'estimate' method
        */
        virtual void estimateParameters(estimators::internal::BaseEstimator &estimator) = 0;

        /**
         * @brief Getter for width of the image
         * @return width of the image
         */
        unsigned int getWidth() const;

        /**
         * @brief Getter for height of the image
         * @return height of the image
         */
        unsigned int getHeight() const;

    };

/**
 * @brief Intrinsic parameters of division model for radial distortion (see A. W. Fitzgibbon "Simultaneous linear estimation of multiple view geometry and lens distortion")
 * \f$ x_u = \frac{x_d}{1 + \lambda_1 ||x_d||^2 + ... \lambda_N ||x_d||^{2N}} \f$
 */
    template<unsigned int N = 1>
    class DivisionModelIntrinsic : public IntrinsicsBase {
        double ppx_;
        double ppy_;
        double f_;
        Eigen::Matrix<double, N, 1> lambdas_;
    public:
        /**
         * @brief Constructor
         * @param ppx_ X-axis coordinate of principal point
         * @param ppy_ Y-axis coordinate of principal point
         * @param f_ Focal length
         * @param lambdas_ Parameters of division model
         */
        DivisionModelIntrinsic(double f_, const Eigen::Matrix<double, N, 1> &lambdas_, double ppx_ = 0, double ppy_ = 0)
                : ppx_(ppx_),
                  ppy_(ppy_),
                  f_(f_),
                  lambdas_(lambdas_) {}

        /**
         @brief Constructor
         * @param ppx_ X-axis coordinate of principal point
         * @param ppy_ Y-axis coordinate of principal point
         * @param f_ Focal length
         */
        explicit DivisionModelIntrinsic(unsigned int w_ = 0, unsigned int h_ = 0, double f_ = 0, double ppx_ = 0,
                                        double ppy_ = 0) : IntrinsicsBase(w_, h_), ppx_(ppx_), ppy_(ppy_),
                                                           f_(f_) {
            lambdas_.setZero();
        }

        /**
         * @brief See base definition above
         * @param Class with 'estimate' principal point, focal length and distortion coefficients
         */

        void estimateParameters(estimators::internal::BaseEstimator &estimator) override {
            try {
                auto &division_model_estimator = dynamic_cast<estimators::internal::DivisionModelIntrinsicsEstimator<N> &>(estimator);
                if (!division_model_estimator.isEstimated())
                    division_model_estimator.estimate();

                ppx_ = division_model_estimator.getPrincipalPointX();
                ppy_ = division_model_estimator.getPrincipalPointY();
                f_ = division_model_estimator.getFocalLength();
                lambdas_ = division_model_estimator.getDistortionCoefficients();
            }
            catch (const std::bad_cast &e) {
                ppx_ = 0;
                ppy_ = 0;
                f_ = 0;
                lambdas_.setZero();
            }

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
    };

/**
 * @brief If there is no distortion coefficients it is pinhole model
 */
    template<>
    class DivisionModelIntrinsic<0> {
        double ppx_;
        double ppy_;
        double f_;

    public:
        /**
         * @brief See definitions above
         */
        explicit DivisionModelIntrinsic(double f_ = 0, double ppx_ = 0, double ppy_ = 0) : ppx_(ppx_), ppy_(ppy_),
                                                                                           f_(f_) {}
    };

    using PinholeIntrinsic = DivisionModelIntrinsic<0>;

    using StandardDivisionModelIntrinsic = DivisionModelIntrinsic<1>;

}
#endif //CAMERA_CALIBRATION_CAMERA_INTRINSICS_H
