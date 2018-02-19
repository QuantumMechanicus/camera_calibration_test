//
// Created by danielbord on 1/22/18.
//

#ifndef CAMERA_CALIBRATION_CAMERA_INTRINSICS_H
#define CAMERA_CALIBRATION_CAMERA_INTRINSICS_H

#include <Eigen/Dense>
#include "../interfaces/AbstractEstimator.h"


namespace intrinsics {

/**
 * @brief Base class to store intrinsic parameters of camera (e. g. width, height, focal length)
 */
    template <typename Derived>
    class IntrinsicsBase {

    protected:
        unsigned int w_;
        unsigned int h_;

        /**
         * @brief Equality operator
         * @param other Instance of intrinsics
         * @return True if type of other and this the same and their fields are equal
         */
        virtual bool isEqualImpl(const IntrinsicsBase<Derived> &other) const = 0;

    public:
        /**
         * @brief Constructor
         * @param w Width of the image
         * @param h Height of the image
         */
        explicit IntrinsicsBase(unsigned int w = 0, unsigned int h = 0) : w_(w), h_(h) {};

        /**
        * @brief Destructor
        */
        virtual ~IntrinsicsBase() = default;

        /**
        * @brief Method for identifying unknown parameters of model
        * @param estimator Class with 'estimate' method
        */
        virtual void estimateParameters(estimators::AbstractEstimator<Derived> &estimator) = 0;


        /**
         * @brief Getter for width of the image
         * @return width of the image
         */
        unsigned int getWidth() const
        {
            return w_;
        }

        /**
         * @brief Getter for height of the image
         * @return height of the image
         */
        unsigned int getHeight() const
        {
            return h_;
        }

        /**
        * @brief Equality operator
        * @param other Instance of intrinsics
        * @return True if type of other and this the same and their fields are equal
        */
        bool operator==(const IntrinsicsBase<Derived>& other) const
        {
            return isEqualImpl(other);
        }

    };

/**
 * @brief Intrinsic parameters of division model for radial distortion (see A. W. Fitzgibbon "Simultaneous linear estimation of multiple view geometry and lens distortion")
 * \f$ x_u = \frac{x_d}{1 + \lambda_1 ||x_d||^2 + ... \lambda_N ||x_d||^{2N}} \f$
 */
    template<int N = 1>
    class DivisionModelIntrinsic : public IntrinsicsBase<DivisionModelIntrinsic<N>> {
        double ppx_{};
        double ppy_{};
        double f_{};
        Eigen::Matrix<double, N, 1> lambdas_;

        /**
         * @brief See definition above
         */
        bool isEqualImpl(const IntrinsicsBase<DivisionModelIntrinsic<N>> &other) const override {
            const auto *other_casted = dynamic_cast<const DivisionModelIntrinsic<N> *>(&other);
            return other_casted != nullptr && ppx_ == other_casted->getPrincipalPointX() &&
                   ppy_ == other_casted->getPrincipalPointY() &&
                   f_ == other_casted->getFocalLength() &&
                   lambdas_ == other_casted->getDistortionCoefficients();
        }

    public:
        /**
         * @brief Constructor
         */
        DivisionModelIntrinsic() = default;

        /**
         * @brief Constructor
         * @param ppx X-axis coordinate of principal point
         * @param ppy Y-axis coordinate of principal point
         * @param f Focal length
         * @param lambdas Parameters of division model
         */
        explicit DivisionModelIntrinsic(const Eigen::Matrix<double, N, 1> &lambdas, unsigned int w = 0, unsigned int h = 0,
                                        double f = 0, double ppx = 0,
                                        double ppy = 0)
                : IntrinsicsBase<DivisionModelIntrinsic>(w,h), ppx_(ppx),
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
        explicit DivisionModelIntrinsic(unsigned int n, const Eigen::Matrix<double, N, 1> &lambdas, unsigned int w = 0, unsigned int h = 0,
                               double f = 0, double ppx = 0,
                               double ppy = 0)
                : IntrinsicsBase<DivisionModelIntrinsic>(w,h), ppx_(ppx),
                  ppy_(ppy),
                  f_(f),
                  lambdas_(lambdas) {}

        /**
         @brief Constructor for non-dynamic lambdas
         * @param ppx X-axis coordinate of principal point
         * @param ppy Y-axis coordinate of principal point
         * @param f Focal length
         */
        DivisionModelIntrinsic(unsigned int w, unsigned int h, double f = 0, double ppx = 0,
                                        double ppy = 0) : IntrinsicsBase<DivisionModelIntrinsic>(w, h), ppx_(ppx), ppy_(ppy),
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
        DivisionModelIntrinsic(unsigned int n, unsigned int w, unsigned int h, double f = 0, double ppx = 0,
                                        double ppy = 0) : IntrinsicsBase<DivisionModelIntrinsic>(w, h), ppx_(ppx), ppy_(ppy),
                                                          f_(f) {
            lambdas_.resize(n, Eigen::NoChange);
            lambdas_.setZero();
        }

        /**
         * @brief See base definition above
         * @param Class with 'estimate' principal point, focal length and distortion coefficients
         */

        void estimateParameters(estimators::AbstractEstimator<DivisionModelIntrinsic<N>> &estimator) override {
            auto res = estimator.getEstimation();
            ppx_= res.getPrincipalPointX();
            ppy_= res.getPrincipalPointY();
            f_ = res.getFocalLength();
            lambdas_ = res.getDistortionCoefficients();
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

        int getNumberOfCofficients() const
        {
            return static_cast<int>(lambdas_.rows());
        }

        /**
         * @param new_size
         */
        void resizeDistortionCoefficients(long new_size)
        {
            assert(N == Eigen::Dynamic && "You can't delete or add coefficients on static model");
            long old_size = lambdas_.rows();
            lambdas_.conservativeResize(new_size, Eigen::NoChange);
            if (new_size > old_size)
                lambdas_.bottomRows(new_size-old_size).setZero();
        }
    };

    using PinholeIntrinsic = DivisionModelIntrinsic<0>;

    using StandardDivisionModelIntrinsic = DivisionModelIntrinsic<1>;

}
#endif //CAMERA_CALIBRATION_CAMERA_INTRINSICS_H
