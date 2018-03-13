//
// Created by danielbord on 2/28/18.
//

#ifndef CAMERA_CALIBRATION_IINTRINSICS_H
#define CAMERA_CALIBRATION_IINTRINSICS_H
namespace intrinsics {
    /**
     * @brief Base class to store intrinsic parameters of camera (e. g. width, height, focal length)
     * @tparam TDerived --- CRTP
     */
    template<typename TDerived>
    class AbstractIntrinsics {

    protected:
        unsigned int w_;
        unsigned int h_;


    public:
        /**
         * @brief Constructor
         * @param w Width of the image
         * @param h Height of the image
         */
        explicit AbstractIntrinsics(unsigned int w = 0, unsigned int h = 0) : w_(w), h_(h) {};


        /**
        * @brief Method for identifying unknown parameters of model
        * @param estimator Class with 'estimate' method or concrete estimation (see different implementations for details)
        */
        template<typename TEstimator>
        void estimateParameter(TEstimator &estimator) {
            static_cast<TDerived *>(this)->estimateParameterImpl(estimator);
        }

        Eigen::Vector2d undistort(const Eigen::Vector2d &p)
        {
            static_cast<TDerived *>(this)->undistortImpl(p);
        }

        /**
         * @brief Getter for width of the image
         * @return width of the image
         */
        unsigned int getWidth() const {
            return w_;
        }

        /**
         * @brief Getter for height of the image
         * @return height of the image
         */
        unsigned int getHeight() const {
            return h_;
        }

        /**
        * @brief Equality operator
        * @param other Instance of intrinsics
        * @return Expected to be true if type of other and this the same and their fields are equal, otherwise false
        */
        bool operator==(const AbstractIntrinsics<TDerived> &other) const {

            return static_cast<TDerived *>(this)->isEqualImpl(other);
        }

    };
}
#endif //CAMERA_CALIBRATION_IINTRINSICS_H
