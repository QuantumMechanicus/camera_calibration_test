//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACT_ESTIMATOR_H
#define CAMERA_CALIBRATION_ABSTRACT_ESTIMATOR_H

namespace estimators {

    /**
     * @brief Base class for different estimators
     * @tparam TEstimatedParameters --- class which concrete estimator can evaluate
     */
    template<typename TEstimatedParameters>
    class AbstractEstimator {
    protected:

        /**
         * @brief Implementation of evaluation process, see estimation()
         */
        virtual void estimateImpl() = 0;

        /**
         * @brief Implementation of getting result, see getEstimation()
         * @param result --- computed estimation
         */
        virtual void getEstimationImpl(TEstimatedParameters &result) = 0;

    public:

        virtual ~AbstractEstimator() = default;

        /**
         * @brief Check if estimation completed
         * @return true if estimate() finished job
         */
        virtual bool isEstimated() const = 0;

        /**
         * @brief Run estimation of TEstimatedParameters class
         */
        void estimate() {
            estimateImpl();
        }

        /**
         * @brief Run estimation and get result
         * @return computed estimation
         */
        TEstimatedParameters getEstimation() {
            if (!isEstimated())
                estimate();
            TEstimatedParameters result;
            getEstimationImpl(result);
            return result;
        }
    };
}

#endif //CAMERA_CALIBRATION_ABSTRACT_ESTIMATOR_H
