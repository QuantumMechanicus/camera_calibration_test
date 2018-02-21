//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACT_ESTIMATOR_H
#define CAMERA_CALIBRATION_ABSTRACT_ESTIMATOR_H

namespace estimators {


    template<typename TEstimatedParameters>
    class AbstractEstimator {
    protected:

        virtual void estimateImpl() = 0;

        virtual void getEstimationImpl(TEstimatedParameters &result) = 0;

    public:

        virtual ~AbstractEstimator() = default;

        virtual bool isEstimated() const = 0;

        void estimate() {
            estimateImpl();
        }

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
