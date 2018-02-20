//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACTESTIMATOR_H
#define CAMERA_CALIBRATION_ABSTRACTESTIMATOR_H

namespace estimators {


    template<typename TEstimatedParameters>
    class AbstractEstimator {
    protected:

        virtual void estimateImpl() = 0;

        virtual void getEstimationImpl(TEstimatedParameters &result) = 0;

    public:

        AbstractEstimator() = default;

        AbstractEstimator(const AbstractEstimator &rhs) = default;

        AbstractEstimator(AbstractEstimator &&rhs) noexcept = default;

        virtual ~AbstractEstimator() = default;

        AbstractEstimator &operator=(AbstractEstimator &&rhs) noexcept = default;

        AbstractEstimator &operator=(const AbstractEstimator &rhs) = default;

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

#endif //CAMERA_CALIBRATION_ABSTRACTESTIMATOR_H
