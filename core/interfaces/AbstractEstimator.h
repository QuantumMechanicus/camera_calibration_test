//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACTESTIMATOR_H
#define CAMERA_CALIBRATION_ABSTRACTESTIMATOR_H

namespace estimators {

    namespace internal
    {
        struct IEstimator
        {
            //virtual void estimate() = 0;

            virtual bool isEstimated() const = 0;

            virtual ~IEstimator() = default;
        };
    }

    template<typename TEstimatedParameters>
    class AbstractEstimator : virtual public internal::IEstimator {
    protected:

        virtual void estimateImpl() = 0;

        virtual void getEstimationImpl(TEstimatedParameters &result) = 0;

    public:

        AbstractEstimator() = default;

        AbstractEstimator(const AbstractEstimator &rhs) = default;

        AbstractEstimator(AbstractEstimator &&rhs) noexcept = default;

        ~AbstractEstimator() override = default;

        AbstractEstimator &operator=(AbstractEstimator &&rhs) noexcept = default;

        AbstractEstimator &operator=(const AbstractEstimator &rhs) = default;

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
