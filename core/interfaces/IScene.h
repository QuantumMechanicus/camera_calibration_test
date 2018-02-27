//
// Created by danielbord on 2/21/18.
//

#ifndef CAMERA_CALIBRATION_ISCENE_H
#define CAMERA_CALIBRATION_ISCENE_H


#include "ITwo_View.h"
#include <memory>

namespace scene {
    template<typename TDerived>
    struct IScene {

        template <typename TLabel, typename TEstimator>
        void estimateCamera(const TLabel &label, TEstimator &estimator)
        {
            static_cast<TDerived*>(this)->estimateCameraImpl(label, estimator);
        }

        //void estimateCameras();

        template <typename TEstimator>
        void estimateStereoPair(size_t label, TEstimator &estimator)
        {
            static_cast<TDerived*>(this)->estimateStereoPair(label, estimator);
        }

        template <typename TEstimator>
        void estimateStereoPairs(TEstimator &estimator)
        {
            static_cast<TDerived*>(this)->estimateStereoPairs(estimator);
        }

    };
}
#endif //CAMERA_CALIBRATION_ISCENE_H
