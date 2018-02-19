//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ITWOVIEW_H
#define CAMERA_CALIBRATION_ITWOVIEW_H

#include "Eigen/Dense"
#include "AbstractEstimator.h"

template<typename Intrinsics>
class ITwoView {
public:
    virtual void estimateFundamentalMatrix(estimators::AbstractEstimator<Eigen::Matrix3d> &estimator) = 0;
};


#endif //CAMERA_CALIBRATION_ITWOVIEW_H
