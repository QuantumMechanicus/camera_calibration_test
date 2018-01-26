//
// Created by danielbord on 1/22/18.
//

#include "Camera_Intrinsics.h"

namespace intrinsics {

    IntrinsicsBase::IntrinsicsBase(unsigned int w, unsigned int h) : w_(w), h_(h) {}

    unsigned int IntrinsicsBase::getWidth() const {
        return w_;
    }

    unsigned int IntrinsicsBase::getHeight() const {
        return h_;
    }

    bool IntrinsicsBase::operator==(const IntrinsicsBase& other) const
    {
        return isEqual(other);
    }
}