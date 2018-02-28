//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_INODE_H
#define CAMERA_CALIBRATION_INODE_H

namespace graph {

    /**
     * @brief Base class for labled vertex
     * @tparam TDerived --- CRTP
     * @tparam TLabel --- information stored in vertex
     */
    template<typename TDerived, typename TLabel = int>
    struct INode {

        using Label_t = TLabel;

        TLabel getLabel() const
        {
            return static_cast<const TDerived*>(this)->getLabelImpl();
        }

    };

}


#endif //CAMERA_CALIBRATION_INODE_H
