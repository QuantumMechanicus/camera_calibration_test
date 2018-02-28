//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_INODE_H
#define CAMERA_CALIBRATION_INODE_H

namespace graph {

    template<typename TDerived, typename TLabel = int>
    struct INode {

        typedef TLabel NodeLabel;

        TLabel getLabel() const
        {
            return static_cast<const TDerived*>(this)->getLabelImpl();
        }

        using NodeLabel_t = TLabel;

    };

}


#endif //CAMERA_CALIBRATION_INODE_H
