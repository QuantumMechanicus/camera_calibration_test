//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_INODE_H
#define CAMERA_CALIBRATION_INODE_H

namespace graph {

    template<typename TLabel = int>
    struct INode {

        typedef TLabel NodeLabel;

        virtual ~INode<TLabel>() = default;

        virtual TLabel getLabel() const = 0;

    };

}


#endif //CAMERA_CALIBRATION_INODE_H
