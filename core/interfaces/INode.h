//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_INODE_H
#define CAMERA_CALIBRATION_INODE_H

namespace graph {

    template<typename TLabel = int>
    struct INode {

        INode<TLabel>() = default;

        INode<TLabel>(const INode<TLabel> &rhs) = default;

        INode<TLabel>(INode<TLabel> &&rhs) noexcept = default;

        virtual ~INode<TLabel>() = default;

        INode<TLabel> &operator=(INode<TLabel> &&rhs) noexcept = default;

        INode<TLabel> &operator=(const INode<TLabel> &rhs) = default;

        //virtual TLabel getLabel() = 0;

    };

}


#endif //CAMERA_CALIBRATION_INODE_H
