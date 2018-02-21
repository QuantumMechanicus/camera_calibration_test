//
// Created by danielbord on 2/21/18.
//

#ifndef CAMERA_CALIBRATION_IGRAPH_H
#define CAMERA_CALIBRATION_IGRAPH_H

#include "INode.h"
#include "Abstract_Edge.h"

namespace graph {
    template <typename TNode = INode<>, typename TEdge = AbstractEdge<TNode>>
    struct IGraph
    {
        virtual ~IGraph() = default;

    };

}
#endif //CAMERA_CALIBRATION_IGRAPH_H
