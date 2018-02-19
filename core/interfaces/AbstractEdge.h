//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACTEDGE_H
#define CAMERA_CALIBRATION_ABSTRACTEDGE_H

#include <memory>
#include "INode.h"

namespace graph {
    template<typename TNode = INode<>, typename TWeight = double>
    class AbstractEdge {
    protected:

        //TODO
        //std::weak_ptr<TNode> st_vertex_{};
        //std::weak_ptr<TNode> end_vertex_{};

        TNode st_vertex_{};
        TNode end_vertex_{};

    public:

        AbstractEdge<TNode, TWeight>() = default;

        AbstractEdge<TNode, TWeight>(const AbstractEdge<TNode, TWeight> &rhs) = default;

        AbstractEdge<TNode, TWeight>(AbstractEdge<TNode, TWeight> &&rhs) noexcept = default;

        /*AbstractEdge<TNode, TWeight>(std::shared_ptr<TNode> st_vertex, std::shared_ptr<TNode> end_vertex) :
                st_vertex_(st_vertex),
                end_vertex_(end_vertex) {}*/

        AbstractEdge<TNode, TWeight>(TNode st_vertex, TNode end_vertex) :
                st_vertex_(std::move(st_vertex)),
                end_vertex_(std::move(end_vertex)) {}


        virtual ~AbstractEdge<TNode, TWeight>() = default;

        AbstractEdge<TNode, TWeight> &operator=(AbstractEdge<TNode, TWeight> &&rhs) noexcept = default;

        AbstractEdge<TNode, TWeight> &operator=(const AbstractEdge<TNode, TWeight> &rhs) = default;

        //virtual TWeight getWeight() const = 0;

        /*bool doesExist() const {
            return !st_vertex_.expired() && !end_vertex_.expired();
        }*/
    };
}


#endif //CAMERA_CALIBRATION_ABSTRACTEDGE_H
