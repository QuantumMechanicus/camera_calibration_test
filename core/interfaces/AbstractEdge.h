//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACTEDGE_H
#define CAMERA_CALIBRATION_ABSTRACTEDGE_H

#include <memory>
#include "INode.h"

namespace graph {
    template<typename TVertex, typename TWeight = double>
    class AbstractEdge {
    protected:
        typename TVertex::TNodeLabel start_vertex_label_;
        typename TVertex::TNodeLabel end_vertex_label_;
        std::shared_ptr<std::map<typename TVertex::TNodeLabel, TVertex>> ptr_to_list_of_vertices_{};

    public:

        AbstractEdge() = default;

        AbstractEdge(const AbstractEdge<TVertex, TWeight> &rhs) = default;

        AbstractEdge(AbstractEdge<TVertex, TWeight> &&rhs) noexcept = default;

        AbstractEdge(typename TVertex::TNodeLabel st_vertex, typename TVertex::TNodeLabel end_vertex,
                     std::shared_ptr<std::map<typename TVertex::TNodeLabel, TVertex>> ptr_to_list_of_vertices) :
                ptr_to_list_of_vertices_(std::move(ptr_to_list_of_vertices)),
                start_vertex_label_(std::move(st_vertex)),
                end_vertex_label_(std::move(end_vertex)){}

        virtual ~AbstractEdge() = default;

        AbstractEdge &operator=(AbstractEdge<TVertex, TWeight> &&rhs) noexcept = default;

        AbstractEdge &operator=(const AbstractEdge<TVertex, TWeight> &rhs) = default;

        const TVertex &getStartVertex() const {
            return ptr_to_list_of_vertices_->at(start_vertex_label_);
        }

        const TVertex &getFinishVertex() const{
            return ptr_to_list_of_vertices_->at(end_vertex_label_);
        }

        //virtual TWeight getWeight() const = 0;

    };
}


#endif //CAMERA_CALIBRATION_ABSTRACTEDGE_H
