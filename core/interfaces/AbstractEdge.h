//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACTEDGE_H
#define CAMERA_CALIBRATION_ABSTRACTEDGE_H

#include <memory>
#include "INode.h"

namespace graph {
    template<typename TVertex, typename TNodeLabel = int, typename TWeight = double>
    class AbstractEdge {
    protected:

        TNodeLabel start_vertex_label_;
        TNodeLabel end_vertex_label_;
        std::shared_ptr<std::map<TNodeLabel, TVertex>> ptr_to_list_of_vertices_{};

        //TODO
        //std::weak_ptr<TNodeLabel> st_vertex_{};
        //std::weak_ptr<TNodeLabel> end_vertex_{};

        //TNodeLabel st_vertex_{};
        //TNodeLabel end_vertex_{};

    public:

        AbstractEdge() = default;

        AbstractEdge(const AbstractEdge<TVertex, TNodeLabel, TWeight> &rhs) = default;

        AbstractEdge(AbstractEdge<TVertex, TNodeLabel, TWeight> &&rhs) noexcept = default;

        /*AbstractEdge<TNodeLabel, TWeight>(std::shared_ptr<TNodeLabel> st_vertex, std::shared_ptr<TNodeLabel> end_vertex) :
                st_vertex_(st_vertex),
                end_vertex_(end_vertex) {}*/

        AbstractEdge(TNodeLabel st_vertex, TNodeLabel end_vertex,
                     std::shared_ptr<std::map<TNodeLabel, TVertex>> ptr_to_list_of_vertices) :
                ptr_to_list_of_vertices_(std::move(ptr_to_list_of_vertices)),
                start_vertex_label_(std::move(st_vertex)),
                end_vertex_label_(std::move(end_vertex)){}


        virtual ~AbstractEdge() = default;

        AbstractEdge &operator=(AbstractEdge<TVertex, TNodeLabel, TWeight> &&rhs) noexcept = default;

        AbstractEdge &operator=(const AbstractEdge<TVertex,TNodeLabel, TWeight> &rhs) = default;

        const TVertex &getStartVertex() const {
            return ptr_to_list_of_vertices_->at(start_vertex_label_);
        }

        const TVertex &getFinishVertex() const{
            return ptr_to_list_of_vertices_->at(end_vertex_label_);
        }


        //virtual TWeight getWeight() const = 0;

        /*bool doesExist() const {
            return !st_vertex_.expired() && !end_vertex_.expired();
        }*/
    };
}


#endif //CAMERA_CALIBRATION_ABSTRACTEDGE_H
