//
// Created by danielbord on 2/19/18.
//

#ifndef CAMERA_CALIBRATION_ABSTRACT_EDGE_H
#define CAMERA_CALIBRATION_ABSTRACT_EDGE_H

#include <memory>
#include <map>
#include "INode.h"

namespace graph {
    template<template<typename , typename > class TDerived, typename TVertex, typename TWeight>
    struct IEdge
    {
        using SDerived = TDerived<TVertex, TWeight>;
        const TVertex &getStartVertex() const {
            return static_cast<const SDerived *>(this)->getStartVertexImpl();
        }

        const TVertex &getFinishVertex() const{
            return static_cast<const SDerived *>(this)->getStartVertexImpl();
        }
    };

    template<typename TVertex, typename TWeight = double>
    class AbstractEdge : public IEdge<AbstractEdge, TVertex, TWeight> {
    protected:
        friend class IEdge<AbstractEdge, TVertex, TWeight>;

        typename TVertex::NodeLabel_t start_vertex_label_;
        typename TVertex::NodeLabel_t end_vertex_label_;
        std::shared_ptr<std::map<typename TVertex::NodeLabel_t, TVertex>> ptr_to_list_of_vertices_{};

        const TVertex &getStartVertexImpl() const {
            return ptr_to_list_of_vertices_->at(start_vertex_label_);
        }

        const TVertex &getFinishVertexImpl() const{
            return ptr_to_list_of_vertices_->at(end_vertex_label_);
        }


    public:
        using Vertex_t = TVertex;
        using Weight_t = TWeight;
        using VertexMap_t = std::map<typename TVertex::NodeLabel_t, TVertex>;

        AbstractEdge() : start_vertex_label_{}, end_vertex_label_{}
        {
                ptr_to_list_of_vertices_ = std::make_shared<VertexMap_t >(VertexMap_t());
        }

        AbstractEdge(typename TVertex::NodeLabel_t st_vertex, typename TVertex::NodeLabel_t end_vertex,
                     std::shared_ptr<VertexMap_t > ptr_to_list_of_vertices) :
                ptr_to_list_of_vertices_(std::move(ptr_to_list_of_vertices)),
                start_vertex_label_(std::move(st_vertex)),
                end_vertex_label_(std::move(end_vertex)){}

        const std::shared_ptr<VertexMap_t> &getVertexListPointer() const {
            return ptr_to_list_of_vertices_;
        }

    };
}


#endif //CAMERA_CALIBRATION_ABSTRACT_EDGE_H
