#pragma once
#include <BoostGraphX/Common.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <utility>
#include <vector>

namespace bglx::detail {
template <typename VertexListGraph>
struct EulerDigraghHelper {
    using Vertex = typename boost::graph_traits<VertexListGraph>::vertex_descriptor;
    using Edge = typename boost::graph_traits<VertexListGraph>::edge_descriptor;

    using DirectedCategory = typename boost::graph_traits<VertexListGraph>::directed_category;
    using IsDirected = std::is_convertible<DirectedCategory, boost::directed_tag>;

    //
    // The algorithms of finding a directed euler cycle and a undirected euler cycle are very similar, both could
    // be solved by Hierholz algorithm.
    //
    // However, the issue is that for the directed case, marking the edge as visited could be very effitient
    // it is more tricky for the undirected case.
    //
    static_assert(IsDirected::value,
        "Calculating euler directed cycle and directed tour requires a directed graph");

    //
    // We want to save the vertex descriptor inside the vertex property(for linked-list iterator).
    //
    // However the vertex descriptor type depends on the graph type which depends on the vertex
    // property type.
    //
    // We therefor "hack here" to use the vertex_descriptor not dependent on vertex property
    // Indeed, the vertex descriptor is just a size_t here. Not that a big concern.
    //
    using DagAuxPure = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
    using VAux = typename DagAuxPure::vertex_descriptor;

    using EList = std::list<Edge>;
    using EListIt = typename std::list<Edge>::iterator;
    using VAuxList = std::list<VAux>;
    using VAuxListIt = typename VAuxList::iterator;

    struct EAuxProp {
        Edge eOrigin;
        std::optional<EListIt> eWalkListIt;
    };

    struct VAuxProp {
        Vertex vOrigin;
        size_t firstUnVisitedEdgeId { 0 };
        // we could not use EAuxListIt here, since once it is default constructed
        // it is a "singular" iterator that we could not even assign to it
        std::optional<VAuxListIt> vUnfinishedListIt;
    };

    using GAux = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::directedS,
        VAuxProp,
        EAuxProp>;

    using EAux = typename GAux::edge_descriptor;

    //
    // Find from the unvisited edges a closed walk starting and ending at v
    // Time complexity: O(E)
    //
    // The new walk is put as a "slice" into the doubly linked list of old closed walk.
    // Small constant time.
    //
    // You could possibly specify a mark edge here.
    // The function returns the start edge of the new circle from v,
    // The bool part of the return pair is true if the start edge of the
    // closed walk corresponds to markEdge
    //
    static std::pair<EListIt, bool> discoverNewClosedWalk(
        VAux v,
        GAux& aux,
        VAuxList& unfinishedVList,
        EList& oldClosedWalk,
        std::optional<EAux> markEdge = {})
    {
        EList newClosedWalk;
        auto vFirst = v;
        if (boost::out_degree(v, aux) == aux[v].firstUnVisitedEdgeId) {
            //we could not find any exit and thus no new closed walk from this v
            return { {}, false };
        }
        std::unique_ptr<EListIt> prevVisitedEOut;
        auto outEdges = boost::out_edges(v, aux);
        auto firstUnvisitedE = *(outEdges.first + aux[v].firstUnVisitedEdgeId);
        bool isMarked = markEdge ? firstUnvisitedE == markEdge.value() : false;

        if (aux[v].firstUnVisitedEdgeId != 0) {
            const auto& visitedE = *(outEdges.first + aux[v].firstUnVisitedEdgeId - 1);
            const auto& itVisitedE = aux[visitedE].eWalkListIt;
            prevVisitedEOut = std::make_unique<EListIt>(*itVisitedE);
        }

        do {
            //find the first exit for current vertex
            auto& vProp = aux[v];
            auto outDeg = boost::out_degree(v, aux);
            auto itEOut = boost::out_edges(v, aux).first + vProp.firstUnVisitedEdgeId;
            auto eOut = *itEOut;
            //this exist edge is now visited, move it from unvisited to visited
            ++vProp.firstUnVisitedEdgeId;

            newClosedWalk.emplace_back(aux[eOut].eOrigin);
            aux[eOut].eWalkListIt = { --newClosedWalk.end() };

            //v was partially visited before, if we have finished the last
            bool justBecomeUnfinished = vProp.firstUnVisitedEdgeId == 1 && vProp.firstUnVisitedEdgeId != outDeg;
            if (justBecomeUnfinished) {
                unfinishedVList.emplace_back(v);
                vProp.vUnfinishedListIt = { --unfinishedVList.end() };
            }

            //it was partially visited before and now need to be moved out
            bool justMovedOutOfUnfinishedVList = vProp.firstUnVisitedEdgeId == outDeg && outDeg > 1;
            if (justMovedOutOfUnfinishedVList) {
                unfinishedVList.erase(*vProp.vUnfinishedListIt);
            }
            v = boost::target(eOut, aux);
        } while (v != vFirst);

        auto eListIt = newClosedWalk.begin();
        if (prevVisitedEOut) {
            //the visited out edge, which need to break the previous linked list
            oldClosedWalk.splice(*prevVisitedEOut, newClosedWalk);
        } else {
            //old closed walk must be empty
            std::swap(oldClosedWalk, newClosedWalk);
        }
        return { eListIt, isMarked };
    }

    struct VCopier {
        VCopier(const VertexListGraph& from, GAux& to)
            : from { from }
            , to { to } {};
        const VertexListGraph& from;
        GAux& to;

        void operator()(Vertex input, VAux output) const
        {
            to[output].vOrigin = input;
        }
    };

    struct ECopier {
        ECopier(const VertexListGraph& from, GAux& to)
            : from { from }
            , to { to } {};

        const VertexListGraph& from;
        GAux& to;

        void operator()(Edge input, EAux output) const
        {
            to[output].eOrigin = input;
        }
    };

    //create a auxilliary dag, this is indeed take half of the whole computing time
    static void makeAuxDag(VertexListGraph& gOrigin, GAux& gAux)
    {
        boost::copy_graph(
            gOrigin,
            gAux,
            boost::vertex_copy(VCopier(gOrigin, gAux)).edge_copy(ECopier(gOrigin, gAux)));
    }

    //find euler cycle using Hierholzer algorithm, with out ensuring any edge being the
    //first edge
    static std::list<Edge> __findEulerCycle_Hierholzer(GAux& gAux, Vertex firstV)
    {
        VAuxList unfinishedVList;
        EList closedWalk;
        auto nrEdges = boost::num_edges(gAux);
        discoverNewClosedWalk(firstV, gAux, unfinishedVList, closedWalk, {});
        while (closedWalk.size() < nrEdges) {
            auto v = unfinishedVList.front();
            discoverNewClosedWalk(v, gAux, unfinishedVList, closedWalk, {});
        }
        return closedWalk;
    }

    //
    // find euler cycle start and end at firstV
    // e_ensure_first will be the first edge in the cycle
    // This is useful for algorithm for Eulerian trail (or Eulerian path) with
    // the start vertex and the end vertex different
	//
    static std::list<Edge> __findEulerCycle_Hierholzer_ensure_first_edge(
        GAux& gAux, Vertex firstV, EAux e_ensure_first)
    {
        VAuxList unfinishedVList;
        EList closedWalk;
        auto eMarkedIt = closedWalk.end();
        auto nrEdges = boost::num_edges(gAux);
        while (closedWalk.size() < nrEdges) {
            if (unfinishedVList.empty() || unfinishedVList.front() == firstV) {
                auto [eItAdd, ifMarked] = discoverNewClosedWalk(firstV, gAux, unfinishedVList, closedWalk, { e_ensure_first });
                if (ifMarked) {
                    eMarkedIt = eItAdd;
                }
            } else {
                auto v = unfinishedVList.front();
                discoverNewClosedWalk(v, gAux, unfinishedVList, closedWalk, {});
            }
        }
        EList closedWalkReordered;
        closedWalkReordered.splice(closedWalkReordered.begin(), closedWalk, eMarkedIt, closedWalk.end());
        closedWalkReordered.splice(closedWalkReordered.end(), closedWalk);
        return closedWalkReordered;
    }
};
}

namespace bglx {
template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_euler_cycle(const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_cycle_start_vertex)
{
    using EulerGraph = detail::EulerDigraghHelper<VertexListGraph>;
    typename EulerGraph::GAux gAux;
    typename EulerGraph::VCopier vCopier(g, gAux);
    typename EulerGraph::ECopier eCopier(g, gAux);
    boost::copy_graph(g, gAux,
        boost::vertex_copy(vCopier).edge_copy(eCopier));
    return EulerGraph::__findEulerCycle_Hierholzer(gAux, euler_cycle_start_vertex);
}

template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_euler_cycle(const VertexListGraph& g)
{
    if (boost::num_vertices(g) == 0) {
        return {};
    }

    auto vStart = *boost::vertices(g).first;
    return find_one_euler_cycle(g, vStart);
}

template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_euler_tour(const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_start_vertex,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_end_vertex)
{
    //assert(euler_c)
}

}
