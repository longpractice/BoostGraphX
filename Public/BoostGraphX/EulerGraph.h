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
    using V = typename boost::graph_traits<VertexListGraph>::vertex_descriptor;
    using E = typename boost::graph_traits<VertexListGraph>::edge_descriptor;

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

    using EList = std::list<E>;
    using EListIt = typename std::list<E>::iterator;
    using VAuxList = std::list<VAux>;
    using VAuxListIt = typename VAuxList::iterator;

    struct EAuxProp {
        E eOrigin;
        std::optional<EListIt> eWalkListIt;
    };

    struct VAuxProp {
        V vOrigin;
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

        void operator()(V input, VAux output) const
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

        void operator()(E input, EAux output) const
        {
            to[output].eOrigin = input;
        }
    };


    //find euler cycle using Hierholzer algorithm, with out ensuring any edge being the
    //first edge
    static std::list<E> findEulerCycle_Hierholzer(GAux& gAux, VAux firstV)
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
    static std::list<E> findEulerCycle_Hierholzer_ensure_first_edge(
        GAux& gAux, V firstV, EAux e_ensure_first)
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

    static std::list<E> findEulerTour(GAux& gAux, VAux firstV, VAux lastV)
    {
        //First add a edge going from lastV to firstV and
        //Then find an euler cycle starting and ending at lastV
        //After erasing the first edge, we have a euler tour from firstV -> lastV
        auto [e, ok] = boost::add_edge(lastV, firstV, gAux);
		auto edges = findEulerCycle_Hierholzer_ensure_first_edge(gAux, lastV, e);
		//erase the first edge
		edges.erase(edges.begin());
		return edges;
    }
};
}

namespace bglx {

template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_cycle(
    const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_cycle_start_vertex)
{
    using GHelper = detail::EulerDigraghHelper<VertexListGraph>;
    typename GHelper::GAux gAux;
    typename GHelper::VAux vAuxStart;

    //while coping, record the vAux correspondent to euler_cycle_start_vertex
    auto vMarkedCopier =
        [&](typename GHelper::V input, typename GHelper::VAux output) {
            gAux[output].vOrigin = input;
            if (input == euler_cycle_start_vertex) {
                vAuxStart = output;
            }
        };
    typename GHelper::ECopier eCopier(g, gAux);
    boost::copy_graph(
        g,
        gAux,
        boost::vertex_copy(vMarkedCopier).
		edge_copy(eCopier)
	);
    return GHelper::findEulerCycle_Hierholzer(gAux, vAuxStart);
}

template <typename VertexListGraph, typename VertexIndexMap>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_cycle(
    const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_cycle_start_vertex,
    const VertexIndexMap& i_map)
{
	using GHelper = detail::EulerDigraghHelper<VertexListGraph>;
	typename GHelper::GAux gAux;
	typename GHelper::VAux vAuxStart;

	//while coping, record the vAux correspondent to euler_cycle_start_vertex
	auto vMarkedCopier =
		[&](typename GHelper::V input, typename GHelper::VAux output) {
		gAux[output].vOrigin = input;
		if (input == euler_cycle_start_vertex) {
			vAuxStart = output;
		}
	};
	typename GHelper::ECopier eCopier(g, gAux);
	boost::copy_graph(
		g,
		gAux,
		boost::vertex_copy(vMarkedCopier).edge_copy(eCopier).vertex_index_map(i_map));
	return GHelper::findEulerCycle_Hierholzer(gAux, vAuxStart);
}

template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_tour(const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_start_vertex,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_end_vertex)
{
    //should use euler cycle functions if the start and tour equals
    assert(euler_tour_start_vertex != euler_tour_end_vertex);
    using EulerGraph = detail::EulerDigraghHelper<VertexListGraph>;
    typename EulerGraph::GAux gAux;
    typename EulerGraph::VAux vAuxStart;
    typename EulerGraph::VAux vAuxLast;

    auto vMarkedCopier =
        [&](typename EulerGraph::V input, typename EulerGraph::VAux output) {
            gAux[output].vOrigin = input;
            if (input == euler_tour_start_vertex) {
                vAuxStart = output;
            } else if (input == euler_tour_end_vertex) {
                vAuxLast = output;
            }
        };
    typename EulerGraph::ECopier eCopier(g, gAux);
    boost::copy_graph(
        g, gAux,
        boost::vertex_copy(vMarkedCopier).edge_copy(eCopier));

    return EulerGraph::findEulerTour(gAux, vAuxStart, vAuxLast);;
}

template <typename VertexListGraph, typename VertexIndexMap>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_tour(const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_start_vertex,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_end_vertex,
    const VertexIndexMap& i_map)
{
    //should use euler cycle functions if the start and tour equals
    assert(euler_tour_start_vertex != euler_tour_end_vertex);
    using EulerGraph = detail::EulerDigraghHelper<VertexListGraph>;
    typename EulerGraph::GAux gAux;
    typename EulerGraph::VAux vAuxStart;
    typename EulerGraph::VAux vAuxLast;

    auto vMarkedCopier =
        [&](typename EulerGraph::V input, typename EulerGraph::VAux output) {
            gAux[output].vOrigin = input;
            if (input == euler_tour_start_vertex) {
                vAuxStart = output;
            } else if (input == euler_tour_end_vertex) {
                vAuxLast = output;
            }
        };
    typename EulerGraph::ECopier eCopier(g, gAux);
    boost::copy_graph(
        g, gAux,
        boost::vertex_copy(vMarkedCopier).edge_copy(eCopier).vertex_index_map(i_map));

    return EulerGraph::findEulerTour(gAux, vAuxStart, vAuxLast);
}
}
