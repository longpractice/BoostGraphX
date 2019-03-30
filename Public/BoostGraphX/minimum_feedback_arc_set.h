#pragma once
#include <BoostGraphX/directed_bidirectional_graph_concept.h>

#include <Boost/concept/requires.hpp>
#include <Boost/graph/adjacency_list.hpp>
/*
	 Works on boost directed graph

	 The graph must be directed and model a concept of BidirectionalGraph,

	 A fast and effective heuristic implementation
	 for solving the minimum feedback arc set problem

	 Based on

	 Eades, Peter, Xuemin Lin, and William F. Smyth.
	 "A fast and effective heuristic for the feedback arc set problem."
	 Information Processing Letters 47.6 (1993): 319-323.

	 But this is a version with the weight of the edge
*/

//Defines the graph concept the method accepts, directed and Bidirectional
template <typename G>
struct Directed_Bidirectional_Graph_Concept {
	BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<G>));
	static_assert(std::is_convertible<
		boost::graph_traits<G>::directed_category, boost::directed_tag>::value);
};

namespace detail {
	template <typename G>
	struct Graph_Info {
		BOOST_CONCEPT_ASSERT((Directed_Bidirectional_Graph_Concept<G>));

		using Graph_V = G::vertex_descriptor;
		using Graph_E = G::edge_descriptor;

		//vertex descriptors and edge descriptors are normally very small
		std::unordered_map<Graph_V, std::set<Graph_E>> in_edges_by_v;
		std::unordered_map<Graph_V, std::set<Graph_E>> out_edges_by_v;

		int in_degree(Graph_V v)
		{
			return in_edges_by_v.find(v)->second.size();
		}

		int out_degree(Graph_V v)
		{
			return out_edges_by_v.find(v)->second.size();
		}

		//add one vertice from graph to this graph info
		void add(Graph_V v, const G& g)
		{
			auto in_edges_it_pair = boost::in_edges(v, g);
			std::set<Graph_E> in_edges(in_edges_it_pair.first, in_edges_it_pair.second);

			auto out_edges_it_pair = boost::out_edges(v, g);
			std::set<Graph_E> out_edges(out_edges_it_pair.first, out_edges_it_pair.second);

			in_edges_by_v.emplace(v, std::move(in_edges));
			out_edges_by_v.emplace(v, std::move(out_edges));
		}

		/*
				Remove one vertex from this data structure
				Return the tuple consists of:
				edges with its v-source unremoved,
				new vertices emerge as sinks, and
				new vertices emerge as sources.
				if false == if_find_non_orphan_edges, we will ignore the term edges with its v-source unremoved
				since when we are removing sinks, this information does not give a feedback arc
			*/
		std::tuple<std::set<Graph_E>, std::unordered_set<Graph_V>, std::unordered_set<Graph_V>>
			remove(
				Graph_V v,
				const G& g,
				bool if_find_non_orphan_edges)
		{
			std::unordered_set<Graph_V> emerged_sinks;
			std::unordered_set<Graph_V> emerged_sources;
			std::set<Graph_E> edges_with_unremoved_vsource;

			/*
				 * First remove the in-edges
				 * in this process, we may generate new sinks
				 */
			const auto& in_edges = in_edges_by_v.find(v)->second;
			for (const auto& in_edge : in_edges) {
				auto src = boost::source(in_edge, g);

				auto it_src_out_edges = out_edges_by_v.find(src);
				if (it_src_out_edges == out_edges_by_v.end()) {
					//should not happen though
					continue;
				}

				if (it_src_out_edges->second.size() == 1) {
					//if the edge to be erased is the only out edge of
					emerged_sinks.insert(it_src_out_edges->first);
				}

				it_src_out_edges->second.erase(in_edge);

				if (if_find_non_orphan_edges) {
					edges_with_unremoved_vsource.insert(in_edge);
				}
			}

			/*
			 * Remove the out-edges
			 * In this process, we may generate new sources
			 */
			const auto& out_edges = out_edges_by_v.find(v)->second;
			for (const auto& out_edge : out_edges) {
				//the target of this out edge, who will be affected
				auto tgt = boost::target(out_edge, g);
				//the in-edges of the tgt vertex
				auto it_tgt_in_edges = in_edges_by_v.find(tgt);
				if (it_tgt_in_edges == in_edges_by_v.end()) {
					//should not happen though
					continue;
				}

				if (it_tgt_in_edges->second.size() == 1) {
					//if the edge to be erased is the only in edge of
					emerged_sources.insert(it_tgt_in_edges->first);
				}

				it_tgt_in_edges->second.erase(out_edge);
			}

			in_edges_by_v.erase(v);
			out_edges_by_v.erase(v);

			return {
				std::move(edges_with_unremoved_vsource),
				std::move(emerged_sinks),
				std::move(emerged_sources)
			};
		}
	};

	/*
	 * The paper use the (out-degree - in-degree) for the criteria to remove vertices
	 * who are not a source nor a sink.
	 * We make a weighed version here
	 * We remove the vertex whose "(in-edge weight sum) - (out-edge weight sum)" is minimum
	 * The biase here is kept in a multimap so that updating and finding min is fast
	 */

	template <typename G>
	class Weighed_Degree_Info {
	public:
		BOOST_CONCEPT_ASSERT((Directed_Bidirectional_Graph_Concept<G>));
		using Graph_V = G::vertex_descriptor;
		using Graph_E = G::edge_descriptor;

		std::multimap<int, Graph_V> weighed_in_out_degree_biases;
		std::unordered_map<Graph_V, std::multimap<int, Graph_V>::iterator> degree_biase_iter_by_v;

		//the one with least
		Graph_V v_with_least_in_out_weighed_degree_biase()
		{
			return { weighed_in_out_degree_biases.begin()->second };
		}

		//W is a weight returning callable, it should have signature of
		//int(const G& g, const G::edge_descriptor& e), that is,
		//it should return a weight of an edge given the graph and the edge descriptor
		template <typename W>
		void add(Graph_V v, const Graph_Info<G>& graph_info, const G& g, W& w)
		{
			auto weighed_sum = int{};

			const auto& in_edges = graph_info.in_edges_by_v.find(v)->second;
			for (auto in_edge : in_edges) {
				weighed_sum += w(g, in_edge);
			}

			const auto& out_edges = graph_info.out_edges_by_v.find(v)->second;
			for (auto out_edge : out_edges) {
				weighed_sum -= w(g, out_edge);
			}
			auto biase_rank_iter = weighed_in_out_degree_biases.emplace(weighed_sum, v);
			degree_biase_iter_by_v.emplace(v, biase_rank_iter);
		}

		template <typename W>
		void remove(DgV v, const Graph_Info<G>& graph_info, const G& g, W& w)
		{
			const auto& edges_in = graph_info.in_edges_by_v.find(v)->second;
			/*
			 * We are going to remove all the in-edges,
			 * therefore the source nodes are affected.
			 */
			for (const auto& edge_in : edges_in) {
				auto v_src = boost::source(edge_in, dg);
				auto weight = w(g, edge_in);
				auto it = degree_biase_iter_by_v.find(v_src);
				if (it != degree_biase_iter_by_v()) {
					auto w_old = get_biase(v_src);
					auto wNew = w_old + w;
					update_biase(v_src, wNew);
				}
			}

			/*
			 * set the out edge's weight to 1
			 */
			const auto& out_edges = graph_info.out_edges_by_v.find(v)->second;
			for (auto e : out_edges) {
				auto v_dst = boost::target(e, g);
				auto weight = w(g, edge_in);
				auto it = degree_biase_iter_by_v.find(v_dst)
					if (it != degree_biase_iter_by_v.end())
					{
						auto w_old = get_biase(v_dst);
						auto w_new = wOld - w;
						update_biase(v_dst, wNew);
					}
			}

			// remove the vertices in related tables
			auto it_biase_rank_iter = degree_biase_iter_by_v.find(v);
			//note you cannot swap the order of the following two steps
			weighed_in_out_degree_biases.erase(it_biase_rank_iter->second);
			degree_biase_iter_by_v.erase(it_biase_rank_iter);
		}

	private:
		float get_biase(DgV v)
		{
			auto biase_iter = degree_biase_iter_by_v.find(v)->second;
			return biase_iter->first;
		}

		bool update_biase(DgV v, float biase)
		{
			auto iter_biase_rank_by_v = degree_biase_iter_by_v.find(v);
			if (iter_biase_rank_by_v == degree_biase_iter_by_v.end()) {
				return false;
			}

			auto iter_biase_rank = iter_biase_rank_by_v->second;

			//update the key
			auto node = weighed_in_out_degree_biases.extract(iter_biase_rank);
			node.key() = biase;
			auto new_iter_v_to_biase_rank_iter = weighed_in_out_degree_biases.insert(std::move(node));

			//update the iterator table
			degree_biase_iter_by_v.erase(iter_biase_rank_by_v);
			degree_biase_iter_by_v.emplace(v, new_iter_v_to_biase_rank_iter);
			return true;
		}
		std::pair<std::unordered_set<DgV>, std::unordered_set<DgV>>
			getInitSourceAndSinks(const Dg& dg)
		{
			std::unordered_set<DgV> sources;
			std::unordered_set<DgV> sinks;
			for (auto v : Range{ boost::vertices(dg) }) {
				if (boost::in_degree(v, dg) == 0) {
					sources.insert(v);
				}
				else if (boost::out_degree(v, dg) == 0) {
					sinks.insert(v);
				}
			}
			return { std::move(sources), std::move(sinks) };
		}
	}

	static std::pair<std::set<DgE>, std::vector<DgV>>
		minFeedbackArcSet(Dg& dg)
	{
		std::set<DgE> feedbackEdges;

		//sink sequence, so called "s2" in the paper, in the reverse order than the paper
		std::vector<DgV> sSinks;
		sSinks.reserve(boost::num_vertices(dg));

		//anything other than s2, so called "s1" in the paper
		std::vector<DgV> sNotSinks;
		sNotSinks.reserve(boost::num_vertices(dg));

		// Some Preparations
		auto boostVertices = boost::vertices(dg);
		std::set<DgV> vLeft(boostVertices.first, boostVertices.second);
		for (auto v : vLeft) {
			std::cout << dg[v].name << ": " << v << std::endl;
		}
		std::cout << "\n\n"
			<< std::endl;
		InOutInfo inoutInfo;
		WeighedDegreeInfo wDegInfo;
		for (auto v : vLeft) {
			inoutInfo.add(v, dg);
		}

		for (auto v : vLeft) {
			wDegInfo.add(v, inoutInfo, dg);
		}

		auto[sources, sinks] = getInitSourceAndSinks(dg);

		int iLoop = 0;
		while (!vLeft.empty()) {
			for (auto src : sources) {
				sNotSinks.emplace_back(src);
				wDegInfo.remove(src, inoutInfo, dg);
				inoutInfo.remove(src, dg, false);
				vLeft.erase(src);
				auto name = dg[src].name;
				std::cout << "Remove source: " << name << ". " << std::endl;
			}
			sources.clear();

			while (!sinks.empty()) {
				auto sink = *sinks.begin();
				sSinks.emplace_back(sink);
				wDegInfo.remove(sink, inoutInfo, dg);
				auto[nonOrphanEdges, newSinks] = inoutInfo.remove(sink, dg, false);
				sinks.insert(newSinks.begin(), newSinks.end());
				vLeft.erase(sink);
				sinks.erase(sink);
				auto name = dg[sink].name;
				std::cout << "Remove sink: " << name << ". " << std::endl;
			}

			auto topBiasedVCond = wDegInfo.mostBiasedV();
			if (!topBiasedVCond) {
				continue;
			}
			auto topBiasedV = topBiasedVCond.value();

			wDegInfo.remove(topBiasedV, inoutInfo, dg);
			auto name = dg[topBiasedV].name;
			std::cout << "Remove top biase: " << name << ". " << std::endl;
			auto[nonOrphanEdges, newSinks] = inoutInfo.remove(topBiasedV, dg, true);
			feedbackEdges.insert(nonOrphanEdges.begin(), nonOrphanEdges.end());
			sinks.insert(newSinks.begin(), newSinks.end());
			vLeft.erase(topBiasedV);
			sNotSinks.emplace_back(topBiasedV);
			iLoop++;
		}

		std::copy(sSinks.rbegin(), sSinks.rend(), std::back_inserter(sNotSinks));
		return { feedbackEdges, sNotSinks };
	}
}