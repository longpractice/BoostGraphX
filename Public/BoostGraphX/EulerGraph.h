#pragma once
#include <BoostGraphX/Common.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <utility>
#include <vector>

namespace bglx::detail
{
	template <typename VertexListGraph>
	struct Euler_Digraph_Helper
	{
		using V = typename boost::graph_traits<VertexListGraph >::vertex_descriptor;
		using E = typename boost::graph_traits<VertexListGraph>::edge_descriptor;

		using Directed_Category = typename boost::graph_traits<VertexListGraph>::directed_category;
		using Is_Directed = std::is_convertible<Directed_Category, boost::directed_tag>;

		//
		// The algorithms of finding a directed euler cycle and a undirected euler cycle are very similar, both could
		// be solved by Hierholz algorithm.
		//
		// However, the issue is that for the directed case, marking the edge as visited could be very effitient
		// it is more tricky for the undirected case.
		//
		static_assert(Is_Directed::value,
			"Calculating euler directed cycle and directed tour requires a directed graph.");

		//
		// We want to save the vertex descriptor inside the vertex property(a linked-list iterator).
		//
		// However, the vertex descriptor type depends on the graph type which depends on the vertex
		// property type.
		//
		// We therefore use this Dag_Aux_Pure::vertex_descriptor independent on vertex property
		// Indeed, the vertex descriptor is just a size_t here. Not that a big concern.
		//
		using Dag_Aux_Pure = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
		using V_Aux = typename Dag_Aux_Pure::vertex_descriptor;

		using V_Aux_List = std::list<V_Aux>;
		using V_Aux_List_It = typename V_Aux_List::iterator;

		using E_List = std::list<E>;
		using E_List_It = typename std::list<E>::iterator;


		struct E_Aux_Prop
		{
			E e_origin;
			// we could not use raw iterator here, since once it is default constructed
			// it is a "singular" iterator that we could not even assign to it
			std::optional<E_List_It> e_walk_list_it;
		};

		struct V_Aux_Prop {
			V v_origin;
			size_t first_unfinished_edge_id{ 0 };
			// we could not use raw iterator here, since once it is default constructed
			// it is a "singular" iterator that we could not even assign to it
			std::optional<V_Aux_List_It> v_unfinished_list_it;
		};

		using G_Aux = boost::adjacency_list<
			boost::vecS,
			boost::vecS,
			boost::directedS,
			V_Aux_Prop,
			E_Aux_Prop>;

		using E_Aux = typename G_Aux::edge_descriptor;

		//
		// Find from the unvisited edges a closed walk starting and ending at v
		//
		// You could possibly specify a mark edge here.
		// The function returns the start edge of the new circle from v,
		// The bool part of the return pair is true if the start edge of the
		// closed walk corresponds to mark_edge
		//
		static std::pair<E_List_It, bool> discover_new_close_walk(
			V_Aux v,
			G_Aux& aux,
			V_Aux_List& unfinished_v_list,
			E_List& old_closed_walk,
			std::optional<E_Aux> mark_edge = {})
		{
			E_List new_closed_walk;
			auto v_first = v;
			if (boost::out_degree(v, aux) == aux[v].first_unfinished_edge_id)
			{
				//we could not find any exit and thus no new closed walk from this v
				return { {}, false };
			}
			std::unique_ptr<E_List_It> prev_visited_out_edge;
			auto out_edges = boost::out_edges(v, aux);
			auto e_first_unvisited = *(out_edges.first + aux[v].first_unfinished_edge_id);
			bool is_marked = mark_edge ? e_first_unvisited == mark_edge.value() : false;

			if (aux[v].first_unfinished_edge_id != 0) {
				const auto& e_visited = *(out_edges.first + aux[v].first_unfinished_edge_id - 1);
				const auto& it_e_visited = aux[e_visited].e_walk_list_it;
				prev_visited_out_edge = std::make_unique<E_List_It>(*it_e_visited);
			}

			do {
				//find the first exit for current vertex
				auto& v_prop = aux[v];
				auto out_degree = boost::out_degree(v, aux);
				auto it_e_out = boost::out_edges(v, aux).first + v_prop.first_unfinished_edge_id;
				auto e_out = *it_e_out;
				//this exist edge is now visited, move it from unvisited to visited
				++v_prop.first_unfinished_edge_id;

				new_closed_walk.emplace_back(aux[e_out].e_origin);
				aux[e_out].e_walk_list_it = { --new_closed_walk.end() };

				//v was partially visited before, if we have finished the last
				bool just_become_unfinished = v_prop.first_unfinished_edge_id == 1 && v_prop.first_unfinished_edge_id != out_degree;
				if (just_become_unfinished) {
					unfinished_v_list.emplace_back(v);
					v_prop.v_unfinished_list_it = { --unfinished_v_list.end() };
				}

				//it was partially visited before and now need to be moved out
				bool justMovedOutOfUnfinishedVList = v_prop.first_unfinished_edge_id == out_degree && out_degree > 1;
				if (justMovedOutOfUnfinishedVList) {
					unfinished_v_list.erase(*v_prop.v_unfinished_list_it);
				}
				v = boost::target(e_out, aux);
			} while (v != v_first);

			auto e_list_it = new_closed_walk.begin();
			if (prev_visited_out_edge) {
				//the visited out edge, which need to break the previous linked list
				old_closed_walk.splice(*prev_visited_out_edge, new_closed_walk);
			}
			else {
				//old closed walk must be empty
				std::swap(old_closed_walk, new_closed_walk);
			}
			return { e_list_it, is_marked };
		}

		struct V_Copier {
			V_Copier(const VertexListGraph& from, G_Aux& to)
				: from{ from }, to{ to } {};
			const VertexListGraph& from;
			G_Aux& to;

			void operator()(V input, V_Aux output) const
			{
				to[output].v_origin = input;
			}
		};

		struct E_Copier {
			E_Copier(const VertexListGraph& from, G_Aux& to)
				: from{ from }, to{ to } {};

			const VertexListGraph& from;
			G_Aux& to;

			void operator()(E input, E_Aux output) const
			{
				to[output].e_origin = input;
			}
		};


		//find euler cycle using Hierholzer algorithm, with out ensuring any edge being the
		//first edge
		static std::list<E> find_euler_cycle_hierholzer(G_Aux& g_aux, V_Aux first_v)
		{
			V_Aux_List unfinished_v_list;
			E_List closed_walk;
			auto nr_edges = boost::num_edges(g_aux);
			discover_new_close_walk(first_v, g_aux, unfinished_v_list, closed_walk, {});
			while (closed_walk.size() < nr_edges) {
				auto v = unfinished_v_list.front();
				discover_new_close_walk(v, g_aux, unfinished_v_list, closed_walk, {});
			}
			return closed_walk;
		}

		//
		// find euler cycle start and end at firstV
		// e_ensure_first will be the first edge in the cycle
		// This is useful for algorithm for Eulerian trail (or Eulerian path) with
		// the start vertex and the end vertex different
		//
		static std::list<E> find_euler_cycle_hierholzer__ensure_first_edge(
			G_Aux& g_aux, V_Aux v_first, E_Aux e_ensure_first)
		{
			V_Aux_List unfinished_v_list;
			E_List closed_walk;
			auto e_marked_it = closed_walk.end();
			auto nr_edges = boost::num_edges(g_aux);
			while (closed_walk.size() < nr_edges) {
				if (unfinished_v_list.empty() || unfinished_v_list.front() == v_first) {
					auto[e_it_first_added, if_marked] = discover_new_close_walk(v_first, g_aux, unfinished_v_list, closed_walk, { e_ensure_first });
					if (if_marked) {
						e_marked_it = e_it_first_added;
					}
				}
				else {
					auto v = unfinished_v_list.front();
					discover_new_close_walk(v, g_aux, unfinished_v_list, closed_walk, {});
				}
			}
			E_List closed_walk_reordered;
			closed_walk_reordered.splice(closed_walk_reordered.begin(), closed_walk, e_marked_it, closed_walk.end());
			closed_walk_reordered.splice(closed_walk_reordered.end(), closed_walk);
			return closed_walk_reordered;
		}

		static std::list<E> find_euler_tour(G_Aux& gAux, V_Aux firstV, V_Aux lastV)
		{
			//First add a edge going from lastV to firstV and
			//Then find an euler cycle starting and ending at lastV
			//After erasing the first edge, we have a euler tour from firstV -> lastV
			auto[e, ok] = boost::add_edge(lastV, firstV, gAux);
			auto edges = find_euler_cycle_hierholzer__ensure_first_edge(gAux, lastV, e);
			//erase the first edge
			edges.erase(edges.begin());
			return edges;
		}
	};
}

namespace bglx {

	template <typename VertexListGraph>
	std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
		find_one_directed_euler_cycle_hierholzer(
			const VertexListGraph& g,
			typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_cycle_start_vertex)
	{
		using GHelper = detail::Euler_Digraph_Helper<VertexListGraph>;
		typename GHelper::G_Aux g_aux;
		typename GHelper::V_Aux v_aux_start;

		//while coping, record the vAux correspondent to euler_cycle_start_vertex
		auto v_marked__copier =
			[&](typename GHelper::V input, typename GHelper::V_Aux output) {
			g_aux[output].v_origin = input;
			if (input == euler_cycle_start_vertex) {
				v_aux_start = output;
			}
		};
		typename GHelper::E_Copier e_copier(g, g_aux);
		boost::copy_graph(
			g,
			g_aux,
			boost::vertex_copy(v_marked__copier).
			edge_copy(e_copier)
		);
		return GHelper::find_euler_cycle_hierholzer(g_aux, v_aux_start);
	}

	template <typename VertexListGraph, typename VertexIndexMap>
	std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
		find_one_directed_euler_cycle_hierholzer(
			const VertexListGraph& g,
			typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_cycle_start_vertex,
			const VertexIndexMap& i_map)
	{
		using GHelper = detail::Euler_Digraph_Helper<VertexListGraph>;
		typename GHelper::G_Aux g_aux;
		typename GHelper::V_Aux v_aux_start;

		//while coping, record the vAux correspondent to euler_cycle_start_vertex
		auto v_marked__copier =
			[&](typename GHelper::V input, typename GHelper::V_Aux output) {
			g_aux[output].v_origin = input;
			if (input == euler_cycle_start_vertex) {
				v_aux_start = output;
			}
		};
		typename GHelper::E_Copier e_copier(g, g_aux);
		boost::copy_graph(
			g,
			g_aux,
			boost::vertex_copy(v_marked__copier).edge_copy(e_copier).vertex_index_map(i_map));
		return GHelper::find_euler_cycle_hierholzer(g_aux, v_aux_start);
	}

	template <typename VertexListGraph>
	std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
		find_one_directed_euler_tour_hierholzer(const VertexListGraph& g,
			typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_start_vertex,
			typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_end_vertex)
	{
		//should use euler cycle functions if the start and tour equals
		assert(euler_tour_start_vertex != euler_tour_end_vertex);
		using GHelper = detail::Euler_Digraph_Helper<VertexListGraph>;
		typename GHelper::G_Aux g_aux;
		typename GHelper::V_Aux v_aux_start;
		typename GHelper::V_Aux v_aux_last;

		auto v_marked__copier =
			[&](typename GHelper::V input, typename GHelper::V_Aux output) {
			g_aux[output].v_origin = input;
			if (input == euler_tour_start_vertex) {
				v_aux_start = output;
			}
			else if (input == euler_tour_end_vertex) {
				v_aux_last = output;
			}
		};
		typename GHelper::E_Copier e_copier(g, g_aux);
		boost::copy_graph(
			g, g_aux,
			boost::vertex_copy(v_marked__copier).edge_copy(e_copier));

		return GHelper::find_euler_tour(g_aux, v_aux_start, v_aux_last);
	}

	template <typename VertexListGraph, typename VertexIndexMap>
	std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
		find_one_directed_euler_tour_hierholzer(const VertexListGraph& g,
			typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_start_vertex,
			typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_tour_end_vertex,
			const VertexIndexMap& i_map)
	{
		//should use euler cycle functions if the start and tour equals
		assert(euler_tour_start_vertex != euler_tour_end_vertex);
		using GHelper = detail::Euler_Digraph_Helper<VertexListGraph>;
		typename GHelper::G_Aux g_aux;
		typename GHelper::V_Aux v_aux_start;
		typename GHelper::V_Aux v_aux_last;

		auto v_marked__copier =
			[&](typename GHelper::V input, typename GHelper::V_Aux output) {
			g_aux[output].v_origin = input;
			if (input == euler_tour_start_vertex)
			{
				v_aux_start = output;
			}
			else if (input == euler_tour_end_vertex)
			{
				v_aux_last = output;
			}
		};
		typename GHelper::E_Copier e_copier(g, g_aux);
		boost::copy_graph
		(
			g, g_aux,
			boost::vertex_copy(v_marked__copier).edge_copy(e_copier).vertex_index_map(i_map)
		);

		return GHelper::find_euler_tour(g_aux, v_aux_start, v_aux_last);
	}
}
