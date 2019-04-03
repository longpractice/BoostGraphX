#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <utility>
#include <vector>
#include <queue>
#include <stack>

namespace bglx
{
	/*
	 * Implements:
	 *
	 * Eades P, Lin X, Smyth WF.
	 * A fast and effective heuristic for the feedback arc set problem.
	 * Information Processing Letters. 1993 Oct 18;47(6):319-23.
	 *
	 * This algorithm works on non-weighted directed graph. It also gives
	 * out the topologically sorted vertices (for a directed acyclic graph without the
	 * computed feedback arc set).
	 *
	 * Usage example:
	 *
	 * ```cpp
	 * 	using Dag = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>;
	 *	typedef boost::erdos_renyi_iterator<boost::minstd_rand, Dag> ERGen;
	 *	boost::minstd_rand gen;
	 *	//Create graph with 100 nodes and edges with probability 0.05
	 *	Dag g(ERGen(gen, 100, 0.05), ERGen(), 100);
	 *	//all the computing happens during construction
	 *	Minimum_Feedback_Arc_Set_Solver_Eades_Lin<Dag> solver(g);
	 *	const auto& feedback_arc_set = solver.feedback_arc_set();
	 *	const auto& topologically_sorted_vertices = solver.topological_sorted_vertices();
	 * ```
	 *
	 * Linear time and space complexity guaranteed.
	 *
	 * Number of edges E and number of vertices V, the feedback edges count is no greater than
	 * E/2 - V/6.
	 *
	 * This therefore works well on sparse graphs.
	 *
	 * The input graph must be directed and model
	 * a VertexListGraphConcept and an IncidenceGraphConcept
	 *
	 * Internally, there are two overloads of the constructor of the solver
	 * The two overloads serve the purpose of:
	 * The first one above is suitable for a graph type with an internal vertex_index property.
	 * For example,  adjacency_list with VertexList=vecS has this property.
	 *
	 * The second one is suitable for a graph type without an internal vertex_index property.
	 * For example, adjacency_list with VertexList=listS or setS.
	 * The vertex index map type of VertexIndexMap must be a model of Readable Property Map
	 * and must map the vertex descriptors of G to the integers in the half-open range [0,num_vertices(G)).
	 * This is normally adapted from a key-value container like std::map,
	 * std::unordered_map using boost::associative_property_map.
	 *
	 */

	template <typename Graph>
	class Minimum_Feedback_Arc_Set_Solver_Eades_Lin
	{
	public:
		BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));

		BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<Graph>));

		static_assert(
			std::is_convertible
			<
			typename boost::graph_traits<Graph>::directed_category,
			boost::directed_tag
			>
			::value
			);

		Minimum_Feedback_Arc_Set_Solver_Eades_Lin(const Graph& g_origin)
		{
			copy_graph(g_origin);
			init(g_origin);
			process();
		}

		template <typename VertexIndexMap>
		Minimum_Feedback_Arc_Set_Solver_Eades_Lin(const Graph& g_origin, const VertexIndexMap& i_map)
		{
			copy_graph(g_origin, i_map);
			init(g_origin);
			process();
		}

		/// read out a list of pairs of vertices, each pair is the source and the target of the feedback arc set
		std::list<std::pair<typename Graph::vertex_descriptor, typename Graph::vertex_descriptor>>&
			feedback_arc_set()
		{
			return feedback_edges_res;
		}

		/// read out a list of vertices topologically sorted without the feedback edges
		std::list<typename Graph::vertex_descriptor>&
			topological_sorted_vertices()
		{
			return s1;
		}


	private:
		using E_Origin = typename Graph::edge_descriptor;
		using V_Origin = typename Graph::vertex_descriptor;


		//just to break the dependency cycle
		struct E_Prop
		{
			size_t edge_id;
		};

		using G_Pure = boost::adjacency_list<
			boost::listS,
			boost::listS,
			boost::directedS,
			boost::no_property,
			E_Prop>;

		using V_Pure = typename G_Pure::vertex_descriptor; //actually just void*
		using E_Pure = typename G_Pure::edge_descriptor;

		using V_List_Pure_It = typename std::list<V_Pure>::iterator;
		//
		//the in-edges are listed in each vertex as a vector of edge ids(index of vector of member edges)
		//a fast and dirty solution would be to use the bidirectionalS for the auxiliary graph(member g)
		//however, removing a perticular in-edge from the in-edge list of a vertex is not constant time
		//(it can be O(E/V) or O(log(E/V)) depending on edge-list container type)
		//
		// For very dense graphs, the number of edges per vertex can be large, we still want this to be 
		// constant time, the solution is to maintain it ourselvse by saving edge-ids in vectors.
		// 
		// Then how could we linear time erase in-edge? We do a swap and pop on the vector. 
		// As it turns out, it is even faster than maintaining a std::list for in-edges and 
		// pass in the iterator for erasing, due to the cache friendliness of std::vector compared to std::list
		//
		using In_Edges = std::vector<size_t>;

		struct V_Prop
		{
			V_Origin v_origin;
			//iterator into the vertices list in the bins
			//this facilitates constant-time moving vertex from one bin to another
			//by splicing the iterator
			V_List_Pure_It v_list_it;
			int64_t deg_diff{};
			In_Edges in_edges;
		};


		///////////////////////////////////////
		////  The graph //////////////////////
		///////////////////////////////////////
		using G = boost::adjacency_list<
			boost::listS,
			boost::listS,
			boost::directedS,
			V_Prop,
			E_Prop
		>;

		G g;

		using V = typename G::vertex_descriptor;
		using E = typename G::edge_descriptor;

		/////////////////Edge list/////////////////////////////////////////
		//extra edge info for quickly looking up the out-edge-iterator when we are 
		//clearing the in-edges for a certain vertex, or the other way around
		struct Edge_Info
		{
			size_t in_edge_id;
			//keep this one, when passed to boost::remove_edge, it is constant time
			typename G::out_edge_iterator out_edge_it_boost;
		};

		std::vector<Edge_Info> edges;

		size_t in_degree(V v)
		{
			const auto& v_prop = g[v];
			return v_prop.in_edges.size();
		}

		void remove_in_edge(V v, size_t in_edge_id)
		{
			auto& in_edges = g[v].in_edges;
			//swap and pop to erase an element in a vector
			in_edges[in_edge_id] = in_edges.back();
			in_edges.pop_back();
			if (in_edge_id != in_edges.size())
			{
				auto end_e_id = in_edges[in_edge_id];
				edges[end_e_id].in_edge_id = in_edge_id;
			}
		}

		void clear_vertex(V v)
		{
			auto out_edges = boost::out_edges(v, g);

			for (auto e_it = out_edges.first; e_it != out_edges.second; ++e_it)
			{
				const auto& e = *e_it;
				const auto& e_prop = g[e];
				auto e_id = e_prop.edge_id;
				const auto& e_info = edges[e_id];
				auto tgt = boost::target(e, g);
				remove_in_edge(tgt, e_info.in_edge_id);
			}
			boost::clear_out_edges(v, g);

			//handle in edges differently
			for (auto e_id : g[v].in_edges)
			{
				boost::remove_edge(edges[e_id].out_edge_it_boost, g);
			}
		}

		/////////////////////////////////////////////////////////////////////////////


		////////////////////////////// Bins related //////////////////////////
		//the bin as stated in the paper
		//since we need to be able to find the largest bin that has vertices in it
		//we also embed a doubly linked list to maintain all the bins with vertices
		struct Bin
		{
			std::list<V> v_list;
			Bin* p_prev_non_empty_bin{ nullptr };
			Bin* p_next_non_empty_bin{ nullptr };
		};

		//points to the last element of the linked list of non-empty bins
		Bin* p_last_non_empty_bin{ nullptr };
		//all the bins, no matter have vertices or not
		std::vector<Bin> bins;

		//max_out_deg + max_in_deg would be maximum number of bins
		//out_deg - in_deg is within range of -max_in_deg and max_out_deg
		int64_t max_out_deg{};
		int64_t max_in_deg{};

		void init_v_to_bin(V v)
		{
			const auto deg_in = in_degree(v);
			const auto deg_out = boost::out_degree(v, g);

			if (deg_in == 0)
			{
				unprocessed_sources.push(v);
			}
			else if (deg_out == 0)
			{
				unprocessed_sinks.push(v);
			}
			else
			{
				int64_t deg_diff = static_cast<int64_t>(deg_out) - static_cast<int64_t>(deg_in);
				auto& bin = get_bin(deg_diff);
				bin.v_list.emplace_back(v);

				auto& vProp = g[v];
				vProp.deg_diff = deg_diff;
				vProp.v_list_it = std::prev(bin.v_list.end());
			}
		}

		void init_occupied_bin_list_prev()
		{
			Bin* current_last_occupied_bin = nullptr;
			for (auto& bin : bins)
			{
				if (bin.v_list.empty())
				{
					continue;
				}
				bin.p_prev_non_empty_bin = current_last_occupied_bin;
				current_last_occupied_bin = &bin;
			}
			p_last_non_empty_bin = current_last_occupied_bin;
		}

		void init_occupied_bin_list_next()
		{
			Bin* current_first_occupied_bin = nullptr;
			for (auto itBin = bins.rbegin(); itBin != bins.rend(); ++itBin)
			{
				auto& bin = *itBin;
				if (bin.v_list.empty())
				{
					continue;
				}
				bin.p_next_non_empty_bin = current_first_occupied_bin;
				current_first_occupied_bin = &bin;
			}
		}


		std::queue<V> unprocessed_sources;
		std::stack<V> unprocessed_sinks;

		std::list<V_Origin> s1;
		std::list<V_Origin> s2;
		using Feedback_Arc_List = std::list<std::pair<V_Origin, V_Origin>>;
		Feedback_Arc_List feedback_edges_res;

		void copy_graph(const Graph& g_origin)
		{
			auto v_copier = [&](const V_Origin v_origin, const V v)
			{
				g[v].v_origin = v_origin;
			};

			auto e_copier = [&](const E_Origin& e_origin, const E& e)
			{};

			boost::copy_graph(g_origin, g, boost::vertex_copy(v_copier).edge_copy(e_copier));
		}

		template <typename VertexIndexMap>
		void copy_graph(const Graph& g_origin, const VertexIndexMap& i_map)
		{
			auto v_copier = [&](const V_Origin v_origin, const V v)
			{
				g[v].v_origin = v_origin;
			};

			auto e_copier = [&](const E_Origin& e_origin, const E& e)
			{};
			boost::copy_graph(g_origin, g, boost::vertex_copy(v_copier).edge_copy(e_copier).vertex_index_map(i_map));
		}

		void init(const Graph& g_origin)
		{
			const auto nr_edges = boost::num_edges(g_origin);
			edges.reserve(nr_edges);

			max_out_deg = std::numeric_limits<int64_t>::min();
			max_in_deg = std::numeric_limits<int64_t>::min();
			for (auto v : boost::make_iterator_range(boost::vertices(g)))
			{
				const auto in_deg = static_cast<int64_t>(in_degree(v));
				const auto out_deg = static_cast<int64_t>(boost::out_degree(v, g));
				if (in_deg > max_in_deg)
				{
					max_in_deg = in_deg;
				}
				if (out_deg > max_out_deg)
				{
					max_out_deg = out_deg;
				}

				auto out_edges = boost::out_edges(v, g);
				for (auto boost_e_it = out_edges.first; boost_e_it != out_edges.second; ++boost_e_it)
				{
					const auto& e = *boost_e_it;
					auto e_id = edges.size();
					auto tgt = boost::target(e, g);
					auto& tgt_in_edges = g[tgt].in_edges;
					tgt_in_edges.emplace_back(e_id);
					const auto e_in_edge_id = tgt_in_edges.size() - 1;
					edges.push_back(Edge_Info{ e_in_edge_id, boost_e_it });
					g[e].edge_id = e_id;
				}
			}

			bins.resize(max_in_deg + max_out_deg + 1);

			for (auto v : boost::make_iterator_range(boost::vertices(g)))
			{
				init_v_to_bin(v);
			}

			init_occupied_bin_list_prev();
			init_occupied_bin_list_next();
		}

		void process_sources()
		{
			while (!unprocessed_sources.empty())
			{
				const auto src = unprocessed_sources.front();
				unprocessed_sources.pop();
				for (auto tgt : boost::make_iterator_range(boost::adjacent_vertices(src, g)))
				{
					decrease_in_degree_on_v(tgt);
				}
				add_to_s1(src);
				clear_vertex(src);
			}
		}

		void process_sinks()
		{
			while (!unprocessed_sinks.empty())
			{
				const auto tgt = unprocessed_sinks.top();
				unprocessed_sinks.pop();

				for (auto e_id : g[tgt].in_edges)
				{
					decrease_out_degree_on_v(boost::source(*edges[e_id].out_edge_it_boost, g));
				}
				add_to_s2(tgt);
				clear_vertex(tgt);
			}
		}

		void process_max_delta_vertex(V v)
		{
			auto deg_diff = g[v].deg_diff;
			auto& bin = get_bin(deg_diff);
			bin.v_list.erase(g[v].v_list_it);
			for (auto tgt : boost::make_iterator_range(boost::adjacent_vertices(v, g)))
			{
				decrease_in_degree_on_v(tgt);
			}

			for (auto e_id : g[v].in_edges)
			{
				auto src = boost::source(*edges[e_id].out_edge_it_boost, g);
				feedback_edges_res.emplace_back(std::pair<V_Origin, V_Origin>{g[src].v_origin, g[v].v_origin});
				decrease_out_degree_on_v(src);
			}

			add_to_s1(v);

			if (bin.v_list.empty())
			{
				handle_emptied_bin(bin);
			}

			clear_vertex(v);
		}

		void process()
		{
			if (!p_last_non_empty_bin)
			{
				//empty graph
				return;
			}
			auto nrVertices = boost::num_vertices(g);
			while (s1.size() + s2.size() != nrVertices)
			{
				while (!(unprocessed_sources.empty() && unprocessed_sinks.empty()))
				{
					process_sources();
					process_sinks();
				}
				if (p_last_non_empty_bin)
				{
					auto v = p_last_non_empty_bin->v_list.back();
					process_max_delta_vertex(v);
				}
			}
			//s1 is the result
			s1.splice(s1.end(), s2);
		}

		int64_t get_bin_id(int64_t deg_diff) const
		{
			return deg_diff + max_in_deg;
		}

		Bin& get_bin(int64_t deg_diff)
		{
			auto id = get_bin_id(deg_diff);
			return bins[id];
		}


		//non-sinks: emplace to the back of s1
		void add_to_s1(V v)
		{
			//new sources are attached to the back
			s1.emplace_back(g[v].v_origin);
		}

		//sinks: emplace to the front of s2
		void add_to_s2(V v)
		{
			//new sinks are attached to the front
			s2.emplace_front(g[v].v_origin);
		}


		void handle_emptied_bin(Bin& bin)
		{
			//this bin has become empty, link the two neighbors of the linked list
			if (bin.p_prev_non_empty_bin)
			{
				bin.p_prev_non_empty_bin->p_next_non_empty_bin = bin.p_next_non_empty_bin;
			}

			if (bin.p_next_non_empty_bin)
			{
				bin.p_next_non_empty_bin->p_prev_non_empty_bin = bin.p_prev_non_empty_bin;
			}
			else
			{
				p_last_non_empty_bin = bin.p_prev_non_empty_bin;
			}
		}

		//this is due to the one of the in-edges is removed
		void handle_newly_occupied_bin_downgrade_v(Bin& bin, Bin& next_occupied_bin)
		{
			//
			// From linked list A<--->B to A <---> C <---> B
			// Two tasks:
			//
			// 1) set the next and prev of C
			// 2) set the next of A and the prev of B
			//
			// the new one was in chain since it might be just emptied but not yet removed from the chain
			if (next_occupied_bin.p_prev_non_empty_bin == &bin)
				return;
			bin.p_next_non_empty_bin = &next_occupied_bin;
			bin.p_prev_non_empty_bin = next_occupied_bin.p_prev_non_empty_bin;

			if (next_occupied_bin.p_prev_non_empty_bin)
			{
				next_occupied_bin.p_prev_non_empty_bin->p_next_non_empty_bin = &bin;
			}
			next_occupied_bin.p_prev_non_empty_bin = &bin;
		}

		void handle_newly_occupied_bin_upgrade_v(Bin& bin, Bin& prev_occupied_bin)
		{
			//the new one was in chain since it might be just emptied but not yet removed from the chain
			if (prev_occupied_bin.p_next_non_empty_bin == &bin)
				return;
			bin.p_next_non_empty_bin = prev_occupied_bin.p_next_non_empty_bin;
			bin.p_prev_non_empty_bin = &prev_occupied_bin;
			if (prev_occupied_bin.p_next_non_empty_bin)
			{
				prev_occupied_bin.p_next_non_empty_bin->p_prev_non_empty_bin = &bin;
			}
			else
			{
				p_last_non_empty_bin = &bin;
			}
			prev_occupied_bin.p_next_non_empty_bin = &bin;
		}

		//if return true, the caller needs to remove the vertex from graph also
		void decrease_in_degree_on_v(V v)
		{
			auto& v_prop = g[v];
			auto old_deg_diff = v_prop.deg_diff;
			auto old_in_deg = in_degree(v);
			auto old_out_deg = boost::out_degree(v, g);

			if (old_out_deg == 0)
			{
				//treat as sinks
				return;
			}

			auto& old_bin = get_bin(old_deg_diff);
			//transfer the vertex to his new home, s1 or new bin
			if (old_in_deg == 1)
			{
				//becomes a source
				if (v_prop.v_list_it != old_bin.v_list.end())
				{
					old_bin.v_list.erase(v_prop.v_list_it);
					v_prop.v_list_it = old_bin.v_list.end();
					unprocessed_sources.push(v);
				}
			}
			else
			{
				const auto new_deg_diff = old_deg_diff + 1;
				v_prop.deg_diff = new_deg_diff;
				//upgrade to another bin
				auto& new_bin = get_bin(new_deg_diff);
				const bool was_new_bin_empty = new_bin.v_list.empty();
				//splice from old bin to new bin, which maintains the iterator validation
				new_bin.v_list.splice(new_bin.v_list.end(), old_bin.v_list, v_prop.v_list_it);
				if (was_new_bin_empty)
				{
					handle_newly_occupied_bin_upgrade_v(new_bin, old_bin);
				}
			}
			if (old_bin.v_list.empty())
			{
				handle_emptied_bin(old_bin);
			}
		}

		//if return true, the caller need to remove the vertex from the graph
		void decrease_out_degree_on_v(V v)
		{
			auto& v_prop = g[v];
			auto old_deg_diff = v_prop.deg_diff;
			auto old_out_deg = boost::out_degree(v, g);
			auto old_in_deg = in_degree(v);

			if (old_out_deg > 0 && old_in_deg == 0)
			{
				//must be classified as a source before
				return;
			}

			auto& old_bin = get_bin(old_deg_diff);

			//transfer the vertex to his new home, s2 or new bin
			if (old_out_deg == 1)
			{
				//becomes a sink
				if (v_prop.v_list_it != old_bin.v_list.end())
				{
					old_bin.v_list.erase(v_prop.v_list_it);
					v_prop.v_list_it = old_bin.v_list.end();
					unprocessed_sinks.push(v);
				}
			}
			else
			{
				auto new_deg_diff = old_deg_diff - 1;
				v_prop.deg_diff = new_deg_diff;
				//downgrade the vertex to lower bin
				auto& new_bin = get_bin(new_deg_diff);
				bool was_new_bin_empty = new_bin.v_list.empty();
				new_bin.v_list.splice(new_bin.v_list.end(), old_bin.v_list, v_prop.v_list_it);
				if (was_new_bin_empty)
				{
					handle_newly_occupied_bin_downgrade_v(new_bin, old_bin);
				}
			}
			if (old_bin.v_list.empty())
			{
				handle_emptied_bin(old_bin);
			}
		}
	};
}
