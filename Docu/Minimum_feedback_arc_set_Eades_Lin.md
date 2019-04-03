This method is based on:

Eades P, Lin X, Smyth WF.
A fast and effective heuristic for the feedback arc set problem.
Information Processing Letters. 1993 Oct 18;47(6):319-23.

This algorithm works on non-weighted directed graph. It also gives out the topologically sorted vertices (for the directed acyclic graph removing the computed feedback arc set).

Usage example:

```cpp
using Dag = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>;

/*
 * Create a random graph with with 300,000 vertices and 2,700,000 edges
 * Using Erdos-Renyi model
 */
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Dag> ERGen;
boost::minstd_rand gen;
//Create graph with 100 nodes and edges with probability 0.05
Dag g(ERGen(gen, 300000, 0.00003), ERGen(), 100);


//all the computing happens during construction
Minimum_Feedback_Arc_Set_Solver_Eades_Lin<Dag> solver(g);
const auto& feedback_arc_set = solver.feedback_arc_set();
const auto& topologically_sorted_vertices = solver.topological_sorted_vertices();
```

Linear time and space complexity guaranteed. In a Erdos-Renyi model graph as shown in the example with 300,000 vertices and 2,700,000 edges. It takes 3.2s with Intel Core i7-7700K 4.2GHz and 16Gb ram.

Number of edges E and number of vertices V, the feedback edges count is no greater than E/2 - V/6.

This therefore works quite well on sparse graphs.

The input graph must be directed and model VertexListGraphConcept and IncidenceGraphConcept.

The computation is done during the construction of the solver.

---

There are two overloads of the constructor of the solver.


```cpp
Minimum_Feedback_Arc_Set_Solver_Eades_Lin::
	Minimum_Feedback_Arc_Set_Solver_Eades_Lin(const Graph& g_origin)
```
The first one above is suitable for a graph type with an internal vertex_index property. For example, adjacency_list with VertexList=vecS has this property.



```cpp
Minimum_Feedback_Arc_Set_Solver_Eades_Lin::
	template <typename VertexIndexMap>
	Minimum_Feedback_Arc_Set_Solver_Eades_Lin(
		const Graph& g_origin, 
		const VertexIndexMap& i_map
	)
```
The second one above is suitable for a graph type without an internal vertex_index property. For example, adjacency_list with VertexList=listS or setS.
The vertex index map type of VertexIndexMap must be a model of Readable Property Map
and must map the vertex descriptors of G to the integers in the half-open range [0,num_vertices(G)).
This is normally adapted from a key-value container like std::map and std::unordered_map using boost::associative_property_map.

For example:

```cpp
using Dag_VSet = boost::adjacency_list<boost::setS, boost::setS, boost::directedS>;
Dag_VSet dag;
auto v0 = boost::add_vertex(dag);
auto v1 = boost::add_vertex(dag);
auto v2 = boost::add_vertex(dag);
boost::add_edge(v1, v2, dag);
boost::add_edge(v2, v0, dag);

std::unordered_map<Dag_VSet::vertex_descriptor, size_t> vMap =
{
	{v0, 1},
	{v1, 2},
	{v2, 0}
};

auto iMap = boost::associative_property_map<
	std::unordered_map<Dag_VSet::vertex_descriptor, size_t>
>(vMap);

Minimum_Feedback_Arc_Set_Solver_Eades_Lin<Dag_VSet> solver(dag, iMap);
REQUIRE(solver.feedback_arc_set().empty());
REQUIRE(
	solver.topological_sorted_vertices()
	== std::list<Dag_VSet::vertex_descriptor>{v1, v2, v0}
);
```

