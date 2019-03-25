# Find Directed Euler Cycle/Circuit/Path/Tour with Hierholzer's algorithm

The methods are in `/Public/BoostGraphX/directed_bidirectional_graph_concept.h`.

For the notion of Euler path, see: https://en.wikipedia.org/wiki/Eulerian_path.

We implement the Hierholzer's algorith to find an Euler cycle or an Euler path. With Linear Time complexity O(N + E). With linear space O(N + E).

With a graph of 1 million vertices and 2 million edges, finding one euler cycle takes around 1.2s in Intel Core i7-7700K 4.2GHz and 16Gb ram.

There are two sets of methods available: 


## Find a directed Eulerian cycle
staring and ending at `euler_cycle_start_vertex`.
```cpp
template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_cycle_hierholzer
(
    const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor 
    euler_cycle_start_vertex
)
```

```cpp
template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_cycle_hierholzer
(
    const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor 
    euler_cycle_start_vertex,
    const VertexIndexMap& i_map)
)
```

Internally, the methods uses boost::copy_graph to create a seperate auxilliary graph.

The two overloads serve the purpose of:

The first one above is suitable for a graph type with an internal vertex_index property. For example,  adjacency_list with VertexList=vecS has this property. The vertex index map is acquired by default using `boost::get(vertex_index, g)` when calling boost::copy_graph.

The second one is suitable for a graph type without an internal vertex_index property. For exmaple, adjacency_list with VertexList=listS or setS. The vertex index map type of VertexIndexMap  must be a model of Readable Property Map and must map the vertex descriptors of G to the integers in the half-open range [0,num_vertices(G)). This is normally adapted from a key-value container like std::map, std::unordered_map using boost::associative_property_map.

## Find a directed Eulerian trail
starts from euler_tour_start_vertex and ends at euler_tour_end_vertex. `euler_trail_start_vertex` must not equal to `euler_trail_end_vertex`.
```cpp
template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_trail_hierholzer
(
    const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_trail_start_vertex,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_trail_end_vertex
)
```

```cpp
template <typename VertexListGraph>
std::list<typename boost::graph_traits<VertexListGraph>::edge_descriptor>
find_one_directed_euler_trail_hierholzer
(
    const VertexListGraph& g,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_trail_start_vertex,
    typename boost::graph_traits<VertexListGraph>::vertex_descriptor euler_trail_end_vertex,
    const VertexIndexMap& i_map
)
```

Internally, the methods uses boost::copy_graph to create a seperate auxilliary graph.

The two overloads serve the purpose of:

The first one above is suitable for a graph type with an internal vertex_index property. For example,  adjacency_list with VertexList=vecS has this property. The vertex index map is acquired by default using `boost::get(vertex_index, g)` when calling boost::copy_graph.

The second one is suitable for a graph type without an internal vertex_index property. For exmaple, adjacency_list with VertexList=listS or setS. The vertex index map type of VertexIndexMap  must be a model of Readable Property Map and must map the vertex descriptors of G to the integers in the half-open range [0,num_vertices(G)). This is normally adapted from a key-value container like std::map, std::unordered_map using boost::associative_property_map.

eg:

```cpp
//the vertex list type is vecS so that we have internal vertex_index property
using Dag = boost::adjacency_list<
    boost::setS,
    boost::vecS,
    boost::directedS
>;

auto dag = Dag{};
auto v1 = boost::add_vertex(dag);
auto cycle = bglx::find_one_directed_euler_cycle_hierholzer(
    dag, v1
);
//the cycle should contain no edges
assert(cycle.empty());


auto v2 = boost::add_vertex(dag);
auto[e, ok] = boost::add_edge(v1, v2, dag);
auto tour = bglx::find_one_directed_euler_trail_hierholzer(
    dag, v1, v2
);
//the tour should contain only one edge
assert（tour == std::list<Edge>{e});

```

```cpp
//here the vertex list type is setS, which means that we do not have internal vertex_index
//property and we need to pass this in manually
using Dag_VSet = 
        boost::adjacency_list<boost::setS, boost::setS, boost::directedS>;
Dag_VSet dag;
auto v0 = boost::add_vertex(dag);
auto v1 = boost::add_vertex(dag);
auto v2 = boost::add_vertex(dag);
auto[e12, ok12] = boost::add_edge(v1, v2, dag);
auto[e20, ok20] = boost::add_edge(v2, v0, dag);

//arbitrary index mapping
std::unordered_map<Dag_VSet::vertex_descriptor, size_t> vMap =
{
    {v0, 1}, 
    {v1, 2},
    {v2, 0}
};

//adapt std::unordered_map to property_map
auto iMap = boost::associative_property_map<
    std::unordered_map<Dag_VSet::vertex_descriptor, size_t>
>(vMap);

auto cycle = bglx::find_one_directed_euler_trail_hierholzer(dag, v1, v0, iMap);
auto cycleRight = cycle == std::list<Dag_VSet::edge_descriptor>{e12, e20};
REQUIRE(cycleRight);
```
