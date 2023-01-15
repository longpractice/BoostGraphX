2023 Note:
Unfortunately, This is discontinued. 
The reason is that I found boost graph representation in many cases a failed abstraction. It has a common mistake of C++ data structure: mixing the data storage problem with the data relations problem. Its defined several ways (std vectors, unordered_sets, or lists) of storage of vertices and edges are not efficient enough. Even though you can make the Vertex in boost graph very small and hold pointer to your own data structure, for large graph, the allocations are still very expensive). 

Mostly now I tend to decouple the problem of data storage from data relations: vertices/edges are intrusive linked list nodes and we find the data we  
want through pointer offsets. For example, if I know some well defined unit always holds 3 edges(particularly, I had the mesh having a triangle holding 3 directly edges), I can well put 3 edges inside the triangle data structure, and I can still have library functions for purely "graph" structure using purely edges and vertices pointers. My point is that a graph library should really offer the ability of split out the concern of data storage to let the author to do it flexibly. 


## BoostGraphX

This repository is aimed to add more methods for boost graph library. I have been a professional C++ programmer for several years. I frequently find many algorithms missing from BGL which I would love to add by myself.

---

#### The methods now include:

* Euler trial/cycle/tour/circuit (Hierholzer Algorithm)

* Minimum feedback arc set (heuristic method based on the paper by Eades, Lin et al.: https://doi.org/10.1016/0020-0190(93)90079-O)


---

#### Keypoints:

1. __Unit testing__. Reliability of C++ code is very important. All methods are unit-tested. The unit tests are in a seperate repo:
   https://github.com/longpractice/TestBoostGraphX


2. __Genericity__. Just as boost graph library, my generic template methods adapts to your type system in the most convenient way. 

3. __Algorithm/Data Structure Efficiency__. The algorithms and data structures are carefully designed and profiled to guarantee a good memory and time efficiency. STL is frequently used a lot, but carefully not to slow things down. 
C++ is also evolving. I keep an eye on new changes such the polymorphic allocator and polymorphic containers. It is not ideal now, but it would give a huge boost of performance(maybe in C++20 or later).

---

#### Limitations:

Note that I use C++17 features like structured binding and extract method of associated container. Make sure you have the latest compiler.

Due to the complicated nature of many algorithms, a copy of graph is made to store associated auxiliary data. These copy guarantees a linear O(N+E) time/memory complexity. These costs are normally unavoidable. But for some very large graphs(eg, a Facebook worldwide social network graph), the user might wish to allocate the memory themselves(for example, on a memory-mapped file). This genericity would be improved later. 


