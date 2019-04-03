## BoostGraphX

I have been a professional C++ programmer for several years. I frequently find many algorithms missing from BGL which I would love to add by myself.

---

#### The methods now include:

* Euler trial/cycle/tour/circuit (Hierholzer Algorithm)

* Minimum feedback arc set (heuristic method based on the paper by Eades, Lin et al.: https://doi.org/10.1016/0020-0190(93)90079-O)


---

#### Keypoints

1. Unit testing. Reliability of C++ code is very important. All methods are unit-tested. The unit tests are in a seperate repo:
   https://github.com/longpractice/TestBoostGraphX


2. Genericity. Just as boost graph library, my generic template methods adapts to your type system in the most convenient way. 

3. Algorithm/Data Structure Efficiency. The algorithms and data structures are carefully designed and profiled to guarantee a good memory and time efficiency. STL is frequently used a lot, but carefully not to slow things down. 
C++ is also evolving. I keep an eye on new changes such the polymorphic allocator and polymorphic containers. It is not ideal now, but it would give a huge boost of performance(maybe in C++20 or later).

---

#### Limitations:

Note that I use C++17 features like structured binding and extract method of associated container. Make sure you have the latest compiler.

Due to the complicated nature of many algorithms, a copy of graph is made to store associated auxiliary data. These copy guarantees a linear O(N+E) time/memory complexity. These costs are normally unavoidable. But for some very large graphs(eg, a Facebook worldwide social network graph), the user might wish to allocate the memory themselves(for example, on a memory-mapped file). This genericity would be improved later. 


