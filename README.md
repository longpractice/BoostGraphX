# BoostGraphX

This repository adds some additional methods to boost graph library.

In order to let my C++ code reliably reusable, I invest large amount of time on:

1. Unit testing. Reliability of C++ code is very important. All methods are unit-tested. The unit tests are in a seperate repo:
   https://github.com/longpractice/TestBoostGraphX


2. Genericity. Just as boost graph library, my generic template methods adapts to your type system in the most convenient way. 

3. Algorithm/Data Structure Efficiency. The algorithms and data structures are carefully designed and profiled to guarantee a good memory and time efficiency. STL is frequently used a lot, but carefully not to slow things down. 
C++ is also evolving. I keep an eye on new changes such the polymorphic allocator and polymorphic containers. It is not ideal now, but it would give a huge boost of performance(maybe in C++20 or later).


Note that I use C++17 features here like structured binding. Make sure you have the latest compiler. I am just used to writing C++17 code and if it is really a high demand to downgrade to only C++11, I may consider rewrite C++17 parts.


