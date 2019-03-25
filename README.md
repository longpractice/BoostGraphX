# BoostGraphX

This repository adds some additional methods to boost graph library.

In order to let my C++ code reliably reusable, I invest large amount to time on:

1. Unit testing. Reliability of C++ code is very important. All methods are unit-tested. The unit tests are in a seperate repo:
   https://github.com/longpractice/TestBoostGraphX


2. Genericity. Just as boost graph library, my generic template methods adapts to your type system in the most convenient way. 

3. Algorithm/Data Structure Efficiency. Our algorithms and data structures are carefully designed and profiled to guarantee a good memory and time effitiency. STL is used a lot, but carefully not to slow things down. 
C++ is also evolving. I keep an eye on new changes like the polymorphic allocator and polymorphic containers. It is not ideal now, but it would give a huge boost of performance(maybe in C++20 or later).



