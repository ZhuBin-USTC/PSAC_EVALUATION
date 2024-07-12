We implement our HE-based U-Test protocol and a k-percentile protocol from NDSS'17 using Seal 4.1.1.

W. J. Lu, S. Kawasaki, and J. Sakuma, “Using fully homomorphic encryption for statistical analysis of categorical, ordinal and numerical data,” in Proceedings of the 24th Annual Network and Distributed System Security Symposium (NDSS), 2017.

``` bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=./install .. 
cmake --build . --target install --parallel
```