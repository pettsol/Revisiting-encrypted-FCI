Compile the decentralized encrypted fast covariance intersection method with g++ using the following command:

g++ main.cpp de\_efci.cpp privacy\_preserving\_aggregation/ppa.cpp privacy\_preserving\_aggregation/SHA-256/sha-256.cpp -o main -lgmp -O3

The number of sensors whose output we fuse and the number of datasets to iterate through can be configured in main.cpp.
