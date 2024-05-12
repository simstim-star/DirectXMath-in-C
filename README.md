## DirectXMath in C
This is an adaptation of the official [DirectXMath C++ library](https://github.com/microsoft/DirectXMath/tree/main) using C.

This was done for educational purposes, in order to learn better about SIMD and DirectX.

## How to build (MSVC)
To run the tests, build with the following CMake command: `cmake -B build-test -S . -DXM_USE_TEST=ON -DCMAKE_BUILD_TYPE=Debug`

To do the same disabling intrinsics: `cmake -B build-test -S . -DXM_NO_INTRINSICS=ON -DXM_USE_TEST=ON -DCMAKE_BUILD_TYPE=Debug`

To run the tests, run the following: `ctest --test-dir build-test -C Debug -V`
