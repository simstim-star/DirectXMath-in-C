## DirectXMath in C
This is an adaptation of the official [DirectXMath C++ library](https://github.com/microsoft/DirectXMath/tree/main) using C.

This was done for educational purposes, in order to learn better about SIMD and DirectX.

## How to build (MSVC)
To run the tests, build with the following CMake command: `cmake -B build-test -S . -DXM_USE_TEST=1 -DUSE_VECTORCALL=1 -DCMAKE_BUILD_TYPE=Debug`
To not use __vectorcall, just change the `-DUSE_VECTORCALL=0`

To do the same disabling intrinsics: `cmake -B build-test -S . -DXM_NO_INTRINSICS=ON -DXM_USE_TEST=ON -DCMAKE_BUILD_TYPE=Debug`

To run the tests, run the following: 

```
cmake --build build-test
ctest --test-dir build-test -C Debug -V
```


## Installing (to use in other projects)

```
cmake -S . -B build-install
cmake --build build-install --target INSTALL
```

This will generate the folder `xmathc`, probably in `C:/`. This will allow us to use `xmathc` in other projects.
