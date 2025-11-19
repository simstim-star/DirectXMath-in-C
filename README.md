## DirectXMath in C
This is an adaptation of the official [DirectXMath C++ library](https://github.com/microsoft/DirectXMath/tree/main) using C.

This was done for educational purposes, in order to learn better about SIMD and DirectX.

## How to build (MSVC)
To run the tests, build with the following CMake command: 

`cmake -B build-test-vectorcall -S . -DXM_USE_TEST=1 -DUSE_VECTORCALL=1 -DCMAKE_BUILD_TYPE=Debug`

To not use __vectorcall, just change the `-DUSE_VECTORCALL=0`

To do the same disabling intrinsics: 

`cmake -B build-test-nointrin -S . -DXM_NO_INTRINSICS=ON -DXM_USE_TEST=ON -DCMAKE_BUILD_TYPE=Debug`

To run the tests, run the following: 

```
cmake --build build-test-vectorcall
ctest --test-dir build-test-vectorcall -C Debug -V
```

## How to build (GCC)
Just add `-G "MinGW Makefiles"`:

```
cmake -B build-test-gcc -S . -DXM_USE_TEST=1 -DUSE_VECTORCALL=1 -DCMAKE_BUILD_TYPE=Debug -G "MinGW Makefiles"
cmake --build build-test-gcc
ctest --test-dir build-test-gcc -C Debug -V
```

## Installing (to use in other projects)

MSVC:

```
cmake -S . -B build-install
cmake --build build-install --target INSTALL --config Release
```

GCC:

```
cmake -S . -B build-install-gcc -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build-install-gcc --target install
```

This will generate the folder `xmathc`, probably in `C:/`. This will allow us to use `xmathc` in other projects.
