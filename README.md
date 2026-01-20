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


## Usage examples

To allow it to use callconv like the original DirectXMath, I had to create macros that hide the implementation details. In the original implementation, this is cleaner because C++ can pass arguments by reference or by value seamlessly to the caller.

```
XMMATRIX A_XMMATRIX = XMLoadFloat3x3(&A); // A is XMFLOAT3X3 
XMMATRIX B_XMMATRIX = XMLoadFloat3x3(&B); // B is XMFLOAT3X3

// XM_MAT_MULT macro will handle if it should pass by value or ref depending on adequate callconv 
XMMATRIX C_XMMATRIX = XM_MAT_MULT(A_XMMATRIX, B_XMMATRIX);

XMFLOAT4X4 result;
// XM_STORE_FLOAT4X4 is a macro for the same reason XM_MAT_MULT 
XM_STORE_FLOAT4X4(&result, C_XMMATRIX);
```

