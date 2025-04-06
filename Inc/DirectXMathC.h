#pragma once


#if _XM_VECTORCALL_
#define XM_CALLCONV __vectorcall
#elif defined(__GNUC__)
#define XM_CALLCONV
#else
#define XM_CALLCONV __fastcall
#endif

#if !defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
#if (defined(_M_IX86) || defined(_M_X64) || __i386__ || __x86_64__)
#define _XM_SSE_INTRINSICS_
#elif !defined(_XM_NO_INTRINSICS_)
#error DirectX Math does not support this target
#endif // (defined(_M_IX86) || defined(_M_X64) || __i386__ || __x86_64__)
#endif // !_XM_SSE_INTRINSICS_ && !_XM_NO_INTRINSICS_

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4514 4820)
// C4514/4820: Off by default noise
#endif
#include <math.h>
#include <float.h>
#include <stdbool.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4987)
// C4987: Off by default noise
#include <intrin.h>
#pragma warning(pop)
#endif

#if (defined(__clang__) || defined(__GNUC__)) && (__x86_64__ || __i386__)
#include <cpuid.h>
#endif
#ifdef _XM_SSE_INTRINSICS_
#include <xmmintrin.h>
#include <emmintrin.h>
#endif // !_XM_NO_INTRINSICS_

#include "sal.h"
#include <assert.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4005 4668)
// C4005/4668: Old header issue
#endif
#include <stdint.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#if defined(__GNUC__)
#define XM_ALIGNED_DATA(x) __attribute__ ((aligned(x)))
#define XM_ALIGNED_STRUCT(x) struct __attribute__ ((aligned(x)))
#else
#define XM_ALIGNED_DATA(x) __declspec(align(x))
#define XM_ALIGNED_STRUCT(x) __declspec(align(x)) struct
#endif

/****************************************************************************
 Conditional intrinsics
 ****************************************************************************/

#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)

#if defined(_XM_NO_MOVNT_)
#define XM_STREAM_PS( p, a ) _mm_store_ps((p), (a))
#define XM256_STREAM_PS( p, a ) _mm256_store_ps((p), (a))
#define XM_SFENCE()
#else
#define XM_STREAM_PS( p, a ) _mm_stream_ps((p), (a))
#define XM256_STREAM_PS( p, a ) _mm256_stream_ps((p), (a))
#define XM_SFENCE() _mm_sfence()
#endif

#if defined(_XM_FMA3_INTRINSICS_)
#define XM_FMADD_PS( a, b, c ) _mm_fmadd_ps((a), (b), (c))
#define XM_FNMADD_PS( a, b, c ) _mm_fnmadd_ps((a), (b), (c))
#else
#define XM_FMADD_PS( a, b, c ) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define XM_FNMADD_PS( a, b, c ) _mm_sub_ps((c), _mm_mul_ps((a), (b)))
#endif

#if defined(_XM_AVX_INTRINSICS_) && defined(_XM_FAVOR_INTEL_)
#define XM_PERMUTE_PS( v, c ) _mm_permute_ps((v), c )
#else
#define XM_PERMUTE_PS( v, c ) _mm_shuffle_ps((v), (v), c )
#endif


#define XM_LOADU_SI16( p ) _mm_cvtsi32_si128(*(unsigned short const*)(p))

#endif // _XM_SSE_INTRINSICS_ && !_XM_NO_INTRINSICS_

/****************************************************************************
       
 Data types
       
****************************************************************************/


#if defined(_XM_NO_INTRINSICS_)
typedef struct __vector4
{
    union
    {
        float       vector4_f32[4];
        uint32_t    vector4_u32[4];
    };
} __vector4;
#endif // _XM_NO_INTRINSICS_

//------------------------------------------------------------------------------
// Vector intrinsic: Four 32 bit floating point components aligned on a 16 byte
// boundary and mapped to hardware vector registers. 
// Works as a proxy for a SIMD hardware register.
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
typedef __m128 XMVECTOR;
#elif defined(_XM_NO_INTRINSICS_)
typedef __vector4 XMVECTOR;
#endif

// Fix-up for (1st-3rd) XMVECTOR parameters that are pass-in-register for x86 and vector call; by reference otherwise
#if ( defined(_M_IX86) || _XM_VECTORCALL_ || __i386__ ) && !defined(_XM_NO_INTRINSICS_)
typedef const XMVECTOR FXMVECTOR;
#define XM_PARAM_F(p) p
#define XM_1V(v1) v1
#define XM_2V(v1,v2) v1,v2
#define XM_3V(v1,v2,v3) v1,v2,v3
#else
typedef const XMVECTOR* FXMVECTOR;
#define XM_PARAM_F(p) *p
#define XM_1V(v1) &v1
#define XM_2V(v1,v2) &v1,&v2
#define XM_3V(v1,v2,v3) &v1,&v2,&v3
#endif

// Fix-up for (4th) XMVECTOR parameter to pass in-register for vector call; by reference otherwise
#if (_XM_VECTORCALL_) && !defined(_XM_NO_INTRINSICS_)
typedef const XMVECTOR GXMVECTOR;
#define XM_PARAM_H(p) p
#define XM_4V(v1,v2,v3,v4) v1,v2,v3,v4
#else
typedef const XMVECTOR* GXMVECTOR;
#define XM_PARAM_G(p) *p
#define XM_4V(v1,v2,v3,v4) &v1,&v2,&v3,&v4
#endif

// Fix-up for (5th & 6th) XMVECTOR parameter to pass in-register for vector call; by reference otherwise
#if ( _XM_VECTORCALL_ ) && !defined(_XM_NO_INTRINSICS_)
typedef const XMVECTOR HXMVECTOR;
#define XM_PARAM_H(p) p
#define XM_5V(v1,v2,v3,v4,v5) v1,v2,v3,v4,v5
#define XM_6V(v1,v2,v3,v4,v5,v6) v1,v2,v3,v4,v5,v6
#else
typedef const XMVECTOR* HXMVECTOR;
#define XM_PARAM_H(p) *p
#define XM_5V(v1,v2,v3,v4,v5) &v1,&v2,&v3,&v4,&v5
#define XM_6V(v1,v2,v3,v4,v5,v6) &v1,&v2,&v3,&v4,&v5,&v6
#endif

// Fix-up for (7th+) XMVECTOR parameters to pass by reference
typedef const XMVECTOR* CXMVECTOR;
#if ( _XM_VECTORCALL_ ) && !defined(_XM_NO_INTRINSICS_)
#define XM_7V(v1,v2,v3,v4,v5,v6,v7) v1,v2,v3,v4,v5,v6,&v7
#else
#define XM_7V(v1,v2,v3,v4,v5,v6,v7) &v1,&v2,&v3,&v4,&v5,&v6,&v7
#endif

//------------------------------------------------------------------------------
// Matrix type: Sixteen 32 bit floating point components aligned on a
// 16 byte boundary and mapped to four hardware vector registers

typedef struct XMMATRIX XMMATRIX;

// Fix-up for (1st) XMMATRIX parameter to pass in-register for vector call; by reference otherwise
#if ( _XM_VECTORCALL_ ) && !defined(_XM_NO_INTRINSICS_)
typedef const XMMATRIX FXMMATRIX;
#define XM_PARAM_MATRIX(m) m
#define XM_1M(m1) m1
#define XM_MATRIX_GET(m, i) m.r[i]
#else
typedef const XMMATRIX* FXMMATRIX;
#define XM_PARAM_MATRIX(m) *m
#define XM_1M(m1) &m1
#define XM_MATRIX_GET(m, i) m->r[i]
#endif

// Fix-up for (2nd+) XMMATRIX parameters to pass by reference
typedef const XMMATRIX* CXMMATRIX;
#if ( _XM_VECTORCALL_ ) && !defined(_XM_NO_INTRINSICS_)
#define XM_2M(m1,m2) m1,&m2
#else
#define XM_2M(m1,m2) &m1,&m2
#endif

#ifdef _XM_NO_INTRINSICS_
struct XMMATRIX
#else
XM_ALIGNED_STRUCT(16) XMMATRIX
#endif
{
#ifdef _XM_NO_INTRINSICS_
    union
    {
        XMVECTOR r[4];
        struct
        {
            float _11, _12, _13, _14;
            float _21, _22, _23, _24;
            float _31, _32, _33, _34;
            float _41, _42, _43, _44;
        };
        float m[4][4];
    };
#else
    XMVECTOR r[4];
#endif
};

// 2D Vector; 32 bit floating point components
typedef struct XMFLOAT2 { float x, y; } XMFLOAT2;
// 2D Vector; 32 bit floating point components aligned on a 16 byte boundary
typedef XM_ALIGNED_STRUCT(16) XMFLOAT2A { float x, y; } XMFLOAT2A;
// 2D Vector; 32 bit signed integer components
typedef struct XMINT2 { int32_t x, y; } XMINT2;
// 2D Vector; 32 bit unsigned integer components
typedef struct XMUINT2 { uint32_t x, y; } XMUINT2;

// 3D Vector; 32 bit floating point components
typedef struct XMFLOAT3 { float x, y, z; } XMFLOAT3;
// 3D Vector; 32 bit floating point components aligned on a 16 byte boundary
typedef XM_ALIGNED_STRUCT(16) XMFLOAT3A { float x, y, z; } XMFLOAT3A;
// 3D Vector; 32 bit signed integer components
typedef struct XMINT3 { int32_t x, y, z; } XMINT3;
// 3D Vector; 32 bit unsigned integer components
typedef struct XMUINT3 { uint32_t x, y, z; } XMUINT3;

// 4D Vector; 32 bit floating point components
typedef struct XMFLOAT4 { float x, y, z, w; } XMFLOAT4;
// 4D Vector; 32 bit floating point components aligned on a 16 byte boundary
typedef XM_ALIGNED_STRUCT(16) XMFLOAT4A { float x, y, z, w; } XMFLOAT4A;
// 4D Vector; 32 bit signed integer components
typedef struct XMINT4 { float x, y, z, w; } XMINT4;
// 4D Vector; 32 bit unsigned integer components
typedef struct XMUINT4 { uint32_t x, y, z, w; } XMUINT4;

// 3x3 Matrix: 32 bit floating point components
typedef struct XMFLOAT3X3
{
    union
    {
        struct
        {
            float _11, _12, _13;
            float _21, _22, _23;
            float _31, _32, _33;
        };
        float m[3][3];
    };
} XMFLOAT3X3;

// 4x3 Row-major Matrix: 32 bit floating point components
typedef struct XMFLOAT4X3
{
    union
    {
        struct
        {
            float _11, _12, _13;
            float _21, _22, _23;
            float _31, _32, _33;
            float _41, _42, _43;
        };
        float m[4][3];
        float f[12];
    };
} XMFLOAT4X3;

// 4x3 Row-major Matrix: 32 bit floating point components aligned on a 16 byte boundary
typedef XM_ALIGNED_STRUCT(16) XMFLOAT4X3A
{
    union
    {
        struct
        {
            float _11, _12, _13;
            float _21, _22, _23;
            float _31, _32, _33;
            float _41, _42, _43;
        };
        float m[4][3];
        float f[12];
    };
} XMFLOAT4X3A;

// 3x4 Column-major Matrix: 32 bit floating point components
typedef struct XMFLOAT3X4
{
    union
    {
        struct
        {
            float _11, _12, _13, _14;
            float _21, _22, _23, _24;
            float _31, _32, _33, _34;
        };
        float m[3][4];
        float f[12];
    };
} XMFLOAT3X4;

// 3x4 Column-major Matrix: 32 bit floating point components aligned on a 16 byte boundary
typedef XM_ALIGNED_STRUCT(16) XMFLOAT3X4A
{
    union
    {
        struct
        {
            float _11, _12, _13, _14;
            float _21, _22, _23, _24;
            float _31, _32, _33, _34;
        };
        float m[3][4];
        float f[12];
    };
} XMFLOAT3X4A;

// 4x4 Matrix: 32 bit floating point components
typedef struct XMFLOAT4X4
{
    union
    {
        struct
        {
            float _11, _12, _13, _14;
            float _21, _22, _23, _24;
            float _31, _32, _33, _34;
            float _41, _42, _43, _44;
        };
        float m[4][4];
    };
} XMFLOAT4X4;

// 4x4 Matrix: 32 bit floating point components aligned on a 16 byte boundary
typedef XM_ALIGNED_STRUCT(16) XMFLOAT4X4A
{
    union
    {
        struct
        {
            float _11, _12, _13, _14;
            float _21, _22, _23, _24;
            float _31, _32, _33, _34;
            float _41, _42, _43, _44;
        };
        float m[4][4];
    };
} XMFLOAT4X4A;

// Conversion types for constants

typedef XM_ALIGNED_STRUCT(16) XMVECTORF32
{
    union
    {
        float f[4];
        XMVECTOR v;
    };
} XMVECTORF32;

typedef XM_ALIGNED_STRUCT(16) XMVECTORI32
{
    union
    {
        int32_t i[4];
        XMVECTOR v;
    };
} XMVECTORI32;

typedef XM_ALIGNED_STRUCT(16) XMVECTORU32
{
    union
    {
        uint32_t u[4];
        XMVECTOR v;
    };
} XMVECTORU32;

#ifdef _MSC_VER
#pragma warning(pop)
#endif

/****************************************************************************
 Load operations
****************************************************************************/

inline XMVECTOR XM_CALLCONV XMLoadFloat3(const XMFLOAT3* pSource);
inline XMVECTOR XM_CALLCONV XMLoadInt(_In_ const uint32_t* pSource);
inline XMVECTOR XM_CALLCONV XMLoadFloat(_In_ const float* pSource);
inline XMMATRIX XM_CALLCONV XMLoadFloat3x3(const XMFLOAT3X3* pSource);
inline void XM_CALLCONV XMStoreFloat4x4(XMFLOAT4X4* pDestination, FXMMATRIX M);
inline void XM_CALLCONV XMStoreFloat4x4A(XMFLOAT4X4A* pDestination, FXMMATRIX M);
inline XMMATRIX XM_CALLCONV XMLoadFloat4x4(const XMFLOAT4X4* pSource);
inline XMMATRIX XM_CALLCONV XMLoadFloat4x4A(const XMFLOAT4X4A* pSource);

/****************************************************************************
 3D vector operations
****************************************************************************/

bool XM_CALLCONV XMVector3Equal(FXMVECTOR V1, FXMVECTOR V2);
bool XM_CALLCONV XMVector3IsInfinite(FXMVECTOR V);
XMVECTOR XM_CALLCONV XMVectorZero();
XMVECTOR XM_CALLCONV XMVectorSet(float x, float y, float z, float w);
float XM_CALLCONV XMVectorGetX(FXMVECTOR V);
float XM_CALLCONV XMVectorGetY(FXMVECTOR V);
float XM_CALLCONV XMVectorGetZ(FXMVECTOR V);
XMVECTOR XM_CALLCONV XMVectorMultiply(FXMVECTOR V1, FXMVECTOR V2);
inline XMVECTOR XM_CALLCONV XMVectorDivide(FXMVECTOR V1, FXMVECTOR V2);
inline XMVECTOR XM_CALLCONV XMVectorReplicate(float Value);
XMVECTOR XM_CALLCONV XMVector3Normalize(FXMVECTOR V);
XMVECTOR XM_CALLCONV XMVector3NormalizeEst(FXMVECTOR V);
XMVECTOR XM_CALLCONV XMVectorNegate(FXMVECTOR V);
XMVECTOR XM_CALLCONV XMVector3Cross(FXMVECTOR V1, FXMVECTOR V2);
XMVECTOR XM_CALLCONV XMVector3Dot(FXMVECTOR V1, FXMVECTOR V2);
XMVECTOR XM_CALLCONV XMVectorSelect(FXMVECTOR V1, FXMVECTOR V2, FXMVECTOR Control);
inline XMVECTOR XM_CALLCONV XMVectorSubtract(FXMVECTOR V1, FXMVECTOR V2);
inline XMVECTOR XM_CALLCONV XMVector3LengthSq(FXMVECTOR V);
inline XMVECTOR XM_CALLCONV XMVectorSqrt(FXMVECTOR V);
inline XMVECTOR XM_CALLCONV XMVector3Length(FXMVECTOR V);
inline bool XM_CALLCONV XMVector3Greater(FXMVECTOR V1, FXMVECTOR V2);
inline XMVECTOR XM_CALLCONV XMVectorLerp(FXMVECTOR V0, FXMVECTOR V1, float t);

float XM_CALLCONV XMVectorGetByIndex(FXMVECTOR V, size_t i);

/****************************************************************************
 Matrix operations
****************************************************************************/

inline bool XM_CALLCONV XMMatrixIsIdentity(FXMMATRIX M);
inline XMMATRIX XM_CALLCONV XMMatrixMultiply(FXMMATRIX M1, CXMMATRIX M2);
inline XMMATRIX XM_CALLCONV XMMatrixTranspose(FXMMATRIX M);
inline XMMATRIX XM_CALLCONV XMMatrixTranslation(float OffsetX, float OffsetY, float OffsetZ);
inline XMMATRIX XM_CALLCONV XMMatrixLookToLH(FXMVECTOR EyePosition, FXMVECTOR EyeDirection, FXMVECTOR UpDirection);
inline XMMATRIX XM_CALLCONV XMMatrixLookToRH(FXMVECTOR EyePosition, FXMVECTOR EyeDirection, FXMVECTOR UpDirection);


/****************************************************************************
 Globals
****************************************************************************/

// The purpose of the following global constants is to prevent redundant
// reloading of the constants when they are referenced by more than one
// separate inline math routine called within the same function.  Declaring
// a constant locally within a routine is sufficient to prevent redundant
// reloads of that constant when that single routine is called multiple
// times in a function, but if the constant is used (and declared) in a
// separate math routine it would be reloaded.

#ifndef XMGLOBALCONST
#if defined(__GNUC__) && !defined(__MINGW32__)
#define XMGLOBALCONST extern const __attribute__((weak))
#else
#define XMGLOBALCONST extern const __declspec(selectany)
#endif
#endif

/****************************************************************************
 Constant definitions
****************************************************************************/

#if defined(__XNAMATH_H__) && defined(XM_PI)
#undef XM_PI
#undef XM_2PI
#undef XM_1DIVPI
#undef XM_1DIV2PI
#undef XM_PIDIV2
#undef XM_PIDIV4
#undef XM_SELECT_0
#undef XM_SELECT_1
#undef XM_PERMUTE_0X
#undef XM_PERMUTE_0Y
#undef XM_PERMUTE_0Z
#undef XM_PERMUTE_0W
#undef XM_PERMUTE_1X
#undef XM_PERMUTE_1Y
#undef XM_PERMUTE_1Z
#undef XM_PERMUTE_1W
#undef XM_CRMASK_CR6
#undef XM_CRMASK_CR6TRUE
#undef XM_CRMASK_CR6FALSE
#undef XM_CRMASK_CR6BOUNDS
#undef XM_CACHE_LINE_SIZE
#endif

#define XM_PI  3.141592654f
#define XM_2PI  6.283185307f
#define XM_1DIVPI  0.318309886f
#define XM_1DIV2PI  0.159154943f
#define XM_PIDIV2  1.570796327f
#define XM_PIDIV4  0.785398163f
#define XM_SELECT_0  0x00000000
#define XM_SELECT_1  0xFFFFFFFF
#define XM_PERMUTE_0X  0
#define XM_PERMUTE_0Y  1
#define XM_PERMUTE_0Z  2
#define XM_PERMUTE_0W  3
#define XM_PERMUTE_1X  4
#define XM_PERMUTE_1Y  5
#define XM_PERMUTE_1Z  6
#define XM_PERMUTE_1W  7
#define XM_SWIZZLE_X  0
#define XM_SWIZZLE_Y  1
#define XM_SWIZZLE_Z  2
#define XM_SWIZZLE_W  3
#define XM_CRMASK_CR6  0x000000F0
#define XM_CRMASK_CR6TRUE  0x00000080
#define XM_CRMASK_CR6FALSE  0x00000020
#define XM_CRMASK_CR6BOUNDS  XM_CRMASK_CR6FALSE
#define XM_CACHE_LINE_SIZE  64

XMGLOBALCONST XMVECTORF32 g_XMIdentityR0 = { { 1.0f, 0.0f, 0.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMIdentityR1 = { { 0.0f, 1.0f, 0.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMIdentityR2 = { { 0.0f, 0.0f, 1.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMIdentityR3 = { { 0.0f, 0.0f, 0.0f, 1.0f } };
XMGLOBALCONST XMVECTORF32 g_XMNegIdentityR0 = { { -1.0f, 0.0f, 0.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMNegIdentityR1 = { { 0.0f, -1.0f, 0.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMNegIdentityR2 = { { 0.0f, 0.0f, -1.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMNegIdentityR3 = { { 0.0f, 0.0f, 0.0f, -1.0f } };
XMGLOBALCONST XMVECTORI32 g_XMInfinity = { { 0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000 } };
XMGLOBALCONST XMVECTORI32 g_XMQNaN = { { 0x7FC00000, 0x7FC00000, 0x7FC00000, 0x7FC00000 } };
XMGLOBALCONST XMVECTORI32 g_XMAbsMask = { { 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF } };
XMGLOBALCONST XMVECTORU32 g_XMMaskXY = { { { 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000 } } };
XMGLOBALCONST XMVECTORU32 g_XMMask3 = { { { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 } } };
XMGLOBALCONST XMVECTORU32 g_XMMaskX = { { { 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 } } };
XMGLOBALCONST XMVECTORU32 g_XMMaskY = { { { 0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000 } } };
XMGLOBALCONST XMVECTORU32 g_XMMaskZ = { { { 0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000 } } };
XMGLOBALCONST XMVECTORU32 g_XMMaskW = { { { 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF } } };
XMGLOBALCONST XMVECTORU32 g_XMSelect0101 = { { { XM_SELECT_0, XM_SELECT_1, XM_SELECT_0, XM_SELECT_1 } } };
XMGLOBALCONST XMVECTORU32 g_XMSelect1010 = { { { XM_SELECT_1, XM_SELECT_0, XM_SELECT_1, XM_SELECT_0 } } };
XMGLOBALCONST XMVECTORI32 g_XMOneHalfMinusEpsilon = { { { 0x3EFFFFFD, 0x3EFFFFFD, 0x3EFFFFFD, 0x3EFFFFFD } } };
XMGLOBALCONST XMVECTORU32 g_XMSelect1000 = { { { XM_SELECT_1, XM_SELECT_0, XM_SELECT_0, XM_SELECT_0 } } };
XMGLOBALCONST XMVECTORU32 g_XMSelect1100 = { { { XM_SELECT_1, XM_SELECT_1, XM_SELECT_0, XM_SELECT_0 } } };
XMGLOBALCONST XMVECTORU32 g_XMSelect1110 = { { { XM_SELECT_1, XM_SELECT_1, XM_SELECT_1, XM_SELECT_0 } } };
XMGLOBALCONST XMVECTORU32 g_XMSelect1011 = { { { XM_SELECT_1, XM_SELECT_0, XM_SELECT_1, XM_SELECT_1 } } };
XMGLOBALCONST XMVECTORF32 g_XMOne = { { 1.0f, 1.0f, 1.0f, 1.0f } };
XMGLOBALCONST XMVECTORF32 g_XMOne3 = { { 1.0f, 1.0f, 1.0f, 0.0f } };
XMGLOBALCONST XMVECTORF32 g_XMZero = { { 0.0f, 0.0f, 0.0f, 0.0f } };


#define XMFLOAT4X4_ZERO (XMFLOAT4X4) {                  \
        ._11 = 0.f, ._12 = 0.f, ._13 = 0.f, ._14 = 0.f, \
        ._21 = 0.f, ._22 = 0.f, ._23 = 0.f, ._24 = 0.f, \
        ._31 = 0.f, ._32 = 0.f, ._33 = 0.f, ._34 = 0.f, \
        ._41 = 0.f, ._42 = 0.f, ._43 = 0.f, ._44 = 0.f, \
    }

#define XMFLOAT3X3_IDENTITY (XMFLOAT3X3) {  \
        ._11 = 1.f, ._12 = .0f, ._13 = .0f, \
        ._21 = .0f, ._22 = 1.f, ._23 = .0f, \
        ._31 = .0f, ._32 = .0f, ._33 = 1.f, \
};

/****************************************************************************
 Miscellaneous operations
****************************************************************************/

inline bool XMScalarNearEqual(float S1,	float S2, float Epsilon);
inline void XMScalarSinCos(float* pSin,float* pCos,float  Value);

/****************************************************************************
 Add inline implementations
****************************************************************************/


#include "DirectXMathCMatrix.inl"
#include "DirectXMathCConvert.inl"
#include "DirectXMathCVector.inl"
#include "DirectXMathCMisc.inl"
