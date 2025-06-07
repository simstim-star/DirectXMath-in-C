#include "DirectXMathC.h"
#pragma once

// Initialize a vector with four floating point values
inline XMVECTOR XM_CALLCONV XMVectorSet
(
    float x,
    float y,
    float z,
    float w
) 
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = { { { x, y, z, w } } };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_set_ps(w, z, y, x);
#endif
}

// Return the X component in an FPU register.
inline float XM_CALLCONV XMVectorGetX(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[0];
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_cvtss_f32(XM_DEREF_F(V));
#endif
}

// Return the Y component in an FPU register.
inline float XM_CALLCONV XMVectorGetY(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[1];
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(1, 1, 1, 1));
    return _mm_cvtss_f32(vTemp);
#endif
}

// Return the Z component in an FPU register.
inline float XM_CALLCONV XMVectorGetZ(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[2];
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(2, 2, 2, 2));
    return _mm_cvtss_f32(vTemp);
#endif
}

// Return the W component in an FPU register.
inline float XM_CALLCONV XMVectorGetW(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[3];
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(3, 3, 3, 3));
    return _mm_cvtss_f32(vTemp);
#endif
}

// Store a component indexed by i into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XM_CALLCONV XMVectorGetByIndexPtr(float* f, FXMVECTOR V, size_t i)
{
    assert(f != NULL);
    assert(i < 4);
    _Analysis_assume_(i < 4);
#if defined(_XM_NO_INTRINSICS_)
    *f = V->vector4_f32[i];
#else
    XMVECTORF32 U;
    U.v = XM_DEREF_F(V);
    *f = U.f[i];
#endif
}

//------------------------------------------------------------------------------

// Store the X component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XM_CALLCONV XMVectorGetXPtr(float* x, FXMVECTOR V)
{
    assert(x != NULL);
#if defined(_XM_NO_INTRINSICS_)
    *x = V->vector4_f32[0];
#elif defined(_XM_SSE_INTRINSICS_)
    _mm_store_ss(x, XM_DEREF_F(V));
#endif
}

// Store the Y component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XM_CALLCONV XMVectorGetYPtr(float* y, FXMVECTOR V)
{
    assert(y != NULL);
#if defined(_XM_NO_INTRINSICS_)
    *y = V->vector4_f32[1];
#elif defined(_XM_SSE4_INTRINSICS_)
    * ((int*)(y)) = _mm_extract_ps(V, 1);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(1, 1, 1, 1));
    _mm_store_ss(y, vResult);
#endif
}

// Store the Z component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XM_CALLCONV XMVectorGetZPtr(float* z, FXMVECTOR V) 
{
    assert(z != NULL);
#if defined(_XM_NO_INTRINSICS_)
    *z = V->vector4_f32[2];
#elif defined(_XM_SSE4_INTRINSICS_)
    * ((int*)(z)) = _mm_extract_ps(V, 2);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(2, 2, 2, 2));
    _mm_store_ss(z, vResult);
#endif
}

// Store the W component into a 32 bit float location in memory.
_Use_decl_annotations_
inline void XM_CALLCONV XMVectorGetWPtr(float* w, FXMVECTOR V) 
{
    assert(w != NULL);
#if defined(_XM_NO_INTRINSICS_)
    *w = V->vector4_f32[3];
#elif defined(_XM_SSE4_INTRINSICS_)
    * ((int*)(w)) = _mm_extract_ps(V, 3);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(3, 3, 3, 3));
    _mm_store_ss(w, vResult);
#endif
}

inline bool XM_CALLCONV XMVector3Equal
(
    FXMVECTOR V1,
    FXMVECTOR V2
) 
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1->vector4_f32[0] == V2->vector4_f32[0]) && (V1->vector4_f32[1] == V2->vector4_f32[1]) && (V1->vector4_f32[2] == V2->vector4_f32[2])) != 0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpeq_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
    return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
}

inline bool XM_CALLCONV XMVector3IsInfinite(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return (XMISINF(V.vector4_f32[0]) ||
        XMISINF(V.vector4_f32[1]) ||
        XMISINF(V.vector4_f32[2]));
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bit
    __m128 vTemp = _mm_and_ps(XM_DEREF_F(V), g_XMAbsMask.v);
    // Compare to infinity
    vTemp = _mm_cmpeq_ps(vTemp, g_XMInfinity.v);
    // If x,y or z are infinity, the signs are true.
    return ((_mm_movemask_ps(vTemp) & 7) != 0);
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorZero()
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = { { { 0.0f, 0.0f, 0.0f, 0.0f } } };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_setzero_ps();
#endif
}


inline XMVECTOR XM_CALLCONV XMVectorMultiply
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V1->vector4_f32[0] * V2->vector4_f32[0],
            V1->vector4_f32[1] * V2->vector4_f32[1],
            V1->vector4_f32[2] * V2->vector4_f32[2],
            V1->vector4_f32[3] * V2->vector4_f32[3]
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_mul_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorDivide
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V1->vector4_f32[0] / V2->vector4_f32[0],
            V1->vector4_f32[1] / V2->vector4_f32[1],
            V1->vector4_f32[2] / V2->vector4_f32[2],
            V1->vector4_f32[3] / V2->vector4_f32[3]
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_div_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
#endif
}

// Initialize a vector with a replicated floating point value
inline XMVECTOR XM_CALLCONV XMVectorReplicate(float Value)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = Value;
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_set_ps1(Value);
#endif
}

// Initialize a vector with a replicated floating point value passed by pointer
_Use_decl_annotations_
inline XMVECTOR XM_CALLCONV XMVectorReplicatePtr(const float* pValue)
{
#if defined(_XM_NO_INTRINSICS_)
    float Value = pValue[0];
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = Value;
    return vResult.v;
#elif defined(_XM_AVX_INTRINSICS_)
    return _mm_broadcast_ss(pValue);
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_load_ps1(pValue);
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3ReciprocalLengthEst(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;
    Result = XMVector3LengthSq(V);
    Result = XMVectorReciprocalSqrtEst(Result);
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    // vTemp has z and y
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 2, 1, 2));
    // x+z, y
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    // y,y,y,y
    vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
    // x+z+y,??,??,??
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    // Splat the length squared
    vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
    // Get the reciprocal
    vLengthSq = _mm_rsqrt_ps(vLengthSq);
    return vLengthSq;
#endif
}

// XMVector3NormalizeEst uses a reciprocal estimate and
// returns QNaN on zero and infinite vectors.

inline XMVECTOR XM_CALLCONV XMVector3NormalizeEst(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XMVector3ReciprocalLength(V);
    Result = XMVectorMultiply(V, Result);
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product
    XMVECTOR vDot = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    // x=Dot.y, y=Dot.z
    XMVECTOR vTemp = XM_PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
    // Result.x = x+y
    vDot = _mm_add_ss(vDot, vTemp);
    // x=Dot.z
    vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
    // Result.x = (x+y)+z
    vDot = _mm_add_ss(vDot, vTemp);
    // Splat x
    vDot = XM_PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
    // Get the reciprocal
    vDot = _mm_rsqrt_ps(vDot);
    // Perform the normalization
    vDot = _mm_mul_ps(vDot, XM_DEREF_F(V));
    return vDot;
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3Normalize(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    float fLength;
    XMVECTOR vResult;

    vResult = XMVector3Length(V);
    fLength = vResult.vector4_f32[0];

    // Prevent divide by zero
    if (fLength > 0)
    {
        fLength = 1.0f / fLength;
    }

    vResult.vector4_f32[0] = V.vector4_f32[0] * fLength;
    vResult.vector4_f32[1] = V.vector4_f32[1] * fLength;
    vResult.vector4_f32[2] = V.vector4_f32[2] * fLength;
    vResult.vector4_f32[3] = V.vector4_f32[3] * fLength;
    return vResult;

#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z only
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 1, 2, 1));
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
    // Prepare for the division
    XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
    // Create zero with a single instruction
    XMVECTOR vZeroMask = _mm_setzero_ps();
    // Test for a divide by zero (Must be FP to detect -0.0)
    vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
    // Failsafe on zero (Or epsilon) length planes
    // If the length is infinity, set the elements to zero
    vLengthSq = _mm_cmpneq_ps(vLengthSq, g_XMInfinity.v);
    // Divide to perform the normalization
    vResult = _mm_div_ps(XM_DEREF_F(V), vResult);
    // Any that are infinity, set to zero
    vResult = _mm_and_ps(vResult, vZeroMask);
    // Select qnan or result based on infinite length
    XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_XMQNaN.v);
    XMVECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
    vResult = _mm_or_ps(vTemp1, vTemp2);
    return vResult;
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorNegate(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTORF32 Result = { { {
            -V->vector4_f32[0],
            -V->vector4_f32[1],
            -V->vector4_f32[2],
            -V->vector4_f32[3]
        } } };
    return Result.v;

#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR Z;
    Z = _mm_setzero_ps();
    return _mm_sub_ps(Z, XM_DEREF_F(V));
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3Cross
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    // [ V1.yXM_ARG_F(V)2.z - V1.zXM_ARG_F(V)2.y, V1.zXM_ARG_F(V)2.x - V1.xXM_ARG_F(V)2.z, V1.xXM_ARG_F(V)2.y - V1.yXM_ARG_F(V)2.x ]
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = { { {
            (V1->vector4_f32[1] * V2->vector4_f32[2]) - (V1->vector4_f32[2] * V2->vector4_f32[1]),
            (V1->vector4_f32[2] * V2->vector4_f32[0]) - (V1->vector4_f32[0] * V2->vector4_f32[2]),
            (V1->vector4_f32[0] * V2->vector4_f32[1]) - (V1->vector4_f32[1] * V2->vector4_f32[0]),
            0.0f
        } } };
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // y1,z1,x1,w1
    XMVECTOR vTemp1 = XM_PERMUTE_PS(XM_DEREF_F(V1), _MM_SHUFFLE(3, 0, 2, 1));
    // z2,x2,y2,w2
    XMVECTOR vTemp2 = XM_PERMUTE_PS(XM_DEREF_F(V2), _MM_SHUFFLE(3, 1, 0, 2));
    // Perform the left operation
    XMVECTOR vResult = _mm_mul_ps(vTemp1, vTemp2);
    // z1,x1,y1,w1
    vTemp1 = XM_PERMUTE_PS(vTemp1, _MM_SHUFFLE(3, 0, 2, 1));
    // y2,z2,x2,w2
    vTemp2 = XM_PERMUTE_PS(vTemp2, _MM_SHUFFLE(3, 1, 0, 2));
    // Perform the right operation
    vResult = XM_FNMADD_PS(vTemp1, vTemp2, vResult);
    // Set w to zero
    return _mm_and_ps(vResult, g_XMMask3.v);
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3Dot
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    float fValue = V1->vector4_f32[0] * V2->vector4_f32[0] + V1->vector4_f32[1] * V2->vector4_f32[1] + V1->vector4_f32[2] * V2->vector4_f32[2];
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = fValue;
    return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product
    XMVECTOR vDot = _mm_mul_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
    // x=Dot.vector4_f32[1], y=Dot.vector4_f32[2]
    XMVECTOR vTemp = XM_PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
    // Result.vector4_f32[0] = x+y
    vDot = _mm_add_ss(vDot, vTemp);
    // x=Dot.vector4_f32[2]
    vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
    // Result.vector4_f32[0] = (x+y)+z
    vDot = _mm_add_ss(vDot, vTemp);
    // Splat x
    return XM_PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorSelect
(
    FXMVECTOR V1,
    FXMVECTOR V2,
    FXMVECTOR Control
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTORU32 Result = { { {
            (V1->vector4_u32[0] & ~Control->vector4_u32[0]) | (V2->vector4_u32[0] & Control->vector4_u32[0]),
            (V1->vector4_u32[1] & ~Control->vector4_u32[1]) | (V2->vector4_u32[1] & Control->vector4_u32[1]),
            (V1->vector4_u32[2] & ~Control->vector4_u32[2]) | (V2->vector4_u32[2] & Control->vector4_u32[2]),
            (V1->vector4_u32[3] & ~Control->vector4_u32[3]) | (V2->vector4_u32[3] & Control->vector4_u32[3]),
        } } };
    return Result.v;

#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp1 = _mm_andnot_ps(XM_DEREF_F(Control), XM_DEREF_F(V1));
    XMVECTOR vTemp2 = _mm_and_ps(XM_DEREF_F(V2), XM_DEREF_F(Control));
    return _mm_or_ps(vTemp1, vTemp2);
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorSubtract
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V1->vector4_f32[0] - V2->vector4_f32[0],
            V1->vector4_f32[1] - V2->vector4_f32[1],
            V1->vector4_f32[2] - V2->vector4_f32[2],
            V1->vector4_f32[3] - V2->vector4_f32[3]
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_sub_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3LengthSq(FXMVECTOR V)
{
    return XMVector3Dot(V, V);
}

inline XMVECTOR XM_CALLCONV XMVectorSqrt(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            sqrtf(V->vector4_f32[0]),
            sqrtf(V->vector4_f32[1]),
            sqrtf(V->vector4_f32[2]),
            sqrtf(V->vector4_f32[3])
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_sqrt_ps(XM_DEREF_F(V));
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3Length(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;

    Result = XMVector3LengthSq(V);
    Result = XMVectorSqrt(&Result);

    return Result;
#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vTemp = _mm_dp_ps(XM_ARG_F(V), XM_ARG_F(V), 0x7f);
    return _mm_sqrt_ps(vTemp);
#elif defined(_XM_SSE3_INTRINSICS_)
    XMVECTOR vLengthSq = _mm_mul_ps(XM_ARG_F(V), XM_ARG_F(V));
    vLengthSq = _mm_and_ps(vLengthSq, g_XMMask3);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    // vTemp has z and y
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 2, 1, 2));
    // x+z, y
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    // y,y,y,y
    vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
    // x+z+y,??,??,??
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    // Splat the length squared
    vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
    // Get the length
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#endif
}

inline bool XM_CALLCONV XMVector3Greater
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1->vector4_f32[0] > V2->vector4_f32[0]) && (V1->vector4_f32[1] > V2->vector4_f32[1]) && (V1->vector4_f32[2] > V2->vector4_f32[2])) != 0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmpgt_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
    return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorLerp
(
    FXMVECTOR V0,
    FXMVECTOR V1,
    float    t
) 
{
    // V0 + t * (V1 - V0)

#if defined(_XM_NO_INTRINSICS_)
    // TODO
    XMVECTOR Scale = XMVectorReplicate(t);
    XMVECTOR Length = XMVectorSubtract(V1, V0);
    return XMVectorMultiplyAdd(Length, Scale, V0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR L = _mm_sub_ps(XM_DEREF_F(V1), XM_DEREF_F(V0));
    XMVECTOR S = _mm_set_ps1(t);
    return XM_FMADD_PS(L, S, XM_DEREF_F(V0));
#endif
}


//------------------------------------------------------------------------------
// Return a floating point value via an index. This is not a recommended
// function to use due to performance loss.
inline float XM_CALLCONV XMVectorGetByIndex(FXMVECTOR V, size_t i)
{
    assert(i < 4);
    _Analysis_assume_(i < 4);
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[i];
#else
    XMVECTORF32 U;
    U.v = XM_DEREF_F(V);
    return U.f[i];
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorScale
(
    FXMVECTOR V,
    float    ScaleFactor
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V->vector4_f32[0] * ScaleFactor,
            V->vector4_f32[1] * ScaleFactor,
            V->vector4_f32[2] * ScaleFactor,
            V->vector4_f32[3] * ScaleFactor
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_set_ps1(ScaleFactor);
    return _mm_mul_ps(vResult, XM_DEREF_F(V));
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorAdd
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V1->vector4_f32[0] + V2->vector4_f32[0],
            V1->vector4_f32[1] + V2->vector4_f32[1],
            V1->vector4_f32[2] + V2->vector4_f32[2],
            V1->vector4_f32[3] + V2->vector4_f32[3]
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_add_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
#endif
}

inline bool XM_CALLCONV XMVector3LessOrEqual
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)
    return (((V1->vector4_f32[0] <= V2->vector4_f32[0]) && (V1->vector4_f32[1] <= V2->vector4_f32[1]) && (V1->vector4_f32[2] <= V2->vector4_f32[2])) != 0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = _mm_cmple_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
    return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorSwizzle
(
    FXMVECTOR V,
    uint32_t E0,
    uint32_t E1,
    uint32_t E2,
    uint32_t E3
)
{
    assert((E0 < 4) && (E1 < 4) && (E2 < 4) && (E3 < 4));
    _Analysis_assume_((E0 < 4) && (E1 < 4) && (E2 < 4) && (E3 < 4));
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V->vector4_f32[E0],
            V->vector4_f32[E1],
            V->vector4_f32[E2],
            V->vector4_f32[E3]
        } } };
    return Result.v;
#elif defined(_XM_AVX_INTRINSICS_)
    unsigned int elem[4] = { E0, E1, E2, E3 };
    __m128i vControl = _mm_loadu_si128((const __m128i*)(&elem[0]));
    return _mm_permutevar_ps(V, vControl);
#else
    const uint32_t* aPtr = (const uint32_t*)(&V);

    XMVECTOR Result;
    uint32_t* pWork = (uint32_t*)(&Result);

    pWork[0] = aPtr[E0];
    pWork[1] = aPtr[E1];
    pWork[2] = aPtr[E2];
    pWork[3] = aPtr[E3];

    return Result;
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorNegativeMultiplySubtract
(
    FXMVECTOR V1,
    FXMVECTOR V2,
    FXMVECTOR V3
) 
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V3->vector4_f32[0] - (V1->vector4_f32[0] * V2->vector4_f32[0]),
            V3->vector4_f32[1] - (V1->vector4_f32[1] * V2->vector4_f32[1]),
            V3->vector4_f32[2] - (V1->vector4_f32[2] * V2->vector4_f32[2]),
            V3->vector4_f32[3] - (V1->vector4_f32[3] * V2->vector4_f32[3])
        } } };
    return Result;
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_FNMADD_PS(XM_DEREF_F(V1), XM_DEREF_F(V2), XM_DEREF_F(V3));
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorMultiplyAdd
(
    FXMVECTOR V1,
    FXMVECTOR V2,
    FXMVECTOR V3
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            V1->vector4_f32[0] * V2.vector4_f32[0] + V3->vector4_f32[0],
            V1->vector4_f32[1] * V2.vector4_f32[1] + V3->vector4_f32[1],
            V1->vector4_f32[2] * V2.vector4_f32[2] + V3->vector4_f32[2],
            V1->vector4_f32[3] * V2.vector4_f32[3] + V3->vector4_f32[3]
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_FMADD_PS(XM_DEREF_F(V1), XM_DEREF_F(V2), XM_DEREF_F(V3));
#endif
}

inline XMVECTOR XM_CALLCONV XMVector4Dot
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTORF32 Result;
    Result.f[0] =
        Result.f[1] =
        Result.f[2] =
        Result.f[3] = V1->vector4_f32[0] * V2->vector4_f32[0] + V1->vector4_f32[1] * V2->vector4_f32[1] + V1->vector4_f32[2] * V2->vector4_f32[2] + V1->vector4_f32[3] * V2->vector4_f32[3];
    return Result.v;

#elif defined(_XM_SSE4_INTRINSICS_)
    return _mm_dp_ps(XM_DEREF_F(V1), XM_DEREF_F(V2), 0xff);
#elif defined(_XM_SSE3_INTRINSICS_)
    XMVECTOR vTemp = _mm_mul_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
    vTemp = _mm_hadd_ps(vTemp, vTemp);
    return _mm_hadd_ps(vTemp, vTemp);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp2 = XM_DEREF_F(V2);
    XMVECTOR vTemp = _mm_mul_ps(XM_DEREF_F(V1), vTemp2);
    vTemp2 = _mm_shuffle_ps(vTemp2, vTemp, _MM_SHUFFLE(1, 0, 0, 0)); // Copy X to the Z position and Y to the W position
    vTemp2 = _mm_add_ps(vTemp2, vTemp);          // Add Z = X+Z; W = Y+W;
    vTemp = _mm_shuffle_ps(vTemp, vTemp2, _MM_SHUFFLE(0, 3, 0, 0));  // Copy W to the Z position
    vTemp = _mm_add_ps(vTemp, vTemp2);           // Add Z and W together
    return XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(2, 2, 2, 2));    // Splat Z and return
#endif
}

inline XMVECTOR XM_CALLCONV XMVector4Length(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XMVector4LengthSq(XM_DEREF_F(V));
    Result = XMVectorSqrt(Result);

    return Result;

#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vTemp = _mm_dp_ps(XM_DEREF_F(V), XM_DEREF_F(V), 0xff);
    return _mm_sqrt_ps(vTemp);
#elif defined(_XM_SSE3_INTRINSICS_)
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y,z and w
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    // vTemp has z and w
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(3, 2, 3, 2));
    // x+z, y+w
    vLengthSq = _mm_add_ps(vLengthSq, vTemp);
    // x+z,x+z,x+z,y+w
    vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 0, 0, 0));
    // ??,??,y+w,y+w
    vTemp = _mm_shuffle_ps(vTemp, vLengthSq, _MM_SHUFFLE(3, 3, 0, 0));
    // ??,??,x+z+y+w,??
    vLengthSq = _mm_add_ps(vLengthSq, vTemp);
    // Splat the length
    vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 2, 2, 2));
    // Get the length
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#endif
}


