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
    * ((int*)(y)) = _mm_extract_ps(XM_DEREF_F(V), 1);
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
    * ((int*)(z)) = _mm_extract_ps(XM_DEREF_F(V), 2);
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
    return (XMISINF(V->vector4_f32[0]) ||
        XMISINF(V->vector4_f32[1]) ||
        XMISINF(V->vector4_f32[2]));
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

// Sets the X component of a vector to a passed floating point value
inline XMVECTOR XM_CALLCONV XMVectorSetX(FXMVECTOR V, float x)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 U = { { {
            x,
            V->vector4_f32[1],
            V->vector4_f32[2],
            V->vector4_f32[3]
        } } };
    return U.v;
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = _mm_set_ss(x);
    vResult = _mm_move_ss(XM_DEREF_F(V), vResult);
    return vResult;
#endif
}

// Sets the Y component of a vector to a passed floating point value
inline XMVECTOR XM_CALLCONV XMVectorSetY(FXMVECTOR V, float y)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 U = { { {
            V->vector4_f32[0],
            y,
            V->vector4_f32[2],
            V->vector4_f32[3]
        } } };
    return U.v;
#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vResult = _mm_set_ss(y);
    vResult = _mm_insert_ps(XM_DEREF_F(V), vResult, 0x10);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap y and x
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(3, 2, 0, 1));
    // Convert input to vector
    XMVECTOR vTemp = _mm_set_ss(y);
    // Replace the x component
    vResult = _mm_move_ss(vResult, vTemp);
    // Swap y and x again
    vResult = XM_PERMUTE_PS(vResult, _MM_SHUFFLE(3, 2, 0, 1));
    return vResult;
#endif
}
// Sets the Z component of a vector to a passed floating point value
inline XMVECTOR XM_CALLCONV XMVectorSetZ(FXMVECTOR V, float z)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 U = { { {
            V->vector4_f32[0],
            V->vector4_f32[1],
            z,
            V->vector4_f32[3]
        } } };
    return U.v;
#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vResult = _mm_set_ss(z);
    vResult = _mm_insert_ps(XM_DEREF_F(V), vResult, 0x20);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap z and x
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(3, 0, 1, 2));
    // Convert input to vector
    XMVECTOR vTemp = _mm_set_ss(z);
    // Replace the x component
    vResult = _mm_move_ss(vResult, vTemp);
    // Swap z and x again
    vResult = XM_PERMUTE_PS(vResult, _MM_SHUFFLE(3, 0, 1, 2));
    return vResult;
#endif
}

// Sets the W component of a vector to a passed floating point value
inline XMVECTOR XM_CALLCONV XMVectorSetW(FXMVECTOR V, float w)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 U = { { {
            V->vector4_f32[0],
            V->vector4_f32[1],
            V->vector4_f32[2],
            w
        } } };
    return U.v;
#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vResult = _mm_set_ss(w);
    vResult = _mm_insert_ps(XM_DEREF_F(V), vResult, 0x30);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Swap w and x
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(0, 2, 1, 3));
    // Convert input to vector
    XMVECTOR vTemp = _mm_set_ss(w);
    // Replace the x component
    vResult = _mm_move_ss(vResult, vTemp);
    // Swap w and x again
    vResult = XM_PERMUTE_PS(vResult, _MM_SHUFFLE(0, 2, 1, 3));
    return vResult;
#endif
}

//------------------------------------------------------------------------------
// Replicate the x component of the vector
inline XMVECTOR XM_CALLCONV XMVectorSplatX(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = V->vector4_f32[0];
    return vResult.v;
#elif defined(_XM_AVX2_INTRINSICS_) && defined(_XM_FAVOR_INTEL_)
    return _mm_broadcastss_ps(XM_DEREF_F(V));
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(0, 0, 0, 0));
#endif
}

//------------------------------------------------------------------------------
// Replicate the y component of the vector
inline XMVECTOR XM_CALLCONV XMVectorSplatY(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = V->vector4_f32[1];
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32(vget_low_f32(XM_DEREF_F(V)), 1);
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(1, 1, 1, 1));
#endif
}

//------------------------------------------------------------------------------
// Replicate the z component of the vector
inline XMVECTOR XM_CALLCONV XMVectorSplatZ(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = V->vector4_f32[2];
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32(vget_high_f32(XM_DEREF_F(V)), 0);
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(2, 2, 2, 2));
#endif
}

//------------------------------------------------------------------------------
// Replicate the w component of the vector
inline XMVECTOR XM_CALLCONV XMVectorSplatW(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult;
    vResult.f[0] =
        vResult.f[1] =
        vResult.f[2] =
        vResult.f[3] = V->vector4_f32[3];
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    return vdupq_lane_f32(vget_high_f32(XM_DEREF_F(V)), 1);
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(3, 3, 3, 3));
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

// XMVector3NormalizeEst uses a reciprocal estimate and
// returns QNaN on zero and infinite vectors.

inline XMVECTOR XM_CALLCONV XMVector3NormalizeEst(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;
    Result = XM_VEC3_RECIPROCAL_LEN(V);
    Result = XM_VEC_MULT(V, Result);
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

    vResult.vector4_f32[0] = V->vector4_f32[0] * fLength;
    vResult.vector4_f32[1] = V->vector4_f32[1] * fLength;
    vResult.vector4_f32[2] = V->vector4_f32[2] * fLength;
    vResult.vector4_f32[3] = V->vector4_f32[3] * fLength;
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

inline XMVECTOR XM_CALLCONV XMVectorMergeXY
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTORU32 Result = { { {
            V1->vector4_u32[0],
            V2->vector4_u32[0],
            V1->vector4_u32[1],
            V2->vector4_u32[1],
        } } };
    return Result.v;

#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_unpacklo_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
#endif
}

//------------------------------------------------------------------------------

inline XMVECTOR XM_CALLCONV XMVectorMergeZW
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTORU32 Result = { { {
            V1->vector4_u32[2],
            V2->vector4_u32[2],
            V1->vector4_u32[3],
            V2->vector4_u32[3]
        } } };
    return Result.v;

#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_unpackhi_ps(XM_DEREF_F(V1), XM_DEREF_F(V2));
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

inline XMVECTOR XM_CALLCONV XMVector3ReciprocalLengthEst(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTOR Result;

    Result = XM_VEC3_LEN_SQ(V);
    Result = XM_VEC_RECIPROCAL_SQRT_EST(Result);

    return Result;

#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vTemp = _mm_dp_ps(XM_DEREF_F(V), XM_DEREF_F(V), 0x7f);
    return _mm_rsqrt_ps(vTemp);
#elif defined(_XM_SSE3_INTRINSICS_)
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(V), XM_DEREF_F(V));
    vLengthSq = _mm_and_ps(vLengthSq, g_XMMask3.v);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_rsqrt_ps(vLengthSq);
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
    // Get the reciprocal
    vLengthSq = _mm_rsqrt_ps(vLengthSq);
    return vLengthSq;
#endif
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

inline XMVECTOR XM_CALLCONV XMVectorReciprocalSqrtEst(FXMVECTOR V) 
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            1.f / sqrtf(V->vector4_f32[0]),
            1.f / sqrtf(V->vector4_f32[1]),
            1.f / sqrtf(V->vector4_f32[2]),
            1.f / sqrtf(V->vector4_f32[3])
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_rsqrt_ps(XM_DEREF_F(V));
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3Length(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Result;

    Result = XM_VEC3_LEN_SQ(V);
    Result = XM_VEC_SQRT(Result);

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
    XMVECTOR Scale = XMVectorReplicate(t);
    XMVECTOR Length = XMVectorSubtract(V1, V0);
    return XM_VEC_MULT_ADD(Length, Scale, V0);

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
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return XM_FNMADD_PS(XM_DEREF_F(V1), XM_DEREF_F(V2), XM_DEREF_F(V3));
#endif
}

inline XMVECTOR XM_CALLCONV XMVectorReciprocal(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 Result = { { {
            1.f / V->vector4_f32[0],
            1.f / V->vector4_f32[1],
            1.f / V->vector4_f32[2],
            1.f / V->vector4_f32[3]
        } } };
    return Result.v;
#elif defined(_XM_SSE_INTRINSICS_)
    return _mm_div_ps(g_XMOne.v, XM_DEREF_F(V));
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
            V1->vector4_f32[0] * V2->vector4_f32[0] + V3->vector4_f32[0],
            V1->vector4_f32[1] * V2->vector4_f32[1] + V3->vector4_f32[1],
            V1->vector4_f32[2] * V2->vector4_f32[2] + V3->vector4_f32[2],
            V1->vector4_f32[3] * V2->vector4_f32[3] + V3->vector4_f32[3]
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

    Result = XM_VEC4_LEN_SQ(V);
    Result = XM_VEC_SQRT(Result);

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

inline XMVECTOR XM_CALLCONV XMVector4LengthSq(FXMVECTOR V)
{
    return XMVector4Dot(V, V);
}

inline XMVECTOR XM_CALLCONV XMVectorPermute
(
    FXMVECTOR V1,
    FXMVECTOR V2,
    uint32_t PermuteX,
    uint32_t PermuteY,
    uint32_t PermuteZ,
    uint32_t PermuteW
)
{
    assert(PermuteX <= 7 && PermuteY <= 7 && PermuteZ <= 7 && PermuteW <= 7);
    _Analysis_assume_(PermuteX <= 7 && PermuteY <= 7 && PermuteZ <= 7 && PermuteW <= 7);

    #if defined(_XM_AVX_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    static const XMVECTORU32 three = { { { 3, 3, 3, 3 } } };

    XM_ALIGNED_DATA(16) unsigned int elem[4] = { PermuteX, PermuteY, PermuteZ, PermuteW };
    __m128i vControl = _mm_load_si128((const __m128i*)(&elem[0]));

    __m128i vSelect = _mm_cmpgt_epi32(vControl, _mm_castps_si128(three.v));
    vControl = _mm_castps_si128(_mm_and_ps(_mm_castsi128_ps(vControl), three.v));

    __m128 shuffled1 = _mm_permutevar_ps(XM_DEREF_F(V1), vControl);
    __m128 shuffled2 = _mm_permutevar_ps(XM_DEREF_F(V2), vControl);

    __m128 masked1 = _mm_andnot_ps(_mm_castsi128_ps(vSelect), shuffled1);
    __m128 masked2 = _mm_and_ps(_mm_castsi128_ps(vSelect), shuffled2);

    return _mm_or_ps(masked1, masked2);
    #else

    const uint32_t* aPtr[2];
    aPtr[0] = (const uint32_t*)(&V1);
    aPtr[1] = (const uint32_t*)(&V2);

    XMVECTOR Result;
    uint32_t* pWork = (uint32_t*)(&Result);

    const uint32_t i0 = PermuteX & 3;
    const uint32_t vi0 = PermuteX >> 2;
    pWork[0] = aPtr[vi0][i0];

    const uint32_t i1 = PermuteY & 3;
    const uint32_t vi1 = PermuteY >> 2;
    pWork[1] = aPtr[vi1][i1];

    const uint32_t i2 = PermuteZ & 3;
    const uint32_t vi2 = PermuteZ >> 2;
    pWork[2] = aPtr[vi2][i2];

    const uint32_t i3 = PermuteW & 3;
    const uint32_t vi3 = PermuteW >> 2;
    pWork[3] = aPtr[vi3][i3];

    return Result;
    #endif
}

inline XMVECTOR XM_CALLCONV XMVector3Transform
(
    FXMVECTOR V,
    FXMMATRIX M
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMVECTOR Z = XMVectorSplatZ(V);
    XMVECTOR Y = XMVectorSplatY(V);
    XMVECTOR X = XMVectorSplatX(V);

    XMVECTOR Result = XM_VEC_MULT_ADD(Z, M->r[2], M->r[3]);
    Result = XM_VEC_MULT_ADD(Y, M->r[1], Result);
    Result = XM_VEC_MULT_ADD(X, M->r[0], Result);

    return Result;

#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vResult = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(2, 2, 2, 2)); // Z
    vResult = XM_FMADD_PS(vResult, XM_DEREF_MATRIX(M).r[2], XM_DEREF_MATRIX(M).r[3]);
    XMVECTOR vTemp = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(1, 1, 1, 1)); // Y
    vResult = XM_FMADD_PS(vTemp, XM_DEREF_MATRIX(M).r[1], vResult);
    vTemp = XM_PERMUTE_PS(XM_DEREF_F(V), _MM_SHUFFLE(0, 0, 0, 0)); // X
    vResult = XM_FMADD_PS(vTemp, XM_DEREF_MATRIX(M).r[0], vResult);
    return vResult;
#endif
}


