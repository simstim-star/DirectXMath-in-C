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
    return _mm_cvtss_f32(*V);
#endif
}

// Return the Y component in an FPU register.
inline float XM_CALLCONV XMVectorGetY(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[1];
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(*V, _MM_SHUFFLE(1, 1, 1, 1));
    return _mm_cvtss_f32(vTemp);
#endif
}

// Return the Z component in an FPU register.
inline float XM_CALLCONV XMVectorGetZ(FXMVECTOR V)
{
#if defined(_XM_NO_INTRINSICS_)
    return V->vector4_f32[2];
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp = XM_PERMUTE_PS(*V, _MM_SHUFFLE(2, 2, 2, 2));
    return _mm_cvtss_f32(vTemp);
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
    XMVECTOR vTemp = _mm_cmpeq_ps(*V1, *V2);
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
    __m128 vTemp = _mm_and_ps(*V, g_XMAbsMask.v);
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
    return _mm_mul_ps(*V1, *V2);
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
    return _mm_div_ps(*V1, *V2);
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
    XMVECTOR vLengthSq = _mm_mul_ps(*V, *V);
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
    XMVECTOR vDot = _mm_mul_ps(*V, *V);
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
    vDot = _mm_mul_ps(vDot, *V);
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
    XMVECTOR vLengthSq = _mm_mul_ps(*V, *V);
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
    vResult = _mm_div_ps(*V, vResult);
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
    return _mm_sub_ps(Z, *V);
#endif
}

inline XMVECTOR XM_CALLCONV XMVector3Cross
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    // [ V1.y*V2.z - V1.z*V2.y, V1.z*V2.x - V1.x*V2.z, V1.x*V2.y - V1.y*V2.x ]
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
    XMVECTOR vTemp1 = XM_PERMUTE_PS(*V1, _MM_SHUFFLE(3, 0, 2, 1));
    // z2,x2,y2,w2
    XMVECTOR vTemp2 = XM_PERMUTE_PS(*V2, _MM_SHUFFLE(3, 1, 0, 2));
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
    XMVECTOR vDot = _mm_mul_ps(*V1, *V2);
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
    XMVECTOR vTemp1 = _mm_andnot_ps(*Control, *V1);
    XMVECTOR vTemp2 = _mm_and_ps(*V2, *Control);
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
    return _mm_sub_ps(*V1, *V2);
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
    return _mm_sqrt_ps(*V);
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
    XMVECTOR vTemp = _mm_dp_ps(*V, *V, 0x7f);
    return _mm_sqrt_ps(vTemp);
#elif defined(_XM_SSE3_INTRINSICS_)
    XMVECTOR vLengthSq = _mm_mul_ps(*V, *V);
    vLengthSq = _mm_and_ps(vLengthSq, g_XMMask3);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
    vLengthSq = _mm_sqrt_ps(vLengthSq);
    return vLengthSq;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z
    XMVECTOR vLengthSq = _mm_mul_ps(*V, *V);
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
    XMVECTOR vTemp = _mm_cmpgt_ps(*V1, *V2);
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
    XMVECTOR L = _mm_sub_ps(*V1, *V0);
    XMVECTOR S = _mm_set_ps1(t);
    return XM_FMADD_PS(L, S, *V0);
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
    return _mm_mul_ps(vResult, *V);
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
    return _mm_add_ps(*V1, *V2);
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
    XMVECTOR vTemp = _mm_cmple_ps(*V1, *V2);
    return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
}