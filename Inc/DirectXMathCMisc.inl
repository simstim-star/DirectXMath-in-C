#pragma once
#include <math.h>

inline bool XMScalarNearEqual
(
    float S1,
    float S2,
    float Epsilon
)
{
    float Delta = S1 - S2;
    return (fabsf(Delta) <= Epsilon);
}

_Use_decl_annotations_
inline void XMScalarSinCos
(
    float* pSin,
    float* pCos,
    float  Value
)
{
    assert(pSin);
    assert(pCos);

    // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
    float quotient = XM_1DIV2PI * Value;
    if (Value >= 0.0f)
    {
        quotient = (float)((int)(quotient + 0.5f));
    }
    else
    {
        quotient = (float)((int)(quotient - 0.5f));
    }
    float y = Value - XM_2PI * quotient;

    // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
    float sign;
    if (y > XM_PIDIV2)
    {
        y = XM_PI - y;
        sign = -1.0f;
    }
    else if (y < -XM_PIDIV2)
    {
        y = -XM_PI - y;
        sign = -1.0f;
    }
    else
    {
        sign = +1.0f;
    }

    float y2 = y * y;

    // 11-degree minimax approximation
    *pSin = (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;

    // 10-degree minimax approximation
    float p = ((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f;
    *pCos = sign * p;
}


inline XMVECTOR XM_CALLCONV XMQuaternionRotationMatrix(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)

    XMVECTORF32 q;
    float r22 = M->m[2][2];
    if (r22 <= 0.f)  // x^2 + y^2 >= z^2 + w^2
    {
        float dif10 = M->m[1][1] - M->m[0][0];
        float omr22 = 1.f - r22;
        if (dif10 <= 0.f)  // x^2 >= y^2
        {
            float fourXSqr = omr22 - dif10;
            float inv4x = 0.5f / sqrtf(fourXSqr);
            q.f[0] = fourXSqr * inv4x;
            q.f[1] = (M->m[0][1] + M->m[1][0]) * inv4x;
            q.f[2] = (M->m[0][2] + M->m[2][0]) * inv4x;
            q.f[3] = (M->m[1][2] - M->m[2][1]) * inv4x;
        }
        else  // y^2 >= x^2
        {
            float fourYSqr = omr22 + dif10;
            float inv4y = 0.5f / sqrtf(fourYSqr);
            q.f[0] = (M->m[0][1] + M->m[1][0]) * inv4y;
            q.f[1] = fourYSqr * inv4y;
            q.f[2] = (M->m[1][2] + M->m[2][1]) * inv4y;
            q.f[3] = (M->m[2][0] - M->m[0][2]) * inv4y;
        }
    }
    else  // z^2 + w^2 >= x^2 + y^2
    {
        float sum10 = M->m[1][1] + M->m[0][0];
        float opr22 = 1.f + r22;
        if (sum10 <= 0.f)  // z^2 >= w^2
        {
            float fourZSqr = opr22 - sum10;
            float inv4z = 0.5f / sqrtf(fourZSqr);
            q.f[0] = (M->m[0][2] + M->m[2][0]) * inv4z;
            q.f[1] = (M->m[1][2] + M->m[2][1]) * inv4z;
            q.f[2] = fourZSqr * inv4z;
            q.f[3] = (M->m[0][1] - M->m[1][0]) * inv4z;
        }
        else  // w^2 >= z^2
        {
            float fourWSqr = opr22 + sum10;
            float inv4w = 0.5f / sqrtf(fourWSqr);
            q.f[0] = (M->m[1][2] - M->m[2][1]) * inv4w;
            q.f[1] = (M->m[2][0] - M->m[0][2]) * inv4w;
            q.f[2] = (M->m[0][1] - M->m[1][0]) * inv4w;
            q.f[3] = fourWSqr * inv4w;
        }
    }
    return q.v;

#elif defined(_XM_SSE_INTRINSICS_)
    static const XMVECTORF32 XMPMMP = { { { +1.0f, -1.0f, -1.0f, +1.0f } } };
    static const XMVECTORF32 XMMPMP = { { { -1.0f, +1.0f, -1.0f, +1.0f } } };
    static const XMVECTORF32 XMMMPP = { { { -1.0f, -1.0f, +1.0f, +1.0f } } };

    XMVECTOR r0 = XM_MATRIX_GET(M,0);  // (r00, r01, r02, 0)
    XMVECTOR r1 = XM_MATRIX_GET(M,1);  // (r10, r11, r12, 0)
    XMVECTOR r2 = XM_MATRIX_GET(M,2);  // (r20, r21, r22, 0)

    // (r00, r00, r00, r00)
    XMVECTOR r00 = XM_PERMUTE_PS(r0, _MM_SHUFFLE(0, 0, 0, 0));
    // (r11, r11, r11, r11)
    XMVECTOR r11 = XM_PERMUTE_PS(r1, _MM_SHUFFLE(1, 1, 1, 1));
    // (r22, r22, r22, r22)
    XMVECTOR r22 = XM_PERMUTE_PS(r2, _MM_SHUFFLE(2, 2, 2, 2));

    // x^2 >= y^2 equivalent to r11 - r00 <= 0
    // (r11 - r00, r11 - r00, r11 - r00, r11 - r00)
    XMVECTOR r11mr00 = _mm_sub_ps(r11, r00);
    XMVECTOR x2gey2 = _mm_cmple_ps(r11mr00, g_XMZero.v);

    // z^2 >= w^2 equivalent to r11 + r00 <= 0
    // (r11 + r00, r11 + r00, r11 + r00, r11 + r00)
    XMVECTOR r11pr00 = _mm_add_ps(r11, r00);
    XMVECTOR z2gew2 = _mm_cmple_ps(r11pr00, g_XMZero.v);

    // x^2 + y^2 >= z^2 + w^2 equivalent to r22 <= 0
    XMVECTOR x2py2gez2pw2 = _mm_cmple_ps(r22, g_XMZero.v);

    // (4*x^2, 4*y^2, 4*z^2, 4*w^2)
    XMVECTOR t0 = XM_FMADD_PS(XMPMMP.v, r00, g_XMOne.v);
    XMVECTOR t1 = _mm_mul_ps(XMMPMP.v, r11);
    XMVECTOR t2 = XM_FMADD_PS(XMMMPP.v, r22, t0);
    XMVECTOR x2y2z2w2 = _mm_add_ps(t1, t2);

    // (r01, r02, r12, r11)
    t0 = _mm_shuffle_ps(r0, r1, _MM_SHUFFLE(1, 2, 2, 1));
    // (r10, r10, r20, r21)
    t1 = _mm_shuffle_ps(r1, r2, _MM_SHUFFLE(1, 0, 0, 0));
    // (r10, r20, r21, r10)
    t1 = XM_PERMUTE_PS(t1, _MM_SHUFFLE(1, 3, 2, 0));
    // (4*x*y, 4*x*z, 4*y*z, unused)
    XMVECTOR xyxzyz = _mm_add_ps(t0, t1);

    // (r21, r20, r10, r10)
    t0 = _mm_shuffle_ps(r2, r1, _MM_SHUFFLE(0, 0, 0, 1));
    // (r12, r12, r02, r01)
    t1 = _mm_shuffle_ps(r1, r0, _MM_SHUFFLE(1, 2, 2, 2));
    // (r12, r02, r01, r12)
    t1 = XM_PERMUTE_PS(t1, _MM_SHUFFLE(1, 3, 2, 0));
    // (4*x*w, 4*y*w, 4*z*w, unused)
    XMVECTOR xwywzw = _mm_sub_ps(t0, t1);
    xwywzw = _mm_mul_ps(XMMPMP.v, xwywzw);

    // (4*x^2, 4*y^2, 4*x*y, unused)
    t0 = _mm_shuffle_ps(x2y2z2w2, xyxzyz, _MM_SHUFFLE(0, 0, 1, 0));
    // (4*z^2, 4*w^2, 4*z*w, unused)
    t1 = _mm_shuffle_ps(x2y2z2w2, xwywzw, _MM_SHUFFLE(0, 2, 3, 2));
    // (4*x*z, 4*y*z, 4*x*w, 4*y*w)
    t2 = _mm_shuffle_ps(xyxzyz, xwywzw, _MM_SHUFFLE(1, 0, 2, 1));

    // (4*x*x, 4*x*y, 4*x*z, 4*x*w)
    XMVECTOR tensor0 = _mm_shuffle_ps(t0, t2, _MM_SHUFFLE(2, 0, 2, 0));
    // (4*y*x, 4*y*y, 4*y*z, 4*y*w)
    XMVECTOR tensor1 = _mm_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 1, 1, 2));
    // (4*z*x, 4*z*y, 4*z*z, 4*z*w)
    XMVECTOR tensor2 = _mm_shuffle_ps(t2, t1, _MM_SHUFFLE(2, 0, 1, 0));
    // (4*w*x, 4*w*y, 4*w*z, 4*w*w)
    XMVECTOR tensor3 = _mm_shuffle_ps(t2, t1, _MM_SHUFFLE(1, 2, 3, 2));

    // Select the row of the tensor-product matrix that has the largest
    // magnitude.
    t0 = _mm_and_ps(x2gey2, tensor0);
    t1 = _mm_andnot_ps(x2gey2, tensor1);
    t0 = _mm_or_ps(t0, t1);
    t1 = _mm_and_ps(z2gew2, tensor2);
    t2 = _mm_andnot_ps(z2gew2, tensor3);
    t1 = _mm_or_ps(t1, t2);
    t0 = _mm_and_ps(x2py2gez2pw2, t0);
    t1 = _mm_andnot_ps(x2py2gez2pw2, t1);
    t2 = _mm_or_ps(t0, t1);

    // Normalize the row.  No division by zero is possible because the
    // quaternion is unit-length (and the row is a nonzero multiple of
    // the quaternion).
    t0 = XMVector4Length(XM_REF_1V(t2));
    return _mm_div_ps(t2, t0);
#endif
}

inline XMVECTOR XM_CALLCONV XMPlaneNormalize(FXMVECTOR P)
{
#if defined(_XM_NO_INTRINSICS_)
    float fLengthSq = sqrtf((P->vector4_f32[0] * P->vector4_f32[0]) + (P->vector4_f32[1] * P->vector4_f32[1]) + (P->vector4_f32[2] * P->vector4_f32[2]));
    // Prevent divide by zero
    if (fLengthSq > 0)
    {
        fLengthSq = 1.0f / fLengthSq;
    }
    XMVECTORF32 vResult = { { {
            P->vector4_f32[0] * fLengthSq,
            P->vector4_f32[1] * fLengthSq,
            P->vector4_f32[2] * fLengthSq,
            P->vector4_f32[3] * fLengthSq
        } } };
    return vResult.v;
#elif defined(_XM_SSE4_INTRINSICS_)
    XMVECTOR vLengthSq = _mm_dp_ps(XM_DEREF_F(P), XM_DEREF_F(P), 0x7f);
    // Prepare for the division
    XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
    // Failsafe on zero (Or epsilon) length planes
    // If the length is infinity, set the elements to zero
    vLengthSq = _mm_cmpneq_ps(vLengthSq, g_XMInfinity.v);
    // Reciprocal mul to perform the normalization
    vResult = _mm_div_ps(XM_DEREF_F(P), vResult);
    // Any that are infinity, set to zero
    vResult = _mm_and_ps(vResult, vLengthSq);
    return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
    // Perform the dot product on x,y and z only
    XMVECTOR vLengthSq = _mm_mul_ps(XM_DEREF_F(P), XM_DEREF_F(P));
    XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 1, 2, 1));
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
    vLengthSq = _mm_add_ss(vLengthSq, vTemp);
    vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
    // Prepare for the division
    XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
    // Failsafe on zero (Or epsilon) length planes
    // If the length is infinity, set the elements to zero
    vLengthSq = _mm_cmpneq_ps(vLengthSq, g_XMInfinity.v);
    // Reciprocal mul to perform the normalization
    vResult = _mm_div_ps(XM_DEREF_F(P), vResult);
    // Any that are infinity, set to zero
    vResult = _mm_and_ps(vResult, vLengthSq);
    return vResult;
#endif
}