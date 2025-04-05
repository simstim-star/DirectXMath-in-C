#pragma once

//------------------------------------------------------------------------------
 // Comparison operations
 //------------------------------------------------------------------------------

#if !defined(_XM_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(push)
#pragma float_control(precise, on)
#endif

// Return true if any entry in the matrix is NaN
inline bool XM_CALLCONV XMMatrixIsNaN(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)
    size_t i = 16;
    const uint32_t* pWork = (const uint32_t*)(&M->m[0][0]);
    do {
        // Fetch value into integer unit
        uint32_t uTest = pWork[0];
        // Remove sign
        uTest &= 0x7FFFFFFFU;
        // NaN is 0x7F800001 through 0x7FFFFFFF inclusive
        uTest -= 0x7F800001U;
        if (uTest < 0x007FFFFFU)
        {
            break;      // NaN found
        }
        ++pWork;        // Next entry
    } while (--i);
    return (i != 0);      // i == 0 if nothing matched
#elif defined(_XM_SSE_INTRINSICS_)
    // Load in registers
    XMVECTOR vX = XM_MATRIX_GET(M, 0);
    XMVECTOR vY = XM_MATRIX_GET(M, 1);
    XMVECTOR vZ = XM_MATRIX_GET(M, 2);
    XMVECTOR vW = XM_MATRIX_GET(M, 3);
    // Test themselves to check for NaN
    vX = _mm_cmpneq_ps(vX, vX);
    vY = _mm_cmpneq_ps(vY, vY);
    vZ = _mm_cmpneq_ps(vZ, vZ);
    vW = _mm_cmpneq_ps(vW, vW);
    // Or all the results
    vX = _mm_or_ps(vX, vZ);
    vY = _mm_or_ps(vY, vW);
    vX = _mm_or_ps(vX, vY);
    // If any tested true, return true
    return (_mm_movemask_ps(vX) != 0);
#else
#endif
}

// Return true if any entry in the matrix is +/-INF
inline bool XM_CALLCONV XMMatrixIsInfinite(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)
    size_t i = 16;
    const uint32_t* pWork = (const uint32_t*)(&M->m[0][0]);
    do {
        // Fetch value into integer unit
        uint32_t uTest = pWork[0];
        // Remove sign
        uTest &= 0x7FFFFFFFU;
        // INF is 0x7F800000
        if (uTest == 0x7F800000U)
        {
            break;      // INF found
        }
        ++pWork;        // Next entry
    } while (--i);
    return (i != 0);      // i == 0 if nothing matched
#elif defined(_XM_SSE_INTRINSICS_)
    // Mask off the sign bits
    XMVECTOR vTemp1 = _mm_and_ps(XM_MATRIX_GET(M,0), g_XMAbsMask.v);
    XMVECTOR vTemp2 = _mm_and_ps(XM_MATRIX_GET(M,1), g_XMAbsMask.v);
    XMVECTOR vTemp3 = _mm_and_ps(XM_MATRIX_GET(M,2), g_XMAbsMask.v);
    XMVECTOR vTemp4 = _mm_and_ps(XM_MATRIX_GET(M,3), g_XMAbsMask.v);
    // Compare to infinity
    vTemp1 = _mm_cmpeq_ps(vTemp1, g_XMInfinity.v);
    vTemp2 = _mm_cmpeq_ps(vTemp2, g_XMInfinity.v);
    vTemp3 = _mm_cmpeq_ps(vTemp3, g_XMInfinity.v);
    vTemp4 = _mm_cmpeq_ps(vTemp4, g_XMInfinity.v);
    // Or the answers together
    vTemp1 = _mm_or_ps(vTemp1, vTemp2);
    vTemp3 = _mm_or_ps(vTemp3, vTemp4);
    vTemp1 = _mm_or_ps(vTemp1, vTemp3);
    // If any are infinity, the signs are true.
    return (_mm_movemask_ps(vTemp1) != 0);
#endif
}

// Return true if the XMMatrix is equal to identity
inline bool XM_CALLCONV XMMatrixIsIdentity(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)
    // Use the integer pipeline to reduce branching to a minimum
    const uint32_t* pWork = (const uint32_t*)(&M->m[0][0]);
    // Convert 1.0f to zero and or them together
    uint32_t uOne = pWork[0] ^ 0x3F800000U;
    // Or all the 0.0f entries together
    uint32_t uZero = pWork[1];
    uZero |= pWork[2];
    uZero |= pWork[3];
    // 2nd row
    uZero |= pWork[4];
    uOne |= pWork[5] ^ 0x3F800000U;
    uZero |= pWork[6];
    uZero |= pWork[7];
    // 3rd row
    uZero |= pWork[8];
    uZero |= pWork[9];
    uOne |= pWork[10] ^ 0x3F800000U;
    uZero |= pWork[11];
    // 4th row
    uZero |= pWork[12];
    uZero |= pWork[13];
    uZero |= pWork[14];
    uOne |= pWork[15] ^ 0x3F800000U;
    // If all zero entries are zero, the uZero==0
    uZero &= 0x7FFFFFFF;    // Allow -0.0f
    // If all 1.0f entries are 1.0f, then uOne==0
    uOne |= uZero;
    return (uOne == 0);
#elif defined(_XM_SSE_INTRINSICS_)
    XMVECTOR vTemp1 = _mm_cmpeq_ps(XM_MATRIX_GET(M,0), g_XMIdentityR0.v);
    XMVECTOR vTemp2 = _mm_cmpeq_ps(XM_MATRIX_GET(M,1), g_XMIdentityR1.v);
    XMVECTOR vTemp3 = _mm_cmpeq_ps(XM_MATRIX_GET(M,2), g_XMIdentityR2.v);
    XMVECTOR vTemp4 = _mm_cmpeq_ps(XM_MATRIX_GET(M,3), g_XMIdentityR3.v);
    vTemp1 = _mm_and_ps(vTemp1, vTemp2);
    vTemp3 = _mm_and_ps(vTemp3, vTemp4);
    vTemp1 = _mm_and_ps(vTemp1, vTemp3);
    return (_mm_movemask_ps(vTemp1) == 0x0f);
#endif
}

//------------------------------------------------------------------------------
// Computation operations
//------------------------------------------------------------------------------

// Perform a 4x4 matrix multiply by a 4x4 matrix
inline XMMATRIX XM_CALLCONV XMMatrixMultiply
(
    FXMMATRIX A,
    CXMMATRIX B
)
{
#if defined(_XM_NO_INTRINSICS_)
    XMMATRIX mResult;
    // Cache the invariants in registers
    float x = A->m[0][0];
    float y = A->m[0][1];
    float z = A->m[0][2];
    float w = A->m[0][3];
    // Perform the operation on the first row
    mResult.m[0][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[0][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[0][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[0][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    // Repeat for all the other rows
    x = A->m[1][0];
    y = A->m[1][1];
    z = A->m[1][2];
    w = A->m[1][3];
    mResult.m[1][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[1][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[1][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[1][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    x = A->m[2][0];
    y = A->m[2][1];
    z = A->m[2][2];
    w = A->m[2][3];
    mResult.m[2][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[2][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[2][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[2][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    x = A->m[3][0];
    y = A->m[3][1];
    z = A->m[3][2];
    w = A->m[3][3];
    mResult.m[3][0] = (B->m[0][0] * x) + (B->m[1][0] * y) + (B->m[2][0] * z) + (B->m[3][0] * w);
    mResult.m[3][1] = (B->m[0][1] * x) + (B->m[1][1] * y) + (B->m[2][1] * z) + (B->m[3][1] * w);
    mResult.m[3][2] = (B->m[0][2] * x) + (B->m[1][2] * y) + (B->m[2][2] * z) + (B->m[3][2] * w);
    mResult.m[3][3] = (B->m[0][3] * x) + (B->m[1][3] * y) + (B->m[2][3] * z) + (B->m[3][3] * w);
    return mResult;
#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX mResult;
    // If the first row of A is [X,Y,Z,W]
    // Splat it into registers
    // Use vW to hold the original row
    XMVECTOR vW = XM_MATRIX_GET(A,0);
    XMVECTOR vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0)); // X X X X
    XMVECTOR vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1)); // Y Y Y Y
    XMVECTOR vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2)); // Z Z Z Z
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3)); // W W W W
    // Perform the operation on the first row
    vX = _mm_mul_ps(vX, B->r[0]); // X*B00 X*B01 X*B02 X*B03
    vY = _mm_mul_ps(vY, B->r[1]); // Y*B10 Y*B11 Y*B12 Y*B13
    vZ = _mm_mul_ps(vZ, B->r[2]); // Z*B20 Z*B21 Z*B22 Z*B23
    vW = _mm_mul_ps(vW, B->r[3]); // W*B30 W*B31 W*B32 W*B33
    // Sum the partials above
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[0] = vX;
    // Repeat for the other 3 rows
    vW = XM_MATRIX_GET(A,1);
    vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
    vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
    vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
    vX = _mm_mul_ps(vX, B->r[0]);
    vY = _mm_mul_ps(vY, B->r[1]);
    vZ = _mm_mul_ps(vZ, B->r[2]);
    vW = _mm_mul_ps(vW, B->r[3]);
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[1] = vX;
    vW = XM_MATRIX_GET(A,2);
    vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
    vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
    vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
    vX = _mm_mul_ps(vX, B->r[0]);
    vY = _mm_mul_ps(vY, B->r[1]);
    vZ = _mm_mul_ps(vZ, B->r[2]);
    vW = _mm_mul_ps(vW, B->r[3]);
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[2] = vX;
    vW = XM_MATRIX_GET(A,3);
    vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
    vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
    vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
    vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
    vX = _mm_mul_ps(vX, B->r[0]);
    vY = _mm_mul_ps(vY, B->r[1]);
    vZ = _mm_mul_ps(vZ, B->r[2]);
    vW = _mm_mul_ps(vW, B->r[3]);
    vX = _mm_add_ps(vX, vZ);
    vY = _mm_add_ps(vY, vW);
    vX = _mm_add_ps(vX, vY);
    mResult.r[3] = vX;
    return mResult;
#endif
}

inline XMMATRIX XM_CALLCONV XMMatrixTranspose(FXMMATRIX M)
{
#if defined(_XM_NO_INTRINSICS_)

    // Original matrix:
    //
    //     m00m01m02m03
    //     m10m11m12m13
    //     m20m21m22m23
    //     m30m31m32m33

    XMMATRIX P;
    P.r[0] = XMVectorMergeXY(XM_MATRIX_GET(M,0), XM_MATRIX_GET(M,2)); // m00m20m01m21
    P.r[1] = XMVectorMergeXY(XM_MATRIX_GET(M,1), XM_MATRIX_GET(M,3)); // m10m30m11m31
    P.r[2] = XMVectorMergeZW(XM_MATRIX_GET(M,0), XM_MATRIX_GET(M,2)); // m02m22m03m23
    P.r[3] = XMVectorMergeZW(XM_MATRIX_GET(M,1), XM_MATRIX_GET(M,3)); // m12m32m13m33

    XMMATRIX MT;
    MT.r[0] = XMVectorMergeXY(P.r[0], P.r[1]); // m00m10m20m30
    MT.r[1] = XMVectorMergeZW(P.r[0], P.r[1]); // m01m11m21m31
    MT.r[2] = XMVectorMergeXY(P.r[2], P.r[3]); // m02m12m22m32
    MT.r[3] = XMVectorMergeZW(P.r[2], P.r[3]); // m03m13m23m33
    return MT;

#elif defined(_XM_SSE_INTRINSICS_)
    // x.x,x.y,y.x,y.y
    XMVECTOR vTemp1 = _mm_shuffle_ps(XM_MATRIX_GET(M,0), XM_MATRIX_GET(M,1), _MM_SHUFFLE(1, 0, 1, 0));
    // x.z,x.w,y.z,y.w
    XMVECTOR vTemp3 = _mm_shuffle_ps(XM_MATRIX_GET(M,0), XM_MATRIX_GET(M,1), _MM_SHUFFLE(3, 2, 3, 2));
    // z.x,z.y,w.x,w.y
    XMVECTOR vTemp2 = _mm_shuffle_ps(XM_MATRIX_GET(M,2), XM_MATRIX_GET(M,3), _MM_SHUFFLE(1, 0, 1, 0));
    // z.z,z.w,w.z,w.w
    XMVECTOR vTemp4 = _mm_shuffle_ps(XM_MATRIX_GET(M,2), XM_MATRIX_GET(M,3), _MM_SHUFFLE(3, 2, 3, 2));

    XMMATRIX mResult;
    // x.x,y.x,z.x,w.x
    mResult.r[0] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(2, 0, 2, 0));
    // x.y,y.y,z.y,w.y
    mResult.r[1] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(3, 1, 3, 1));
    // x.z,y.z,z.z,w.z
    mResult.r[2] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(2, 0, 2, 0));
    // x.w,y.w,z.w,w.w
    mResult.r[3] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(3, 1, 3, 1));
    return mResult;
#endif
}

inline XMMATRIX XM_CALLCONV XMMatrixTranslation
(
    float OffsetX,
    float OffsetY,
    float OffsetZ
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX M;
    M.m[0][0] = 1.0f;
    M.m[0][1] = 0.0f;
    M.m[0][2] = 0.0f;
    M.m[0][3] = 0.0f;

    M.m[1][0] = 0.0f;
    M.m[1][1] = 1.0f;
    M.m[1][2] = 0.0f;
    M.m[1][3] = 0.0f;

    M.m[2][0] = 0.0f;
    M.m[2][1] = 0.0f;
    M.m[2][2] = 1.0f;
    M.m[2][3] = 0.0f;

    M.m[3][0] = OffsetX;
    M.m[3][1] = OffsetY;
    M.m[3][2] = OffsetZ;
    M.m[3][3] = 1.0f;
    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX M;
    M.r[0] = g_XMIdentityR0.v;
    M.r[1] = g_XMIdentityR1.v;
    M.r[2] = g_XMIdentityR2.v;
    M.r[3] = XMVectorSet(OffsetX, OffsetY, OffsetZ, 1.f);
    return M;
#endif
}

inline XMMATRIX XM_CALLCONV XMMatrixLookToLH
(
    FXMVECTOR EyePosition,
    FXMVECTOR EyeDirection,
    FXMVECTOR UpDirection
)
{
    XMVECTOR Zero = XMVectorZero();
    assert(!XMVector3Equal(EyeDirection, XM_1V(g_XMZero.v)));
    assert(!XMVector3IsInfinite(EyeDirection));
    assert(!XMVector3Equal(UpDirection, XM_1V(g_XMZero.v)));
    assert(!XMVector3IsInfinite(UpDirection));

    XMVECTOR R2 = XMVector3Normalize(EyeDirection);

    XMVECTOR R0 = XMVector3Cross(UpDirection, XM_1V(R2));
    R0 = XMVector3Normalize(XM_1V(R0));

    XMVECTOR R1 = XMVector3Cross(XM_2V(R2, R0));

    XMVECTOR NegEyePosition = XMVectorNegate(EyePosition);

    XMVECTOR D0 = XMVector3Dot(XM_2V(R0, NegEyePosition));
    XMVECTOR D1 = XMVector3Dot(XM_2V(R1, NegEyePosition));
    XMVECTOR D2 = XMVector3Dot(XM_2V(R2, NegEyePosition));

    XMMATRIX M;
    M.r[0] = XMVectorSelect(XM_3V(D0, R0, g_XMSelect1110.v));
    M.r[1] = XMVectorSelect(XM_3V(D1, R1, g_XMSelect1110.v));
    M.r[2] = XMVectorSelect(XM_3V(D2, R2, g_XMSelect1110.v));
    M.r[3] = g_XMIdentityR3.v;

    M = XMMatrixTranspose(XM_1M(M));

    return M;
}

inline XMMATRIX XM_CALLCONV XMMatrixLookToRH
(
    FXMVECTOR EyePosition,
    FXMVECTOR EyeDirection,
    FXMVECTOR UpDirection
)
{
    XMVECTOR NegEyeDirection = XMVectorNegate(EyeDirection);
    return XMMatrixLookToLH(EyePosition, XM_1V(NegEyeDirection), UpDirection);
}

inline XMMATRIX XM_CALLCONV XMMatrixPerspectiveFovRH
(
    float FovAngleY,
    float AspectRatio,
    float NearZ,
    float FarZ
)
{
    assert(NearZ > 0.f && FarZ > 0.f);
    assert(!XMScalarNearEqual(FovAngleY, 0.0f, 0.00001f * 2.0f));
    assert(!XMScalarNearEqual(AspectRatio, 0.0f, 0.00001f));
    assert(!XMScalarNearEqual(FarZ, NearZ, 0.00001f));

#if defined(_XM_NO_INTRINSICS_)

    float    SinFov;
    float    CosFov;
    XMScalarSinCos(&SinFov, &CosFov, 0.5f * FovAngleY);

    float Height = CosFov / SinFov;
    float Width = Height / AspectRatio;
    float fRange = FarZ / (NearZ - FarZ);

    XMMATRIX M;
    M.m[0][0] = Width;
    M.m[0][1] = 0.0f;
    M.m[0][2] = 0.0f;
    M.m[0][3] = 0.0f;

    M.m[1][0] = 0.0f;
    M.m[1][1] = Height;
    M.m[1][2] = 0.0f;
    M.m[1][3] = 0.0f;

    M.m[2][0] = 0.0f;
    M.m[2][1] = 0.0f;
    M.m[2][2] = fRange;
    M.m[2][3] = -1.0f;

    M.m[3][0] = 0.0f;
    M.m[3][1] = 0.0f;
    M.m[3][2] = fRange * NearZ;
    M.m[3][3] = 0.0f;
    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    float    SinFov;
    float    CosFov;
    XMScalarSinCos(&SinFov, &CosFov, 0.5f * FovAngleY);
    float fRange = FarZ / (NearZ - FarZ);
    // Note: This is recorded on the stack
    float Height = CosFov / SinFov;
    XMVECTOR rMem = {
        Height / AspectRatio,
        Height,
        fRange,
        fRange * NearZ
    };
    // Copy from memory to SSE register
    XMVECTOR vValues = rMem;
    XMVECTOR vTemp = _mm_setzero_ps();
    // Copy x only
    vTemp = _mm_move_ss(vTemp, vValues);
    // Height / AspectRatio,0,0,0
    XMMATRIX M;
    M.r[0] = vTemp;
    // 0,Height,0,0
    vTemp = vValues;
    vTemp = _mm_and_ps(vTemp, g_XMMaskY.v);
    M.r[1] = vTemp;
    // x=fRange,y=-fRange * NearZ,0,-1.0f
    vTemp = _mm_setzero_ps();
    vValues = _mm_shuffle_ps(vValues, g_XMNegIdentityR3.v, _MM_SHUFFLE(3, 2, 3, 2));
    // 0,0,fRange,-1.0f
    vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
    M.r[2] = vTemp;
    // 0,0,fRange * NearZ,0.0f
    vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
    M.r[3] = vTemp;
    return M;
#endif
}