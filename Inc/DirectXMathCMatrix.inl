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

//------------------------------------------------------------------------------
// Return the inverse and the determinant of a 4x4 matrix
_Use_decl_annotations_
inline XMMATRIX XM_CALLCONV XMMatrixInverse
(
    XMVECTOR* pDeterminant,
    FXMMATRIX  M
)
{
#if defined(_XM_NO_INTRINSICS_)

    XMMATRIX MT = XMMatrixTranspose(M);

    XMVECTOR V0[4], V1[4];
    V0[0] = XM_VEC_SWIZZLE(MT.r[2], XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_Y, XM_SWIZZLE_Y);
    V1[0] = XM_VEC_SWIZZLE(MT.r[3], XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Z, XM_SWIZZLE_W);
    V0[1] = XM_VEC_SWIZZLE(MT.r[0], XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_Y, XM_SWIZZLE_Y);
    V1[1] = XM_VEC_SWIZZLE(MT.r[1], XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Z, XM_SWIZZLE_W);
    V0[2] = XM_VEC_PERMUTE(MT.r[2], MT.r[0],XM_PERMUTE_0X, XM_PERMUTE_0Z, XM_PERMUTE_1X, XM_PERMUTE_1Z);
    V1[2] = XM_VEC_PERMUTE(MT.r[3], MT.r[1],XM_PERMUTE_0Y, XM_PERMUTE_0W, XM_PERMUTE_1Y, XM_PERMUTE_1W);

    XMVECTOR D0 = XM_VEC_MULT(V0[0], V1[0]);
    XMVECTOR D1 = XM_VEC_MULT(V0[1], V1[1]);
    XMVECTOR D2 = XM_VEC_MULT(V0[2], V1[2]);

    V0[0] = XM_VEC_SWIZZLE(MT.r[2], XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Z, XM_SWIZZLE_W);
    V1[0] = XM_VEC_SWIZZLE(MT.r[3], XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_Y, XM_SWIZZLE_Y);
    V0[1] = XM_VEC_SWIZZLE(MT.r[0], XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Z, XM_SWIZZLE_W);
    V1[1] = XM_VEC_SWIZZLE(MT.r[1], XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_Y, XM_SWIZZLE_Y);
    V0[2] = XM_VEC_PERMUTE(MT.r[2], MT.r[0], XM_PERMUTE_0Y, XM_PERMUTE_0W, XM_PERMUTE_1Y, XM_PERMUTE_1W);
    V1[2] = XM_VEC_PERMUTE(MT.r[3], MT.r[1], XM_PERMUTE_0X, XM_PERMUTE_0Z, XM_PERMUTE_1X, XM_PERMUTE_1Z);

    D0 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[0], V1[0], D0);
    D1 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[1], V1[1], D1);
    D2 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[2], V1[2], D2);

    V0[0] = XM_VEC_SWIZZLE(MT.r[1], XM_SWIZZLE_Y, XM_SWIZZLE_Z, XM_SWIZZLE_X, XM_SWIZZLE_Y);
    V1[0] = XM_VEC_PERMUTE(D0, D2, XM_PERMUTE_1Y, XM_PERMUTE_0Y, XM_PERMUTE_0W, XM_PERMUTE_0X);
    V0[1] = XM_VEC_SWIZZLE(MT.r[0], XM_SWIZZLE_Z, XM_SWIZZLE_X, XM_SWIZZLE_Y, XM_SWIZZLE_X);
    V1[1] = XM_VEC_PERMUTE(D0, D2, XM_PERMUTE_0W, XM_PERMUTE_1Y, XM_PERMUTE_0Y, XM_PERMUTE_0Z);
    V0[2] = XM_VEC_SWIZZLE(MT.r[3], XM_SWIZZLE_Y, XM_SWIZZLE_Z, XM_SWIZZLE_X, XM_SWIZZLE_Y);
    V1[2] = XM_VEC_PERMUTE(D1, D2, XM_PERMUTE_1W, XM_PERMUTE_0Y, XM_PERMUTE_0W, XM_PERMUTE_0X);
    V0[3] = XM_VEC_SWIZZLE(MT.r[2], XM_SWIZZLE_Z, XM_SWIZZLE_X, XM_SWIZZLE_Y, XM_SWIZZLE_X);
    V1[3] = XM_VEC_PERMUTE(D1, D2, XM_PERMUTE_0W, XM_PERMUTE_1W, XM_PERMUTE_0Y, XM_PERMUTE_0Z);

    XMVECTOR C0 = XM_VEC_MULT(V0[0], V1[0]);
    XMVECTOR C2 = XM_VEC_MULT(V0[1], V1[1]);
    XMVECTOR C4 = XM_VEC_MULT(V0[2], V1[2]);
    XMVECTOR C6 = XM_VEC_MULT(V0[3], V1[3]);

    V0[0] = XM_VEC_SWIZZLE(MT.r[1], XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Y, XM_SWIZZLE_Z);
    V1[0] = XM_VEC_PERMUTE(D0, D2, XM_PERMUTE_0W, XM_PERMUTE_0X, XM_PERMUTE_0Y, XM_PERMUTE_1X);
    V0[1] = XM_VEC_SWIZZLE(MT.r[0], XM_SWIZZLE_W, XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Y);
    V1[1] = XM_VEC_PERMUTE(D0, D2, XM_PERMUTE_0Z, XM_PERMUTE_0Y, XM_PERMUTE_1X, XM_PERMUTE_0X);
    V0[2] = XM_VEC_SWIZZLE(MT.r[3], XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Y, XM_SWIZZLE_Z);
    V1[2] = XM_VEC_PERMUTE(D1, D2, XM_PERMUTE_0W, XM_PERMUTE_0X, XM_PERMUTE_0Y, XM_PERMUTE_1Z);
    V0[3] = XM_VEC_SWIZZLE(MT.r[2], XM_SWIZZLE_W, XM_SWIZZLE_Z, XM_SWIZZLE_W, XM_SWIZZLE_Y);
    V1[3] = XM_VEC_PERMUTE(D1, D2, XM_PERMUTE_0Z, XM_PERMUTE_0Y, XM_PERMUTE_1Z, XM_PERMUTE_0X);

    C0 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[0], V1[0], C0);
    C2 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[1], V1[1], C2);
    C4 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[2], V1[2], C4);
    C6 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[3], V1[3], C6);

    V0[0] = XM_VEC_SWIZZLE(MT.r[1], XM_SWIZZLE_W, XM_SWIZZLE_X, XM_SWIZZLE_W, XM_SWIZZLE_X);
    V1[0] = XM_VEC_PERMUTE(D0, D2, XM_PERMUTE_0Z, XM_PERMUTE_1Y, XM_PERMUTE_1X, XM_PERMUTE_0Z);
    V0[1] = XM_VEC_SWIZZLE(MT.r[0], XM_SWIZZLE_Y, XM_SWIZZLE_W, XM_SWIZZLE_X, XM_SWIZZLE_Z);
    V1[1] = XM_VEC_PERMUTE(D0, D2, XM_PERMUTE_1Y, XM_PERMUTE_0X, XM_PERMUTE_0W, XM_PERMUTE_1X);
    V0[2] = XM_VEC_SWIZZLE(MT.r[3], XM_SWIZZLE_W, XM_SWIZZLE_X, XM_SWIZZLE_W, XM_SWIZZLE_X);
    V1[2] = XM_VEC_PERMUTE(D1, D2, XM_PERMUTE_0Z, XM_PERMUTE_1W, XM_PERMUTE_1Z, XM_PERMUTE_0Z);
    V0[3] = XM_VEC_SWIZZLE(MT.r[2], XM_SWIZZLE_Y, XM_SWIZZLE_W, XM_SWIZZLE_X, XM_SWIZZLE_Z);
    V1[3] = XM_VEC_PERMUTE(D1, D2, XM_PERMUTE_1W, XM_PERMUTE_0X, XM_PERMUTE_0W, XM_PERMUTE_1Z);

    XMVECTOR C1 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[0], V1[0], C0);
    C0 = XM_VEC_MULT_ADD(V0[0], V1[0], C0);
    XMVECTOR C3 = XM_VEC_MULT_ADD(V0[1], V1[1], C2);
    C2 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[1], V1[1], C2);
    XMVECTOR C5 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[2], V1[2], C4);
    C4 = XM_VEC_MULT_ADD(V0[2], V1[2], C4);
    XMVECTOR C7 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[3], V1[3], C6);
    C6 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0[3], V1[3], C6);

    XMMATRIX R;
    R.r[0] = XM_VEC_SELECT(C0, C1, g_XMSelect0101.v);
    R.r[1] = XM_VEC_SELECT(C2, C3, g_XMSelect0101.v);
    R.r[2] = XM_VEC_SELECT(C4, C5, g_XMSelect0101.v);
    R.r[3] = XM_VEC_SELECT(C6, C7, g_XMSelect0101.v);

    XMVECTOR Determinant = XM_VEC4_DOT(R.r[0], MT.r[0]);

    if (pDeterminant != nullptr)
        *pDeterminant = Determinant;

    XMVECTOR Reciprocal = XM_VEC_RECIPROCAL(Determinant);

    XMMATRIX Result;
    Result.r[0] = XM_VEC_MULT(R.r[0], Reciprocal);
    Result.r[1] = XM_VEC_MULT(R.r[1], Reciprocal);
    Result.r[2] = XM_VEC_MULT(R.r[2], Reciprocal);
    Result.r[3] = XM_VEC_MULT(R.r[3], Reciprocal);
    return Result;

#elif defined(_XM_SSE_INTRINSICS_)
    // Transpose matrix
    XMVECTOR vTemp1 = _mm_shuffle_ps(XM_DEREF_MATRIX(M).r[0], XM_DEREF_MATRIX(M).r[1], _MM_SHUFFLE(1, 0, 1, 0));
    XMVECTOR vTemp3 = _mm_shuffle_ps(XM_DEREF_MATRIX(M).r[0], XM_DEREF_MATRIX(M).r[1], _MM_SHUFFLE(3, 2, 3, 2));
    XMVECTOR vTemp2 = _mm_shuffle_ps(XM_DEREF_MATRIX(M).r[2], XM_DEREF_MATRIX(M).r[3], _MM_SHUFFLE(1, 0, 1, 0));
    XMVECTOR vTemp4 = _mm_shuffle_ps(XM_DEREF_MATRIX(M).r[2], XM_DEREF_MATRIX(M).r[3], _MM_SHUFFLE(3, 2, 3, 2));

    XMMATRIX MT;
    MT.r[0] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(2, 0, 2, 0));
    MT.r[1] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(3, 1, 3, 1));
    MT.r[2] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(2, 0, 2, 0));
    MT.r[3] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(3, 1, 3, 1));

    XMVECTOR V00 = XM_PERMUTE_PS(MT.r[2], _MM_SHUFFLE(1, 1, 0, 0));
    XMVECTOR V10 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(3, 2, 3, 2));
    XMVECTOR V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(1, 1, 0, 0));
    XMVECTOR V11 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(3, 2, 3, 2));
    XMVECTOR V02 = _mm_shuffle_ps(MT.r[2], MT.r[0], _MM_SHUFFLE(2, 0, 2, 0));
    XMVECTOR V12 = _mm_shuffle_ps(MT.r[3], MT.r[1], _MM_SHUFFLE(3, 1, 3, 1));

    XMVECTOR D0 = _mm_mul_ps(V00, V10);
    XMVECTOR D1 = _mm_mul_ps(V01, V11);
    XMVECTOR D2 = _mm_mul_ps(V02, V12);

    V00 = XM_PERMUTE_PS(MT.r[2], _MM_SHUFFLE(3, 2, 3, 2));
    V10 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(1, 1, 0, 0));
    V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(3, 2, 3, 2));
    V11 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(1, 1, 0, 0));
    V02 = _mm_shuffle_ps(MT.r[2], MT.r[0], _MM_SHUFFLE(3, 1, 3, 1));
    V12 = _mm_shuffle_ps(MT.r[3], MT.r[1], _MM_SHUFFLE(2, 0, 2, 0));

    D0 = XM_FNMADD_PS(V00, V10, D0);
    D1 = XM_FNMADD_PS(V01, V11, D1);
    D2 = XM_FNMADD_PS(V02, V12, D2);
    // V11 = D0Y,D0W,D2Y,D2Y
    V11 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(1, 1, 3, 1));
    V00 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(1, 0, 2, 1));
    V10 = _mm_shuffle_ps(V11, D0, _MM_SHUFFLE(0, 3, 0, 2));
    V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(0, 1, 0, 2));
    V11 = _mm_shuffle_ps(V11, D0, _MM_SHUFFLE(2, 1, 2, 1));
    // V13 = D1Y,D1W,D2W,D2W
    XMVECTOR V13 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(3, 3, 3, 1));
    V02 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(1, 0, 2, 1));
    V12 = _mm_shuffle_ps(V13, D1, _MM_SHUFFLE(0, 3, 0, 2));
    XMVECTOR V03 = XM_PERMUTE_PS(MT.r[2], _MM_SHUFFLE(0, 1, 0, 2));
    V13 = _mm_shuffle_ps(V13, D1, _MM_SHUFFLE(2, 1, 2, 1));

    XMVECTOR C0 = _mm_mul_ps(V00, V10);
    XMVECTOR C2 = _mm_mul_ps(V01, V11);
    XMVECTOR C4 = _mm_mul_ps(V02, V12);
    XMVECTOR C6 = _mm_mul_ps(V03, V13);

    // V11 = D0X,D0Y,D2X,D2X
    V11 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(0, 0, 1, 0));
    V00 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(2, 1, 3, 2));
    V10 = _mm_shuffle_ps(D0, V11, _MM_SHUFFLE(2, 1, 0, 3));
    V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(1, 3, 2, 3));
    V11 = _mm_shuffle_ps(D0, V11, _MM_SHUFFLE(0, 2, 1, 2));
    // V13 = D1X,D1Y,D2Z,D2Z
    V13 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(2, 2, 1, 0));
    V02 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(2, 1, 3, 2));
    V12 = _mm_shuffle_ps(D1, V13, _MM_SHUFFLE(2, 1, 0, 3));
    V03 = XM_PERMUTE_PS(MT.r[2], _MM_SHUFFLE(1, 3, 2, 3));
    V13 = _mm_shuffle_ps(D1, V13, _MM_SHUFFLE(0, 2, 1, 2));

    C0 = XM_FNMADD_PS(V00, V10, C0);
    C2 = XM_FNMADD_PS(V01, V11, C2);
    C4 = XM_FNMADD_PS(V02, V12, C4);
    C6 = XM_FNMADD_PS(V03, V13, C6);

    V00 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(0, 3, 0, 3));
    // V10 = D0Z,D0Z,D2X,D2Y
    V10 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(1, 0, 2, 2));
    V10 = XM_PERMUTE_PS(V10, _MM_SHUFFLE(0, 2, 3, 0));
    V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(2, 0, 3, 1));
    // V11 = D0X,D0W,D2X,D2Y
    V11 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(1, 0, 3, 0));
    V11 = XM_PERMUTE_PS(V11, _MM_SHUFFLE(2, 1, 0, 3));
    V02 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(0, 3, 0, 3));
    // V12 = D1Z,D1Z,D2Z,D2W
    V12 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(3, 2, 2, 2));
    V12 = XM_PERMUTE_PS(V12, _MM_SHUFFLE(0, 2, 3, 0));
    V03 = XM_PERMUTE_PS(MT.r[2], _MM_SHUFFLE(2, 0, 3, 1));
    // V13 = D1X,D1W,D2Z,D2W
    V13 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(3, 2, 3, 0));
    V13 = XM_PERMUTE_PS(V13, _MM_SHUFFLE(2, 1, 0, 3));

    V00 = _mm_mul_ps(V00, V10);
    V01 = _mm_mul_ps(V01, V11);
    V02 = _mm_mul_ps(V02, V12);
    V03 = _mm_mul_ps(V03, V13);
    XMVECTOR C1 = _mm_sub_ps(C0, V00);
    C0 = _mm_add_ps(C0, V00);
    XMVECTOR C3 = _mm_add_ps(C2, V01);
    C2 = _mm_sub_ps(C2, V01);
    XMVECTOR C5 = _mm_sub_ps(C4, V02);
    C4 = _mm_add_ps(C4, V02);
    XMVECTOR C7 = _mm_add_ps(C6, V03);
    C6 = _mm_sub_ps(C6, V03);

    C0 = _mm_shuffle_ps(C0, C1, _MM_SHUFFLE(3, 1, 2, 0));
    C2 = _mm_shuffle_ps(C2, C3, _MM_SHUFFLE(3, 1, 2, 0));
    C4 = _mm_shuffle_ps(C4, C5, _MM_SHUFFLE(3, 1, 2, 0));
    C6 = _mm_shuffle_ps(C6, C7, _MM_SHUFFLE(3, 1, 2, 0));
    C0 = XM_PERMUTE_PS(C0, _MM_SHUFFLE(3, 1, 2, 0));
    C2 = XM_PERMUTE_PS(C2, _MM_SHUFFLE(3, 1, 2, 0));
    C4 = XM_PERMUTE_PS(C4, _MM_SHUFFLE(3, 1, 2, 0));
    C6 = XM_PERMUTE_PS(C6, _MM_SHUFFLE(3, 1, 2, 0));
    // Get the determinant
    XMVECTOR vTemp = XM_VEC4_DOT(C0, MT.r[0]);
    if (pDeterminant != NULL)
        *pDeterminant = vTemp;
    vTemp = _mm_div_ps(g_XMOne.v, vTemp);
    XMMATRIX mResult;
    mResult.r[0] = _mm_mul_ps(C0, vTemp);
    mResult.r[1] = _mm_mul_ps(C2, vTemp);
    mResult.r[2] = _mm_mul_ps(C4, vTemp);
    mResult.r[3] = _mm_mul_ps(C6, vTemp);
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
    assert(!XMVector3Equal(EyeDirection, XM_REF_1V(g_XMZero.v)));
    assert(!XMVector3IsInfinite(EyeDirection));
    assert(!XMVector3Equal(UpDirection, XM_REF_1V(g_XMZero.v)));
    assert(!XMVector3IsInfinite(UpDirection));

    XMVECTOR R2 = XMVector3Normalize(EyeDirection);

    XMVECTOR R0 = XMVector3Cross(UpDirection, XM_REF_1V(R2));
    R0 = XM_VEC3_NORM(R0);

    XMVECTOR R1 = XM_VEC3_CROSS(R2, R0);

    XMVECTOR NegEyePosition = XMVectorNegate(EyePosition);

    XMVECTOR D0 = XM_VEC3_DOT(R0, NegEyePosition);
    XMVECTOR D1 = XM_VEC3_DOT(R1, NegEyePosition);
    XMVECTOR D2 = XM_VEC3_DOT(R2, NegEyePosition);

    XMMATRIX M;
    M.r[0] = XM_VEC_SELECT(D0, R0, g_XMSelect1110.v);
    M.r[1] = XM_VEC_SELECT(D1, R1, g_XMSelect1110.v);
    M.r[2] = XM_VEC_SELECT(D2, R2, g_XMSelect1110.v);
    M.r[3] = g_XMIdentityR3.v;

    return XM_MAT_TRANSP(M);
}

inline XMMATRIX XM_CALLCONV XMMatrixLookToRH
(
    FXMVECTOR EyePosition,
    FXMVECTOR EyeDirection,
    FXMVECTOR UpDirection
)
{
    XMVECTOR NegEyeDirection = XMVectorNegate(EyeDirection);
    return XMMatrixLookToLH(EyePosition, XM_REF_1V(NegEyeDirection), UpDirection);
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

#define XM3RANKDECOMPOSE(a, b, c, x, y, z)      \
    if((x) < (y))                   \
    {                               \
        if((y) < (z))               \
        {                           \
            (a) = 2;                \
            (b) = 1;                \
            (c) = 0;                \
        }                           \
        else                        \
        {                           \
            (a) = 1;                \
                                    \
            if((x) < (z))           \
            {                       \
                (b) = 2;            \
                (c) = 0;            \
            }                       \
            else                    \
            {                       \
                (b) = 0;            \
                (c) = 2;            \
            }                       \
        }                           \
    }                               \
    else                            \
    {                               \
        if((x) < (z))               \
        {                           \
            (a) = 2;                \
            (b) = 0;                \
            (c) = 1;                \
        }                           \
        else                        \
        {                           \
            (a) = 0;                \
                                    \
            if((y) < (z))           \
            {                       \
                (b) = 2;            \
                (c) = 1;            \
            }                       \
            else                    \
            {                       \
                (b) = 1;            \
                (c) = 2;            \
            }                       \
        }                           \
    }

#define XM3_DECOMP_EPSILON 0.0001f

_Use_decl_annotations_
inline bool XM_CALLCONV XMMatrixDecompose
(
    XMVECTOR* outScale,
    XMVECTOR* outRotQuat,
    XMVECTOR* outTrans,
    FXMMATRIX M
)
{
    static const XMVECTOR* pvCanonicalBasis[3] = {
        &g_XMIdentityR0.v,
        &g_XMIdentityR1.v,
        &g_XMIdentityR2.v
    };

    assert(outScale != NULL);
    assert(outRotQuat != NULL);
    assert(outTrans != NULL);

    // Get the translation
    outTrans[0] = XM_DEREF_MATRIX(M).r[3];

    XMVECTOR* ppvBasis[3];
    XMMATRIX matTemp;
    ppvBasis[0] = &matTemp.r[0];
    ppvBasis[1] = &matTemp.r[1];
    ppvBasis[2] = &matTemp.r[2];

    matTemp.r[0] = XM_DEREF_MATRIX(M).r[0];
    matTemp.r[1] = XM_DEREF_MATRIX(M).r[1];
    matTemp.r[2] = XM_DEREF_MATRIX(M).r[2];
    matTemp.r[3] = g_XMIdentityR3.v;

    float* pfScales = (float*)(outScale);

    size_t a, b, c;
    XMVECTOR len1 = XM_VEC3_LEN(ppvBasis[0][0]);
    XMVECTOR len2 = XM_VEC3_LEN(ppvBasis[1][0]);
    XMVECTOR len3 = XM_VEC3_LEN(ppvBasis[2][0]);
    XM_VECX_PTR(&pfScales[0], len1);
    XM_VECX_PTR(&pfScales[1], len2);
    XM_VECX_PTR(&pfScales[2], len3);
    pfScales[3] = 0.f;

    XM3RANKDECOMPOSE(a, b, c, pfScales[0], pfScales[1], pfScales[2])

        if (pfScales[a] < XM3_DECOMP_EPSILON)
        {
            ppvBasis[a][0] = pvCanonicalBasis[a][0];
        }
    ppvBasis[a][0] = XM_VEC3_NORM(ppvBasis[a][0]);

    if (pfScales[b] < XM3_DECOMP_EPSILON)
    {
        size_t aa, bb, cc;
        float fAbsX, fAbsY, fAbsZ;

        fAbsX = fabsf(XM_VECX(ppvBasis[a][0]));
        fAbsY = fabsf(XM_VECY(ppvBasis[a][0]));
        fAbsZ = fabsf(XM_VECZ(ppvBasis[a][0]));

        XM3RANKDECOMPOSE(aa, bb, cc, fAbsX, fAbsY, fAbsZ)

            ppvBasis[b][0] = XM_VEC3_CROSS(ppvBasis[a][0], pvCanonicalBasis[cc][0]);
    }

    ppvBasis[b][0] = XM_VEC3_NORM(ppvBasis[b][0]);

    if (pfScales[c] < XM3_DECOMP_EPSILON)
    {
        ppvBasis[c][0] = XM_VEC3_CROSS(ppvBasis[a][0], ppvBasis[b][0]);
    }

    ppvBasis[c][0] = XM_VEC3_NORM(ppvBasis[c][0]);

    XMVECTOR det = XM_DETERMINANT(matTemp);
    float fDet = XM_VECX(det);

    // use Kramer's rule to check for handedness of coordinate system
    if (fDet < 0.0f)
    {
        // switch coordinate system by negating the scale and inverting the basis vector on the x-axis
        pfScales[a] = -pfScales[a];
        ppvBasis[a][0] = XM_VEC_NEG(ppvBasis[a][0]);

        fDet = -fDet;
    }

    fDet -= 1.0f;
    fDet *= fDet;

    if (XM3_DECOMP_EPSILON < fDet)
    {
        // Non-SRT matrix encountered
        return false;
    }

    // generate the quaternion from the matrix
    outRotQuat[0] = XM_QUAT_ROT_MAT(matTemp);
    return true;
}

/***************************
 * M = 
 *     m00, m01, m02, m03
 *     m10, m11, m12, m13
 *     m20, m21, m22, m23
 *     m30, m31, m32, m33
 ***************************/
inline XMVECTOR XM_CALLCONV XMMatrixDeterminant(FXMMATRIX M) {
    static const XMVECTORF32 Sign = { { { 1.0f, -1.0f, 1.0f, -1.0f } } };
    
    XMVECTOR V0 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[2], XM_SWIZZLE_Y, XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_X); // [m21, m20, m20, m20]
    XMVECTOR V1 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[3], XM_SWIZZLE_Z, XM_SWIZZLE_Z, XM_SWIZZLE_Y, XM_SWIZZLE_Y); // [m32, m32, m31, m31]
    XMVECTOR V2 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[2], XM_SWIZZLE_Y, XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_X); // [m21, m20, m20, m20]
    XMVECTOR V3 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[3], XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_Z); // [m33, m33, m33, m32]
    XMVECTOR V4 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[2], XM_SWIZZLE_Z, XM_SWIZZLE_Z, XM_SWIZZLE_Y, XM_SWIZZLE_Y); // [m22, m22, m21, m21]
    XMVECTOR V5 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[3], XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_Z); // [m33, m33, m33, m32]


    XMVECTOR P0 = XM_VEC_MULT(V0, V1); // [m21 * m32, m20 * m32, m20 * m31, m20 * m31]
    XMVECTOR P1 = XM_VEC_MULT(V2, V3); // [m21 * m33, m20 * m33, m20 * m33, m20 * m32]
    XMVECTOR P2 = XM_VEC_MULT(V4, V5); // [m22 * m33, m22 * m33, m21 * m33, m21 * m32]

    V0 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[2], XM_SWIZZLE_Z, XM_SWIZZLE_Z, XM_SWIZZLE_Y, XM_SWIZZLE_Y); // [m22, m22, m21, m21]
    V1 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[3], XM_SWIZZLE_Y, XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_X); // [m31, m30, m30, m30]
    V2 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[2], XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_Z); // [m23, m23, m23, m22]
    V3 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[3], XM_SWIZZLE_Y, XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_X); // [m31, m30, m30, m30]
    V4 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[2], XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_Z); // [m23, m23, m23, m22]
    V5 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[3], XM_SWIZZLE_Z, XM_SWIZZLE_Z, XM_SWIZZLE_Y, XM_SWIZZLE_Y); // [m32, m32, m31, m31]

    P0 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V0, V1, P0);
    P1 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V2, V3, P1);
    P2 = XM_VEC_NEGATIVE_MULT_SUBTRACT(V4, V5, P2);
    
    V0 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[1], XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_W, XM_SWIZZLE_Z); // [m13, m13, m13, m12]
    V1 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[1], XM_SWIZZLE_Z, XM_SWIZZLE_Z, XM_SWIZZLE_Y, XM_SWIZZLE_Y); // [m12, m12, m11, m11]
    V2 = XM_VEC_SWIZZLE(XM_DEREF_MATRIX(M).r[1], XM_SWIZZLE_Y, XM_SWIZZLE_X, XM_SWIZZLE_X, XM_SWIZZLE_X); // [m11, m10, m10, m10]

    XMVECTOR S = XM_VEC_MULT(XM_DEREF_MATRIX(M).r[0], Sign.v);
    XMVECTOR R = XM_VEC_MULT(V0, P0);
    R = XM_VEC_NEGATIVE_MULT_SUBTRACT(V1, P1, R);
    R = XM_VEC_MULT_ADD(V2, P2, R);

    return XM_VEC4_DOT(S, R);
}


inline XMMATRIX XMMatrixScale(FXMMATRIX M, float S)
{
    XMMATRIX mResult;
    mResult.r[0] = XMVectorScale(XM_REF_1V(XM_DEREF_MATRIX(M).r[0]), S);
    mResult.r[1] = XMVectorScale(XM_REF_1V(XM_DEREF_MATRIX(M).r[1]), S);
    mResult.r[2] = XMVectorScale(XM_REF_1V(XM_DEREF_MATRIX(M).r[2]), S);
    mResult.r[3] = XMVectorScale(XM_REF_1V(XM_DEREF_MATRIX(M).r[3]), S);
    return mResult;
}

inline XMMATRIX XM_CALLCONV XMMatrixTranslationFromVector(FXMVECTOR Offset)
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

    M.m[3][0] = Offset->vector4_f32[0];
    M.m[3][1] = Offset->vector4_f32[1];
    M.m[3][2] = Offset->vector4_f32[2];
    M.m[3][3] = 1.0f;
    return M;

#elif defined(_XM_SSE_INTRINSICS_)
    XMMATRIX M;
    M.r[0] = g_XMIdentityR0.v;
    M.r[1] = g_XMIdentityR1.v;
    M.r[2] = g_XMIdentityR2.v;
    M.r[3] = XMVectorSelect(XM_REF_1V(g_XMIdentityR3.v), Offset, XM_REF_1V(g_XMSelect1110.v));
    return M;
#endif
}