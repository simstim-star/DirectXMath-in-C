#include <stdlib.h>
#include <math.h>

#include "DirectXMathC.h"
#include "TestCommons.h"

TEST_DECLARE(mat3_is_identity);
TEST_DECLARE(mat3_mul);

inline void TESTS_Matrix() {
    TEST(mat3_is_identity);
    TEST(mat3_mul);
}

#define EPSILON 0.0000009
#define drand48()  ((float)(rand() / (RAND_MAX + 1.0)))

XMFLOAT3X3 GenerateRandomMatrix();

TEST_DECLARE(mat3_is_identity) {
    XMFLOAT3X3 I = XMFLOAT3X3_IDENTITY;
    XMMATRIX I_XMMATRIX = XMLoadFloat3x3(&I);
    TEST_ASSERT(XMMatrixIsIdentity(&I_XMMATRIX))

    XMFLOAT3X3 notI = (XMFLOAT3X3){
        ._11 = 1.1f, ._12 = .0f, ._13 = .0f,
        ._21 = .0f, ._22 = 1.f, ._23 = .0f,
        ._31 = .0f, ._32 = .0f, ._33 = 1.f,
    };
    XMMATRIX notI_XMMATRIX = XMLoadFloat3x3(&notI);
    TEST_ASSERT(!XMMatrixIsIdentity(&notI_XMMATRIX))

    TEST_SUCCESS(mat3_is_identity);
}

TEST_DECLARE(mat3_mul) {
    XMFLOAT3X3 A = GenerateRandomMatrix();
    XMFLOAT3X3 B = GenerateRandomMatrix();

    XMFLOAT4X4 C = XMFLOAT4X4_ZERO;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; k++) {
                C.m[i][j] += A.m[i][k] * B.m[k][j];
            }
        }
    }
    XMMATRIX A_XMMATRIX = XMLoadFloat3x3(&A);
    XMMATRIX B_XMMATRIX = XMLoadFloat3x3(&B);
    XMMATRIX C_XMMATRIX = XMMatrixMultiply(&A_XMMATRIX, &B_XMMATRIX);
    XMFLOAT4X4 result;
    XMStoreFloat4x4(&result, &C_XMMATRIX);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            TEST_ASSERT(fabsf(result.m[i][j] - C.m[i][j]) < EPSILON)
        }
    }
    TEST_SUCCESS(mat3_mul);
}

XMFLOAT3X3 GenerateRandomMatrix() {
    return (XMFLOAT3X3) {
        ._11 = drand48(), ._12 = drand48(), ._13 = drand48(),
        ._21 = drand48(), ._22 = drand48(), ._23 = drand48(),
        ._31 = drand48(), ._32 = drand48(), ._33 = drand48(),
    };
}