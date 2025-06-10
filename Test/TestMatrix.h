#include <stdlib.h>
#include <math.h>

#include "DirectXMathC.h"
#include "TestCommons.h"

#define MATRIX_TEST_EPSILON 0.0000009
#define drand48()  ((float)(rand() / (RAND_MAX + 1.0)))

/*******************************************
* Tests forward declarations               *
********************************************/

TEST_DECLARE(mat3_is_nan);
TEST_DECLARE(mat3_is_inf);
TEST_DECLARE(mat3_is_identity);
TEST_DECLARE(mat3_mul);
TEST_DECLARE(alloc_4x4);
TEST_DECLARE(determinant);

/*******************************************
* Utilitary functions forward declarations *
********************************************/

XMFLOAT3X3 GenerateRandomMatrix();

// Write a floating point value into byte aligned memory
void WriteFloat(float fInput, char* pOutput);

// Read a floating point value from byte aligned memory
float ReadFloat(const char* pInput);

// Print a matrix
inline void printmatrixe(const XMMATRIX* xmmatrix);

enum COMPARISON CompareXMVECTOR(XMVECTOR a, XMVECTOR b, int NumElements);

float getdet(int size, const float* a);

/*******************************************
* Tests entry point                        *
********************************************/

inline void TESTS_Matrix() {
    TEST(mat3_is_nan);
    TEST(mat3_is_inf);
    TEST(mat3_is_identity);
    TEST(mat3_mul);
    TEST(mat3_mul);
    TEST(alloc_4x4);
    TEST(determinant);
}

/*******************************************
* Tests implementations                    *
********************************************/

TEST_DECLARE(mat3_is_nan) {
    XMFLOAT3X3 NaNMatrix = (XMFLOAT3X3){
        ._11 = NAN, ._12 = .0f, ._13 = .0f,
        ._21 = .0f, ._22 = 1.f, ._23 = .0f,
        ._31 = .0f, ._32 = .0f, ._33 = 1.f,
    };
    XMMATRIX NaNMatrix_XMMATRIX = XMLoadFloat3x3(&NaNMatrix);
    TEST_ASSERT(XM_MAT_IS_NAN(NaNMatrix_XMMATRIX));

    XMFLOAT3X3 I = XMFLOAT3X3_IDENTITY;
    XMMATRIX I_XMMATRIX = XMLoadFloat3x3(&I);
    TEST_ASSERT(!XM_MAT_IS_NAN(I_XMMATRIX));

    TEST_SUCCESS(mat3_is_nan);
}

TEST_DECLARE(mat3_is_inf) {
    XMFLOAT3X3 InfMatrix = (XMFLOAT3X3){
        ._11 = INFINITY, ._12 = .0f, ._13 = .0f,
        ._21 = .0f, ._22 = 1.f, ._23 = .0f,
        ._31 = .0f, ._32 = .0f, ._33 = 1.f,
    };
    XMMATRIX InfMatrix_XMMATRIX = XMLoadFloat3x3(&InfMatrix);
    TEST_ASSERT(XM_MAT_IS_INF(InfMatrix_XMMATRIX));

    XMFLOAT3X3 I = XMFLOAT3X3_IDENTITY;
    XMMATRIX I_XMMATRIX = XMLoadFloat3x3(&I);
    TEST_ASSERT(!XM_MAT_IS_INF(I_XMMATRIX));

    TEST_SUCCESS(mat3_is_inf);
}

TEST_DECLARE(mat3_is_identity) {
    XMFLOAT3X3 I = XMFLOAT3X3_IDENTITY;
    XMMATRIX I_XMMATRIX = XMLoadFloat3x3(&I);
    TEST_ASSERT(XM_MAT_IS_ID(I_XMMATRIX));

    XMFLOAT3X3 notI = (XMFLOAT3X3){
        ._11 = 1.1f, ._12 = .0f, ._13 = .0f,
        ._21 = .0f, ._22 = 1.f, ._23 = .0f,
        ._31 = .0f, ._32 = .0f, ._33 = 1.f,
    };
    XMMATRIX notI_XMMATRIX = XMLoadFloat3x3(&notI);
    TEST_ASSERT(!XM_MAT_IS_ID(notI_XMMATRIX));

    TEST_SUCCESS(mat3_is_identity);
}

TEST_DECLARE(mat3_mul) {
    for (int i = 0; i < 100000; ++i) {
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
        XMMATRIX C_XMMATRIX = XM_MAT_MULT(A_XMMATRIX, B_XMMATRIX);
        XMFLOAT4X4 result;
        XM_STORE_FLOAT4X4(&result, C_XMMATRIX);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                TEST_ASSERT(fabsf(result.m[i][j] - C.m[i][j]) < MATRIX_TEST_EPSILON)
            }
        }
    }
    TEST_SUCCESS(mat3_mul);
}

TEST_DECLARE(alloc_4x4) {
    TEST_ASSERT(sizeof(XMFLOAT4X4) == 4 * 16);
    const int csize = 256 + 65536;
    char *c = _aligned_malloc(csize, 65536);
    if (!c) { 
        TEST_LOG_FAILED("Alloc 4x4 - malloc failed");
        TEST_ASSERT(false);
    }
    const intptr_t pc64k = 65536;

    const int offset = 0;
    const int floatcount = 16;
    for (intptr_t j = pc64k - 16; j < pc64k + 16; j++) {
        intptr_t i;
        for (i = 0; i < csize; i++) {
            c[i] = (char)(~i & 0xff);
        }
        intptr_t first = offset + j;
        intptr_t last = offset + j + 4 * floatcount;

        for (i = 0; i < floatcount; i++) {
            WriteFloat((float)i, (char*)&c[first + (i * 4)]);
        }

        XMMATRIX m = XMLoadFloat4x4((const XMFLOAT4X4*)&c[offset + j]);
        for (i = 0; i < first; i++) {
            if (c[i] != (char)(~i & 0xff)) {
                TEST_LOG_FAILED("Alloc 4x4 - corrupted byte");
                fprintf(stderr, "%p corrupted byte %p: %x ... %x\n", (void*)(j), (void*)(i), c[i], (unsigned char)(~i));
                TEST_ASSERT(false);
            }
        }
        for (i = last; i < csize; i++) {
            if (c[i] != (char)(~i & 0xff)) {
                TEST_LOG_FAILED("Alloc 4x4 - corrupted byte");
                fprintf(stderr, "%p corrupted byte %p: %x ... %x\n", (void*)(j), (void*)(i), c[i], (unsigned char)(~i));
                TEST_ASSERT(false);
            }
        }
        for (i = 0; i < floatcount; i++) {
            float f = ReadFloat((char*)&(c[first + (i * 4)]));
            float check = (float)i;
            float mf = XM_VEC_IDX(m.r[i / 4], i % 4);
            if (f != check) {
                TEST_LOG_FAILED("Alloc 4x4 - corrupted source float");
                fprintf(stderr, "%p corrupted source float %p: %x ... %x\n", (void*)(j), (void*)(i), f, check);
                TEST_ASSERT(false);
            }
            if (mf != check) {
                TEST_LOG_FAILED("Alloc 4x4 - improperly read float");
                fprintf(stderr, "%p corrupted read float %p: %x ... %x\n", (void*)(j), (void*)(i), mf, check);
                TEST_ASSERT(false);
            }
        }
    }
    _aligned_free(c);
    TEST_SUCCESS(alloc_4x4);
}

TEST_DECLARE(determinant) {
    //XMMatrixDeterminant
    int i, j, k;
    for (k = 0; k < 10; k++) {
        float tmp[4][4];
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                tmp[i][j] = ((float)drand48()) / 100.f;
            }
        }
#ifdef _MSC_VER
#pragma warning( suppress : 6385 )
#endif
        XMMATRIX m = (XMMATRIX){
            .r = tmp[0][0],
        };
        enum COMPARISON c;
        XMVECTOR check = XMVectorReplicate(getdet(4, &tmp[0][0]));
        XMVECTOR r = XM_DETERMINANT(m);
        c = CompareXMVECTOR(r, check, 4);
        if (c > WITHIN4096) {
            fprintf(stderr, "%s: \n", "determinant");
            printmatrixe(&m);
            fprintf(stderr, "%f %f %f %f ... %f %f %f %f (%d)\n",
                XM_VECX(r), XM_VECY(r), XM_VECZ(r), XM_VECW(r),
                XM_VECX(check), XM_VECY(check), XM_VECZ(check), XM_VECW(check), c);
            TEST_ASSERT(false);
        }
    }
    TEST_SUCCESS(determinant);
}

/**************************************
* Utilitary functions implementations *
***************************************/

XMFLOAT3X3 GenerateRandomMatrix() {
    return (XMFLOAT3X3) {
        ._11 = drand48(), ._12 = drand48(), ._13 = drand48(),
        ._21 = drand48(), ._22 = drand48(), ._23 = drand48(),
        ._31 = drand48(), ._32 = drand48(), ._33 = drand48(),
    };
}

void WriteFloat(float fInput, char* pOutput)
{
    // Ensure the data is float aligned
    union {
        char m_cArray[4];
        float m_fTemp;
    } Temp;
    // Store the float to aligned memory
    Temp.m_fTemp = fInput;
    // Do a "memcpy" to unaligned memory
    // Note: This is an unavoidable Load/Hit/Store
    pOutput[0] = Temp.m_cArray[0];
    pOutput[1] = Temp.m_cArray[1];
    pOutput[2] = Temp.m_cArray[2];
    pOutput[3] = Temp.m_cArray[3];
}

float ReadFloat(const char* pInput)
{
    // Ensure the data is float aligned
    union {
        char m_cArray[4];
        float m_fTemp;
    } Temp;
    // Copy the data to an aligned buffer
    Temp.m_cArray[0] = pInput[0];
    Temp.m_cArray[1] = pInput[1];
    Temp.m_cArray[2] = pInput[2];
    Temp.m_cArray[3] = pInput[3];
    // Fetch the floating point number
    // from float aligned memory

    // Note: This is an unavoidable Load/Hit/Store
    return Temp.m_fTemp;
}


inline void printmatrixe(const XMMATRIX *xmmatrix)
{
    float m[4][4];
    memcpy_s(&m, sizeof(float) * 16, &xmmatrix, sizeof(XMMATRIX));
    fprintf(stdout, "%f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n",
        m[0][0], m[0][1], m[0][2], m[0][3],
        m[1][0], m[1][1], m[1][2], m[1][3],
        m[2][0], m[2][1], m[2][2], m[2][3],
        m[3][0], m[3][1], m[3][2], m[3][3]);
}

enum COMPARISON Compare(float a, float b)
{
    if (isnan(a) && isnan(b)) return EXACT;
    if (isnan(a) || isnan(b)) return WAYOFF;
    if (!isfinite(a) && !isfinite(b)) {
        if (copysign(1.0f, a) == copysign(1.0f, b)) return EXACT;
        else return WAYOFF;
    }
    if (!isfinite(a) || !isfinite(b)) return WAYOFF;

    if (a == b) return EXACT;
    float f = fabsf(b - a);
    if (f <= TESTEPSILON) return WITHINEPSILON;
    if (f <= TEST2EPSILON) return WITHIN2EPSILON;
    if (f <= TEST10EPSILON) return WITHIN10EPSILON;
    if (f <= TEST100EPSILON) return WITHIN100EPSILON;
    if (f <= fabsf(a / 4096)) return WITHIN4096;
    if (f <= TESTBIGEPSILON) return WITHINBIGEPSILON;
    if (f <= TESTBIGGEREPSILON) return WITHINBIGGEREPSILON;
    if (f <= TESTHUGEEPSILON) return WITHINHUGEEPSILON;
    if (f <= fabsf(a / 256)) return WITHIN1_256;
    if (f <= fabsf(a / 64)) return WITHIN1_64;
    if (f <= fabsf(a / 16)) return WITHIN1_16;
    if (f <= fabsf(a / 8)) return WITHIN1_8;
    if (f <= fabsf(a / 4)) return CLOSE;

    return WAYOFF;
}

enum COMPARISON CompareXMVECTOR(XMVECTOR a, XMVECTOR b, int NumElements) {

    assert((NumElements <= 4) && (NumElements >= 0));
    int i;
    enum COMPARISON ret = EXACT;
    enum COMPARISON r;

    for (i = 0; i < NumElements; i++) {
        r = Compare(XM_VEC_IDX(a, i), XM_VEC_IDX(b, i));
        ret = (r > ret) ? r : ret;
    }

    return ret;
}


float getdet(int size, const float* a)
{
    float t[16];
    int it = 0;
    float res = 0;
    float sign = 1;
    int i, x, y;
    if (size == 1) return *a;
    for (i = 0; i < size; i++) {
        for (x = 0; x < size; x++) {
            if (x == i) x++;
            if (x == size) break;
            for (y = 1; y < size; y++) {
                if (it > 15) break;
                t[it++] = a[y * size + x];
            }
        }
        //float gd = getdet(size-1,t);
        //printi ("size: %d, i: %d, gd: %f\n",size,i,gd);
        res += sign * a[i] * getdet(size - 1, t);
        sign = -sign;
        it = 0;
    }
    return res;
}