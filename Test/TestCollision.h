#include <stdlib.h>
#include <math.h>

#include "DirectXCollisionC.h"
#include "TestCommons.h"

#define printsh(v) fprintf(stderr, "%s: center=%f,%f,%f  radius=%f\n", #v, v.center.x, v.center.y, v.center.z, v.r)

#define printct(t) switch(t) { case DISJOINT: fprintf(stderr, "%s: DISJOINT", #t); break; \
                               case INTERSECTS: fprintf(stderr, "%s: INTERSECTS", #t); break; \
                               case CONTAINS: fprintf(stderr, "%s: CONTAINS", #t); break; } fprintf(stderr, "\n")


#define XM_RAND_MAX (0x7fff)
#define XM_RAND() (rand()%(XM_RAND_MAX+1))

XMVECTOR rand_vec16(void);

TEST_DECLARE(bounding_sphere_merged);
TEST_DECLARE(bounding_sphere_from_points);

inline void TESTS_Collision() {
    TEST(bounding_sphere_merged);
    TEST(bounding_sphere_from_points);
}

#define COLLISION_TEST_EPSILON 0.00001

inline bool SpheresEqual(const DirectX_BoundingSphere *const s1, const DirectX_BoundingSphere *const s2)
{
    return ((fabs(s1->center.x - s2->center.x) < COLLISION_TEST_EPSILON)
        && (fabs(s1->center.y - s2->center.y) < COLLISION_TEST_EPSILON)
        && (fabs(s1->center.z - s2->center.z) < COLLISION_TEST_EPSILON)
        && (fabs(s1->r - s2->r) < COLLISION_TEST_EPSILON));
}

TEST_DECLARE(bounding_sphere_merged) {
    bool success = true;
    DirectX_BoundingSphere sht;

    const DirectX_BoundingSphere unit = { {0,0,0}, 1 };
    DirectX_BoundingSphere_Merged(&sht, &unit, &unit);
    if (!SpheresEqual(&sht, &unit))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - empty merge");
        printsh(unit);
        printsh(sht);
        success = false;
    }


    const DirectX_BoundingSphere _small = { {0.1f, 0.0f, 0.0f }, 0.5f };
    DirectX_BoundingSphere_Merged(&sht, &unit, &_small);
    if (!SpheresEqual(&unit, &sht))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - small test");
        printsh(unit);
        printsh(_small);
        printsh(sht);
        success = false;
    }

    const DirectX_BoundingSphere big = { {1.f, 2.f, 3.f}, 10.0f };
    DirectX_BoundingSphere_Merged(&sht, &unit, &big);
    if (!SpheresEqual(&sht, &big))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - big test");
        printsh(unit);
        printsh(big);
        printsh(sht);
        success = false;
    }

    const DirectX_BoundingSphere _far = { {10.f, -5.f, 4.f}, 2.0f };
    DirectX_BoundingSphere_Merged(&sht, &_far, &big);

    const DirectX_BoundingSphere result = { {2.354666f, 0.946371f, 3.150518f}, 11.722761f };
    if (!SpheresEqual(&sht, &result))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - _far-big test ");
        printsh(_far);
        printsh(big);
        printsh(result);
        printsh(sht);
        success = false;
    }

    {
        const DirectX_BoundingSphere sph2 = { {2.f, 0, 0}, 1.0f };

        DirectX_BoundingSphere_Merged(&sht, &unit, &sph2);

        const DirectX_BoundingSphere result2 = { {1.f, 0, 0}, 2.0f };
        if (!SpheresEqual(&sht, &result2))
        {
            TEST_LOG_FAILED("bounding_sphere_merged - create merge test1");
            printsh(unit);
            printsh(sph2);
            printsh(sht);
            printsh(result2);
            success = false;
        }
    }

    {
        const DirectX_BoundingSphere sph1 = { {0, 0, 0}, 5.f };
        const DirectX_BoundingSphere sph2 = { {2.f, 0, 0}, 1.f };

        DirectX_BoundingSphere_Merged(&sht, &sph1, &sph2);

        if (!SpheresEqual(&sht, &sph1))
        {
            TEST_LOG_FAILED("bounding_sphere_merged - create merge test2");
            printsh(sph1);
            printsh(sph2);
            printsh(sht);
            success = false;
        }
    }

    {
        const DirectX_BoundingSphere sph2 = { {2.f, 0, 0}, 5.f };

        DirectX_BoundingSphere_Merged(&sht, &unit, &sph2);

        if (!SpheresEqual(&sht, &sph2))
        {
            TEST_LOG_FAILED("bounding_sphere_merged - create merge test3");
            printsh(unit);
            printsh(sph2);
            printsh(sht);
            success = false;
        }
    }

    {
        const DirectX_BoundingSphere sph2 = { {0, 0, 0}, 2.f };

        DirectX_BoundingSphere_Merged(&sht, &unit, &sph2);

        if (!SpheresEqual(&sht, &sph2))
        {
            TEST_LOG_FAILED("bounding_sphere_merged - create merge test4");
            printsh(unit);
            printsh(sph2);
            printsh(sht);
            success = false;
        }
    }

    if (success) TEST_SUCCESS(bounding_sphere_merged);
    else         TEST_ASSERT(false);
}


TEST_DECLARE(bounding_sphere_from_points) {
    bool success = true;

    {
        const float points[] = { 0, 0, 0, 1.f, 0, 0 };

        DirectX_BoundingSphere sht;
        DirectX_BoundingSphere_FromPoints(&sht, 2, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const DirectX_BoundingSphere result = { {0.5f, 0, 0}, 0.5f };
        if (!SpheresEqual(&sht, &result))
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - creating from points test1");
            printsh(sht);
            printsh(result);
            success = false;
        }
    }

    {
        const float points[] = { 0, 0, 0, 0.5f, 0, 1.0f };

        DirectX_BoundingSphere sht;
        DirectX_BoundingSphere_FromPoints(&sht, 2, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const DirectX_BoundingSphere result = { {0.25f, 0, 0.5f}, 0.559017f };
        if (!SpheresEqual(&sht, &result))
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - creating from points test2");
            printsh(sht);
            printsh(result);
            success = false;
        }
    }

    {
        const float points[] = { 0, 0, 0, 0.0f, 0.5f, 1.0f };

        DirectX_BoundingSphere sht;
        DirectX_BoundingSphere_FromPoints(&sht, 2, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const DirectX_BoundingSphere result = { {0, 0.25f, 0.5f}, 0.559017f };
        if (!SpheresEqual(&sht, &result))
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - creating from points test3");
            printsh(sht);
            printsh(result);
            success = false;
        }
    }

    {
        float irt2 = 1.f / sqrtf(2.f);
        float irt3 = 1.f / sqrtf(3.f);
        const float points[] = { 0, 1.f, 0,
                                 0, -irt2, -irt2,
                                 -irt3, -irt3, irt3,
                                 irt3, -irt3, irt3 };

        DirectX_BoundingSphere sht;
        DirectX_BoundingSphere_FromPoints(&sht, 4, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const DirectX_BoundingSphere result = { { -0.0621097758f, 0.01741325f, -0.187598482f}, 1.16094756f };
        if (!SpheresEqual(&sht, &result))
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - creating from points test4");
            printsh(sht);
            printsh(result);
            success = false;
        }

        /* TODO: implement when we have bounding box
        // See that the bound above is tighter than a sphere created from a bounding box of the same points
        dx_bounding_sphere box;
        dx_bounding_sphere_from_points(&box, 4, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        dx_bounding_sphere shv;
        BoundingSphere::CreateFromBoundingBox(shv, box);

        if (sht.r > shv.r)
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - expected tighter bounds");
            printsh(sht);
            printsh(shv);
            success = false;
        }
        */
    }

    {
        XMFLOAT3 points[32];

        for (size_t i = 0; i < 32; ++i)
        {
            XMVECTOR rand_vec = rand_vec16();
            XMStoreFloat3(&points[i], XM_REF_1V(rand_vec));
        }

        DirectX_BoundingSphere sht;
        DirectX_BoundingSphere_FromPoints(&sht, 32, points, sizeof(XMFLOAT3));

        // Expand a bit to ensure Contains works for all input points on all platforms
        sht.r += COLLISION_TEST_EPSILON;

        for (size_t i = 0; i < 32; ++i)
        {
            XMVECTOR p = XMLoadFloat3(&points[i]);
            DirectX_ContainmentType ct = DirectX_BoundingSphere_ContainsPoint(&sht, XM_REF_1V(p));
            if (ct != CONTAINS)
            {
                TEST_LOG_FAILED("bounding_sphere_from_points - Sphere-Point verification test");
                printsh(sht);
                printct(ct);
                success = false;
            }
        }
    }

    if (success) TEST_SUCCESS(bounding_sphere_from_points);
    else         TEST_ASSERT(false);
}



// UTILS

#define XM_RAND_MAX (0x7fff)
#define XM_RAND() (rand()%(XM_RAND_MAX+1))

XMVECTOR rand_vec16(void)
{
    // Use scalar math, not vector math, to
    // generate the values. This will prevent
    // accidental false positives if the vector
    // operations being used are the ones that
    // are failing.

    // Grab the raw integer value and convert to a float
    float fx = (float)(XM_RAND());
    float fy = (float)(XM_RAND());
    float fz = (float)(XM_RAND());
    float fw = (float)(XM_RAND());
    // The float is 0 to XM_RAND_MAX inclusive
    fx = fx / ((float)(XM_RAND_MAX) / 32.767f);   // Normalize to 0-32.767f
    fy = fy / ((float)(XM_RAND_MAX) / 32.767f);
    fz = fz / ((float)(XM_RAND_MAX) / 32.767f);
    fw = fw / ((float)(XM_RAND_MAX) / 32.767f);
    fx = fx - 16.0f;          // Convert 0-32.767f to -16f - 16.767f
    fy = fy - 16.0f;
    fz = fz - 16.0f;
    fw = fw - 16.0f;
    // This data type ensures there is no code being generated.
    XMVECTORF32 vResult = { { { fx, fy, fz, fw } } };
    // Return the vector value
    return vResult.v;
}