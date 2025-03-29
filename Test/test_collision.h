#include <stdlib.h>
#include <math.h>

#include "DirectXCCollision.h"
#include "TestCommons.h"
#include <DirectXCCollision.h>

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

#define EPSILON 0.0000009

inline bool spheres_equal(const dx_bounding_sphere *const s1, const dx_bounding_sphere *const s2)
{
    return ((fabs(s1->center.x - s2->center.x) < EPSILON)
        && (fabs(s1->center.y - s2->center.y) < EPSILON)
        && (fabs(s1->center.z - s2->center.z) < EPSILON)
        && (fabs(s1->r - s2->r) < EPSILON));
}

TEST_DECLARE(bounding_sphere_merged) {
    bool success = true;
    dx_bounding_sphere sht;

    const dx_bounding_sphere unit = { {0,0,0}, 1 };
    dx_bounding_sphere_merged(&sht, &unit, &unit);
    if (!spheres_equal(&sht, &unit))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - empty merge");
        printsh(unit);
        printsh(sht);
        success = false;
    }


    const dx_bounding_sphere _small = { {0.1f, 0.0f, 0.0f }, 0.5f };
    dx_bounding_sphere_merged(&sht, &unit, &_small);
    if (!spheres_equal(&unit, &sht))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - small test");
        printsh(unit);
        printsh(_small);
        printsh(sht);
        success = false;
    }

    const dx_bounding_sphere big = { {1.f, 2.f, 3.f}, 10.0f };
    dx_bounding_sphere_merged(&sht, &unit, &big);
    if (!spheres_equal(&sht, &big))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - big test");
        printsh(unit);
        printsh(big);
        printsh(sht);
        success = false;
    }

    const dx_bounding_sphere _far = { {10.f, -5.f, 4.f}, 2.0f };
    dx_bounding_sphere_merged(&sht, &_far, &big);

    const dx_bounding_sphere result = { {2.354666f, 0.946371f, 3.150518f}, 11.722761f };
    if (!spheres_equal(&sht, &result))
    {
        TEST_LOG_FAILED("bounding_sphere_merged - _far-big test ");
        printsh(_far);
        printsh(big);
        printsh(result);
        printsh(sht);
        success = false;
    }

    {
        const dx_bounding_sphere sph2 = { {2.f, 0, 0}, 1.0f };

        dx_bounding_sphere_merged(&sht, &unit, &sph2);

        const dx_bounding_sphere result2 = { {1.f, 0, 0}, 2.0f };
        if (!spheres_equal(&sht, &result2))
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
        const dx_bounding_sphere sph1 = { {0, 0, 0}, 5.f };
        const dx_bounding_sphere sph2 = { {2.f, 0, 0}, 1.f };

        dx_bounding_sphere_merged(&sht, &sph1, &sph2);

        if (!spheres_equal(&sht, &sph1))
        {
            TEST_LOG_FAILED("bounding_sphere_merged - create merge test2");
            printsh(sph1);
            printsh(sph2);
            printsh(sht);
            success = false;
        }
    }

    {
        const dx_bounding_sphere sph2 = { {2.f, 0, 0}, 5.f };

        dx_bounding_sphere_merged(&sht, &unit, &sph2);

        if (!spheres_equal(&sht, &sph2))
        {
            TEST_LOG_FAILED("bounding_sphere_merged - create merge test3");
            printsh(unit);
            printsh(sph2);
            printsh(sht);
            success = false;
        }
    }

    {
        const dx_bounding_sphere sph2 = { {0, 0, 0}, 2.f };

        dx_bounding_sphere_merged(&sht, &unit, &sph2);

        if (!spheres_equal(&sht, &sph2))
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

        dx_bounding_sphere sht;
        dx_bounding_sphere_from_points(&sht, 2, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const dx_bounding_sphere result = { {0.5f, 0, 0}, 0.5f };
        if (!spheres_equal(&sht, &result))
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - creating from points test1");
            printsh(sht);
            printsh(result);
            success = false;
        }
    }

    {
        const float points[] = { 0, 0, 0, 0.5f, 0, 1.0f };

        dx_bounding_sphere sht;
        dx_bounding_sphere_from_points(&sht, 2, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const dx_bounding_sphere result = { {0.25f, 0, 0.5f}, 0.559017f };
        if (!spheres_equal(&sht, &result))
        {
            TEST_LOG_FAILED("bounding_sphere_from_points - creating from points test2");
            printsh(sht);
            printsh(result);
            success = false;
        }
    }

    {
        const float points[] = { 0, 0, 0, 0.0f, 0.5f, 1.0f };

        dx_bounding_sphere sht;
        dx_bounding_sphere_from_points(&sht, 2, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const dx_bounding_sphere result = { {0, 0.25f, 0.5f}, 0.559017f };
        if (!spheres_equal(&sht, &result))
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

        dx_bounding_sphere sht;
        dx_bounding_sphere_from_points(&sht, 4, (const XMFLOAT3*)(points), sizeof(XMFLOAT3));

        const dx_bounding_sphere result = { { -0.0621097758f, 0.01741325f, -0.187598482f}, 1.16094756f };
        if (!spheres_equal(&sht, &result))
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
            XMStoreFloat3(&points[i], &rand_vec);
        }

        dx_bounding_sphere sht;
        dx_bounding_sphere_from_points(&sht, 32, points, sizeof(XMFLOAT3));

        // Expand a bit to ensure Contains works for all input points on all platforms
        sht.r += EPSILON;

        for (size_t i = 0; i < 32; ++i)
        {
            XMVECTOR p = XMLoadFloat3(&points[i]);
            dx_containment_type ct = dx_bounding_sphere_contains_point(&sht, &p);
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