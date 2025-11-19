// Based on https://github.com/recp/cglm/blob/master/test/include/common.h

#pragma once

#include <stdio.h>

#define TEST(FUN) test_ ## FUN()

#ifdef _MSC_VER
#  define TEST_DECLARE(FUN) inline void test_ ## FUN(void)
#else
#  ifdef __GNUC__
#    define TEST_DECLARE(FUN) static inline void test_ ## FUN(void)
#  endif
#endif

#define TEST_SUCCESS(FUN) fprintf(stdout, GREEN  "  " OK_TEXT RESET " %-*s  \n", 2, # FUN)
#define TEST_LOG_FAILED(TEST_NAME) fprintf(stderr, RED "%s" RESET ": unit test failed\n", TEST_NAME)

#define TEST_ASSERT(expr)                          \
  if (!(expr)) {                                   \
    fprintf(stderr,                                \
            RED "  assert fail" RESET              \
            " in " BOLDCYAN "%s " RESET            \
            "on " BOLDMAGENTA "line %d" RESET      \
            " : " BOLDWHITE " ASSERT(%s)\n" RESET, \
            __FILE__,                              \
            __LINE__,                              \
            #expr);                                \
    return;                                        \
  }

# define OK_TEXT    "ok:"
# define FAIL_TEXT  "fail:"

#define RESET       "\033[0m"
#define RED         "\033[31m"             /* Red */
#define GREEN       "\033[32m"             /* Green */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

//comparisons that are equal return 0 (EXACT)
//comparisons that are within TESTEPSILON return 1 (WITHINEPSILON)
//comparisons that are within TEST2EPSILON return 2 (WITHIN2EPSILON)
//...
//comparisons that are within 100%/4096 return 5 (WITHIN4096)
//...
//comparisons that are within 25% return 10 (CLOSE)
//comparisons that are beyond that, or involve a mismatched NAN or INF return WAYOFF.
enum COMPARISON {
    EXACT, WITHINEPSILON, WITHIN2EPSILON,
    WITHIN10EPSILON, WITHIN100EPSILON, WITHIN4096,
    WITHINBIGEPSILON, WITHINBIGGEREPSILON, WITHINHUGEEPSILON, WITHIN1_256,
    WITHIN1_64, WITHIN1_16, WITHIN1_8, CLOSE, WAYOFF
};

#define TESTEPSILON  1.192092896e-7f
#define TEST10EPSILON  1.192092896e-6f
#define TEST100EPSILON  1.192092896e-6f
#define TEST2EPSILON .00000023841859f
#define TESTBIGEPSILON .001f
#define TESTBIGGEREPSILON .0025f
#define TESTHUGEEPSILON .01f