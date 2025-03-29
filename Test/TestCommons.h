// Based on https://github.com/recp/cglm/blob/master/test/include/common.h

#pragma once

#include <stdio.h>

#define TEST(FUN) test_ ## FUN()
#define TEST_DECLARE(FUN) inline void test_ ## FUN(void)
#define TEST_SUCCESS(FUN) fprintf(stderr, GREEN  "  " OK_TEXT RESET " %-*s  \n", 2, # FUN)
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