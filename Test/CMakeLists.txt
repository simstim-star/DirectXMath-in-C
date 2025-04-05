cmake_minimum_required(VERSION 3.8.2)

# List all files containing tests. (Change as needed)
set(TESTFILES Runner.c TestCommons.h TestMatrix.h TestCollision.h)
set(TEST_MAIN Tests)

add_executable(${TEST_MAIN} ${TESTFILES})

target_link_libraries(${TEST_MAIN} PRIVATE XMathC)

add_test(
    NAME ${PROJECT_NAME}.${TEST_MAIN}
    COMMAND ${TEST_MAIN})

add_custom_target(check
  make
  COMMAND ${CMAKE_CTEST_COMMAND} -V
  DEPENDS ${PROJECT_NAME})