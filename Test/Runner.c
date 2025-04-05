#include <stdio.h>

#include "TestMatrix.h"
#include "TestCollision.h"

void PrintCallConvention() {
#if defined(_XM_VECTORCALL_)
    fprintf(stdout, "[[Call Convention]] :: __vectorcall \n");
#elif defined(__GNUC__)
    fprintf(stdout, "[[Call Convention]] :: none (GNUC) \n");
#else
    fprintf(stdout, "[[Call Convention]] :: __fastcall \n");
#endif
}

int main()
{
    PrintCallConvention();
    TESTS_Matrix();
    TESTS_Collision();
}
