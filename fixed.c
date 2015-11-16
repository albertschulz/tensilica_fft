#include "fixed.h"

/*
        fix_mpy() - fixed-point multiplication
*/
fixed fix_mpy(fixed a, fixed b)
{
    FIX_MPY(a,a,b);
    return a;
}
