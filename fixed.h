#ifndef FIXED_H
#define FIXED_H

// Definition of fixed - Data Type
#ifndef fixed
#define fixed short
#endif

/* FIX_MPY() - fixed-point multiplication macro.
   This macro is a statement, not an expression (uses asm).
   BEWARE: make sure _DX is not clobbered by evaluating (A) or DEST.
   args are all of type fixed.
   Scaling ensures that 32767*32767 = 32767. */
#define FIX_MPY(DEST,A,B)       DEST = ((long)(A) * (long)(B))>>15

// Multiply two fixed-point numbers
fixed fix_mpy(fixed a, fixed b);

#endif // FIXED_H
