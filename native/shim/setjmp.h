/* A fatal C error stops the current WebAssembly call. JavaScript discards
   the instance on a trap. The numeric computation uses upstream code. */
#ifndef TRACEPICK_SETJMP_H
#define TRACEPICK_SETJMP_H
typedef int jmp_buf[1];
#define setjmp(buf) (0)
#define longjmp(buf, value) __builtin_trap()
#endif
