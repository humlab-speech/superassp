/*  smileConsole.h -- console output shim for the superassp build of openSMILE.
 *
 *  CRAN requires that compiled code not write to stdout/stderr itself: output
 *  has to reach R's console.  openSMILE cannot simply call Rprintf, because the
 *  static library is linked both by the R package and by the SMILExtract
 *  executable, which is a standalone process without an R runtime.
 *
 *  The library therefore formats messages itself and hands the finished text to
 *  writers installed by the host:
 *
 *    - the R package installs writers backed by Rprintf/REvprintf
 *      (src/opensmile_wrapper.cpp),
 *    - SMILExtract installs writers backed by fwrite(stdout/stderr)
 *      (progsrc/smilextract/SMILExtract.cpp).
 *
 *  With no writer installed output is discarded, which keeps the library free of
 *  any host dependency.
 */

#ifndef SMILE_CONSOLE_H_
#define SMILE_CONSOLE_H_

#include <stddef.h>

/*  Format-string checking, mirroring R's own R_PRINTF_FORMAT (R_ext/Print.h).
 *  GCC on Windows validates against MSVCRT's printf by default, which would
 *  turn every C99 conversion (e.g. %zu) into a spurious warning even though the
 *  shim formats with the C library's own vsnprintf; "gnu_printf" describes what
 *  that call accepts.  Clang needs no such distinction, and on pre-UCRT msvcrt
 *  there is nothing to validate against, so the attribute is dropped there.  */
#if defined(__GNUC__)
# ifdef _WIN32
#  if defined(_UCRT) || ((__MSVCRT_VERSION__ >= 0x1400) || \
                        (__MSVCRT_VERSION__ >= 0xE00 && __MSVCRT_VERSION__ < 0x1000))
#   if defined(__clang__)
#    define SMILE_CONSOLE_FORMAT(M, N) __attribute__ ((format (printf, M, N)))
#   else
#    define SMILE_CONSOLE_FORMAT(M, N) __attribute__ ((format (gnu_printf, M, N)))
#   endif
#  else
#   define SMILE_CONSOLE_FORMAT(M, N)
#  endif
# else
#  define SMILE_CONSOLE_FORMAT(M, N) __attribute__ ((format (printf, M, N)))
# endif
#else
# define SMILE_CONSOLE_FORMAT(M, N)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*  A writer receives NUL-free text of the given length; hosts are responsible
 *  for flushing.  Passing NULL discards that stream.  */
typedef void (*smile_console_writer)(const char *text, size_t len);

void smile_console_set_writer(smile_console_writer out, smile_console_writer err);

/*  Whether the error stream is a terminal (drives colour escapes in the
 *  logger).  Hosts that do not know should leave the default of 0.  */
void smile_console_set_tty(int err_is_tty);
int  smile_console_is_tty(void);

int smile_console_printf(const char *fmt, ...) SMILE_CONSOLE_FORMAT(1, 2);
int smile_console_error(const char *fmt, ...) SMILE_CONSOLE_FORMAT(1, 2);
int smile_console_puts(const char *s);      /*  appends a newline, like puts()  */
int smile_console_putchar(int c);

#ifdef __cplusplus
}  /*  extern "C"  */

#include <ostream>

/*  std::cout / std::cerr replacements; stream formatting is preserved.  */
std::ostream &smile_console_out();
std::ostream &smile_console_err();

#endif

#endif /* SMILE_CONSOLE_H_ */
