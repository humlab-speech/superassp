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

#if defined(__GNUC__)
# define SMILE_CONSOLE_FORMAT(M, N) __attribute__ ((format (printf, M, N)))
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
