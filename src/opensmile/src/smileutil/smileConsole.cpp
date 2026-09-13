/*  smileConsole.cpp -- see smileConsole.h for the rationale.  */

#include <smileutil/smileConsole.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*  MSVC's <stdarg.h> has no va_copy; the shim only needs a forward copy.  */
#ifndef va_copy
#define va_copy(dst, src) ((dst) = (src))
#endif

#ifdef __cplusplus
#include <ostream>
#include <streambuf>
#endif

static smile_console_writer console_out_writer = NULL;
static smile_console_writer console_err_writer = NULL;
static int console_err_is_tty = 0;

void smile_console_set_writer(smile_console_writer out, smile_console_writer err)
{
  console_out_writer = out;
  console_err_writer = err;
}

void smile_console_set_tty(int err_is_tty)
{
  console_err_is_tty = err_is_tty;
}

int smile_console_is_tty(void)
{
  return console_err_is_tty;
}

static void console_emit_text(smile_console_writer writer, const char *text, size_t len)
{
  if (writer != NULL && text != NULL && len > 0) {
    writer(text, len);
  }
}

/*  Format into a stack buffer, falling back to the heap for long messages, and
 *  hand the result to the writer.  The library formats rather than the host, so
 *  hosts never have to deal with format strings or va_list.  */
static void console_emit(smile_console_writer writer, const char *fmt, va_list ap)
{
  char stackbuf[1024];
  char *heap;
  va_list ap_size;
  int n;

  if (writer == NULL || fmt == NULL) {
    return;
  }

  va_copy(ap_size, ap);
  n = vsnprintf(NULL, 0, fmt, ap_size);
  va_end(ap_size);
  if (n < 0) {
    return;
  }

  if ((size_t)n < sizeof(stackbuf)) {
    vsnprintf(stackbuf, sizeof(stackbuf), fmt, ap);
    writer(stackbuf, (size_t)n);
    return;
  }

  heap = (char *)malloc((size_t)n + 1);
  if (heap == NULL) {
    return;
  }
  vsnprintf(heap, (size_t)n + 1, fmt, ap);
  writer(heap, (size_t)n);
  free(heap);
}

int smile_console_printf(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  console_emit(console_out_writer, fmt, ap);
  va_end(ap);
  return 0;
}

int smile_console_error(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  console_emit(console_err_writer, fmt, ap);
  va_end(ap);
  return 0;
}

int smile_console_puts(const char *s)
{
  if (s == NULL) {
    return 0;
  }
  console_emit_text(console_out_writer, s, strlen(s));
  console_emit_text(console_out_writer, "\n", 1);
  return 0;
}

int smile_console_putchar(int c)
{
  const char ch = (char)c;
  console_emit_text(console_out_writer, &ch, 1);
  return c;
}

#ifdef __cplusplus

namespace {

/*  Forwards stream writes to the installed writer, so newmat and any other
 *  C++ code that prints through std::cout keeps its formatting and its
 *  semantics.  */
class SmileConsoleBuf : public std::streambuf {
 public:
  explicit SmileConsoleBuf(int is_err) : is_err_(is_err) {}

 protected:
  virtual int overflow(int c)
  {
    if (c != EOF) {
      const char ch = (char)c;
      emit(&ch, 1);
    }
    return c;
  }

  virtual std::streamsize xsputn(const char *s, std::streamsize n)
  {
    emit(s, (size_t)n);
    return n;
  }

 private:
  void emit(const char *s, size_t n)
  {
    console_emit_text(is_err_ ? console_err_writer : console_out_writer, s, n);
  }

  int is_err_;
};

}  /*  namespace  */

std::ostream &smile_console_out()
{
  static SmileConsoleBuf buf(0);
  static std::ostream out(&buf);
  return out;
}

std::ostream &smile_console_err()
{
  static SmileConsoleBuf buf(1);
  static std::ostream err(&buf);
  return err;
}

#endif /* __cplusplus */
