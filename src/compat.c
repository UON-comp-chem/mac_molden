#include <strings.h>
#include <string.h>
#include <glob.h>
#include <unistd.h>
#include <fnmatch.h>

#if defined(HAVE_SYS_SELECT_H)
#include <sys/select.h>
#endif
#undef glob

extern "C" int glob_old(const char * pattern, int flags, int (*errfunc) (const char *epath, int eerrno), glob_t *pglob);

#ifdef __i386__
    __asm__(".symver glob_old,glob@GLIBC_2.1");
#elif defined(__amd64__)
    __asm__(".symver glob_old,glob@GLIBC_2.2.5");
#elif defined(__arm__)
    __asm(".symver glob_old,glob@GLIBC_2.4");
#elif defined(__aarch64__)
    __asm__(".symver glob_old,glob@GLIBC_2.17");
#endif

extern "C" int __wrap_glob(const char * pattern, int flags, int (*errfunc) (const char *epath, int eerrno), glob_t *pglob)
{
    return glob_old(pattern, flags, errfunc, pglob);
}


#undef explicit_bzero
void explicit_bzero (void *s, size_t len)
{
    memset (s, '\0', len);
    /* Compiler barrier.  */
    asm volatile ("" ::: "memory");
}

void __explicit_bzero_chk (void *dst, size_t len, size_t dstlen)
{
    if (__glibc_unlikely (dstlen < len))
        __chk_fail ();
    memset (dst, '\0', len);
    asm volatile ("" ::: "memory");
}
#define strong_alias (__explicit_bzero_chk, __explicit_bzero_chk_internal)


