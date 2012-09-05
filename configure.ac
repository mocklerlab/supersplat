#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Supersplat], [2.0], [DBryant@ddpsc.org])
AM_INIT_AUTOMAKE([Supersplat], [2.0])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries. 
AC_ARG_WITH([seqan-include-path],
  [AS_HELP_STRING([--with-seqan-include-path],
    [location of SeqAn, defaults to /usr/include/seqan])],
  [SEQAN_CFLAGS="-I${withval}/core/include -I${withval}/extras/include"],
  [SEQAN_CFLAGS='-I/usr/include/seqan/core/include -I/usr/include/seqan/extras/include'])
AC_SUBST([SEQAN_CFLAGS])

# Checks for header files.
AC_CHECK_HEADERS([float.h stddef.h stdint.h stdlib.h string.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove pow sqrt strchr strstr strtol strtoul])

AC_OUTPUT([Makefile])
