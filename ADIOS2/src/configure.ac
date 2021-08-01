
AC_PREREQ([2.63])
AC_INIT([TCLB_ADIOS], 0.01)
#AC_CONFIG_HEADERS([config.h:./src/config.h.in])
AC_CONFIG_FILES([config.mk:src/config.mk.in])

AC_LANG([C++])

AC_ARG_WITH([adios_dir],
	AC_HELP_STRING([--with-adios_dir=adios_dir],
		[specify the directory containing adios_dir /usr/local/adios2]),
	[ADIOS_DIR="$withval"])



test -z "$CXX" && CXX=mpic++
test -z "$CC" && CC=mpicc

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_HEADERS([float.h malloc.h stddef.h stdint.h stdlib.h string.h wchar.h],[],[AC_MSG_ERROR([Cannot find standart headers])])


AC_ARG_WITH([mpi-include],
	AC_HELP_STRING([--with-mpi-include=mpi_include],
		[specify the full path to your mpi headers, e.g. /usr/include/mpi]),
	[MPI_INCLUDE="$withval"])

AC_ARG_WITH([mpi-lib],
	AC_HELP_STRING([--with-mpi-lib=mpi_lib],
		[specify the full path to your mpi shared libraries, e.g. /usr/lib/openmpi/]),
	[MPI_LIB="$withval"])



AC_MSG_CHECKING([MPI include path])
if test -z "${MPI_INCLUDE}"; then
	if test -z "${MPI}"; then
		MPI_INCLUDE="/usr/include/mpi"
	else
		MPI_INCLUDE="${MPI}/include"
	fi
	AC_MSG_RESULT([${MPI_INCLUDE}])
else
	AC_MSG_RESULT([${MPI_INCLUDE}])
fi
CPPFLAGS="${CPPFLAGS} -I${MPI_INCLUDE}"

AC_MSG_CHECKING([MPI library path])
if test -z "${MPI_LIB}"; then
	if test -z "${MPI}"; then
		MPI_LIB="/usr/lib/openmpi/lib"
	else
		MPI_LIB="${MPI}/lib"
	fi
	AC_MSG_RESULT([${MPI_LIB}])
else
	AC_MSG_RESULT([${MPI_LIB}])
fi
LDFLAGS="${LDFLAGS} -L${MPI_LIB}"
AC_CHECK_HEADERS([mpi.h],[],[AC_MSG_ERROR([Didn't find mpi.h])])
AC_CHECK_LIB([mpi], [MPI_Recv],[],[
	AC_CHECK_FUNC([MPI_Recv],[],[
		AC_MSG_ERROR([Didn't find MPI Library])
	])
])




adios2_config="${ADIOS_DIR}/bin/adios2-config"
AC_CHECK_PROG(adios2_config,${adios2_config},${adios2_config},adios2-config)

AC_MSG_NOTICE([ADIOS2 config: ${adios2_config}])

AC_MSG_CHECKING([ADIOS2 embandinng: linker flags])
LD=$(${adios2_config} --cxx-libs)
LDFLAGS="${LDFLAGS} ${LD}"
AC_MSG_RESULT(["${LD}"])


AC_MSG_CHECKING([ADIOS2 embandinng: linker flags])
LD=$(${adios2_config} --cxx-flags)

CPPFLAGS="${CPPFLAGS} ${LD}"
AC_MSG_RESULT(["${LD}"])


AC_OUTPUT