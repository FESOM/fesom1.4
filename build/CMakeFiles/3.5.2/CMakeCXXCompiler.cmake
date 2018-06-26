set(CMAKE_CXX_COMPILER "/opt/mpi/bullxmpi_mlx/1.2.8.3/bin/mpicxx")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Intel")
set(CMAKE_CXX_COMPILER_VERSION "17.0.2.20170213")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "98")
set(CMAKE_CXX_COMPILE_FEATURES "")
set(CMAKE_CXX98_COMPILE_FEATURES "")
set(CMAKE_CXX11_COMPILE_FEATURES "")
set(CMAKE_CXX14_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_SIMULATE_VERSION "")

set(CMAKE_AR "/sw/rhel6-x64/gcc/binutils-2.24-gccsys/bin/ar")
set(CMAKE_RANLIB "/sw/rhel6-x64/gcc/binutils-2.24-gccsys/bin/ranlib")
set(CMAKE_LINKER "/sw/rhel6-x64/gcc/binutils-2.24-gccsys/bin/ld")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;mm;CPP)
set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()




set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "mpi_cxx;mpi;dl;imf;m;numa;rt;nsl;util;imf;m;dl;imf;svml;irng;stdc++;m;ipgo;decimal;cilkrts;stdc++;irc;pthread;svml;c;irc_s;dl;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/opt/mpi/bullxmpi_mlx/1.2.8.3/lib;/sw/rhel6-x64/intel/intel-17.0.2/mkl/lib/intel64;/sw/rhel6-x64/intel/intel-17.0.2/compilers_and_libraries_2017.2.174/linux/compiler/lib/intel64_lin;/sw/rhel6-x64/gcc/gcc-4.8.2/lib/gcc/x86_64-unknown-linux-gnu/4.8.2;/sw/rhel6-x64/gcc/gcc-4.8.2/lib64;/lib64;/usr/lib64;/sw/rhel6-x64/gcc/gcc-4.8.2/lib;/lib;/usr/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
