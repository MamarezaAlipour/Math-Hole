// Copyright (c) 2022 Parisa Khaleghi
// All rights reserved

#ifndef MATHHOLE_CXX17_STAT_HXX
#define MATHHOLE_CXX17_STAT_HXX


// If <mathhole/cxx17/exec.hxx> has already been included,
// pull in the parallel/simd backend declarations
#if defined (MATHHOLE_CXX17_EXEC_HXX)
  #include <mathhole/detail/stat_exec_impl.hxx>
// Otherwise, only pull in serial backend
#else
  #include <mathhole/detail/stat_impl.hxx>
#endif


#endif