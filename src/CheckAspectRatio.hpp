
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

#ifndef CHECKASPECTRATIO_HPP
#define CHECKASPECTRATIO_HPP
#ifndef USE_LAIK
#define USE_LAIK
#endif
extern int CheckAspectRatio(double smallest_ratio, int x, int y, int z, const char *what, bool DoIo);
#endif // CHECKASPECTRATIO_HPP

