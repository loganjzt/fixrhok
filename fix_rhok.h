/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing authors:
		Zhitong Jiang
   Contact:

------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(rhok,FixRhok)

#else

#ifndef FIX_RHOK_H
#define FIX_RHOK_H

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#include <complex>
#include "fix.h"

namespace LAMMPS_NS {

class FixRhok : public Fix {
 public:
  FixRhok(class LAMMPS *, int, char **);
  ~FixRhok();

  int setmask();
  void init();
  void post_force(int);
  double memory_usage();
  void end_of_step();
  void compute_vector();

 private:
  int nloop;				// number of wave vector that controled
  int nevery;

  int *alpha;
  double *kappa_re;
  double *kappa_im;

  double *rhokStar_re;
  double *rhokStar_im;

  int *wavenumber;

  double *globalRhok;		// current rhok value
  double *localRhok;	
  double lbox[3];

};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix phonon command...

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No atom found for fix phonon!

Self-explanatory. Number of atoms in the group that was passed to
fix-phonon is less than 1.

E: Can not open output file %s"

Self-explanatory.

E: Illegal fix_modify command

Self-explanatory.

E: Could not find fix_modify temp ID

Self-explanatory.

E: Fix_modify temp ID does not compute temperature

Self-explanatory.

*/
