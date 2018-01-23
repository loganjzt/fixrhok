/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Zhitong Jiang

   Contact:

------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_rhok.h"
#include "atom.h"
#include "atom_masks.h"
#include "compute.h"
#include "domain.h"
#include "force.h"
#include "region.h"
#include "variable.h"
#include "group.h"
#include "lattice.h"
#include "modify.h"
#include "update.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"

#include "math.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRhok::FixRhok(LAMMPS *lmp,  int narg, char **arg) : Fix(lmp, narg, arg)
{
  virial_flag = 1;

  //if (narg < 5 ) error->all(FLERR,"Illegal fix vector");
  if( (narg - 4) % 6 != 0) error->all(FLERR,"Illegal input");
  nloop = (narg-4)/6;

  alpha = new int[nloop];

  kappa_re = new double[nloop];
  kappa_im = new double[nloop];

  rhokStar_re = new double[nloop];
  rhokStar_im = new double[nloop];
  wavenumber = new int[nloop];

  localRhok = new double[2*nloop]; 
  globalRhok = new double[2*nloop]; 

  for(int i =0;i<nloop;i++){
    if(strcmp(arg[ i*6 + 3 ],"x") == 0 ) alpha[i] = 0;
    if(strcmp(arg[ i*6 + 3 ],"y") == 0 ) alpha[i] = 1;
    if(strcmp(arg[ i*6 + 3 ],"z") == 0 ) alpha[i] = 2;

    kappa_re[i] = atof(arg[ i*6 + 4]);
    rhokStar_re[i] = atof(arg[ i*6 + 5]);

    kappa_im[i] = atof(arg[ i*6 + 6 ]);
    rhokStar_im[i] = atof(arg[ i*6 + 7 ]);

    wavenumber[i] = atoi(arg[ i*6 + 8 ]);
  }
  nevery = atoi(arg[narg-1]);
  if(nevery <= 0) error->all(FLERR,"Illegal fix print step");
} // end of constructor

/* ---------------------------------------------------------------------- */

FixRhok::~FixRhok()
{
  // delete locally stored array
}

/* ---------------------------------------------------------------------- */

int FixRhok::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixRhok::init(){
  if(comm->me==0)
  for(int li = 0 ; li < nloop; li++)
  printf("#wavenumber = %d\tkappa = %e\trhok*=%e\n#",wavenumber[li],kappa_re[li],rhokStar_re[li]);
}

/* ---------------------------------------------------------------------- */
void FixRhok::post_force(int vflag)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double v[6];
  double **f=atom->f;

  lbox[0] = domain->boxhi[0] - domain->boxlo[0] ;
  lbox[1] = domain->boxhi[1] - domain->boxlo[1] ;
  lbox[2] = domain->boxhi[2] - domain->boxlo[2] ;

  // energy and virial setup
  if(vflag) v_setup(vflag);
  else evflag = 0;

  compute_vector();

  for (int i = 0; i < nlocal; i++){
    for(int li = 0; li < nloop; li++){
      if (mask[i] & groupbit) {
        f[i][alpha[li]] += kappa_re[li]*(globalRhok[2*li]-rhokStar_re[li]) * sin( 2.0 * M_PI * x[i][alpha[li]] * double(wavenumber[li]) / lbox[alpha[li]] ) * 2.0 * M_PI * double(wavenumber[li]) / lbox[alpha[li]];
        f[i][alpha[li]] += kappa_im[li]*(globalRhok[2*li+1]-rhokStar_im[li]) * cos( 2.0 * M_PI * x[i][alpha[li]] * double(wavenumber[li]) / lbox[alpha[li]] ) * 2.0 * M_PI * double(wavenumber[li]) / lbox[alpha[li]];
      }
    }   
  } 

/*  if (evflag) {
    for(int i = 0 ; i < nlocal){
      v[0] = f[i][0] * x[i][0];
      v[1] = f[i][1] * x[i][1];
      v[2] = f[i][2] * x[i][2];
      v[3] = f[i][0] * x[i][1];
      v[4] = f[i][0] * x[i][2];
      v[5] = f[i][1] * x[i][2];
      v_tally(i,v);  
    }
  }
*/

}   // end of post_force()

/* ---------------------------------------------------------------------- */

double FixRhok::memory_usage()
{
  double bytes;

  return bytes;
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 * to output the ftrho
 * --------------------------------------------------------------------*/
void FixRhok::end_of_step()
{
  double step = update->ntimestep - update->beginstep;
  if(comm->me == 0) {
    printf("%.2f",step);
    for(int i=0;i<nloop;i++){
	  printf("\t %e \t %e",
        globalRhok[2*i],globalRhok[2*i+1]);
    }
    printf("\n");
  }
}   // end of postprocess

/* --------------------------------------------------------------------*/
void FixRhok::compute_vector(){

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for(int li=0;li<nloop;li++){
    localRhok[2*li] = 0.0;
    localRhok[2*li+1] = 0.0;
    for(int i = 0; i < nlocal; i++){
	  if(mask[i] & groupbit){
        localRhok[2*li] +=  cos( 2.0 * M_PI * x[i][alpha[li]] * double(wavenumber[li]) / lbox[alpha[li]] );	// rhok[0] is real, 1 is imaginary
        localRhok[2*li+1] += - sin( 2.0 * M_PI * x[i][alpha[li]] * double(wavenumber[li]) / lbox[alpha[li]]);	
      }
    }
  }

  MPI_Allreduce(localRhok,globalRhok,2*nloop,MPI_DOUBLE,MPI_SUM,world);	// sum over all processor
}
