#include "mpMap.h"

void pf(int *obs, int*id, int *mother, int *father, int *nfinals, int *nfounders, int *nped, int *out)
{
 int j, k;
 int *funnel;

 funnel = (int*) R_alloc(*nfounders, sizeof(int));
 
 for (j=0; j<*nfinals; j++)
 {
   pedfunnel(obs[j], id, mother, father, funnel, *nfounders, *nped);

   for (k=0; k<*nfounders; k++)
	out[j*(*nfounders)+k] = funnel[k];
 } // end of loop over individuals

 return;
}
