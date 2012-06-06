#include "mpMap.h"

// set this up so that i can run a scan over a grid of positions, calculate
// likelihood at each point. 
void lkhd3pt(int *finalg, int *founderg, int *id, int *mother, int *father, int *nfinals, int *nfounders, int *nmrk, int *nped, double *out, double *r12, double *r23, double *r13, int *left, int *mid, int *right, int *npos)
{
 int i, j, k, l, m;
 int *funnel, *geni, *hapi;
 double *hprob, *p, *x;
 int *cf1, *cf2, *cf3;
 int *fou;
 double hp;

 fou = (int*) R_alloc(*nfounders*3, sizeof(int));
 p = (double*) R_alloc(10, sizeof(double));
 x = (double*) R_alloc(10, sizeof(double));
 hprob = (double*) R_alloc(20, sizeof(double));
 funnel = (int*) R_alloc(*nfounders, sizeof(int));
 geni = (int*) R_alloc(*nmrk, sizeof(int));
 hapi = (int*) R_alloc(*nmrk, sizeof(int));
 cf1 = (int*) R_alloc(*nfounders, sizeof(int));
 cf2 = (int*) R_alloc(*nfounders, sizeof(int));
 cf3 = (int*) R_alloc(*nfounders, sizeof(int));
 
 for (i=0; i<*npos; i++)
	out[i] = 0;

 for (i=0; i<*npos; i++)
 for (j=0; j<*nfinals; j++) 
 {
   // calculate haplotype probabilities, then sum over all possibilities?
   if (*nfounders==4) {
	hp4way(r12[i], r23[i], r13[i], p, x);
	for (k=0; k<10; k++)	hprob[k] = p[k]*x[k]; }

   if (*nfounders==8)
	hp8way(r12[i], r23[i], r13[i], hprob);

   geni[0] = finalg[j*(*nmrk+1)+left[i]];
   geni[1] = finalg[j*(*nmrk+1)+mid[i]]; 
   geni[2] = finalg[j*(*nmrk+1)+right[i]];
   for (l=0; l<*nfounders; l++) {
 	fou[l] = founderg[l*(*nmrk)+left[i]-1];
   	fou[(*nfounders)+l] = founderg[l*(*nmrk)+mid[i]-1];
	fou[2*(*nfounders)+l] = founderg[l*(*nmrk)+right[i]-1]; }

   pedfunnel(finalg[j*(*nmrk+1)], id, mother, father, funnel, *nfounders, *nped);
   for (k=0; k<*nfounders; k++)
   {
 	cf1[k] = (fou[(funnel[k]-1)]==geni[0]);
	cf2[k] = (fou[*nfounders+(funnel[k]-1)]==geni[1]);
	cf3[k] = (fou[2*(*nfounders)+(funnel[k]-1)]==geni[2]);
	if (geni[0]<0) cf1[k] = 1;
	if (geni[1]<0) cf2[k] = 1;
 	if (geni[2]<0) cf3[k] = 1;
   }

  hp = 0;
  for (k=0; k<*nfounders; k++)
  for (l=0; l<*nfounders; l++)
  for (m=0; m<*nfounders; m++)
  if (cf1[k]*cf2[l]*cf3[m]) 
  {
	hapi[0] = k;
	hapi[1] = l;
	hapi[2] = m;
	hp += hprob[hclass(hapi, *nfounders)];
  }

  out[i] += log10(hp);

  // end of loop over individuals
 }// end of loop over positions

 return;
}
