#include "mpMap.h"
// Also need to write R wrapper

void pr2pt4way(double r, double *prob)
{
  prob[0] = (1-r)/(4+8*r);
  prob[1] = r/(4+8*r);
  prob[2] = r/(4+8*r);
}

void pr2pt8way(double r, double *prob)
{
  prob[0] = (1-r)*(1-r)/(8+16*r);
  prob[1] = r*(1-r)/(8+16*r);
  prob[2] = r/2/(8+16*r);
}

void pr2ptirip(double r, int s, int n, double *prob)
{
  int x=pow(2,n);
  prob[0]=(pow(1-r, n+s-1)/x+(2*r+1-pow(1-r, s-1))/x/x)/(1+2*r); 
  prob[1]=prob[2]=(1-x*prob[0])/(x*(x-1));
}

int hclass2pt(int h1, int h2)
{
  int hsmall = (h1<=h2)*h1+(h2<h1)*h2;
  int hbig = (h1>h2)*h1+(h2>=h1)*h2;
  return(2-2*(h1==h2)-((hsmall%2)==1)*((hbig-hsmall)==1));
}

// Note that this requires that the id, mother, father portions of the pedigree are a) integers and b) the founders are in the first 8 positions

void rfhaps(int *finalg, int *founderg, int *id, int *mother, int *father, int *pair1, int *pair2, int *nfinals, int *nfounders, int *nmrk, int *npairs, int *nped, int *ngen, double *thvec, int *nr, double *pairth)
{
 int i, j, k, l, r;
 int p1, p2;
 int *cf1, *cf2, *funnel, *geni;
 double *probclass;
 double theta; 
 double hp;
 int n = (int) log(*nfounders)/log(2);

 probclass = (double*) R_alloc(3, sizeof(double));

 cf1 = (int*) R_alloc(*nfounders, sizeof(int));
 cf2 = (int*) R_alloc(*nfounders, sizeof(int));
 funnel = (int*) R_alloc(*nfounders, sizeof(int));
 geni = (int*) R_alloc(*nmrk, sizeof(int));

 for (i=0; i<(*npairs)*(*nr); i++)
 	pairth[i] = 0;

 for (i=0; i<*nfinals; i++)
 {
   for (j=0; j<*nmrk; j++)
	geni[j] = finalg[i*(*nmrk+1)+j+1];

   pedfunnel(finalg[i*(*nmrk+1)], id, mother, father, funnel, *nfounders, *nped);
   for (k=0; k<*npairs; k++)
   {
     p1 = pair1[k]-1;
     p2 = pair2[k]-1;

     if ((geni[p1]>-1)*(geni[p2]>-1))
     {
      for (j=0; j<*nfounders; j++)
      {
	cf1[j] = (founderg[(funnel[j]-1)*(*nmrk)+p1]==geni[p1]);
	cf2[j] = (founderg[(funnel[j]-1)*(*nmrk)+p2]==geni[p2]);
      }

      for (r=0; r<*nr; r++)
      {
   	theta = thvec[r];
 	if (*ngen==0) {
	  if (*nfounders==4)  pr2pt4way(theta, probclass);
	  if (*nfounders==8)  pr2pt8way(theta, probclass);
	}
	if (*ngen > 0)
	  pr2ptirip(theta, *ngen, n, probclass);


	// need to construct the combinations of cf1 and cf2

	hp = 0;
	for (j=0; j<*nfounders; j++)
	for (l=0; l<*nfounders; l++)
	if (cf1[j]*cf2[l]) 
	  hp += probclass[hclass2pt(j, l)];

	pairth[k*(*nr)+r] += log10(hp);

       } // end of theta loop

   } // end of NA check

   } // end of npairs loop

 } // end of nfinals loop

 return;
}
