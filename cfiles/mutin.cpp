#include <math.h> 
#include <stdlib.h> 
#include "mex.h" 
  
  
double sqr(double x) 
{ 
 return x*x; 
} 
  
struct prvek{ 
              double dat; 
              int ind; 
            }; 
  
struct prvekz{ 
              int x1,y1,x2,y2; 
              void set(int nx1,int ny1,int nx2,int ny2){ 
                        x1=nx1; x2=nx2; y1=ny1; y2=ny2;}; 
            }; 
  
struct zasobnik{ 
       int index; 
       prvekz data[100]; 
       zasobnik(){index=0;}; 
       void add(int x1, int y1, int x2, int y2) 
                   { data[++index].x1=x1; 
                     data[index].x2=x2; 
                     data[index].y1=y1; 
                     data[index].y2=y2;}; 
       void add(prvekz X){data[++index]=X;}; 
       void uber(){index--;}; 
       prvekz &aktiv(){return data[index];}; 
               }; 
  
int sort_function( const void *a, const void *b) 
{ 
 int k=-1; 
 prvek *A,*B; 
 A=(prvek*)a;B=(prvek*)b; 
 if (A->dat > B->dat) k=1; 
 if (A->dat == B->dat) k=0; 
 return k; 
} 
  
double MI1(double *X,double *Y,int Pocet) 
{ 
 prvek *x=new prvek[Pocet]; 
 prvek *y=new prvek[Pocet]; 
 //prvek x[Pocet];prvek y[Pocet]; 
 int *xx=new int[Pocet]; 
 int *yy=new int[Pocet]; 
// int xx[Pocet];int yy[Pocet]; 
 for(int i=0;i<Pocet;i++) { 
                            x[i].dat=X[i]; 
                            x[i].ind=i; 
                            y[i].dat=Y[i]; 
                            y[i].ind=i; 
                          } 
 qsort((void *)x, Pocet, sizeof(prvek), sort_function); 
 qsort((void *)y, Pocet, sizeof(prvek), sort_function); 
 for(int i=0;i<Pocet;i++) { 
                           xx[x[i].ind]=i; 
                           yy[y[i].ind]=i; 
                          } 
 zasobnik marg; 
 prvekz amarg[4]; 
 int poc[100],kon[100]; 
 double m=0; 
 double tst; 
 int Nex,avex,avey; 
// int poradi[Pocet],apor[Pocet]; 
 int *poradi=new int[Pocet]; 
 int *apor=new int[Pocet]; 
 for(int i=0;i<Pocet;i++) poradi[i]=i; //indexy se pocitaji od nuly 
 int run=0; 
 int NN[4]={0,0,0,0}; 
 int apoc,akon; 
 int Nx,Ny; 
 poc[1]=0; 
 kon[1]=Pocet-1; 
 marg.add(0,0,Pocet-1,Pocet-1); 
 while (marg.index>0) 
 { 
  run++; 
  NN[0]=0;NN[1]=0;NN[2]=0;NN[3]=0; 
  apoc=poc[marg.index]; 
  akon=kon[marg.index]; 
  for(int i=apoc;i<=akon;i++) apor[i]=poradi[i]; 
  Nex=akon-apoc+1; 
  avex=(marg.aktiv().x1+marg.aktiv().x2)/2; 
  avey=(marg.aktiv().y1+marg.aktiv().y2)/2; 
  for (int i=apoc;i<=akon;i++) 
    { 
     if (xx[apor[i]]<=avex) if (yy[apor[i]]<=avey) NN[0]++;else NN[1]++; 
     else if (yy[apor[i]]<=avey) NN[2]++;else NN[3]++; 
    }//for  NN=sum(I) 
  amarg[0].set(marg.aktiv().x1, marg.aktiv().y1, avex, avey); 
  amarg[1].set(marg.aktiv().x1, avey+1, avex, marg.aktiv().y2); 
  amarg[2].set(avex+1, marg.aktiv().y1, marg.aktiv().x2, avey); 
  amarg[3].set(avex+1, avey+1, marg.aktiv().x2, marg.aktiv().y2); 
  tst=(double)4*(sqr(NN[0]-(double)Nex/4)+sqr(NN[1]-(double)Nex/4) 
                +sqr(NN[2]-(double)Nex/4)+sqr(NN[3]-(double)Nex/4))/Nex; 
  if ((tst>7.8)||(run==1)) 
    { 
     marg.uber(); 
     int ap=apoc, ak=akon; 
     for(int k=0;k<4;k++) 
       { 
        if (NN[k]>2){ 
                     marg.add(amarg[k]); 
                     akon=apoc+NN[k]-1; 
                     poc[marg.index]=apoc; 
                     kon[marg.index]=akon; 
                     int j=apoc; 
                     switch (k) { 
                       case 0:for(int i=ap;i<=ak;i++) 
                               if ((xx[apor[i]]<=avex)&&(yy[apor[i]]<=avey)) 
                                              poradi[j++]=apor[i]; 
                              break; 
                       case 1:for(int i=ap;i<=ak;i++) 
                               if ((xx[apor[i]]<=avex)&&(yy[apor[i]]>avey)) 
                                              poradi[j++]=apor[i]; 
                              break; 
                       case 2:for(int i=ap;i<=ak;i++) 
                               if ((xx[apor[i]]>avex)&&(yy[apor[i]]<=avey)) 
                                              poradi[j++]=apor[i]; 
                              break; 
                       case 3:for(int i=ap;i<=ak;i++) 
                               if ((xx[apor[i]]>avex)&&(yy[apor[i]]>avey)) 
                                              poradi[j++]=apor[i]; 
                                };//switch 
                     apoc=akon+1; 
                    } 
        else        { 
                     if (NN[k]>0) { 
                        Nx=amarg[k].x2-amarg[k].x1+1; 
                        Ny=amarg[k].y2-amarg[k].y1+1; 
                        m+=(double)NN[k]*log((double)NN[k]/(Nx*Ny)); 
                                  }//if 
                    } //else 
       }//for 
    } //if test 
   else 
   { 
    Nx=marg.aktiv().x2-marg.aktiv().x1+1; 
    Ny=marg.aktiv().y2-marg.aktiv().y1+1; 
    m+=(double)Nex*log((double)Nex/(Nx*Ny)); 
    marg.uber(); 
   }//else test 
 } //while 
 delete []x; 
 delete []xx; 
 delete []y; 
 delete []yy; 
 delete []poradi; 
 delete []apor; 
 return m/Pocet+log(Pocet); 
} 
  
  
void mexFunction( 
                 int nlhs,       mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[] 
   ) 
{ 
  double *yp; 
  double *Zdroj1,*Zdroj2,*vystup; 
  int           Pocet; 
  unsigned int rows,cols; 
  
  /* Check for proper number of arguments */ 
  
  if (nrhs > 1) { 
    mexErrMsgTxt("MUTIN requires one input argument."); 
  } else if (nlhs > 1) { 
    mexErrMsgTxt("MUTIN requires one output argument."); 
  } 
  
  
  /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
  
  rows = mxGetM(prhs[0]); 
  cols = mxGetN(prhs[0]); 
  if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || 
      mxIsSparse(prhs[0])  || !mxIsDouble(prhs[0]) || 
      (cols != 2)) { 
    mexErrMsgTxt("MUTIN requires that Y be a N x 2 vector."); 
  } 
  
  
  /* Create a matrix scalar) for the return argument */ 
  
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); 
  
  
  /* Assign pointers to the various parameters */ 
  
  Zdroj1 = mxGetPr(prhs[0]); 
  Zdroj2 = Zdroj1+rows; 
  vystup = mxGetPr(plhs[0]); 
  /* Do the actual computations in a subroutine */ 
  
  *vystup = MI1(Zdroj1,Zdroj2,rows); 
  return; 
} 
    