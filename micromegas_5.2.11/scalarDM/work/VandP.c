#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=14;
static ModelPrtclsStr ModelPrtcls_[14]=
{
  {"ve","ve~",0, 12, "0","0",1,1,1,0}
, {"vm","vm~",0, 14, "0","0",1,1,1,0}
, {"vt","vt~",0, 16, "0","0",1,1,1,0}
, {"e-","e+",0, 11, "Me","0",1,1,2,-3}
, {"mu-","mu+",0, 13, "Mmu","0",1,1,2,-3}
, {"ta-","ta+",0, 15, "Ml","0",1,1,2,-3}
, {"u","u~",0, 2, "Mu","0",1,3,6,2}
, {"c","c~",0, 4, "Mc","0",1,3,6,2}
, {"t","t~",0, 6, "Mt","WT",1,3,6,2}
, {"d","d~",0, 1, "Md","0",1,3,6,-1}
, {"s","s~",0, 3, "Ms","0",1,3,6,-1}
, {"b","b~",0, 5, "Mb","0",1,3,6,-1}
, {"~dm","~DM",0, 100100, "mdm","0",0,1,1,0}
, {"Ap","Ap",1, 100023, "mAp","wAp",2,1,3,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=52;
int nModelFunc=7;
static int nCurrentVars=51;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[59]={
 "aEWM1","Mh","Gf","aS","ymdo","ymup","yms","ymc","ymb","ymt"
,"yme","ymm","ymtau","QS","yHEs0lR","lams0h","lams0","MZ","Me","Mmu"
,"Ml","Mu","Mc","Mt","Md","Ms","Mb","mh","mHE","E"
,"Pi","mdm","mAp","gD","eps","eps_e","eps_mu","eps_tau","eps_nu","eps_d"
,"eps_u","eps_s","eps_c","eps_t","eps_b","MW","alfSMZ","McMc","MbMb","Mtp"
,"EE","Q","CW","SW","LamQCD","GF","Mb","Mt","Mc"};
char**varNames=varNames_;
static REAL varValues_[59]={
   1.279000E+02,  1.250000E+02,  1.166370E-05,  1.184000E-01,  5.040000E-03,  2.550000E-03,  1.010000E-01,  1.270000E+00,  4.700000E+00,  1.720000E+02
,  5.110000E-04,  1.056600E-01,  1.777000E+00,  1.000000E+02,  2.000000E-05,  2.000000E-02,  1.000000E-01,  9.118760E+01,  5.110000E-04,  1.056600E-01
,  1.777000E+00,  2.550000E-03,  1.270000E+00,  1.720000E+02,  5.040000E-03,  1.010000E-01,  4.700000E+00,  1.250000E+02,  1.000000E+02,  2.718282E+00
,  3.141593E+00,  1.000000E+00,  3.000000E+00,  5.000000E-01,  1.000000E-11,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  8.038500E+01,  1.184000E-01,  1.200000E+00,  4.230000E+00,  1.725000E+02
,  3.134300E-01,  1.000000E+02};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   nCurrentVars=52;
   V[52]=V[45]/(V[17]);
   if(!isfinite(V[52]) || FError) return 52;
   nCurrentVars=53;
   V[53]=Sqrt(1-Pow(V[52],2));
   if(!isfinite(V[53]) || FError) return 53;
   nCurrentVars=54;
   V[54]=initQCD5(V[46],V[47],V[48],V[49]);
   if(!isfinite(V[54]) || FError) return 54;
   nCurrentVars=55;
   V[55]=Pow(V[50],2)/(Pow(2*V[53]*V[45],2))/(M_SQRT2);
   if(!isfinite(V[55]) || FError) return 55;
 FirstQ:
 cErr=1;
   nCurrentVars=56;
   V[56]=MbEff(V[51]);
   if(!isfinite(V[56]) || FError) return 56;
   nCurrentVars=57;
   V[57]=MtEff(V[51]);
   if(!isfinite(V[57]) || FError) return 57;
   nCurrentVars=58;
   V[58]=McEff(V[51]);
   if(!isfinite(V[58]) || FError) return 58;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
