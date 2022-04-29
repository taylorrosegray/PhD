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
, {"~dm","~DM",0, 5000531, "mdm","0",2,1,3,0}
, {"Ap","Ap",1, 2001, "mAp","wAp",2,1,3,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=47;
int nModelFunc=3;
static int nCurrentVars=46;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[50]={
 "aEWM1","Mh","Gf","aS","ymdo","ymup","yms","ymc","ymb","ymt"
,"yme","ymm","ymtau","QS","yHEs0lR","lams0h","lams0","MZ","Me","Mmu"
,"Ml","Mu","Mc","Mt","Md","Ms","Mb","mh","mHE","E"
,"Pi","mdm","mAp","gD","gSM","MW","alfSMZ","McMc","MbMb","Mtp"
,"EE","Q","b5","b6R","b6I","b7R","b7I","CW","SW","LamQCD"
};
char**varNames=varNames_;
static REAL varValues_[50]={
   1.279000E+02,  1.250000E+02,  1.166370E-05,  1.184000E-01,  5.040000E-03,  2.550000E-03,  1.010000E-01,  1.270000E+00,  4.700000E+00,  1.720000E+02
,  5.110000E-04,  1.056600E-01,  1.777000E+00,  1.000000E+02,  2.000000E-05,  2.000000E-02,  1.000000E-01,  9.118760E+01,  5.110000E-04,  1.056600E-01
,  1.777000E+00,  2.550000E-03,  1.270000E+00,  1.720000E+02,  5.040000E-03,  1.010000E-01,  4.700000E+00,  1.250000E+02,  1.000000E+02,  2.718282E+00
,  3.141593E+00,  1.000000E+00,  3.000000E+00,  5.000000E-01,  1.000000E+00,  8.038500E+01,  1.184000E-01,  1.200000E+00,  4.230000E+00,  1.725000E+02
,  3.134300E-01,  1.000000E+02,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00};
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
     if(i==nModelVars)     return 0;
   }
  cErr=1;
   nCurrentVars=47;
   V[47]=V[35]/(V[17]);
   if(!isfinite(V[47]) || FError) return 47;
   nCurrentVars=48;
   V[48]=Sqrt(1-Pow(V[47],2));
   if(!isfinite(V[48]) || FError) return 48;
   nCurrentVars=49;
   V[49]=initQCD5(V[36],V[37],V[38],V[39]);
   if(!isfinite(V[49]) || FError) return 49;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
