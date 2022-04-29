#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=19;
static ModelPrtclsStr ModelPrtcls_[19]=
{
  {"G","G",1, 21, "0","0",2,8,16,0}
, {"A","A",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","wZ",2,1,3,0}
, {"W+","W-",0, 24, "MW","wW",2,1,3,3}
, {"h","h",1, 25, "Mh","wh",0,1,1,0}
, {"e","E",0, 11, "0","0",1,1,2,-3}
, {"ne","Ne",0, 12, "0","0",1,1,1,0}
, {"m","M",0, 13, "Mm","0",1,1,2,-3}
, {"nm","Nm",0, 14, "0","0",1,1,1,0}
, {"l","L",0, 15, "Ml","0",1,1,2,-3}
, {"nl","Nl",0, 16, "0","0",1,1,1,0}
, {"d","D",0, 1, "Mm","0",1,3,6,-1}
, {"u","U",0, 2, "Mm","0",1,3,6,2}
, {"s","S",0, 3, "Mm","0",1,3,6,-1}
, {"c","C",0, 4, "Mc","0",1,3,6,2}
, {"b","B",0, 5, "Mb","0",1,3,6,-1}
, {"t","T",0, 6, "Mt","wt",1,3,6,2}
, {"Zp","Zp",1, 100023, "MZp","wZp",2,1,3,0}
, {"~dm","~Dm",0, 100100, "Mdm","0",1,1,2,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=28;
int nModelFunc=7;
static int nCurrentVars=27;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[35]={
 "EE","alfSMZ","Mm","Ml","Mq","Q","McMc","MbMb","Mtp","MZ"
,"MW","Mh","Mdm","MZp","gZp","YpU","YpD","YpS","YpC","YpB"
,"YpT","YpE","YpM","YpL","YpNe","YpNm","YpNl","gZpDm","CW","SW"
,"LamQCD","GF","Mb","Mt","Mc"};
char**varNames=varNames_;
static REAL varValues_[35]={
   3.134300E-01,  1.184000E-01,  1.057000E-01,  1.777000E+00,  5.000000E-02,  1.000000E+02,  1.200000E+00,  4.230000E+00,  1.725000E+02,  9.118800E+01
,  8.038500E+01,  1.250000E+02,  1.000000E+02,  3.000000E+02,  1.000000E-05,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00
,  1.000000E+00, -3.000000E+00, -3.000000E+00, -3.000000E+00, -3.000000E+00, -3.000000E+00, -3.000000E+00,  1.000000E-10};
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
   nCurrentVars=28;
   V[28]=V[10]/(V[9]);
   if(!isfinite(V[28]) || FError) return 28;
   nCurrentVars=29;
   V[29]=Sqrt(1-Pow(V[28],2));
   if(!isfinite(V[29]) || FError) return 29;
   nCurrentVars=30;
   V[30]=initQCD5(V[1],V[6],V[7],V[8]);
   if(!isfinite(V[30]) || FError) return 30;
   nCurrentVars=31;
   V[31]=Pow(V[0],2)/(Pow(2*V[29]*V[10],2))/(M_SQRT2);
   if(!isfinite(V[31]) || FError) return 31;
 FirstQ:
 cErr=1;
   nCurrentVars=32;
   V[32]=MbEff(V[5]);
   if(!isfinite(V[32]) || FError) return 32;
   nCurrentVars=33;
   V[33]=MtEff(V[5]);
   if(!isfinite(V[33]) || FError) return 33;
   nCurrentVars=34;
   V[34]=McEff(V[5]);
   if(!isfinite(V[34]) || FError) return 34;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
