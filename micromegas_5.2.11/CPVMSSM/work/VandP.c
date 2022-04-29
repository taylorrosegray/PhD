#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=48;
static ModelPrtclsStr ModelPrtcls_[48]=
{
  {"A","A",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","wZ",2,1,3,0}
, {"W+","W-",0, 24, "MW","wW",2,1,3,3}
, {"G","G",1, 21, "0","0",2,8,16,0}
, {"e","E",0, 11, "0","0",1,1,2,-3}
, {"ne","Ne",0, 12, "0","0",1,1,1,0}
, {"m","M",0, 13, "0","0",1,1,2,-3}
, {"nm","Nm",0, 14, "0","0",1,1,1,0}
, {"l","L",0, 15, "Ml","0",1,1,2,-3}
, {"nl","Nl",0, 16, "0","0",1,1,1,0}
, {"u","U",0, 2, "Mq","0",1,3,6,2}
, {"d","D",0, 1, "Mq","0",1,3,6,-1}
, {"s","S",0, 3, "Mq","0",1,3,6,-1}
, {"c","C",0, 4, "Mc","0",1,3,6,2}
, {"b","B",0, 5, "Mb","0",1,3,6,-1}
, {"t","T",0, 6, "Mt","wt",1,3,6,2}
, {"h1","h1",1, 25, "Mh1","wh1",0,1,1,0}
, {"h2","h2",1, 35, "Mh2","wh2",0,1,1,0}
, {"h3","h3",1, 36, "Mh3","wh3",0,1,1,0}
, {"H+","H-",0, 37, "MHc","wHc",0,1,1,3}
, {"~1+","~1-",0, 1000024, "MC1","wC1",1,1,2,3}
, {"~2+","~2-",0, 1000037, "MC2","wC2",1,1,2,3}
, {"~o1","~o1",1, 1000022, "MNE1","0",1,1,2,0}
, {"~o2","~o2",1, 1000023, "MNE2","wNE2",1,1,2,0}
, {"~o3","~o3",1, 1000025, "MNE3","wNE3",1,1,2,0}
, {"~o4","~o4",1, 1000035, "MNE4","wNE4",1,1,2,0}
, {"~g","~g",1, 1000021, "aM3","wSG",1,8,16,0}
, {"~eL","~EL",0, 1000011, "MSeL","wSeL",0,1,1,-3}
, {"~eR","~ER",0, 2000011, "MSeR","wSeR",0,1,1,-3}
, {"~mL","~ML",0, 1000013, "MSmL","wSmL",0,1,1,-3}
, {"~mR","~MR",0, 2000013, "MSmR","wSmR",0,1,1,-3}
, {"~l1","~L1",0, 1000015, "MSl1","wSl1",0,1,1,-3}
, {"~l2","~L2",0, 2000015, "MSl2","wSl2",0,1,1,-3}
, {"~ne","~Ne",0, 1000012, "MSne","wSne",0,1,1,0}
, {"~nm","~Nm",0, 1000014, "MSnm","wSnm",0,1,1,0}
, {"~nl","~Nl",0, 1000016, "MSnl","wSnl",0,1,1,0}
, {"~dL","~DL",0, 1000001, "MSdL","wSdL",0,3,3,-1}
, {"~dR","~DR",0, 2000001, "MSdR","wSdR",0,3,3,-1}
, {"~uL","~UL",0, 1000002, "MSuL","wSuL",0,3,3,2}
, {"~uR","~UR",0, 2000002, "MSuR","wSuR",0,3,3,2}
, {"~sL","~SL",0, 1000003, "MSsL","wSsL",0,3,3,-1}
, {"~sR","~SR",0, 2000003, "MSsR","wSsR",0,3,3,-1}
, {"~cL","~CL",0, 1000004, "MScL","wScL",0,3,3,2}
, {"~cR","~CR",0, 2000004, "MScR","wScR",0,3,3,2}
, {"~b1","~B1",0, 1000005, "MSb1","wSb1",0,3,3,-1}
, {"~b2","~B2",0, 2000005, "MSb2","wSb2",0,3,3,-1}
, {"~t1","~T1",0, 1000006, "MSt1","wSt1",0,3,3,2}
, {"~t2","~T2",0, 2000006, "MSt2","wSt2",0,3,3,2}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=48;
int nModelFunc=279;
static int nCurrentVars=47;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[327]={
 "alfSMZ","Mtp","MbMb","McMc","Q","EE","SW","Ml","MHc","aMu"
,"fiMu","Mq","Au","Ad","aM1","aM2","aM3","fiM1","fiM2","fiM3"
,"Ml2","Ml3","Mr2","Mr3","aAt","fiAt","aAb","fiAb","aAl","fiAl"
,"aAm","fiAm","aAe","fiAe","aAu","fiAu","aAd","fiAd","Mq2","Mq3"
,"Mu2","Mu3","Md2","Md3","MZ","Maux","tb","MSG","Lqcd","MbMt"
,"McMt","SHok","RdOK","CW","MW","S2W","C2W","sb","cb","s2b"
,"c2b","Mh1","Mh2","Mh3","Zh11","Zh12","Zh13","Zh21","Zh22","Zh23"
,"Zh31","Zh32","Zh33","la1","la2","la3","la4","la5r","la6r","la7r"
,"la5i","la6i","la7i","MNE1","MNE2","MNE3","MNE4","Zn11r","Zn12r","Zn13r"
,"Zn14r","Zn21r","Zn22r","Zn23r","Zn24r","Zn31r","Zn32r","Zn33r","Zn34r","Zn41r"
,"Zn42r","Zn43r","Zn44r","Zn11i","Zn12i","Zn13i","Zn14i","Zn21i","Zn22i","Zn23i"
,"Zn24i","Zn31i","Zn32i","Zn33i","Zn34i","Zn41i","Zn42i","Zn43i","Zn44i","MC1"
,"MC2","Zv11r","Zv12r","Zv21r","Zv22r","Zv11i","Zv12i","Zv21i","Zv22i","Zu11r"
,"Zu12r","Zu21r","Zu22r","Zu11i","Zu12i","Zu21i","Zu22i","Tu3r","Tu3i","Td3r"
,"Td3i","MSl1","MSl2","Zl11r","Zl21r","Zl21i","MSt1","MSt2","Zt11r","Zt21r"
,"Zt21i","MSb1","MSb2","Zb11r","Zb21r","Zb21i","PI","mur","mui","Alr"
,"Ali","Abr","Abi","Atr","Ati","Adr","Adi","Aur","Aui","sfiM3"
,"cfiM3","MG1r","MG1i","MG2r","MG2i","MG3r","MG3i","sfi3","cfi3","MSne"
,"MSnm","MSnl","MSeL","MSeR","MSmL","MSmR","MSuL","MScL","MSuR","MScR"
,"MSdL","MSsL","MSdR","MSsR","Mt","Mb","Mc","Mbp","Mcp","Quq"
,"Qdq","ah1F_c","ah1F_b","ah1F_t","ah1F_l","ah2F_c","ah2F_b","ah2F_t","ah2F_l","ah3F_c"
,"ah3F_b","ah3F_t","ah3F_l","ah1S_eL","ah1S_eR","ah1S_mL","ah1S_mR","ah1S_uL","ah1S_uR","ah1S_cL"
,"ah1S_cR","ah1S_dL","ah1S_dR","ah1S_sL","ah1S_sR","ah1S_l1","ah1S_l2","ah1S_t1","ah1S_t2","ah1S_b1"
,"ah1S_b2","ah2S_eL","ah2S_eR","ah2S_mL","ah2S_mR","ah2S_uL","ah2S_uR","ah2S_cL","ah2S_cR","ah2S_dL"
,"ah2S_dR","ah2S_sL","ah2S_sR","ah2S_l1","ah2S_l2","ah2S_t1","ah2S_t2","ah2S_b1","ah2S_b2","ah3S_eL"
,"ah3S_eR","ah3S_mL","ah3S_mR","ah3S_uL","ah3S_uR","ah3S_cL","ah3S_cR","ah3S_dL","ah3S_dR","ah3S_sL"
,"ah3S_sR","ah3S_l1","ah3S_l2","ah3S_t1","ah3S_t2","ah3S_b1","ah3S_b2","ah1V_W","B00000","B00001"
,"B00002","B00003","ah1S_Hc","ah2V_W","B00004","B00005","B00006","B00007","ah2S_Hc","ah3V_W"
,"B00008","B00009","B00010","B00011","ah3S_Hc","ah1F_c1","ah1F_c2","ah2F_c1","ah2F_c2","ah3F_c1"
,"ah3F_c2","imah1F_c","imah1F_b","imah1F_t","imah1F_l","imah2F_c","imah2F_b","imah2F_t","imah2F_l","imah3F_c"
,"imah3F_b","imah3F_t","imah3F_l","imah1F_c1","imah1F_c2","imah2F_c1","imah2F_c2","imah3F_c1","imah3F_c2","aQCDh1"
,"Rqcdh1","aQCDh2","Rqcdh2","aQCDh3","Rqcdh3","LGGh1","LGGh2","LGGh3","LAAh1","LAAh2"
,"LAAh3","imLGGh1","imLGGh2","imLGGh3","imLAAh1","imLAAh2","imLAAh3"};
char**varNames=varNames_;
static REAL varValues_[327]={
   1.184000E-01,  1.730700E+02,  4.200000E+00,  1.300000E+00,  1.000000E+02,  3.133300E-01,  4.740000E-01,  1.777000E+00,  2.500000E+02,  1.500000E+02
,  0.000000E+00,  5.000000E-02,  0.000000E+00,  0.000000E+00,  1.000000E+02,  2.000000E+02,  1.000000E+03,  0.000000E+00,  0.000000E+00,  7.000000E+01
,  1.000000E+03,  5.000000E+02,  1.000000E+03,  5.000000E+02,  1.000000E+03,  9.000000E+01,  1.000000E+03,  9.000000E+01,  1.000000E+03,  9.000000E+01
,  1.000000E+02,  0.000000E+00,  1.000000E+03,  9.000000E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.000000E+03,  5.000000E+02
,  1.000000E+03,  5.000000E+02,  1.000000E+03,  5.000000E+02,  9.118700E+01,  1.000000E+00,  5.000000E+01,  1.000000E+03};
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
   nCurrentVars=48;
   V[48]=initQCD5(V[0],V[3],V[2],V[1]);
   if(!isfinite(V[48]) || FError) return 48;
   nCurrentVars=49;
   V[49]=MbRun(V[1]);
   if(!isfinite(V[49]) || FError) return 49;
   nCurrentVars=50;
   V[50]=McRun(V[1]);
   if(!isfinite(V[50]) || FError) return 50;
   nCurrentVars=51;
   V[51]=cpHiggs(V[5],V[0],V[1],V[50],V[49],V[7],V[6],V[46],V[8],V[9],V[10],V[14],V[17],V[15],V[18],V[16],V[19],V[39],V[41],V[43],V[21],V[23],V[24],V[25],V[26],V[27],V[28],V[29],V[34],V[35],V[36],V[37],V[32],V[33],V[30],V[31]);
   if(!isfinite(V[51]) || FError) return 51;
   nCurrentVars=52;
   V[52]=slhaRead("cpsuperh2_slha.out",4);
   if(!isfinite(V[52]) || FError) return 52;
   nCurrentVars=53;
   V[53]=Sqrt(1-Pow(V[6],2));
   if(!isfinite(V[53]) || FError) return 53;
   nCurrentVars=54;
   V[54]=V[44]*V[53];
   if(!isfinite(V[54]) || FError) return 54;
   nCurrentVars=55;
   V[55]=2*V[6]*V[53];
   if(!isfinite(V[55]) || FError) return 55;
   nCurrentVars=56;
   V[56]=Pow(V[53],2)-Pow(V[6],2);
   if(!isfinite(V[56]) || FError) return 56;
   nCurrentVars=57;
   V[57]=V[46]/(Sqrt(1+Pow(V[46],2)));
   if(!isfinite(V[57]) || FError) return 57;
   nCurrentVars=58;
   V[58]=Sqrt(1-Pow(V[57],2));
   if(!isfinite(V[58]) || FError) return 58;
   nCurrentVars=59;
   V[59]=2*V[57]*V[58];
   if(!isfinite(V[59]) || FError) return 59;
   nCurrentVars=60;
   V[60]=Pow(V[58],2)-Pow(V[57],2);
   if(!isfinite(V[60]) || FError) return 60;
   nCurrentVars=61;
   V[61]=slhaVal("MASS",V[44],1,25);
   if(!isfinite(V[61]) || FError) return 61;
   nCurrentVars=62;
   V[62]=slhaVal("MASS",V[44],1,35);
   if(!isfinite(V[62]) || FError) return 62;
   nCurrentVars=63;
   V[63]=slhaVal("MASS",V[44],1,36);
   if(!isfinite(V[63]) || FError) return 63;
   nCurrentVars=64;
   V[64]=slhaVal("CVHMIX",V[44],2,1,1);
   if(!isfinite(V[64]) || FError) return 64;
   nCurrentVars=65;
   V[65]=slhaVal("CVHMIX",V[44],2,1,2);
   if(!isfinite(V[65]) || FError) return 65;
   nCurrentVars=66;
   V[66]=slhaVal("CVHMIX",V[44],2,1,4)*V[58]-slhaVal("CVHMIX",V[44],2,1,3)*V[57];
   if(!isfinite(V[66]) || FError) return 66;
   nCurrentVars=67;
   V[67]=slhaVal("CVHMIX",V[44],2,2,1);
   if(!isfinite(V[67]) || FError) return 67;
   nCurrentVars=68;
   V[68]=slhaVal("CVHMIX",V[44],2,2,2);
   if(!isfinite(V[68]) || FError) return 68;
   nCurrentVars=69;
   V[69]=slhaVal("CVHMIX",V[44],2,2,4)*V[58]-slhaVal("CVHMIX",V[44],2,2,3)*V[57];
   if(!isfinite(V[69]) || FError) return 69;
   nCurrentVars=70;
   V[70]=slhaVal("CVHMIX",V[44],2,3,1);
   if(!isfinite(V[70]) || FError) return 70;
   nCurrentVars=71;
   V[71]=slhaVal("CVHMIX",V[44],2,3,2);
   if(!isfinite(V[71]) || FError) return 71;
   nCurrentVars=72;
   V[72]=slhaVal("CVHMIX",V[44],2,3,4)*V[58]-slhaVal("CVHMIX",V[44],2,3,3)*V[57];
   if(!isfinite(V[72]) || FError) return 72;
   nCurrentVars=73;
   V[73]=slhaVal("HCOUPLINGS",V[44],1,1);
   if(!isfinite(V[73]) || FError) return 73;
   nCurrentVars=74;
   V[74]=slhaVal("HCOUPLINGS",V[44],1,2);
   if(!isfinite(V[74]) || FError) return 74;
   nCurrentVars=75;
   V[75]=slhaVal("HCOUPLINGS",V[44],1,3);
   if(!isfinite(V[75]) || FError) return 75;
   nCurrentVars=76;
   V[76]=slhaVal("HCOUPLINGS",V[44],1,4);
   if(!isfinite(V[76]) || FError) return 76;
   nCurrentVars=77;
   V[77]=slhaVal("HCOUPLINGS",V[44],1,5);
   if(!isfinite(V[77]) || FError) return 77;
   nCurrentVars=78;
   V[78]=slhaVal("HCOUPLINGS",V[44],1,6);
   if(!isfinite(V[78]) || FError) return 78;
   nCurrentVars=79;
   V[79]=slhaVal("HCOUPLINGS",V[44],1,7);
   if(!isfinite(V[79]) || FError) return 79;
   nCurrentVars=80;
   V[80]=slhaVal("IMHCOUPLINGS",V[44],1,5);
   if(!isfinite(V[80]) || FError) return 80;
   nCurrentVars=81;
   V[81]=slhaVal("IMHCOUPLINGS",V[44],1,6);
   if(!isfinite(V[81]) || FError) return 81;
   nCurrentVars=82;
   V[82]=slhaVal("IMHCOUPLINGS",V[44],1,7);
   if(!isfinite(V[82]) || FError) return 82;
   nCurrentVars=83;
   V[83]=slhaVal("MASS",V[44],1,1000022);
   if(!isfinite(V[83]) || FError) return 83;
   nCurrentVars=84;
   V[84]=slhaVal("MASS",V[44],1,1000023);
   if(!isfinite(V[84]) || FError) return 84;
   nCurrentVars=85;
   V[85]=slhaVal("MASS",V[44],1,1000025);
   if(!isfinite(V[85]) || FError) return 85;
   nCurrentVars=86;
   V[86]=slhaVal("MASS",V[44],1,1000035);
   if(!isfinite(V[86]) || FError) return 86;
   nCurrentVars=87;
   V[87]=slhaVal("NMIX",V[44],2,1,1);
   if(!isfinite(V[87]) || FError) return 87;
   nCurrentVars=88;
   V[88]=slhaVal("NMIX",V[44],2,1,2);
   if(!isfinite(V[88]) || FError) return 88;
   nCurrentVars=89;
   V[89]=slhaVal("NMIX",V[44],2,1,3);
   if(!isfinite(V[89]) || FError) return 89;
   nCurrentVars=90;
   V[90]=slhaVal("NMIX",V[44],2,1,4);
   if(!isfinite(V[90]) || FError) return 90;
   nCurrentVars=91;
   V[91]=slhaVal("NMIX",V[44],2,2,1);
   if(!isfinite(V[91]) || FError) return 91;
   nCurrentVars=92;
   V[92]=slhaVal("NMIX",V[44],2,2,2);
   if(!isfinite(V[92]) || FError) return 92;
   nCurrentVars=93;
   V[93]=slhaVal("NMIX",V[44],2,2,3);
   if(!isfinite(V[93]) || FError) return 93;
   nCurrentVars=94;
   V[94]=slhaVal("NMIX",V[44],2,2,4);
   if(!isfinite(V[94]) || FError) return 94;
   nCurrentVars=95;
   V[95]=slhaVal("NMIX",V[44],2,3,1);
   if(!isfinite(V[95]) || FError) return 95;
   nCurrentVars=96;
   V[96]=slhaVal("NMIX",V[44],2,3,2);
   if(!isfinite(V[96]) || FError) return 96;
   nCurrentVars=97;
   V[97]=slhaVal("NMIX",V[44],2,3,3);
   if(!isfinite(V[97]) || FError) return 97;
   nCurrentVars=98;
   V[98]=slhaVal("NMIX",V[44],2,3,4);
   if(!isfinite(V[98]) || FError) return 98;
   nCurrentVars=99;
   V[99]=slhaVal("NMIX",V[44],2,4,1);
   if(!isfinite(V[99]) || FError) return 99;
   nCurrentVars=100;
   V[100]=slhaVal("NMIX",V[44],2,4,2);
   if(!isfinite(V[100]) || FError) return 100;
   nCurrentVars=101;
   V[101]=slhaVal("NMIX",V[44],2,4,3);
   if(!isfinite(V[101]) || FError) return 101;
   nCurrentVars=102;
   V[102]=slhaVal("NMIX",V[44],2,4,4);
   if(!isfinite(V[102]) || FError) return 102;
   nCurrentVars=103;
   V[103]=slhaVal("IMNMIX",V[44],2,1,1);
   if(!isfinite(V[103]) || FError) return 103;
   nCurrentVars=104;
   V[104]=slhaVal("IMNMIX",V[44],2,1,2);
   if(!isfinite(V[104]) || FError) return 104;
   nCurrentVars=105;
   V[105]=slhaVal("IMNMIX",V[44],2,1,3);
   if(!isfinite(V[105]) || FError) return 105;
   nCurrentVars=106;
   V[106]=slhaVal("IMNMIX",V[44],2,1,4);
   if(!isfinite(V[106]) || FError) return 106;
   nCurrentVars=107;
   V[107]=slhaVal("IMNMIX",V[44],2,2,1);
   if(!isfinite(V[107]) || FError) return 107;
   nCurrentVars=108;
   V[108]=slhaVal("IMNMIX",V[44],2,2,2);
   if(!isfinite(V[108]) || FError) return 108;
   nCurrentVars=109;
   V[109]=slhaVal("IMNMIX",V[44],2,2,3);
   if(!isfinite(V[109]) || FError) return 109;
   nCurrentVars=110;
   V[110]=slhaVal("IMNMIX",V[44],2,2,4);
   if(!isfinite(V[110]) || FError) return 110;
   nCurrentVars=111;
   V[111]=slhaVal("IMNMIX",V[44],2,3,1);
   if(!isfinite(V[111]) || FError) return 111;
   nCurrentVars=112;
   V[112]=slhaVal("IMNMIX",V[44],2,3,2);
   if(!isfinite(V[112]) || FError) return 112;
   nCurrentVars=113;
   V[113]=slhaVal("IMNMIX",V[44],2,3,3);
   if(!isfinite(V[113]) || FError) return 113;
   nCurrentVars=114;
   V[114]=slhaVal("IMNMIX",V[44],2,3,4);
   if(!isfinite(V[114]) || FError) return 114;
   nCurrentVars=115;
   V[115]=slhaVal("IMNMIX",V[44],2,4,1);
   if(!isfinite(V[115]) || FError) return 115;
   nCurrentVars=116;
   V[116]=slhaVal("IMNMIX",V[44],2,4,2);
   if(!isfinite(V[116]) || FError) return 116;
   nCurrentVars=117;
   V[117]=slhaVal("IMNMIX",V[44],2,4,3);
   if(!isfinite(V[117]) || FError) return 117;
   nCurrentVars=118;
   V[118]=slhaVal("IMNMIX",V[44],2,4,4);
   if(!isfinite(V[118]) || FError) return 118;
   nCurrentVars=119;
   V[119]=slhaVal("MASS",V[44],1,1000024);
   if(!isfinite(V[119]) || FError) return 119;
   nCurrentVars=120;
   V[120]=slhaVal("MASS",V[44],1,1000037);
   if(!isfinite(V[120]) || FError) return 120;
   nCurrentVars=121;
   V[121]=slhaVal("VMIX",V[44],2,1,1);
   if(!isfinite(V[121]) || FError) return 121;
   nCurrentVars=122;
   V[122]=slhaVal("VMIX",V[44],2,1,2);
   if(!isfinite(V[122]) || FError) return 122;
   nCurrentVars=123;
   V[123]=slhaVal("VMIX",V[44],2,2,1);
   if(!isfinite(V[123]) || FError) return 123;
   nCurrentVars=124;
   V[124]=slhaVal("VMIX",V[44],2,2,2);
   if(!isfinite(V[124]) || FError) return 124;
   nCurrentVars=125;
   V[125]=slhaVal("IMVMIX",V[44],2,1,1);
   if(!isfinite(V[125]) || FError) return 125;
   nCurrentVars=126;
   V[126]=slhaVal("IMVMIX",V[44],2,1,2);
   if(!isfinite(V[126]) || FError) return 126;
   nCurrentVars=127;
   V[127]=slhaVal("IMVMIX",V[44],2,2,1);
   if(!isfinite(V[127]) || FError) return 127;
   nCurrentVars=128;
   V[128]=slhaVal("IMVMIX",V[44],2,2,2);
   if(!isfinite(V[128]) || FError) return 128;
   nCurrentVars=129;
   V[129]=slhaVal("UMIX",V[44],2,1,1);
   if(!isfinite(V[129]) || FError) return 129;
   nCurrentVars=130;
   V[130]=slhaVal("UMIX",V[44],2,1,2);
   if(!isfinite(V[130]) || FError) return 130;
   nCurrentVars=131;
   V[131]=slhaVal("UMIX",V[44],2,2,1);
   if(!isfinite(V[131]) || FError) return 131;
   nCurrentVars=132;
   V[132]=slhaVal("UMIX",V[44],2,2,2);
   if(!isfinite(V[132]) || FError) return 132;
   nCurrentVars=133;
   V[133]=slhaVal("IMUMIX",V[44],2,1,1);
   if(!isfinite(V[133]) || FError) return 133;
   nCurrentVars=134;
   V[134]=slhaVal("IMUMIX",V[44],2,1,2);
   if(!isfinite(V[134]) || FError) return 134;
   nCurrentVars=135;
   V[135]=slhaVal("IMUMIX",V[44],2,2,1);
   if(!isfinite(V[135]) || FError) return 135;
   nCurrentVars=136;
   V[136]=slhaVal("IMUMIX",V[44],2,2,2);
   if(!isfinite(V[136]) || FError) return 136;
   nCurrentVars=137;
   V[137]=slhaVal("THRESHOLD",V[44],2,6,1);
   if(!isfinite(V[137]) || FError) return 137;
   nCurrentVars=138;
   V[138]=slhaVal("THRESHOLD",V[44],2,6,2);
   if(!isfinite(V[138]) || FError) return 138;
   nCurrentVars=139;
   V[139]=slhaVal("THRESHOLD",V[44],2,5,1);
   if(!isfinite(V[139]) || FError) return 139;
   nCurrentVars=140;
   V[140]=slhaVal("THRESHOLD",V[44],2,5,2);
   if(!isfinite(V[140]) || FError) return 140;
   nCurrentVars=141;
   V[141]=slhaVal("MASS",V[44],1,1000015);
   if(!isfinite(V[141]) || FError) return 141;
   nCurrentVars=142;
   V[142]=slhaVal("MASS",V[44],1,2000015);
   if(!isfinite(V[142]) || FError) return 142;
   nCurrentVars=143;
   V[143]=slhaVal("STAUMIX",V[44],2,1,1);
   if(!isfinite(V[143]) || FError) return 143;
   nCurrentVars=144;
   V[144]=slhaVal("STAUMIX",V[44],2,2,1);
   if(!isfinite(V[144]) || FError) return 144;
   nCurrentVars=145;
   V[145]=slhaVal("IMSTAUMIX",V[44],2,2,1);
   if(!isfinite(V[145]) || FError) return 145;
   nCurrentVars=146;
   V[146]=slhaVal("MASS",V[44],1,1000006);
   if(!isfinite(V[146]) || FError) return 146;
   nCurrentVars=147;
   V[147]=slhaVal("MASS",V[44],1,2000006);
   if(!isfinite(V[147]) || FError) return 147;
   nCurrentVars=148;
   V[148]=slhaVal("STOPMIX",V[44],2,1,1);
   if(!isfinite(V[148]) || FError) return 148;
   nCurrentVars=149;
   V[149]=slhaVal("STOPMIX",V[44],2,2,1);
   if(!isfinite(V[149]) || FError) return 149;
   nCurrentVars=150;
   V[150]=slhaVal("IMSTOPMIX",V[44],2,2,1);
   if(!isfinite(V[150]) || FError) return 150;
   nCurrentVars=151;
   V[151]=slhaVal("MASS",V[44],1,1000005);
   if(!isfinite(V[151]) || FError) return 151;
   nCurrentVars=152;
   V[152]=slhaVal("MASS",V[44],1,2000005);
   if(!isfinite(V[152]) || FError) return 152;
   nCurrentVars=153;
   V[153]=slhaVal("SBOTMIX",V[44],2,1,1);
   if(!isfinite(V[153]) || FError) return 153;
   nCurrentVars=154;
   V[154]=slhaVal("SBOTMIX",V[44],2,2,1);
   if(!isfinite(V[154]) || FError) return 154;
   nCurrentVars=155;
   V[155]=slhaVal("IMSBOTMIX",V[44],2,2,1);
   if(!isfinite(V[155]) || FError) return 155;
   nCurrentVars=156;
   V[156]=3.14159*1;
   if(!isfinite(V[156]) || FError) return 156;
   nCurrentVars=157;
   V[157]=V[9]*Cos(V[10]*V[156]/(180));
   if(!isfinite(V[157]) || FError) return 157;
   nCurrentVars=158;
   V[158]=V[9]*Sin(V[10]*V[156]/(180));
   if(!isfinite(V[158]) || FError) return 158;
   nCurrentVars=159;
   V[159]=V[28]*Cos(V[29]*V[156]/(180));
   if(!isfinite(V[159]) || FError) return 159;
   nCurrentVars=160;
   V[160]=V[28]*Sin(V[29]*V[156]/(180));
   if(!isfinite(V[160]) || FError) return 160;
   nCurrentVars=161;
   V[161]=V[26]*Cos(V[27]*V[156]/(180));
   if(!isfinite(V[161]) || FError) return 161;
   nCurrentVars=162;
   V[162]=V[26]*Sin(V[27]*V[156]/(180));
   if(!isfinite(V[162]) || FError) return 162;
   nCurrentVars=163;
   V[163]=V[24]*Cos(V[25]*V[156]/(180));
   if(!isfinite(V[163]) || FError) return 163;
   nCurrentVars=164;
   V[164]=V[24]*Sin(V[25]*V[156]/(180));
   if(!isfinite(V[164]) || FError) return 164;
   nCurrentVars=165;
   V[165]=V[36]*Cos(V[37]*V[156]/(180));
   if(!isfinite(V[165]) || FError) return 165;
   nCurrentVars=166;
   V[166]=V[36]*Sin(V[37]*V[156]/(180));
   if(!isfinite(V[166]) || FError) return 166;
   nCurrentVars=167;
   V[167]=V[34]*Cos(V[35]*V[156]/(180));
   if(!isfinite(V[167]) || FError) return 167;
   nCurrentVars=168;
   V[168]=V[34]*Sin(V[35]*V[156]/(180));
   if(!isfinite(V[168]) || FError) return 168;
   nCurrentVars=169;
   V[169]=Sin(V[19]*V[156]/(180));
   if(!isfinite(V[169]) || FError) return 169;
   nCurrentVars=170;
   V[170]=Cos(V[19]*V[156]/(180));
   if(!isfinite(V[170]) || FError) return 170;
   nCurrentVars=171;
   V[171]=V[14]*Cos(V[17]*V[156]/(180));
   if(!isfinite(V[171]) || FError) return 171;
   nCurrentVars=172;
   V[172]=V[14]*Sin(V[17]*V[156]/(180));
   if(!isfinite(V[172]) || FError) return 172;
   nCurrentVars=173;
   V[173]=V[15]*Cos(V[18]*V[156]/(180));
   if(!isfinite(V[173]) || FError) return 173;
   nCurrentVars=174;
   V[174]=V[15]*Sin(V[18]*V[156]/(180));
   if(!isfinite(V[174]) || FError) return 174;
   nCurrentVars=175;
   V[175]=V[16]*Cos(V[19]*V[156]/(180));
   if(!isfinite(V[175]) || FError) return 175;
   nCurrentVars=176;
   V[176]=V[16]*Sin(V[19]*V[156]/(180));
   if(!isfinite(V[176]) || FError) return 176;
   nCurrentVars=177;
   V[177]=Sin(V[19]*V[156]/(360));
   if(!isfinite(V[177]) || FError) return 177;
   nCurrentVars=178;
   V[178]=Cos(V[19]*V[156]/(360));
   if(!isfinite(V[178]) || FError) return 178;
   nCurrentVars=179;
   V[179]=Sqrt(V[60]*Pow(V[54],2)/(2)+Pow(V[53],2)*Pow(V[20],2))/(V[53]);
   if(!isfinite(V[179]) || FError) return 179;
   nCurrentVars=180;
   V[180]=Sqrt(V[60]*Pow(V[54],2)/(2)+Pow(V[53],2)*Pow(V[20],2))/(V[53]);
   if(!isfinite(V[180]) || FError) return 180;
   nCurrentVars=181;
   V[181]=Sqrt(V[60]*Pow(V[54],2)/(2)+Pow(V[53],2)*Pow(V[21],2))/(V[53]);
   if(!isfinite(V[181]) || FError) return 181;
   nCurrentVars=182;
   V[182]=Sqrt(Pow(V[53],2)*Pow(V[20],2)-V[60]*(Pow(V[53],2)-Pow(V[6],2))*Pow(V[54],2)/(2))/(V[53]);
   if(!isfinite(V[182]) || FError) return 182;
   nCurrentVars=183;
   V[183]=Sqrt(Pow(V[53],2)*Pow(V[22],2)-Pow(V[6],2)*V[60]*Pow(V[54],2))/(V[53]);
   if(!isfinite(V[183]) || FError) return 183;
   nCurrentVars=184;
   V[184]=Sqrt(Pow(V[53],2)*Pow(V[20],2)-V[60]*(Pow(V[53],2)-Pow(V[6],2))*Pow(V[54],2)/(2))/(V[53]);
   if(!isfinite(V[184]) || FError) return 184;
   nCurrentVars=185;
   V[185]=Sqrt(Pow(V[53],2)*Pow(V[22],2)-Pow(V[6],2)*V[60]*Pow(V[54],2))/(V[53]);
   if(!isfinite(V[185]) || FError) return 185;
   nCurrentVars=186;
   V[186]=Sqrt(Pow(V[38],2)+Pow(V[54],2)/(Pow(V[53],2))*(1/(double)((2))-2/(double)((3))*Pow(V[6],2))*V[60]);
   if(!isfinite(V[186]) || FError) return 186;
   nCurrentVars=187;
   V[187]=Sqrt(Pow(V[38],2)+Pow(V[54],2)/(Pow(V[53],2))*(1/(double)((2))-2/(double)((3))*Pow(V[6],2))*V[60]);
   if(!isfinite(V[187]) || FError) return 187;
   nCurrentVars=188;
   V[188]=Sqrt(Pow(V[40],2)+Pow(V[54],2)/(Pow(V[53],2))*2/(double)((3))*Pow(V[6],2)*V[60]);
   if(!isfinite(V[188]) || FError) return 188;
   nCurrentVars=189;
   V[189]=Sqrt(Pow(V[40],2)+Pow(V[54],2)/(Pow(V[53],2))*2/(double)((3))*Pow(V[6],2)*V[60]);
   if(!isfinite(V[189]) || FError) return 189;
   nCurrentVars=190;
   V[190]=Sqrt(Pow(V[38],2)-Pow(V[54],2)/(Pow(V[53],2))*(1/(double)((2))-1/(double)((3))*Pow(V[6],2))*V[60]);
   if(!isfinite(V[190]) || FError) return 190;
   nCurrentVars=191;
   V[191]=Sqrt(Pow(V[38],2)-Pow(V[54],2)/(Pow(V[53],2))*(1/(double)((2))-1/(double)((3))*Pow(V[6],2))*V[60]);
   if(!isfinite(V[191]) || FError) return 191;
   nCurrentVars=192;
   V[192]=Sqrt(Pow(V[42],2)-1/(double)((3))*Pow(V[6],2)*Pow(V[54],2)/(Pow(V[53],2))*V[60]);
   if(!isfinite(V[192]) || FError) return 192;
   nCurrentVars=193;
   V[193]=Sqrt(Pow(V[42],2)-1/(double)((3))*Pow(V[6],2)*Pow(V[54],2)/(Pow(V[53],2))*V[60]);
   if(!isfinite(V[193]) || FError) return 193;
 FirstQ:
 cErr=1;
   nCurrentVars=194;
   V[194]=MtEff(V[4]);
   if(!isfinite(V[194]) || FError) return 194;
   nCurrentVars=195;
   V[195]=MbEff(V[4]);
   if(!isfinite(V[195]) || FError) return 195;
   nCurrentVars=196;
   V[196]=McEff(V[4]);
   if(!isfinite(V[196]) || FError) return 196;
   nCurrentVars=197;
   V[197]=bPoleMass();
   if(!isfinite(V[197]) || FError) return 197;
   nCurrentVars=198;
   V[198]=V[3]*(1+4/(double)((3))*alphaQCD(V[3])/(V[156]));
   if(!isfinite(V[198]) || FError) return 198;
   nCurrentVars=199;
   V[199]=4/(double)((9));
   if(!isfinite(V[199]) || FError) return 199;
   nCurrentVars=200;
   V[200]=1/(double)((9));
   if(!isfinite(V[200]) || FError) return 200;
   nCurrentVars=201;
   V[201]=-V[5]/(V[54])*V[196]/(V[6])*V[65]/(V[57])/(2)/(V[196]);
   if(!isfinite(V[201]) || FError) return 201;
   nCurrentVars=202;
   V[202]=-V[5]/(V[54])*V[195]/(V[6])/(V[59])*(V[57]*V[64]-V[57]*V[139]*V[64]+V[58]*V[139]*V[65]+V[140]*V[66])/(V[195]);
   if(!isfinite(V[202]) || FError) return 202;
   nCurrentVars=203;
   V[203]=-V[5]/(V[54])*V[194]/(V[6])*V[65]/(V[57])/(2)/(V[194]);
   if(!isfinite(V[203]) || FError) return 203;
   nCurrentVars=204;
   V[204]=-V[5]/(V[54])*V[7]/(V[6])*V[64]/(V[58])/(2)/(V[7]);
   if(!isfinite(V[204]) || FError) return 204;
   nCurrentVars=205;
   V[205]=-V[5]/(V[54])*V[196]/(V[6])*V[68]/(V[57])/(2)/(V[196]);
   if(!isfinite(V[205]) || FError) return 205;
   nCurrentVars=206;
   V[206]=-V[5]/(V[54])*V[195]/(V[6])/(V[59])*(V[57]*V[67]-V[57]*V[139]*V[67]+V[58]*V[139]*V[68]+V[140]*V[69])/(V[195]);
   if(!isfinite(V[206]) || FError) return 206;
   nCurrentVars=207;
   V[207]=-V[5]/(V[54])*V[194]/(V[6])*V[68]/(V[57])/(2)/(V[194]);
   if(!isfinite(V[207]) || FError) return 207;
   nCurrentVars=208;
   V[208]=-V[5]/(V[54])*V[7]/(V[6])*V[67]/(V[58])/(2)/(V[7]);
   if(!isfinite(V[208]) || FError) return 208;
   nCurrentVars=209;
   V[209]=-V[5]/(V[54])*V[196]/(V[6])*V[71]/(V[57])/(2)/(V[196]);
   if(!isfinite(V[209]) || FError) return 209;
   nCurrentVars=210;
   V[210]=-V[5]/(V[54])*V[195]/(V[6])/(V[59])*(V[57]*V[70]-V[57]*V[139]*V[70]+V[58]*V[139]*V[71]+V[140]*V[72])/(V[195]);
   if(!isfinite(V[210]) || FError) return 210;
   nCurrentVars=211;
   V[211]=-V[5]/(V[54])*V[194]/(V[6])*V[71]/(V[57])/(2)/(V[194]);
   if(!isfinite(V[211]) || FError) return 211;
   nCurrentVars=212;
   V[212]=-V[5]/(V[54])*V[7]/(V[6])*V[70]/(V[58])/(2)/(V[7]);
   if(!isfinite(V[212]) || FError) return 212;
   nCurrentVars=213;
   V[213]=V[56]/(Pow(V[53],2))*V[5]*V[54]/(V[6])*(V[58]*V[64]-V[57]*V[65])/(2)/(Pow(V[182],2));
   if(!isfinite(V[213]) || FError) return 213;
   nCurrentVars=214;
   V[214]=1/(Pow(V[53],2))*V[5]*V[54]*V[6]*(V[58]*V[64]-V[57]*V[65])/(Pow(V[183],2));
   if(!isfinite(V[214]) || FError) return 214;
   nCurrentVars=215;
   V[215]=V[56]/(Pow(V[53],2))*V[5]*V[54]/(V[6])*(V[58]*V[64]-V[57]*V[65])/(2)/(Pow(V[184],2));
   if(!isfinite(V[215]) || FError) return 215;
   nCurrentVars=216;
   V[216]=1/(Pow(V[53],2))*V[5]*V[54]*V[6]*(V[58]*V[64]-V[57]*V[65])/(Pow(V[185],2));
   if(!isfinite(V[216]) || FError) return 216;
   nCurrentVars=217;
   V[217]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[54],2)*V[64]*V[59]-3*Pow(V[54],2)*V[64]*V[59]-8*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+6*Pow(V[57],2)*Pow(V[54],2)*V[65]-12*Pow(V[53],2)*Pow(V[11],2)*V[65])/(12)/(Pow(V[186],2));
   if(!isfinite(V[217]) || FError) return 217;
   nCurrentVars=218;
   V[218]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(Pow(V[6],2)*Pow(V[54],2)*V[64]*V[59]-2*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+3*Pow(V[53],2)*Pow(V[11],2)*V[65])/(3)/(Pow(V[188],2));
   if(!isfinite(V[218]) || FError) return 218;
   nCurrentVars=219;
   V[219]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[54],2)*V[64]*V[59]-3*Pow(V[54],2)*V[64]*V[59]-8*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+6*Pow(V[57],2)*Pow(V[54],2)*V[65]-12*Pow(V[53],2)*Pow(V[196],2)*V[65])/(12)/(Pow(V[187],2));
   if(!isfinite(V[219]) || FError) return 219;
   nCurrentVars=220;
   V[220]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(Pow(V[6],2)*Pow(V[54],2)*V[64]*V[59]-2*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+3*Pow(V[53],2)*Pow(V[196],2)*V[65])/(3)/(Pow(V[189],2));
   if(!isfinite(V[220]) || FError) return 220;
   nCurrentVars=221;
   V[221]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-6*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[6],2)*Pow(V[54],2)*V[65]*V[59]+3*Pow(V[54],2)*V[65]*V[59]+12*Pow(V[53],2)*Pow(V[11],2)*V[64])/(12)/(Pow(V[190],2));
   if(!isfinite(V[221]) || FError) return 221;
   nCurrentVars=222;
   V[222]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(2*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-Pow(V[6],2)*Pow(V[54],2)*V[65]*V[59]-6*Pow(V[53],2)*Pow(V[11],2)*V[64])/(6)/(Pow(V[192],2));
   if(!isfinite(V[222]) || FError) return 222;
   nCurrentVars=223;
   V[223]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-6*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[6],2)*Pow(V[54],2)*V[65]*V[59]+3*Pow(V[54],2)*V[65]*V[59]+12*Pow(V[53],2)*Pow(V[11],2)*V[64])/(12)/(Pow(V[191],2));
   if(!isfinite(V[223]) || FError) return 223;
   nCurrentVars=224;
   V[224]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(2*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-Pow(V[6],2)*Pow(V[54],2)*V[65]*V[59]-6*Pow(V[53],2)*Pow(V[11],2)*V[64])/(6)/(Pow(V[193],2));
   if(!isfinite(V[224]) || FError) return 224;
   nCurrentVars=225;
   V[225]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]+2*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-4*Pow(V[6],2)*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]+4*Pow(V[6],2)*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[6],2)*Pow(V[144],2)*Pow(V[54],2)*V[65]*V[59]-Pow(V[143],2)*Pow(V[54],2)*V[65]*V[59]+2*Pow(V[6],2)*Pow(V[143],2)*Pow(V[54],2)*V[65]*V[59]-2*Pow(V[6],2)*Pow(V[145],2)*Pow(V[54],2)*V[65]*V[59]-4*Pow(V[53],2)*V[143]*V[144]*V[7]*V[65]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[58]*V[7]*V[66]*V[157]+4*Pow(V[53],2)*V[143]*V[144]*V[58]*V[7]*V[66]*V[158]+4*Pow(V[53],2)*V[143]*V[145]*V[7]*V[65]*V[158]-4*Pow(V[53],2)*Pow(V[7],2)*V[64]+4*Pow(V[53],2)*V[143]*V[144]*V[159]*V[7]*V[64]+4*Pow(V[53],2)*V[143]*V[145]*V[57]*V[159]*V[7]*V[66]-4*Pow(V[53],2)*V[143]*V[144]*V[57]*V[160]*V[7]*V[66]+4*Pow(V[53],2)*V[143]*V[145]*V[160]*V[7]*V[64])/(4)/(Pow(V[141],2));
   if(!isfinite(V[225]) || FError) return 225;
   nCurrentVars=226;
   V[226]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-4*Pow(V[6],2)*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]+4*Pow(V[6],2)*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[6],2)*Pow(V[144],2)*Pow(V[54],2)*V[65]*V[59]+Pow(V[145],2)*Pow(V[54],2)*V[65]*V[59]+Pow(V[144],2)*Pow(V[54],2)*V[65]*V[59]+2*Pow(V[6],2)*Pow(V[143],2)*Pow(V[54],2)*V[65]*V[59]-2*Pow(V[6],2)*Pow(V[145],2)*Pow(V[54],2)*V[65]*V[59]-4*Pow(V[53],2)*V[143]*V[144]*V[7]*V[65]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[58]*V[7]*V[66]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[7]*V[65]*V[158]+4*Pow(V[53],2)*V[143]*V[144]*V[58]*V[7]*V[66]*V[158]+4*Pow(V[53],2)*Pow(V[7],2)*V[64]+4*Pow(V[53],2)*V[143]*V[144]*V[159]*V[7]*V[64]+4*Pow(V[53],2)*V[143]*V[145]*V[57]*V[159]*V[7]*V[66]+4*Pow(V[53],2)*V[143]*V[145]*V[160]*V[7]*V[64]-4*Pow(V[53],2)*V[143]*V[144]*V[57]*V[160]*V[7]*V[66])/(4)/(Pow(V[142],2));
   if(!isfinite(V[226]) || FError) return 226;
   nCurrentVars=227;
   V[227]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[149],2)*Pow(V[54],2)*V[64]*V[59]+3*Pow(V[148],2)*Pow(V[54],2)*V[64]*V[59]-4*Pow(V[6],2)*Pow(V[148],2)*Pow(V[54],2)*V[64]*V[59]+4*Pow(V[6],2)*Pow(V[150],2)*Pow(V[54],2)*V[64]*V[59]-8*Pow(V[6],2)*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]-6*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+8*Pow(V[6],2)*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]-8*Pow(V[6],2)*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+12*Pow(V[53],2)*V[148]*V[149]*V[194]*V[64]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[57]*V[194]*V[66]*V[157]-12*Pow(V[53],2)*V[148]*V[149]*V[57]*V[194]*V[66]*V[158]-12*Pow(V[53],2)*V[148]*V[150]*V[194]*V[64]*V[158]+12*Pow(V[53],2)*Pow(V[194],2)*V[65]-12*Pow(V[53],2)*V[148]*V[149]*V[163]*V[194]*V[65]-12*Pow(V[53],2)*V[148]*V[150]*V[58]*V[163]*V[194]*V[66]+12*Pow(V[53],2)*V[148]*V[149]*V[58]*V[164]*V[194]*V[66]-12*Pow(V[53],2)*V[148]*V[150]*V[164]*V[194]*V[65])/(12)/(Pow(V[146],2));
   if(!isfinite(V[227]) || FError) return 227;
   nCurrentVars=228;
   V[228]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[149],2)*Pow(V[54],2)*V[64]*V[59]-3*Pow(V[150],2)*Pow(V[54],2)*V[64]*V[59]-3*Pow(V[149],2)*Pow(V[54],2)*V[64]*V[59]-4*Pow(V[6],2)*Pow(V[148],2)*Pow(V[54],2)*V[64]*V[59]+4*Pow(V[6],2)*Pow(V[150],2)*Pow(V[54],2)*V[64]*V[59]-8*Pow(V[6],2)*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+6*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+6*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+8*Pow(V[6],2)*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]-8*Pow(V[6],2)*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[65]+12*Pow(V[53],2)*V[148]*V[149]*V[194]*V[64]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[57]*V[194]*V[66]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[194]*V[64]*V[158]-12*Pow(V[53],2)*V[148]*V[149]*V[57]*V[194]*V[66]*V[158]-12*Pow(V[53],2)*Pow(V[194],2)*V[65]-12*Pow(V[53],2)*V[148]*V[149]*V[163]*V[194]*V[65]-12*Pow(V[53],2)*V[148]*V[150]*V[58]*V[163]*V[194]*V[66]-12*Pow(V[53],2)*V[148]*V[150]*V[164]*V[194]*V[65]+12*Pow(V[53],2)*V[148]*V[149]*V[58]*V[164]*V[194]*V[66])/(12)/(Pow(V[147],2));
   if(!isfinite(V[228]) || FError) return 228;
   nCurrentVars=229;
   V[229]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]+6*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-4*Pow(V[6],2)*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]+4*Pow(V[6],2)*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[6],2)*Pow(V[154],2)*Pow(V[54],2)*V[65]*V[59]-3*Pow(V[153],2)*Pow(V[54],2)*V[65]*V[59]+2*Pow(V[6],2)*Pow(V[153],2)*Pow(V[54],2)*V[65]*V[59]-2*Pow(V[6],2)*Pow(V[155],2)*Pow(V[54],2)*V[65]*V[59]-12*Pow(V[53],2)*V[153]*V[154]*V[195]*V[65]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[58]*V[195]*V[66]*V[157]+12*Pow(V[53],2)*V[153]*V[154]*V[58]*V[195]*V[66]*V[158]+12*Pow(V[53],2)*V[153]*V[155]*V[195]*V[65]*V[158]-12*Pow(V[53],2)*Pow(V[195],2)*V[64]+12*Pow(V[53],2)*V[153]*V[154]*V[161]*V[195]*V[64]+12*Pow(V[53],2)*V[153]*V[155]*V[57]*V[161]*V[195]*V[66]-12*Pow(V[53],2)*V[153]*V[154]*V[57]*V[162]*V[195]*V[66]+12*Pow(V[53],2)*V[153]*V[155]*V[162]*V[195]*V[64])/(12)/(Pow(V[151],2));
   if(!isfinite(V[229]) || FError) return 229;
   nCurrentVars=230;
   V[230]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-6*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-6*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-4*Pow(V[6],2)*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]+4*Pow(V[6],2)*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[64]-2*Pow(V[6],2)*Pow(V[154],2)*Pow(V[54],2)*V[65]*V[59]+3*Pow(V[155],2)*Pow(V[54],2)*V[65]*V[59]+3*Pow(V[154],2)*Pow(V[54],2)*V[65]*V[59]+2*Pow(V[6],2)*Pow(V[153],2)*Pow(V[54],2)*V[65]*V[59]-2*Pow(V[6],2)*Pow(V[155],2)*Pow(V[54],2)*V[65]*V[59]-12*Pow(V[53],2)*V[153]*V[154]*V[195]*V[65]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[58]*V[195]*V[66]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[195]*V[65]*V[158]+12*Pow(V[53],2)*V[153]*V[154]*V[58]*V[195]*V[66]*V[158]+12*Pow(V[53],2)*Pow(V[195],2)*V[64]+12*Pow(V[53],2)*V[153]*V[154]*V[161]*V[195]*V[64]+12*Pow(V[53],2)*V[153]*V[155]*V[57]*V[161]*V[195]*V[66]+12*Pow(V[53],2)*V[153]*V[155]*V[162]*V[195]*V[64]-12*Pow(V[53],2)*V[153]*V[154]*V[57]*V[162]*V[195]*V[66])/(12)/(Pow(V[152],2));
   if(!isfinite(V[230]) || FError) return 230;
   nCurrentVars=231;
   V[231]=V[56]/(Pow(V[53],2))*V[5]*V[54]/(V[6])*(V[58]*V[67]-V[57]*V[68])/(2)/(Pow(V[182],2));
   if(!isfinite(V[231]) || FError) return 231;
   nCurrentVars=232;
   V[232]=1/(Pow(V[53],2))*V[5]*V[54]*V[6]*(V[58]*V[67]-V[57]*V[68])/(Pow(V[183],2));
   if(!isfinite(V[232]) || FError) return 232;
   nCurrentVars=233;
   V[233]=V[56]/(Pow(V[53],2))*V[5]*V[54]/(V[6])*(V[58]*V[67]-V[57]*V[68])/(2)/(Pow(V[184],2));
   if(!isfinite(V[233]) || FError) return 233;
   nCurrentVars=234;
   V[234]=1/(Pow(V[53],2))*V[5]*V[54]*V[6]*(V[58]*V[67]-V[57]*V[68])/(Pow(V[185],2));
   if(!isfinite(V[234]) || FError) return 234;
   nCurrentVars=235;
   V[235]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[54],2)*V[67]*V[59]-3*Pow(V[54],2)*V[67]*V[59]-8*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+6*Pow(V[57],2)*Pow(V[54],2)*V[68]-12*Pow(V[53],2)*Pow(V[11],2)*V[68])/(12)/(Pow(V[186],2));
   if(!isfinite(V[235]) || FError) return 235;
   nCurrentVars=236;
   V[236]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(Pow(V[6],2)*Pow(V[54],2)*V[67]*V[59]-2*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+3*Pow(V[53],2)*Pow(V[11],2)*V[68])/(3)/(Pow(V[188],2));
   if(!isfinite(V[236]) || FError) return 236;
   nCurrentVars=237;
   V[237]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[54],2)*V[67]*V[59]-3*Pow(V[54],2)*V[67]*V[59]-8*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+6*Pow(V[57],2)*Pow(V[54],2)*V[68]-12*Pow(V[53],2)*Pow(V[196],2)*V[68])/(12)/(Pow(V[187],2));
   if(!isfinite(V[237]) || FError) return 237;
   nCurrentVars=238;
   V[238]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(Pow(V[6],2)*Pow(V[54],2)*V[67]*V[59]-2*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+3*Pow(V[53],2)*Pow(V[196],2)*V[68])/(3)/(Pow(V[189],2));
   if(!isfinite(V[238]) || FError) return 238;
   nCurrentVars=239;
   V[239]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-6*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[6],2)*Pow(V[54],2)*V[68]*V[59]+3*Pow(V[54],2)*V[68]*V[59]+12*Pow(V[53],2)*Pow(V[11],2)*V[67])/(12)/(Pow(V[190],2));
   if(!isfinite(V[239]) || FError) return 239;
   nCurrentVars=240;
   V[240]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(2*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-Pow(V[6],2)*Pow(V[54],2)*V[68]*V[59]-6*Pow(V[53],2)*Pow(V[11],2)*V[67])/(6)/(Pow(V[192],2));
   if(!isfinite(V[240]) || FError) return 240;
   nCurrentVars=241;
   V[241]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-6*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[6],2)*Pow(V[54],2)*V[68]*V[59]+3*Pow(V[54],2)*V[68]*V[59]+12*Pow(V[53],2)*Pow(V[11],2)*V[67])/(12)/(Pow(V[191],2));
   if(!isfinite(V[241]) || FError) return 241;
   nCurrentVars=242;
   V[242]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(2*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-Pow(V[6],2)*Pow(V[54],2)*V[68]*V[59]-6*Pow(V[53],2)*Pow(V[11],2)*V[67])/(6)/(Pow(V[193],2));
   if(!isfinite(V[242]) || FError) return 242;
   nCurrentVars=243;
   V[243]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]+2*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-4*Pow(V[6],2)*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]+4*Pow(V[6],2)*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[6],2)*Pow(V[144],2)*Pow(V[54],2)*V[68]*V[59]-Pow(V[143],2)*Pow(V[54],2)*V[68]*V[59]+2*Pow(V[6],2)*Pow(V[143],2)*Pow(V[54],2)*V[68]*V[59]-2*Pow(V[6],2)*Pow(V[145],2)*Pow(V[54],2)*V[68]*V[59]-4*Pow(V[53],2)*V[143]*V[144]*V[7]*V[68]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[58]*V[7]*V[69]*V[157]+4*Pow(V[53],2)*V[143]*V[144]*V[58]*V[7]*V[69]*V[158]+4*Pow(V[53],2)*V[143]*V[145]*V[7]*V[68]*V[158]-4*Pow(V[53],2)*Pow(V[7],2)*V[67]+4*Pow(V[53],2)*V[143]*V[144]*V[159]*V[7]*V[67]+4*Pow(V[53],2)*V[143]*V[145]*V[57]*V[159]*V[7]*V[69]-4*Pow(V[53],2)*V[143]*V[144]*V[57]*V[160]*V[7]*V[69]+4*Pow(V[53],2)*V[143]*V[145]*V[160]*V[7]*V[67])/(4)/(Pow(V[141],2));
   if(!isfinite(V[243]) || FError) return 243;
   nCurrentVars=244;
   V[244]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-4*Pow(V[6],2)*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]+4*Pow(V[6],2)*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[6],2)*Pow(V[144],2)*Pow(V[54],2)*V[68]*V[59]+Pow(V[145],2)*Pow(V[54],2)*V[68]*V[59]+Pow(V[144],2)*Pow(V[54],2)*V[68]*V[59]+2*Pow(V[6],2)*Pow(V[143],2)*Pow(V[54],2)*V[68]*V[59]-2*Pow(V[6],2)*Pow(V[145],2)*Pow(V[54],2)*V[68]*V[59]-4*Pow(V[53],2)*V[143]*V[144]*V[7]*V[68]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[58]*V[7]*V[69]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[7]*V[68]*V[158]+4*Pow(V[53],2)*V[143]*V[144]*V[58]*V[7]*V[69]*V[158]+4*Pow(V[53],2)*Pow(V[7],2)*V[67]+4*Pow(V[53],2)*V[143]*V[144]*V[159]*V[7]*V[67]+4*Pow(V[53],2)*V[143]*V[145]*V[57]*V[159]*V[7]*V[69]+4*Pow(V[53],2)*V[143]*V[145]*V[160]*V[7]*V[67]-4*Pow(V[53],2)*V[143]*V[144]*V[57]*V[160]*V[7]*V[69])/(4)/(Pow(V[142],2));
   if(!isfinite(V[244]) || FError) return 244;
   nCurrentVars=245;
   V[245]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[149],2)*Pow(V[54],2)*V[67]*V[59]+3*Pow(V[148],2)*Pow(V[54],2)*V[67]*V[59]-4*Pow(V[6],2)*Pow(V[148],2)*Pow(V[54],2)*V[67]*V[59]+4*Pow(V[6],2)*Pow(V[150],2)*Pow(V[54],2)*V[67]*V[59]-8*Pow(V[6],2)*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]-6*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+8*Pow(V[6],2)*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]-8*Pow(V[6],2)*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+12*Pow(V[53],2)*V[148]*V[149]*V[194]*V[67]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[57]*V[194]*V[69]*V[157]-12*Pow(V[53],2)*V[148]*V[149]*V[57]*V[194]*V[69]*V[158]-12*Pow(V[53],2)*V[148]*V[150]*V[194]*V[67]*V[158]+12*Pow(V[53],2)*Pow(V[194],2)*V[68]-12*Pow(V[53],2)*V[148]*V[149]*V[163]*V[194]*V[68]-12*Pow(V[53],2)*V[148]*V[150]*V[58]*V[163]*V[194]*V[69]+12*Pow(V[53],2)*V[148]*V[149]*V[58]*V[164]*V[194]*V[69]-12*Pow(V[53],2)*V[148]*V[150]*V[164]*V[194]*V[68])/(12)/(Pow(V[146],2));
   if(!isfinite(V[245]) || FError) return 245;
   nCurrentVars=246;
   V[246]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[149],2)*Pow(V[54],2)*V[67]*V[59]-3*Pow(V[150],2)*Pow(V[54],2)*V[67]*V[59]-3*Pow(V[149],2)*Pow(V[54],2)*V[67]*V[59]-4*Pow(V[6],2)*Pow(V[148],2)*Pow(V[54],2)*V[67]*V[59]+4*Pow(V[6],2)*Pow(V[150],2)*Pow(V[54],2)*V[67]*V[59]-8*Pow(V[6],2)*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+6*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+6*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+8*Pow(V[6],2)*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]-8*Pow(V[6],2)*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[68]+12*Pow(V[53],2)*V[148]*V[149]*V[194]*V[67]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[57]*V[194]*V[69]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[194]*V[67]*V[158]-12*Pow(V[53],2)*V[148]*V[149]*V[57]*V[194]*V[69]*V[158]-12*Pow(V[53],2)*Pow(V[194],2)*V[68]-12*Pow(V[53],2)*V[148]*V[149]*V[163]*V[194]*V[68]-12*Pow(V[53],2)*V[148]*V[150]*V[58]*V[163]*V[194]*V[69]-12*Pow(V[53],2)*V[148]*V[150]*V[164]*V[194]*V[68]+12*Pow(V[53],2)*V[148]*V[149]*V[58]*V[164]*V[194]*V[69])/(12)/(Pow(V[147],2));
   if(!isfinite(V[246]) || FError) return 246;
   nCurrentVars=247;
   V[247]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]+6*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-4*Pow(V[6],2)*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]+4*Pow(V[6],2)*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[6],2)*Pow(V[154],2)*Pow(V[54],2)*V[68]*V[59]-3*Pow(V[153],2)*Pow(V[54],2)*V[68]*V[59]+2*Pow(V[6],2)*Pow(V[153],2)*Pow(V[54],2)*V[68]*V[59]-2*Pow(V[6],2)*Pow(V[155],2)*Pow(V[54],2)*V[68]*V[59]-12*Pow(V[53],2)*V[153]*V[154]*V[195]*V[68]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[58]*V[195]*V[69]*V[157]+12*Pow(V[53],2)*V[153]*V[154]*V[58]*V[195]*V[69]*V[158]+12*Pow(V[53],2)*V[153]*V[155]*V[195]*V[68]*V[158]-12*Pow(V[53],2)*Pow(V[195],2)*V[67]+12*Pow(V[53],2)*V[153]*V[154]*V[161]*V[195]*V[67]+12*Pow(V[53],2)*V[153]*V[155]*V[57]*V[161]*V[195]*V[69]-12*Pow(V[53],2)*V[153]*V[154]*V[57]*V[162]*V[195]*V[69]+12*Pow(V[53],2)*V[153]*V[155]*V[162]*V[195]*V[67])/(12)/(Pow(V[151],2));
   if(!isfinite(V[247]) || FError) return 247;
   nCurrentVars=248;
   V[248]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-6*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-6*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-4*Pow(V[6],2)*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]+4*Pow(V[6],2)*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[67]-2*Pow(V[6],2)*Pow(V[154],2)*Pow(V[54],2)*V[68]*V[59]+3*Pow(V[155],2)*Pow(V[54],2)*V[68]*V[59]+3*Pow(V[154],2)*Pow(V[54],2)*V[68]*V[59]+2*Pow(V[6],2)*Pow(V[153],2)*Pow(V[54],2)*V[68]*V[59]-2*Pow(V[6],2)*Pow(V[155],2)*Pow(V[54],2)*V[68]*V[59]-12*Pow(V[53],2)*V[153]*V[154]*V[195]*V[68]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[58]*V[195]*V[69]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[195]*V[68]*V[158]+12*Pow(V[53],2)*V[153]*V[154]*V[58]*V[195]*V[69]*V[158]+12*Pow(V[53],2)*Pow(V[195],2)*V[67]+12*Pow(V[53],2)*V[153]*V[154]*V[161]*V[195]*V[67]+12*Pow(V[53],2)*V[153]*V[155]*V[57]*V[161]*V[195]*V[69]+12*Pow(V[53],2)*V[153]*V[155]*V[162]*V[195]*V[67]-12*Pow(V[53],2)*V[153]*V[154]*V[57]*V[162]*V[195]*V[69])/(12)/(Pow(V[152],2));
   if(!isfinite(V[248]) || FError) return 248;
   nCurrentVars=249;
   V[249]=V[56]/(Pow(V[53],2))*V[5]*V[54]/(V[6])*(V[58]*V[70]-V[57]*V[71])/(2)/(Pow(V[182],2));
   if(!isfinite(V[249]) || FError) return 249;
   nCurrentVars=250;
   V[250]=1/(Pow(V[53],2))*V[5]*V[54]*V[6]*(V[58]*V[70]-V[57]*V[71])/(Pow(V[183],2));
   if(!isfinite(V[250]) || FError) return 250;
   nCurrentVars=251;
   V[251]=V[56]/(Pow(V[53],2))*V[5]*V[54]/(V[6])*(V[58]*V[70]-V[57]*V[71])/(2)/(Pow(V[184],2));
   if(!isfinite(V[251]) || FError) return 251;
   nCurrentVars=252;
   V[252]=1/(Pow(V[53],2))*V[5]*V[54]*V[6]*(V[58]*V[70]-V[57]*V[71])/(Pow(V[185],2));
   if(!isfinite(V[252]) || FError) return 252;
   nCurrentVars=253;
   V[253]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[54],2)*V[70]*V[59]-3*Pow(V[54],2)*V[70]*V[59]-8*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+6*Pow(V[57],2)*Pow(V[54],2)*V[71]-12*Pow(V[53],2)*Pow(V[11],2)*V[71])/(12)/(Pow(V[186],2));
   if(!isfinite(V[253]) || FError) return 253;
   nCurrentVars=254;
   V[254]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(Pow(V[6],2)*Pow(V[54],2)*V[70]*V[59]-2*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+3*Pow(V[53],2)*Pow(V[11],2)*V[71])/(3)/(Pow(V[188],2));
   if(!isfinite(V[254]) || FError) return 254;
   nCurrentVars=255;
   V[255]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[54],2)*V[70]*V[59]-3*Pow(V[54],2)*V[70]*V[59]-8*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+6*Pow(V[57],2)*Pow(V[54],2)*V[71]-12*Pow(V[53],2)*Pow(V[196],2)*V[71])/(12)/(Pow(V[187],2));
   if(!isfinite(V[255]) || FError) return 255;
   nCurrentVars=256;
   V[256]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(Pow(V[6],2)*Pow(V[54],2)*V[70]*V[59]-2*Pow(V[6],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+3*Pow(V[53],2)*Pow(V[196],2)*V[71])/(3)/(Pow(V[189],2));
   if(!isfinite(V[256]) || FError) return 256;
   nCurrentVars=257;
   V[257]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-6*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[6],2)*Pow(V[54],2)*V[71]*V[59]+3*Pow(V[54],2)*V[71]*V[59]+12*Pow(V[53],2)*Pow(V[11],2)*V[70])/(12)/(Pow(V[190],2));
   if(!isfinite(V[257]) || FError) return 257;
   nCurrentVars=258;
   V[258]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(2*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-Pow(V[6],2)*Pow(V[54],2)*V[71]*V[59]-6*Pow(V[53],2)*Pow(V[11],2)*V[70])/(6)/(Pow(V[192],2));
   if(!isfinite(V[258]) || FError) return 258;
   nCurrentVars=259;
   V[259]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-6*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[6],2)*Pow(V[54],2)*V[71]*V[59]+3*Pow(V[54],2)*V[71]*V[59]+12*Pow(V[53],2)*Pow(V[11],2)*V[70])/(12)/(Pow(V[191],2));
   if(!isfinite(V[259]) || FError) return 259;
   nCurrentVars=260;
   V[260]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(2*Pow(V[6],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-Pow(V[6],2)*Pow(V[54],2)*V[71]*V[59]-6*Pow(V[53],2)*Pow(V[11],2)*V[70])/(6)/(Pow(V[193],2));
   if(!isfinite(V[260]) || FError) return 260;
   nCurrentVars=261;
   V[261]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]+2*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-4*Pow(V[6],2)*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]+4*Pow(V[6],2)*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[6],2)*Pow(V[144],2)*Pow(V[54],2)*V[71]*V[59]-Pow(V[143],2)*Pow(V[54],2)*V[71]*V[59]+2*Pow(V[6],2)*Pow(V[143],2)*Pow(V[54],2)*V[71]*V[59]-2*Pow(V[6],2)*Pow(V[145],2)*Pow(V[54],2)*V[71]*V[59]-4*Pow(V[53],2)*V[143]*V[144]*V[7]*V[71]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[58]*V[7]*V[72]*V[157]+4*Pow(V[53],2)*V[143]*V[144]*V[58]*V[7]*V[72]*V[158]+4*Pow(V[53],2)*V[143]*V[145]*V[7]*V[71]*V[158]-4*Pow(V[53],2)*Pow(V[7],2)*V[70]+4*Pow(V[53],2)*V[143]*V[144]*V[159]*V[7]*V[70]+4*Pow(V[53],2)*V[143]*V[145]*V[57]*V[159]*V[7]*V[72]-4*Pow(V[53],2)*V[143]*V[144]*V[57]*V[160]*V[7]*V[72]+4*Pow(V[53],2)*V[143]*V[145]*V[160]*V[7]*V[70])/(4)/(Pow(V[141],2));
   if(!isfinite(V[261]) || FError) return 261;
   nCurrentVars=262;
   V[262]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[144],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-4*Pow(V[6],2)*Pow(V[143],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]+4*Pow(V[6],2)*Pow(V[145],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[6],2)*Pow(V[144],2)*Pow(V[54],2)*V[71]*V[59]+Pow(V[145],2)*Pow(V[54],2)*V[71]*V[59]+Pow(V[144],2)*Pow(V[54],2)*V[71]*V[59]+2*Pow(V[6],2)*Pow(V[143],2)*Pow(V[54],2)*V[71]*V[59]-2*Pow(V[6],2)*Pow(V[145],2)*Pow(V[54],2)*V[71]*V[59]-4*Pow(V[53],2)*V[143]*V[144]*V[7]*V[71]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[58]*V[7]*V[72]*V[157]+4*Pow(V[53],2)*V[143]*V[145]*V[7]*V[71]*V[158]+4*Pow(V[53],2)*V[143]*V[144]*V[58]*V[7]*V[72]*V[158]+4*Pow(V[53],2)*Pow(V[7],2)*V[70]+4*Pow(V[53],2)*V[143]*V[144]*V[159]*V[7]*V[70]+4*Pow(V[53],2)*V[143]*V[145]*V[57]*V[159]*V[7]*V[72]+4*Pow(V[53],2)*V[143]*V[145]*V[160]*V[7]*V[70]-4*Pow(V[53],2)*V[143]*V[144]*V[57]*V[160]*V[7]*V[72])/(4)/(Pow(V[142],2));
   if(!isfinite(V[262]) || FError) return 262;
   nCurrentVars=263;
   V[263]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[149],2)*Pow(V[54],2)*V[70]*V[59]+3*Pow(V[148],2)*Pow(V[54],2)*V[70]*V[59]-4*Pow(V[6],2)*Pow(V[148],2)*Pow(V[54],2)*V[70]*V[59]+4*Pow(V[6],2)*Pow(V[150],2)*Pow(V[54],2)*V[70]*V[59]-8*Pow(V[6],2)*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]-6*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+8*Pow(V[6],2)*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]-8*Pow(V[6],2)*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+12*Pow(V[53],2)*V[148]*V[149]*V[194]*V[70]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[57]*V[194]*V[72]*V[157]-12*Pow(V[53],2)*V[148]*V[149]*V[57]*V[194]*V[72]*V[158]-12*Pow(V[53],2)*V[148]*V[150]*V[194]*V[70]*V[158]+12*Pow(V[53],2)*Pow(V[194],2)*V[71]-12*Pow(V[53],2)*V[148]*V[149]*V[163]*V[194]*V[71]-12*Pow(V[53],2)*V[148]*V[150]*V[58]*V[163]*V[194]*V[72]+12*Pow(V[53],2)*V[148]*V[149]*V[58]*V[164]*V[194]*V[72]-12*Pow(V[53],2)*V[148]*V[150]*V[164]*V[194]*V[71])/(12)/(Pow(V[146],2));
   if(!isfinite(V[263]) || FError) return 263;
   nCurrentVars=264;
   V[264]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[57])*(4*Pow(V[6],2)*Pow(V[149],2)*Pow(V[54],2)*V[70]*V[59]-3*Pow(V[150],2)*Pow(V[54],2)*V[70]*V[59]-3*Pow(V[149],2)*Pow(V[54],2)*V[70]*V[59]-4*Pow(V[6],2)*Pow(V[148],2)*Pow(V[54],2)*V[70]*V[59]+4*Pow(V[6],2)*Pow(V[150],2)*Pow(V[54],2)*V[70]*V[59]-8*Pow(V[6],2)*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+6*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+6*Pow(V[149],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+8*Pow(V[6],2)*Pow(V[148],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]-8*Pow(V[6],2)*Pow(V[150],2)*Pow(V[57],2)*Pow(V[54],2)*V[71]+12*Pow(V[53],2)*V[148]*V[149]*V[194]*V[70]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[57]*V[194]*V[72]*V[157]-12*Pow(V[53],2)*V[148]*V[150]*V[194]*V[70]*V[158]-12*Pow(V[53],2)*V[148]*V[149]*V[57]*V[194]*V[72]*V[158]-12*Pow(V[53],2)*Pow(V[194],2)*V[71]-12*Pow(V[53],2)*V[148]*V[149]*V[163]*V[194]*V[71]-12*Pow(V[53],2)*V[148]*V[150]*V[58]*V[163]*V[194]*V[72]-12*Pow(V[53],2)*V[148]*V[150]*V[164]*V[194]*V[71]+12*Pow(V[53],2)*V[148]*V[149]*V[58]*V[164]*V[194]*V[72])/(12)/(Pow(V[147],2));
   if(!isfinite(V[264]) || FError) return 264;
   nCurrentVars=265;
   V[265]=1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]+6*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-4*Pow(V[6],2)*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]+4*Pow(V[6],2)*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[6],2)*Pow(V[154],2)*Pow(V[54],2)*V[71]*V[59]-3*Pow(V[153],2)*Pow(V[54],2)*V[71]*V[59]+2*Pow(V[6],2)*Pow(V[153],2)*Pow(V[54],2)*V[71]*V[59]-2*Pow(V[6],2)*Pow(V[155],2)*Pow(V[54],2)*V[71]*V[59]-12*Pow(V[53],2)*V[153]*V[154]*V[195]*V[71]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[58]*V[195]*V[72]*V[157]+12*Pow(V[53],2)*V[153]*V[154]*V[58]*V[195]*V[72]*V[158]+12*Pow(V[53],2)*V[153]*V[155]*V[195]*V[71]*V[158]-12*Pow(V[53],2)*Pow(V[195],2)*V[70]+12*Pow(V[53],2)*V[153]*V[154]*V[161]*V[195]*V[70]+12*Pow(V[53],2)*V[153]*V[155]*V[57]*V[161]*V[195]*V[72]-12*Pow(V[53],2)*V[153]*V[154]*V[57]*V[162]*V[195]*V[72]+12*Pow(V[53],2)*V[153]*V[155]*V[162]*V[195]*V[70])/(12)/(Pow(V[151],2));
   if(!isfinite(V[265]) || FError) return 265;
   nCurrentVars=266;
   V[266]=-1/(Pow(V[53],2))*V[5]/(V[54])/(V[6])/(V[58])*(4*Pow(V[6],2)*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-6*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-6*Pow(V[154],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-4*Pow(V[6],2)*Pow(V[153],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]+4*Pow(V[6],2)*Pow(V[155],2)*Pow(V[58],2)*Pow(V[54],2)*V[70]-2*Pow(V[6],2)*Pow(V[154],2)*Pow(V[54],2)*V[71]*V[59]+3*Pow(V[155],2)*Pow(V[54],2)*V[71]*V[59]+3*Pow(V[154],2)*Pow(V[54],2)*V[71]*V[59]+2*Pow(V[6],2)*Pow(V[153],2)*Pow(V[54],2)*V[71]*V[59]-2*Pow(V[6],2)*Pow(V[155],2)*Pow(V[54],2)*V[71]*V[59]-12*Pow(V[53],2)*V[153]*V[154]*V[195]*V[71]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[58]*V[195]*V[72]*V[157]+12*Pow(V[53],2)*V[153]*V[155]*V[195]*V[71]*V[158]+12*Pow(V[53],2)*V[153]*V[154]*V[58]*V[195]*V[72]*V[158]+12*Pow(V[53],2)*Pow(V[195],2)*V[70]+12*Pow(V[53],2)*V[153]*V[154]*V[161]*V[195]*V[70]+12*Pow(V[53],2)*V[153]*V[155]*V[57]*V[161]*V[195]*V[72]+12*Pow(V[53],2)*V[153]*V[155]*V[162]*V[195]*V[70]-12*Pow(V[53],2)*V[153]*V[154]*V[57]*V[162]*V[195]*V[72])/(12)/(Pow(V[152],2));
   if(!isfinite(V[266]) || FError) return 266;
   nCurrentVars=267;
   V[267]=V[5]*V[54]/(V[6])*(V[58]*V[64]+V[57]*V[65])/(Pow(V[54],2));
   if(!isfinite(V[267]) || FError) return 267;
   nCurrentVars=268;
   V[268]=2*V[58]*Pow(V[57],2)*V[64]*V[73]+2*Pow(V[58],2)*V[57]*V[65]*V[74]+Pow(V[58],3)*V[64]*V[75]+Pow(V[57],3)*V[65]*V[75]-Pow(V[58],2)*V[57]*V[65]*V[76];
   if(!isfinite(V[268]) || FError) return 268;
   nCurrentVars=269;
   V[269]=V[268]-V[58]*Pow(V[57],2)*V[64]*V[76]-2*Pow(V[58],2)*V[57]*V[65]*V[77]-2*V[58]*Pow(V[57],2)*V[64]*V[77]+V[66]*V[80]*V[59]+3*Pow(V[57],3)*V[64]*V[78];
   if(!isfinite(V[269]) || FError) return 269;
   nCurrentVars=270;
   V[270]=V[269]-2*V[57]*V[64]*V[78]+V[58]*Pow(V[57],2)*V[65]*V[78]-Pow(V[57],2)*V[66]*V[81]-3*V[58]*Pow(V[57],2)*V[65]*V[79]+V[58]*V[65]*V[79];
   if(!isfinite(V[270]) || FError) return 270;
   nCurrentVars=271;
   V[271]=V[270]+Pow(V[58],2)*V[57]*V[64]*V[79]-Pow(V[58],2)*V[66]*V[82];
   if(!isfinite(V[271]) || FError) return 271;
   nCurrentVars=272;
   V[272]=2/(V[5])*V[54]*V[6]*V[271]/(Pow(V[8],2));
   if(!isfinite(V[272]) || FError) return 272;
   nCurrentVars=273;
   V[273]=V[5]*V[54]/(V[6])*(V[58]*V[67]+V[57]*V[68])/(Pow(V[54],2));
   if(!isfinite(V[273]) || FError) return 273;
   nCurrentVars=274;
   V[274]=2*V[58]*Pow(V[57],2)*V[67]*V[73]+2*Pow(V[58],2)*V[57]*V[68]*V[74]+Pow(V[58],3)*V[67]*V[75]+Pow(V[57],3)*V[68]*V[75]-Pow(V[58],2)*V[57]*V[68]*V[76];
   if(!isfinite(V[274]) || FError) return 274;
   nCurrentVars=275;
   V[275]=V[274]-V[58]*Pow(V[57],2)*V[67]*V[76]-2*Pow(V[58],2)*V[57]*V[68]*V[77]-2*V[58]*Pow(V[57],2)*V[67]*V[77]+V[69]*V[80]*V[59]+3*Pow(V[57],3)*V[67]*V[78];
   if(!isfinite(V[275]) || FError) return 275;
   nCurrentVars=276;
   V[276]=V[275]-2*V[57]*V[67]*V[78]+V[58]*Pow(V[57],2)*V[68]*V[78]-Pow(V[57],2)*V[69]*V[81]-3*V[58]*Pow(V[57],2)*V[68]*V[79]+V[58]*V[68]*V[79];
   if(!isfinite(V[276]) || FError) return 276;
   nCurrentVars=277;
   V[277]=V[276]+Pow(V[58],2)*V[57]*V[67]*V[79]-Pow(V[58],2)*V[69]*V[82];
   if(!isfinite(V[277]) || FError) return 277;
   nCurrentVars=278;
   V[278]=2/(V[5])*V[54]*V[6]*V[277]/(Pow(V[8],2));
   if(!isfinite(V[278]) || FError) return 278;
   nCurrentVars=279;
   V[279]=V[5]*V[54]/(V[6])*(V[58]*V[70]+V[57]*V[71])/(Pow(V[54],2));
   if(!isfinite(V[279]) || FError) return 279;
   nCurrentVars=280;
   V[280]=2*V[58]*Pow(V[57],2)*V[70]*V[73]+2*Pow(V[58],2)*V[57]*V[71]*V[74]+Pow(V[58],3)*V[70]*V[75]+Pow(V[57],3)*V[71]*V[75]-Pow(V[58],2)*V[57]*V[71]*V[76];
   if(!isfinite(V[280]) || FError) return 280;
   nCurrentVars=281;
   V[281]=V[280]-V[58]*Pow(V[57],2)*V[70]*V[76]-2*Pow(V[58],2)*V[57]*V[71]*V[77]-2*V[58]*Pow(V[57],2)*V[70]*V[77]+V[72]*V[80]*V[59]+3*Pow(V[57],3)*V[70]*V[78];
   if(!isfinite(V[281]) || FError) return 281;
   nCurrentVars=282;
   V[282]=V[281]-2*V[57]*V[70]*V[78]+V[58]*Pow(V[57],2)*V[71]*V[78]-Pow(V[57],2)*V[72]*V[81]-3*V[58]*Pow(V[57],2)*V[71]*V[79]+V[58]*V[71]*V[79];
   if(!isfinite(V[282]) || FError) return 282;
   nCurrentVars=283;
   V[283]=V[282]+Pow(V[58],2)*V[57]*V[70]*V[79]-Pow(V[58],2)*V[72]*V[82];
   if(!isfinite(V[283]) || FError) return 283;
   nCurrentVars=284;
   V[284]=2/(V[5])*V[54]*V[6]*V[283]/(Pow(V[8],2));
   if(!isfinite(V[284]) || FError) return 284;
   nCurrentVars=285;
   V[285]=-V[5]/(V[6])*M_SQRT2*(V[64]*V[130]*V[121]+V[57]*V[66]*V[130]*V[125]+V[64]*V[134]*V[125]-V[57]*V[66]*V[134]*V[121]+V[65]*V[129]*V[122]-V[58]*V[66]*V[133]*V[122]+V[65]*V[133]*V[126]+V[58]*V[66]*V[129]*V[126])/(2)/(V[119]);
   if(!isfinite(V[285]) || FError) return 285;
   nCurrentVars=286;
   V[286]=-V[5]/(V[6])*M_SQRT2*(V[64]*V[132]*V[123]+V[57]*V[66]*V[132]*V[127]+V[64]*V[136]*V[127]-V[57]*V[66]*V[136]*V[123]+V[65]*V[131]*V[124]-V[58]*V[66]*V[135]*V[124]+V[65]*V[135]*V[128]+V[58]*V[66]*V[131]*V[128])/(2)/(V[120]);
   if(!isfinite(V[286]) || FError) return 286;
   nCurrentVars=287;
   V[287]=-V[5]/(V[6])*M_SQRT2*(V[67]*V[130]*V[121]+V[57]*V[69]*V[130]*V[125]+V[67]*V[134]*V[125]-V[57]*V[69]*V[134]*V[121]+V[68]*V[129]*V[122]-V[58]*V[69]*V[133]*V[122]+V[68]*V[133]*V[126]+V[58]*V[69]*V[129]*V[126])/(2)/(V[119]);
   if(!isfinite(V[287]) || FError) return 287;
   nCurrentVars=288;
   V[288]=-V[5]/(V[6])*M_SQRT2*(V[67]*V[132]*V[123]+V[57]*V[69]*V[132]*V[127]+V[67]*V[136]*V[127]-V[57]*V[69]*V[136]*V[123]+V[68]*V[131]*V[124]-V[58]*V[69]*V[135]*V[124]+V[68]*V[135]*V[128]+V[58]*V[69]*V[131]*V[128])/(2)/(V[120]);
   if(!isfinite(V[288]) || FError) return 288;
   nCurrentVars=289;
   V[289]=-V[5]/(V[6])*M_SQRT2*(V[70]*V[130]*V[121]+V[57]*V[72]*V[130]*V[125]+V[70]*V[134]*V[125]-V[57]*V[72]*V[134]*V[121]+V[71]*V[129]*V[122]-V[58]*V[72]*V[133]*V[122]+V[71]*V[133]*V[126]+V[58]*V[72]*V[129]*V[126])/(2)/(V[119]);
   if(!isfinite(V[289]) || FError) return 289;
   nCurrentVars=290;
   V[290]=-V[5]/(V[6])*M_SQRT2*(V[70]*V[132]*V[123]+V[57]*V[72]*V[132]*V[127]+V[70]*V[136]*V[127]-V[57]*V[72]*V[136]*V[123]+V[71]*V[131]*V[124]-V[58]*V[72]*V[135]*V[124]+V[71]*V[135]*V[128]+V[58]*V[72]*V[131]*V[128])/(2)/(V[120]);
   if(!isfinite(V[290]) || FError) return 290;
   nCurrentVars=291;
   V[291]=V[5]/(V[54])*V[196]/(V[6])*V[66]/(V[46])/(2)/(V[196])/(2);
   if(!isfinite(V[291]) || FError) return 291;
   nCurrentVars=292;
   V[292]=V[5]/(V[54])*V[195]/(V[6])/(V[59])*(V[58]*V[140]*V[65]-V[139]*V[66]-V[57]*V[140]*V[64]+Pow(V[57],2)*V[66])/(V[195])/(2);
   if(!isfinite(V[292]) || FError) return 292;
   nCurrentVars=293;
   V[293]=V[5]/(V[54])*V[194]/(V[6])*V[66]/(V[46])/(2)/(V[194])/(2);
   if(!isfinite(V[293]) || FError) return 293;
   nCurrentVars=294;
   V[294]=V[5]/(V[54])*V[7]/(V[6])*V[66]*V[46]/(2)/(V[7])/(2);
   if(!isfinite(V[294]) || FError) return 294;
   nCurrentVars=295;
   V[295]=V[5]/(V[54])*V[196]/(V[6])*V[69]/(V[46])/(2)/(V[196])/(2);
   if(!isfinite(V[295]) || FError) return 295;
   nCurrentVars=296;
   V[296]=V[5]/(V[54])*V[195]/(V[6])/(V[59])*(V[58]*V[140]*V[68]-V[139]*V[69]-V[57]*V[140]*V[67]+Pow(V[57],2)*V[69])/(V[195])/(2);
   if(!isfinite(V[296]) || FError) return 296;
   nCurrentVars=297;
   V[297]=V[5]/(V[54])*V[194]/(V[6])*V[69]/(V[46])/(2)/(V[194])/(2);
   if(!isfinite(V[297]) || FError) return 297;
   nCurrentVars=298;
   V[298]=V[5]/(V[54])*V[7]/(V[6])*V[69]*V[46]/(2)/(V[7])/(2);
   if(!isfinite(V[298]) || FError) return 298;
   nCurrentVars=299;
   V[299]=V[5]/(V[54])*V[196]/(V[6])*V[72]/(V[46])/(2)/(V[196])/(2);
   if(!isfinite(V[299]) || FError) return 299;
   nCurrentVars=300;
   V[300]=V[5]/(V[54])*V[195]/(V[6])/(V[59])*(V[58]*V[140]*V[71]-V[139]*V[72]-V[57]*V[140]*V[70]+Pow(V[57],2)*V[72])/(V[195])/(2);
   if(!isfinite(V[300]) || FError) return 300;
   nCurrentVars=301;
   V[301]=V[5]/(V[54])*V[194]/(V[6])*V[72]/(V[46])/(2)/(V[194])/(2);
   if(!isfinite(V[301]) || FError) return 301;
   nCurrentVars=302;
   V[302]=V[5]/(V[54])*V[7]/(V[6])*V[72]*V[46]/(2)/(V[7])/(2);
   if(!isfinite(V[302]) || FError) return 302;
   nCurrentVars=303;
   V[303]=-V[5]/(V[6])*M_SQRT2*(V[58]*V[66]*V[133]*V[126]-V[65]*V[129]*V[126]+V[58]*V[66]*V[129]*V[122]+V[65]*V[133]*V[122]+V[57]*V[66]*V[134]*V[125]+V[64]*V[134]*V[121]+V[57]*V[66]*V[130]*V[121]-V[64]*V[130]*V[125])/(2)/(V[119])/(2);
   if(!isfinite(V[303]) || FError) return 303;
   nCurrentVars=304;
   V[304]=-V[5]/(V[6])*M_SQRT2*(V[58]*V[66]*V[135]*V[128]-V[65]*V[131]*V[128]+V[58]*V[66]*V[131]*V[124]+V[65]*V[135]*V[124]+V[57]*V[66]*V[136]*V[127]+V[64]*V[136]*V[123]+V[57]*V[66]*V[132]*V[123]-V[64]*V[132]*V[127])/(2)/(V[120])/(2);
   if(!isfinite(V[304]) || FError) return 304;
   nCurrentVars=305;
   V[305]=-V[5]/(V[6])*M_SQRT2*(V[58]*V[69]*V[133]*V[126]-V[68]*V[129]*V[126]+V[58]*V[69]*V[129]*V[122]+V[68]*V[133]*V[122]+V[57]*V[69]*V[134]*V[125]+V[67]*V[134]*V[121]+V[57]*V[69]*V[130]*V[121]-V[67]*V[130]*V[125])/(2)/(V[119])/(2);
   if(!isfinite(V[305]) || FError) return 305;
   nCurrentVars=306;
   V[306]=-V[5]/(V[6])*M_SQRT2*(V[58]*V[69]*V[135]*V[128]-V[68]*V[131]*V[128]+V[58]*V[69]*V[131]*V[124]+V[68]*V[135]*V[124]+V[57]*V[69]*V[136]*V[127]+V[67]*V[136]*V[123]+V[57]*V[69]*V[132]*V[123]-V[67]*V[132]*V[127])/(2)/(V[120])/(2);
   if(!isfinite(V[306]) || FError) return 306;
   nCurrentVars=307;
   V[307]=-V[5]/(V[6])*M_SQRT2*(V[58]*V[72]*V[133]*V[126]-V[71]*V[129]*V[126]+V[58]*V[72]*V[129]*V[122]+V[71]*V[133]*V[122]+V[57]*V[72]*V[134]*V[125]+V[70]*V[134]*V[121]+V[57]*V[72]*V[130]*V[121]-V[70]*V[130]*V[125])/(2)/(V[119])/(2);
   if(!isfinite(V[307]) || FError) return 307;
   nCurrentVars=308;
   V[308]=-V[5]/(V[6])*M_SQRT2*(V[58]*V[72]*V[135]*V[128]-V[71]*V[131]*V[128]+V[58]*V[72]*V[131]*V[124]+V[71]*V[135]*V[124]+V[57]*V[72]*V[136]*V[127]+V[70]*V[136]*V[123]+V[57]*V[72]*V[132]*V[123]-V[70]*V[132]*V[127])/(2)/(V[120])/(2);
   if(!isfinite(V[308]) || FError) return 308;
   nCurrentVars=309;
   V[309]=alphaQCD(V[61])/(V[156]);
   if(!isfinite(V[309]) || FError) return 309;
   nCurrentVars=310;
   V[310]=Sqrt(1+V[309]*(149/(double)((12))+V[309]*(68.6482-V[309]*212.447)));
   if(!isfinite(V[310]) || FError) return 310;
   nCurrentVars=311;
   V[311]=alphaQCD(V[62])/(V[156]);
   if(!isfinite(V[311]) || FError) return 311;
   nCurrentVars=312;
   V[312]=Sqrt(1+V[311]*(149/(double)((12))+V[311]*(68.6482-V[311]*212.447)));
   if(!isfinite(V[312]) || FError) return 312;
   nCurrentVars=313;
   V[313]=alphaQCD(V[63])/(V[156]);
   if(!isfinite(V[313]) || FError) return 313;
   nCurrentVars=314;
   V[314]=Sqrt(1+V[313]*(149/(double)((12))+V[313]*(68.6482-V[313]*212.447)));
   if(!isfinite(V[314]) || FError) return 314;
   nCurrentVars=315;
   V[315]=-Cabs(hGGeven(V[61],V[309],15,1,3,V[198],V[201],1,3,V[197],V[202],1,3,V[1],V[203],0,3,V[186],V[217],0,3,V[188],V[218],0,3,V[190],V[221],0,3,V[192],V[222],0,3,V[191],V[223],0,3,V[193],V[224],0,3,V[187],V[219],0,3,V[189],V[220],0,3,V[151],V[229],0,3,V[152],V[230],0,3,V[146],V[227],0,3,V[147],V[228]));
   if(!isfinite(V[315]) || FError) return 315;
   nCurrentVars=316;
   V[316]=-Cabs(hGGeven(V[62],V[311],15,1,3,V[198],V[205],1,3,V[197],V[206],1,3,V[1],V[207],0,3,V[186],V[235],0,3,V[188],V[236],0,3,V[190],V[239],0,3,V[192],V[240],0,3,V[191],V[241],0,3,V[193],V[242],0,3,V[187],V[237],0,3,V[189],V[238],0,3,V[151],V[247],0,3,V[152],V[248],0,3,V[146],V[245],0,3,V[147],V[246]));
   if(!isfinite(V[316]) || FError) return 316;
   nCurrentVars=317;
   V[317]=-Cabs(hGGeven(V[63],V[313],15,1,3,V[198],V[209],1,3,V[197],V[210],1,3,V[1],V[211],0,3,V[186],V[253],0,3,V[188],V[254],0,3,V[190],V[257],0,3,V[192],V[258],0,3,V[191],V[259],0,3,V[193],V[260],0,3,V[187],V[255],0,3,V[189],V[256],0,3,V[151],V[265],0,3,V[152],V[266],0,3,V[146],V[263],0,3,V[147],V[264]));
   if(!isfinite(V[317]) || FError) return 317;
   nCurrentVars=318;
   V[318]=-Cabs(V[199]*hAAeven(V[61],V[309],8,1,3,V[1],V[203],1,3,V[198],V[201],0,3,V[186],V[217],0,3,V[188],V[218],0,3,V[187],V[219],0,3,V[189],V[220],0,3,V[146],V[227],0,3,V[147],V[228])+V[200]*hAAeven(V[61],V[309],7,1,3,V[197],V[202],0,3,V[190],V[221],0,3,V[192],V[222],0,3,V[191],V[223],0,3,V[193],V[224],0,3,V[151],V[229],0,3,V[152],V[230])+hAAeven(V[61],V[309],11,2,1,V[54],V[267],1,1,V[119],V[285],1,1,V[120],V[286],1,1,V[7],V[204],0,1,V[8],V[272],0,1,V[182],V[213],0,1,V[183],V[214],0,1,V[184],V[215],0,1,V[185],V[216],0,1,V[141],V[225],0,1,V[142],V[226]));
   if(!isfinite(V[318]) || FError) return 318;
   nCurrentVars=319;
   V[319]=-Cabs(V[199]*hAAeven(V[62],V[311],8,1,3,V[1],V[207],1,3,V[198],V[205],0,3,V[186],V[235],0,3,V[188],V[236],0,3,V[187],V[237],0,3,V[189],V[238],0,3,V[146],V[245],0,3,V[147],V[246])+V[200]*hAAeven(V[62],V[311],7,1,3,V[197],V[206],0,3,V[190],V[239],0,3,V[192],V[240],0,3,V[191],V[241],0,3,V[193],V[242],0,3,V[151],V[247],0,3,V[152],V[248])+hAAeven(V[62],V[311],11,2,1,V[54],V[273],1,1,V[119],V[287],1,1,V[120],V[288],1,1,V[7],V[208],0,1,V[8],V[278],0,1,V[182],V[231],0,1,V[183],V[232],0,1,V[184],V[233],0,1,V[185],V[234],0,1,V[141],V[243],0,1,V[142],V[244]));
   if(!isfinite(V[319]) || FError) return 319;
   nCurrentVars=320;
   V[320]=-Cabs(V[199]*hAAeven(V[63],V[313],8,1,3,V[1],V[211],1,3,V[198],V[209],0,3,V[186],V[253],0,3,V[188],V[254],0,3,V[187],V[255],0,3,V[189],V[256],0,3,V[146],V[263],0,3,V[147],V[264])+V[200]*hAAeven(V[63],V[313],7,1,3,V[197],V[210],0,3,V[190],V[257],0,3,V[192],V[258],0,3,V[191],V[259],0,3,V[193],V[260],0,3,V[151],V[265],0,3,V[152],V[266])+hAAeven(V[63],V[313],11,2,1,V[54],V[279],1,1,V[119],V[289],1,1,V[120],V[290],1,1,V[7],V[212],0,1,V[8],V[284],0,1,V[182],V[249],0,1,V[183],V[250],0,1,V[184],V[251],0,1,V[185],V[252],0,1,V[141],V[261],0,1,V[142],V[262]));
   if(!isfinite(V[320]) || FError) return 320;
   nCurrentVars=321;
   V[321]=-Cabs(hGGodd(V[61],V[309],3,1,3,V[198],V[291],1,3,V[197],V[292],1,3,V[1],V[293]));
   if(!isfinite(V[321]) || FError) return 321;
   nCurrentVars=322;
   V[322]=-Cabs(hGGodd(V[62],V[311],3,1,3,V[198],V[295],1,3,V[197],V[296],1,3,V[1],V[297]));
   if(!isfinite(V[322]) || FError) return 322;
   nCurrentVars=323;
   V[323]=-Cabs(hGGodd(V[63],V[313],3,1,3,V[198],V[299],1,3,V[197],V[300],1,3,V[1],V[301]));
   if(!isfinite(V[323]) || FError) return 323;
   nCurrentVars=324;
   V[324]=-Cabs(V[200]*hAAodd(V[61],V[309],1,1,3,V[197],V[292])+V[199]*hAAodd(V[61],V[309],2,1,3,V[1],V[293],1,3,V[198],V[291])+hAAodd(V[61],V[309],3,1,1,V[7],V[294],1,1,V[119],V[303],1,1,V[120],V[304]));
   if(!isfinite(V[324]) || FError) return 324;
   nCurrentVars=325;
   V[325]=-Cabs(V[200]*hAAodd(V[62],V[311],1,1,3,V[197],V[296])+V[199]*hAAodd(V[62],V[311],2,1,3,V[1],V[297],1,3,V[198],V[295])+hAAodd(V[62],V[311],3,1,1,V[7],V[298],1,1,V[119],V[305],1,1,V[120],V[306]));
   if(!isfinite(V[325]) || FError) return 325;
   nCurrentVars=326;
   V[326]=-Cabs(V[200]*hAAodd(V[63],V[313],1,1,3,V[197],V[300])+V[199]*hAAodd(V[63],V[313],2,1,3,V[1],V[301],1,3,V[198],V[299])+hAAodd(V[63],V[313],3,1,1,V[7],V[302],1,1,V[119],V[307],1,1,V[120],V[308]));
   if(!isfinite(V[326]) || FError) return 326;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
