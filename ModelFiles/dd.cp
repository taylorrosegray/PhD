/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/

#define cs_vs_mdm //calculates cross sections and couplings vs mdm
#define cs_vs_mzp //calculates cross sections and couplings vs mZp
#define CLEAN
/*===== End of DEFINE  settings ===== */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
static double medMass=findValW("MZp"); // mediator mass
static double myDDfactor( double E, int A) {  return pow( 1/(1+2*(A*0.94)*(E*1E-6)/medMass/medMass),2);}

// Mass Ranges
double mMin=0.1,mMax=1.0E4;   // DM mass range
double mzMin=1.0E-3,mzMax=1.0E4; //Mz' mass range
double Nmass=0.938; /*nucleon mass*/

int main(int argc,char** argv)
{ 
  int err;
  char cdmName[10];
  int spin2, charge3,cdim;
  char*expName="";
  ForceUG=1;  // to Force Unitary Gauge assign 1
  //useSLHAwidth=0;
  VZdecay=0; VWdecay=0;  

  if(argc==1) // If no input file (data.par)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }

  err=readVar(argv[1]);

  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}

  err=sortOddParticles(cdmName);

  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  if(CDM1) 
  { 
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

double Mcdm_og;
Mcdm_og = findValW("Mdm"); // DM mass
printf("mdm = %.1e",Mcdm_og);

#ifdef cs_vs_mdm
{
double m;
medMass=findValW("MZp");
m=medMass;
      
double csDarkSide=0.,csCresst=0.,csPICO=0.,csXENON1T=0.,csCresstn=0.,csXENON1Tn=0.,csXENON1Tp=0.,csPICOp=0.;
double fDarkSide=0.,fCresst=0.,fPICO=0.,fXENON1T=0.,fCresstn=0.,fXENON1Tn=0.,fXENON1Tp=0.,fPICOp=0.;        
double vfvxDarkSide=0.,vfvxCresst=0.,vfvxPICO,vfvxXENON1T=0.,vfvxCresstn=0.,vfvxXENON1Tn=0.,vfvxXENON1Tp=0.,vfvxPICOp=0.;             

// Print to file "limits_mzp.dat"
char temp[100];
sprintf(temp,"/home/taylor/Documents/Share/limits_mzp%.1E.dat",findValW("MZp"));
ofstream outrate (temp);
outrate << "mdm";
outrate <<"\t"<< "cs_DarkSide" <<"\t"<< "vfvx_DarkSide";
outrate <<"\t"<< "cs_Cresst" <<"\t"<<  "vfvx_Cresst";
outrate  <<"\t"<< "cs_PICO" <<"\t"<<  "vfvx_PICO";
outrate <<"\t"<< "cs_XENON1T" <<"\t"<< "vfvx_XENON1T";

outrate <<"\t"<< "cs_SD_Cresst" <<"\t"<< "vfvx_SD_Cresst";
outrate <<"\t"<< "cs_SD_PICO" <<"\t"<<  "vfvx_SD_PICO";
outrate  <<"\t"<< "cs_SD_XENON1Tp" <<"\t"<<  "vfvx_SD_XENON1Tp";
outrate  <<"\t"<< "cs_SD_XENON1Tn" <<"\t"<<  "vfvx_SD_XENON1Tn";
outrate << endl;
 

if(CDM1)
  { 
  dNdEfact=myDDfactor; // DD_factorCS() takes this into account
  double dNdE[RE_DIM];
  double dNdEarr[10][RE_DIM];
   
for(int i=0;i<1000;i++)
  { Mcdm=mMin*pow(mMax/mMin, (i+0.5)/1000.);
    assignValW("Mdm",Mcdm);
    Mcdm1=Mcdm;
    
    // SI limits
    csDarkSide = DS50_90(Mcdm)*1e36;
    fDarkSide = DD_factorCS(DarkSide_2018,0.1,Maxwell, csDarkSide, csDarkSide,0,0,&expName);
    csDarkSide=fDarkSide*csDarkSide*1E-36;
    vfvxDarkSide = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csDarkSide/3.9204E-28); 

    csCresst = CRESST_III_90(Mcdm)*1E36;
    fCresst= DD_factorCS(CRESST_2019,0.1,Maxwell, csCresst, csCresst,0,0,&expName);
    csCresst=fCresst*csCresst*1E-36;  
    vfvxCresst = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csCresst/3.9204E-28); 

    csPICO=PICO60_90(Mcdm)*1E36;
    fPICO = DD_factorCS(PICO_2019,0.1,Maxwell, csPICO, csPICO,0,0,&expName);
    csPICO = csPICO*fPICO*1E-36;
    vfvxPICO = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csPICO/3.9204E-28); 

    csXENON1T = XENON1T_90(Mcdm)*1E36;
    fXENON1T = DD_factorCS(XENON1T_2018,0.1,Maxwell, csXENON1T, csXENON1T,0,0,&expName);
    csXENON1T=csXENON1T*fXENON1T*1E-36;     
    vfvxXENON1T = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csXENON1T/3.9204E-28); 

    // SD limits
    double spin_fractions=0.35; // (-0.42+0.85-0.08)

    setSpinDepFF(XENON1T_2018,SHELL);
    csXENON1Tp = XENON1T_SDp_90(Mcdm)*1e36;
    fXENON1Tp = DD_factorCS(XENON1T_2018,0.1,Maxwell,0,0,csXENON1Tp,0,&expName);
    csXENON1Tp = fXENON1Tp*csXENON1Tp*1E-36;
    vfvxXENON1Tp = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csXENON1Tp/3.9204E-28)/spin_fractions;

    csCresstn = CRESST_III_SDn_90(Mcdm)*1E36;
    fCresstn= DD_factorCS(CRESST_2019,0.1,Maxwell,0,0,0,csCresstn,&expName);
    csCresstn = fCresstn*csCresstn*1E-36;  
    vfvxCresstn = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csCresstn/3.9204E-28)/spin_fractions; 

    setSpinDepFF(PICO_2019,SHELL);
    csPICOp=PICO60_SDp_90(Mcdm)*1E36;
    fPICOp = DD_factorCS(PICO_2019,0.1,Maxwell,0,0,csPICOp,0,&expName);
    csPICOp = csPICOp*fPICOp*1E-36;
    vfvxPICOp = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+Mcdm))*sqrt(csPICOp/3.9204E-28)/spin_fractions; 

    setSpinDepFF(XENON1T_2018,SHELL);
    csXENON1Tn = XENON1T_SDn_90(Mcdm)*1e36;
    fXENON1Tn = DD_factorCS(XENON1T_2018,0.1,Maxwell, 0, 0,0,csXENON1Tn,&expName);
    csXENON1Tn = fXENON1Tn*csXENON1Tn*1E-36;
    vfvxXENON1Tn = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csXENON1Tn/3.9204E-28)/spin_fractions; 

    outrate << scientific;
    outrate << Mcdm;
    outrate <<"\t"<< csDarkSide <<"\t"<<  vfvxDarkSide;
    outrate <<"\t"<< csCresst   <<"\t"<<  vfvxCresst;
    outrate <<"\t"<< csPICO     <<"\t"<<  vfvxPICO;
    outrate <<"\t"<< csXENON1T  <<"\t"<<  vfvxXENON1T;
    outrate <<"\t"<< csCresstn  <<"\t"<<  vfvxCresstn;
    outrate <<"\t"<< csPICOp    <<"\t"<<  vfvxPICOp;
    outrate <<"\t"<< csXENON1Tp <<"\t"<<  vfvxXENON1Tp;
    outrate <<"\t"<< csXENON1Tn <<"\t"<<  vfvxXENON1Tn;
    outrate << endl;
  }
}
}
#endif

#ifdef cs_vs_mzp
{
double m; //medMass
Mcdm = Mcdm_og; // DM mass
assignValW("Mdm",Mcdm);
Mcdm1=Mcdm;

double csDarkSide=0.,csCresst=0.,csPICO=0.,csXENON1T=0.,csCresstn=0.,csXENON1Tn=0.,csXENON1Tp=0.,csPICOp=0.;
double fDarkSide=0.,fCresst=0.,fPICO=0.,fXENON1T=0.,fCresstn=0.,fXENON1Tn=0.,fXENON1Tp=0.,fPICOp=0.;        
double vfvxDarkSide=0.,vfvxCresst=0.,vfvxPICO,vfvxXENON1T=0.,vfvxCresstn=0.,vfvxXENON1Tn=0.,vfvxXENON1Tp=0.,vfvxPICOp=0.;        

// Print to file "limits_mdm.dat"
char temp[100];
sprintf(temp,"/home/taylor/Documents/Share/limits_mdm%.1E.dat",findValW("Mdm"));
ofstream outrate (temp);
outrate << "mZp";
outrate <<"\t"<< "cs_SI_DarkSide" <<"\t"<< "vfvx_SI_DarkSide";
outrate <<"\t"<< "cs_SI_Cresst" <<"\t"<<  "vfvx_SI_Cresst";
outrate  <<"\t"<< "cs_SI_PICO" <<"\t"<<  "vfvx_SI_PICO";
outrate <<"\t"<< "cs_SI_XENON1T" <<"\t"<< "vfvx_SI_XENON1T";

outrate <<"\t"<< "cs_SD_Cresst" <<"\t"<< "vfvx_SD_Cresst";
outrate <<"\t"<< "cs_SD_PICO" <<"\t"<<  "vfvx_SD_PICO";
outrate  <<"\t"<< "cs_SD_XENON1Tp" <<"\t"<<  "vfvx_SD_XENON1Tp";
outrate  <<"\t"<< "cs_SD_XENON1Tn" <<"\t"<<  "vfvx_SD_XENON1Tn";
outrate << endl;
 

if(CDM1)
  { 
  dNdEfact=myDDfactor; // DD_factorCS() takes this into account
  double dNdE[RE_DIM];
  double dNdEarr[10][RE_DIM];
  
   
for(int i=0;i<1000;i++)
  { medMass = mzMin*pow(mzMax/mzMin, (i+0.5)/1000.);
    m=medMass;
    assignValW("MZp",medMass);

    // SI limits
    csDarkSide = DS50_90(Mcdm)*1e36;
    fDarkSide = DD_factorCS(DarkSide_2018,0.1,Maxwell, csDarkSide, csDarkSide,0,0,&expName);
    csDarkSide=fDarkSide*csDarkSide*1E-36;
    vfvxDarkSide = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csDarkSide/3.9204E-28); 

    csCresst = CRESST_III_90(Mcdm)*1E36;
    fCresst= DD_factorCS(CRESST_2019,0.1,Maxwell, csCresst, csCresst,0,0,&expName);
    csCresst=fCresst*csCresst*1E-36;  
    vfvxCresst = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csCresst/3.9204E-28); 

    csPICO=PICO60_90(Mcdm)*1E36;
    fPICO = DD_factorCS(PICO_2019,0.1,Maxwell, csPICO, csPICO,0,0,&expName);
    csPICO = csPICO*fPICO*1E-36;
    vfvxPICO = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csPICO/3.9204E-28); 

    csXENON1T = XENON1T_90(Mcdm)*1E36;
    fXENON1T = DD_factorCS(XENON1T_2018,0.1,Maxwell, csXENON1T, csXENON1T,0,0,&expName);
    csXENON1T=csXENON1T*fXENON1T*1E-36;     
    vfvxXENON1T = m*m* sqrt(M_PI)/3./(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csXENON1T/3.9204E-28); 

    // SD limits
    double spin_fractions=0.35; // (-0.42+0.85-0.08)

    setSpinDepFF(XENON1T_2018,SHELL);
    csXENON1Tp = XENON1T_SDp_90(Mcdm)*1e36;
    fXENON1Tp = DD_factorCS(XENON1T_2018,0.1,Maxwell,0,0,csXENON1Tp,0,&expName);
    csXENON1Tp = fXENON1Tp*csXENON1Tp*1E-36;
    vfvxXENON1Tp = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csXENON1Tp/3.9204E-28)/spin_fractions;

    csCresstn = CRESST_III_SDn_90(Mcdm)*1E36;
    fCresstn= DD_factorCS(CRESST_2019,0.1,Maxwell,0,0,0,csCresstn,&expName);
    csCresstn = fCresstn*csCresstn*1E-36;  
    vfvxCresstn = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csCresstn/3.9204E-28)/spin_fractions; 

    setSpinDepFF(PICO_2019,SHELL);
    csPICOp=PICO60_SDp_90(Mcdm)*1E36;
    fPICOp = DD_factorCS(PICO_2019,0.1,Maxwell,0,0,csPICOp,0,&expName);
    csPICOp = csPICOp*fPICOp*1E-36;
    vfvxPICOp = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+Mcdm))*sqrt(csPICOp/3.9204E-28)/spin_fractions; 

    setSpinDepFF(XENON1T_2018,SHELL);
    csXENON1Tn = XENON1T_SDn_90(Mcdm)*1e36;
    fXENON1Tn = DD_factorCS(XENON1T_2018,0.1,Maxwell, 0, 0,0,csXENON1Tn,&expName);
    csXENON1Tn = fXENON1Tn*csXENON1Tn*1E-36;
    vfvxXENON1Tn = m*m* sqrt(M_PI/3.)/(Nmass*Mcdm/(Nmass+ Mcdm))*sqrt(csXENON1Tn/3.9204E-28)/spin_fractions; 
    
    // Print values to file
    outrate << scientific;
    outrate << medMass;
    outrate <<"\t"<< csDarkSide <<"\t"<<  vfvxDarkSide;
    outrate <<"\t"<< csCresst   <<"\t"<<  vfvxCresst;
    outrate <<"\t"<< csPICO     <<"\t"<<  vfvxPICO;
    outrate <<"\t"<< csXENON1T  <<"\t"<<  vfvxXENON1T;
    outrate <<"\t"<< csCresstn  <<"\t"<<  vfvxCresstn;
    outrate <<"\t"<< csPICOp    <<"\t"<<  vfvxPICOp;
    outrate <<"\t"<< csXENON1Tp <<"\t"<<  vfvxXENON1Tp;
    outrate <<"\t"<< csXENON1Tn <<"\t"<<  vfvxXENON1Tn;
    outrate << endl;
  } 
}
}
#endif

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.in  smodels.log  smodels.out  summary.*");  
#endif 

  return 0;
}
