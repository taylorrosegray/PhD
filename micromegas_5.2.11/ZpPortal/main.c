/*====== Modules ===============
   Keys to switch on
   various modules of micrOMEGAs
================================*/
//#define MASSES_INFO
  /* Display information about mass spectrum  */

//#define CONSTRAINTS

//#define MONOJET
//#define HIGGSBOUNDS
//#define HIGGSSIGNALS
//#define LILITH
//#define SMODELS

#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN      /*  Calculate relic density in Freeze-in scenario  */

//#define INDIRECT_DETECTION
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation;
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line
     of sight;
     Calculate galactic propagation of positrons and antiprotons.
  */

//#define RESET_FORMFACTORS
  /* Modify default nucleus form factors,
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */
//#define CDM_NUCLEON
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */

//#define CDM_NUCLEUS
   // Calculate  exclusion rate for direct detection experiments Xenon1T, DarkSide50, CRESST-III, and PICO-60
            

//#define NEUTRINO
 /*  Neutrino signal of DM annihilation in Sun and Earth */

//#define DECAYS

//#define CROSS_SECTIONS

/*===== end of Modules  ======*/

/*===== Options ========*/
#define SHOWPLOTS
     /* Display  graphical plots on the screen */
//#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

   double Y[50],dY[50];
   for(int i=0; i<50;i++){ Y[i]=1+i; dY[i]=sqrt(i+1)/2;}

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */

  VZdecay=1; VWdecay=1;

  if(argc==1)
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
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2)
  {
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2);
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

     //// parameters
  //double MDM = 1.0;
  //double GZP = 1E-2;
 // assignValW("Mdm",MDM);

  //assignValW("gZp",GZP);
  // char *mZp = "MZp", *Mdm = "Mdm", *gZp = "gZp";
  // double MZP;
  // double logMZi = -3;
  // double logMZf = 3;
  // double logMXi = -3;
  // double logMXf = 3;
  // double GZP;
  // //double MZPl[] = {1.0,2.0,5.0,10.0,100.0};
  // //double GZPl[] = {1.0E-10,1.0E-3,1.0E-2,1.0E+1};
  // double MDM;
  // int i;
  // int k = 10;
  FILE *file;
  char filename[25];
  snprintf(filename, sizeof(filename), "Omega.dat");
  file = fopen(filename,"w");
//for (i = 0; i < k; ++i){
  //   //MZP = pow(10,logMZi+i*(logMZf-logMZi)/k);
  //   MDM = pow(10,logMXi+i*(logMXf-logMXi)/k);
  //   //MZP = MZPl[i];
  //  // assignValW(Mdm,MDM);
  //   //printf("mZp = %.2f \n",findValW("MZp"));
  //   //GZP = GZPl[i];
  //   //assignValW(gZp,GZP);
  //   //printf("MZp = %.2e \n",findValW("MZp"));


#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.2E pb\n",csLim);
}
#endif

#ifdef MONOJET
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);
}
#endif

#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=3, NHch=1; // number of neutral and charged Higgs particles.
   int HB_id[3]={0,0,0},HB_result[3];
   double  HB_obsratio[3],HS_observ=-1,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch); 
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) if(HB_id[i]) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif 
#ifdef HIGGSSIGNALS
   if(HS_observ>=0)
   {
     printf("HiggsSignals(%s)\n",HS_version); 
     printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
   }
#endif   
}
#endif


#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char Lilith_version[50];
   if(LilithMO("Lilith_in.xml"))
   {
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n",
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }
   } else printf("LILITH: there is no Higgs candidate\n");
}
#endif


#ifdef SMODELS
{ int status=0, smodelsOK=0; 
  double Rvalue, Rexpected, SmoLsig, SmoLmax, SmoLSM;
  char analysis[50]={},topology[100]={},smodelsInfo[100];
  int LHCrun=LHC8|LHC13;  //  LHC8  - 8TeV; LHC13  - 13TeV;   

  printf("\n\n=====  LHC constraints with SModelS  =====\n\n");

#include "../include/SMODELS.inc" // SLHA interface with SModelS

  printf("SModelS %s \n",smodelsInfo);
  if(smodelsOK) 
  { printf(" highest r-value = %.2E",Rvalue); 
    if(Rvalue>0) 
    { printf(" from %s, topology: %s ",analysis,topology);
      if(Rexpected>0) 
      { printf("\n expected r = %.2E ",Rexpected);
        if(SmoLsig>0) 
        { printf("\n -2log (L_signal, L_max, L_SM) = %.2E %.2E %.2E", 
                  -2*log(SmoLsig),-2*log(SmoLmax),-2*log(SmoLSM)); }
      }
    }  
    if(status==1) { printf("\n excluded by SMS results"); }
    else if(status==0) printf("\n not excluded"); 
    else if(status==-1) printf("\n not not tested by results in SModelS database"); 
    printf("\n");
  } else system("cat smodels.err"); // problem: see smodels.err
}   
#endif

// FILE *file;
// char filename[25];
// snprintf(filename, sizeof(filename), "channels_%d.dat",i);
// file = fopen(filename,"w");
#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Omega;
  int i,err;
  printf("\n==== Calculation of relic density =====\n");
  VZdecay=0; VWdecay=0; cleanDecayTable();
  if( CDM1 && CDM2)
  {

    Omega= darkOmega2(fast,Beps);
    printf("Omega_1h^2=%.2E\n", Omega*(1-fracCDM2));
    printf("Omega_2h^2=%.2E\n", Omega*fracCDM2);

  } else
  {  double Xf;
    // float T;
    // float Ti = -2.0; //10^Ti: min temperature
    // float Tf = 2.0; //10^Tf: max temperature
    // int k = 10.0;
    // double YEQ;
    // double Y;
    Omega=darkOmegaFO(&Xf,fast,Beps);
    //Omega=darkOmega(&Xf,fast,Beps,&err);
    //fprintf(file,"%.3E \t %.3E \n ",MZP,Omega);
    fprintf(file,"%.3E \t %.3E \n ",findValW("Mdm"),Omega);


    //Y = YF(1.0);
    //printChannels(Xf,cut,Beps,0,file);
    printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
    //printf("YF = %.2e",Y);

  //   FILE *file;
  //  char filename[25];
   //snprintf(filename, sizeof(filename), "Y.dat");
  // file = fopen(filename,"w");
  // for (i = 0; i < k+1; ++i)
  //   {T = pow(10,Ti+i*(Tf-Ti)/k);
  //   Y = YF(T);
  //   YEQ = Yeq(T);
  //   printf("Y = %.3E\n",Y);
  //   //print to file
  //   fprintf(file,"%.3E \t %.3E \t %.3E \n ",T,Y,YEQ);
  //   };
  //   fclose(file);
     //if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  }
}
#endif

#ifdef FREEZEIN
{
  double TR=1e10;
  double omegaFi;
  toFeebleList(CDM1);
//  toFeebleList("Zp");

  VWdecay=0; VZdecay=0;

#ifdef SHOWPLOTS
  omegaFi=darkOmegaFiDecay(TR,"Zp",0,1);
#else
   omegaFi=darkOmegaFiDecay(TR,"Zp",0,0);
#endif
  printf("OmegaFiDecay(Zp)= %E\n", omegaFi);

  omegaFi=darkOmegaFi(TR,&err);
  printf("Omega freeze-in=%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);

  if(CDM1 && CDM2)
  { double omg1,omg2;
    sort2FiDm(&omg1,&omg2);
    printf("Omega1=%.3E Omega2=%.3E\n", omg1,omg2);
  }
}
#endif


#ifdef INDIRECT_DETECTION
{
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;

printf("\n==== Indirect detection =======\n");

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.

       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections
    */



  {
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);
  }

  {
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);
  }

  {
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  );
#ifdef SHOWPLOTS
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);
  }
}
#endif

#ifdef RESET_FORMFACTORS
{
/*
   The user has approach to form factors  which specifies quark contents
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma";
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");

  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
//                    To restore default form factors of  version 2  call
     calcScalarQuarkFF(0.553,18.9,55.,243.5);


  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call
//  calcScalarQuarkFF(0.56,20.2,34,42);


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;
  double csSIp,csSIn,csSDp,csSDn, csSIp_,csSIn_,csSDp_,csSDn_;

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");

    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM1);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] );

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    csSIp=  SCcoeff*pA0[0]*pA0[0];  csSIp_=  SCcoeff*pA0[1]*pA0[1];
    csSDp=3*SCcoeff*pA5[0]*pA5[0];  csSDp_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn=  SCcoeff*nA0[0]*nA0[0];  csSIn_=  SCcoeff*nA0[1]*nA0[1];
    csSDn=3*SCcoeff*nA5[0]*nA5[0];  csSDn_=3*SCcoeff*nA5[1]*nA5[1];
                    
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIp,csSIp_,csSDp,csSDp_);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIn,csSIn_,csSDn,csSDn_);      
}
#endif

#ifdef CDM_NUCLEUS
{ char* expName; 
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
       if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n");  
}


#endif

#ifdef NEUTRINO
if(!CDM1 || !CDM2)
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;

 printf("\n===============Neutrino Telescope=======  for  ");
 if(forSun) printf("Sun\n"); else printf("Earth\n");

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
#endif
{
    printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));
}

/* Upward events */

  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);
#endif
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu);
#endif
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}
#endif


#ifdef DECAYS
{ char*  pname = "Zp";
  txtList L;
  double width;
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
  }
}
#endif

#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;

  next=nextOdd(1,&nextM);
  if(next && nextM<1000)
  {
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);

     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  }
}

#endif

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f   smodels.*");
#endif

fclose(file);


  killPlots();
  return 0;
}