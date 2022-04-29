#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */
//#define CROSS_SECTIONS 
//#define THERMALLY_AV_CROSS_SECTIONS

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"

int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
  VZdecay=0; VWdecay=0;  

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
     //printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1); 
     //if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     //if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2) 
  { 
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }


// #ifdef OMEGA
// { int fast=1;
//   double Beps=1.E-4, cut=0.01;
//   double Omega;  
//   int i,err; 
//   printf("\n==== Calculation of relic density =====\n");   

//   if(CDM1 && CDM2) 
//   {
  
//     Omega= darkOmega2(fast,Beps);



//   displayPlot("vs11","T",Tend,Tstart,0,5
//       ,"vs1100",0,vs1100F,NULL
//       ,"vs1120",0,vs1120F,NULL
//       ,"vs1122",0,vs1122F,NULL
//       ,"vs1110",0,vs1110F,NULL
//       ,"vs1112",0,vs1112F,NULL
//       );
      
//   displayPlot("vs12","T",Tend,Tstart,0,4
//              ,"vs1210",0,vs1210F,NULL
//              ,"vs1222",0,vs1222F,NULL
//              ,"vs1120",0,vs1120F,NULL
//              ,"vs1211",0,vs1211F,NULL);
                                
//   displayPlot("vs22","T",Tend,Tstart,0,5
//              ,"vs2200",0,vs2200F,NULL
//              ,"vs2210",0,vs2210F,NULL
//              ,"vs2211",0,vs2211F,NULL
//              ,"vs2220",0,vs2220F,NULL
//              ,"vs2221",0,vs2221F,NULL
//              );

//   displayPlot("dY","T",  Tend,Tstart,0,2,"dY1",0,dY1F,NULL,"dY2",0,dY2F,NULL);
//   displayPlot("Y","T",   Tend,Tstart,0,2,"Y1" ,0,Y1F,NULL,"Y2",0,Y2F,NULL);

                                
      
      
    
  
//     printf("Omega_1h^2=%.2E\n", Omega*(1-fracCDM2));
//     printf("Omega_2h^2=%.2E\n", Omega*fracCDM2);
//   } else
//   {  double Xf;
//      Omega=darkOmega(&Xf,fast,Beps,&err);
//      printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
//      if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
//   }
// }

// #endif

#ifdef OMEGA
	{ int fast=1;
	double Beps=1.E-6, cut=0.01;
	int prcnt=0;
	double Omega;  
	double MDM;
	double MAP;
	int i, err; 
	//printf("\n==== Calculation of relic density =====\n");
	double Xf;
	double Y;
	FILE *file;
  char filename[25];
  snprintf(filename, sizeof(filename), "Omega.dat");
  file = fopen(filename,"w");
  	//double Tl[] = {1E2,1E1,1.0,0.1,0.01};
	// double log_mdmi = -1.0;
	// double log_mdmf = 5.0;
	//int k = 10;
	//assignValW("g4",1.0E+2);

	 //    printf("\t Xf=%.2e Omega=%.2e\n",Xf,Omega);};	
	//for (i = 0; i < k+1; ++i){
		// MDM = pow(10,log_mdmi+i*(log_mdmf-log_mdmi)/k);
		// //MDM = 0.1;
		// //MAP = MDM*3.0;
		// MAP = 1.0E1;
		// //printf("%.2e \n",mdm);
		// assignValW("mdm",MDM);
		// assignValW("mAp",MAP);
		// printf("m_chi = %.2e, mAp = %.2e", findValW("mdm"),findValW("mAp"));
	 //    Omega=darkOmegaFO(&Xf,fast,Beps);
	 //    //Omega=darkOmega(&Xf,fast,Beps,&err);
	//printf("m_chi = %.2e, mAp = %.2e", findValW("mdm"),findValW("mAp"));
	//Omega=darkOmega(&Xf,fast,Beps,&err);
	Omega=darkOmegaFO(&Xf,fast,Beps);

	FILE *channels;
	char filename_chan[40];
	snprintf(filename_chan, sizeof(filename_chan), "channels/channels_mdm%.2e.dat",findValW("Mdm"));
  channels = fopen(filename_chan,"w");
	printChannels(Xf,cut,Beps,prcnt,channels);
	fclose(channels);
	// for (i = 0; i < 5; ++i){
	// Y = YF(Tl[i]);
	// printf("%.2e",Y);};
	//printChannels(Xf,cut,Beps,prcnt,file);
	printf("\t Xf=%.2e Omega=%.2e\n",Xf,Omega);
	fprintf(file,"%.3E \t %.3E \n",Xf,Omega);
	fclose(file);


    //if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
    }  
#endif


#ifdef FREEZEIN
{
  double TR=1E10;
  double omegaFi;  
  toFeebleList(CDM1);
  VWdecay=0; VZdecay=0;
  
  omegaFi=darkOmegaFi(TR,&err);
  printf("omega freeze-in=%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);
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

#ifdef THERMALLY_AV_CROSS_SECTIONS
	{//printf("\n==== Calculation of thermally averaged cross section =====\n");
	double Beps=1.E-4;
	int fast=1;
	int i;
	float T;
	float Ti = -5.0; //10^Ti: min temperature
	float Tf = 5.0; //10^Tf: max temperature
	float sigmav;
	int k = 100;
	FILE *file;
	char filename[25];
	snprintf(filename, sizeof(filename), "sigmav.dat");
	file = fopen(filename,"w");
	for (i = 0; i < k+1; ++i)
		{T = pow(10,Ti+i*(Tf-Ti)/k);
		//printf("i = %d \t T=%.2E \t", i, T);
		sigmav = vSigma(T,Beps,fast);
		//printf("sigmav = %.3E\n",sigmav);

		//print to file
		fprintf(file,"%.3E \t %.3E \n",T,sigmav);
		};
		fclose(file);
}
#endif

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
//  system("rm -f   smodels.*");  
#endif 

  killPlots();
  return 0;
}
