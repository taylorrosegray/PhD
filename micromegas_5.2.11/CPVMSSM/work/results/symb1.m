(*
     ==============================
     *  CalcHEP  3.8.9 *
     ==============================
  process  G(p1)+G(p2)->G(p3)+G(p4)
*)

parameters={
 Mtp -> 1.73070000000*10^(2)
,McMc -> 1.30000000000*10^(0)
,Q -> 1.00000000000*10^(2)
,EE -> 3.13330000000*10^(-1)
,SW -> 4.74000000000*10^(-1)
,aMu -> 1.50000000000*10^(2)
,fiMu -> 0.00000000000*10^(0)
,Mq -> 5.00000000000*10^(-2)
,aAt -> 1.00000000000*10^(3)
,fiAt -> 9.00000000000*10^(1)
,aAb -> 1.00000000000*10^(3)
,fiAb -> 9.00000000000*10^(1)
,Mq2 -> 1.00000000000*10^(3)
,Mu2 -> 1.00000000000*10^(3)
,Md2 -> 1.00000000000*10^(3)
,MZ -> 9.11870000000*10^(1)
,tb -> 5.00000000000*10^(1)
,Mh1 -> slhaVal["MASS",MZ,1,25]
,Mh2 -> slhaVal["MASS",MZ,1,35]
,Mh3 -> slhaVal["MASS",MZ,1,36]
,Zh11 -> slhaVal["CVHMIX",MZ,2,1,1]
,Zh12 -> slhaVal["CVHMIX",MZ,2,1,2]
,Zh13 -> -slhaVal["CVHMIX",MZ,2,1,3]*sb+slhaVal["CVHMIX",MZ,2,1,4]*cb
,Zh21 -> slhaVal["CVHMIX",MZ,2,2,1]
,Zh22 -> slhaVal["CVHMIX",MZ,2,2,2]
,Zh23 -> -slhaVal["CVHMIX",MZ,2,2,3]*sb+slhaVal["CVHMIX",MZ,2,2,4]*cb
,Zh31 -> slhaVal["CVHMIX",MZ,2,3,1]
,Zh32 -> slhaVal["CVHMIX",MZ,2,3,2]
,Zh33 -> -slhaVal["CVHMIX",MZ,2,3,3]*sb+slhaVal["CVHMIX",MZ,2,3,4]*cb
,Td3r -> slhaVal["THRESHOLD",MZ,2,5,1]
,Td3i -> slhaVal["THRESHOLD",MZ,2,5,2]
,MSt1 -> slhaVal["MASS",MZ,1,1000006]
,MSt2 -> slhaVal["MASS",MZ,1,2000006]
,Zt11r -> slhaVal["STOPMIX",MZ,2,1,1]
,Zt21r -> slhaVal["STOPMIX",MZ,2,2,1]
,Zt21i -> slhaVal["IMSTOPMIX",MZ,2,2,1]
,MSb1 -> slhaVal["MASS",MZ,1,1000005]
,MSb2 -> slhaVal["MASS",MZ,1,2000005]
,Zb11r -> slhaVal["SBOTMIX",MZ,2,1,1]
,Zb21r -> slhaVal["SBOTMIX",MZ,2,2,1]
,Zb21i -> slhaVal["IMSBOTMIX",MZ,2,2,1]
,Mt -> MtEff[Q]
,Mb -> MbEff[Q]
,Mc -> McEff[Q]
,Mbp -> bPoleMass[]
,Mcp -> McMc*(1+4/3*alphaQCD[McMc]/PI)
,aQCDh1 -> alphaQCD[Mh1]/PI
,aQCDh2 -> alphaQCD[Mh2]/PI
,aQCDh3 -> alphaQCD[Mh3]/PI
,LGGh1 -> -
 cabs[hGGeven[Mh1,aQCDh1,15,1,3,Mcp,ah1F_c,1,3,Mbp,ah1F_b,1,3,Mtp,ah1F_t,0,3,MSuL,ah1S_uL,0,3,MSuR,ah1S_uR,0,3,MSdL,ah1S_dL,0,3,MSdR,ah1S_dR,0,3,MSsL,ah1S_sL,0,3,MSsR,ah1S_sR,0,3,MScL,ah1S_cL,0,3,MScR,ah1S_cR,0,3,MSb1,ah1S_b1,0,3,MSb2,ah1S_b2,0,3,MSt1,ah1S_t1,0,3,MSt2,ah1S_t2]]
 
,LGGh2 -> -
 cabs[hGGeven[Mh2,aQCDh2,15,1,3,Mcp,ah2F_c,1,3,Mbp,ah2F_b,1,3,Mtp,ah2F_t,0,3,MSuL,ah2S_uL,0,3,MSuR,ah2S_uR,0,3,MSdL,ah2S_dL,0,3,MSdR,ah2S_dR,0,3,MSsL,ah2S_sL,0,3,MSsR,ah2S_sR,0,3,MScL,ah2S_cL,0,3,MScR,ah2S_cR,0,3,MSb1,ah2S_b1,0,3,MSb2,ah2S_b2,0,3,MSt1,ah2S_t1,0,3,MSt2,ah2S_t2]]
 
,LGGh3 -> -
 cabs[hGGeven[Mh3,aQCDh3,15,1,3,Mcp,ah3F_c,1,3,Mbp,ah3F_b,1,3,Mtp,ah3F_t,0,3,MSuL,ah3S_uL,0,3,MSuR,ah3S_uR,0,3,MSdL,ah3S_dL,0,3,MSdR,ah3S_dR,0,3,MSsL,ah3S_sL,0,3,MSsR,ah3S_sR,0,3,MScL,ah3S_cL,0,3,MScR,ah3S_cR,0,3,MSb1,ah3S_b1,0,3,MSb2,ah3S_b2,0,3,MSt1,ah3S_t1,0,3,MSt2,ah3S_t2]]
 
,imLGGh1 -> -
 cabs[hGGodd[Mh1,aQCDh1,3,1,3,Mcp,imah1F_c,1,3,Mbp,imah1F_b,1,3,Mtp,imah1F_t]]
 
,imLGGh2 -> -
 cabs[hGGodd[Mh2,aQCDh2,3,1,3,Mcp,imah2F_c,1,3,Mbp,imah2F_b,1,3,Mtp,imah2F_t]]
 
,imLGGh3 -> -
 cabs[hGGodd[Mh3,aQCDh3,3,1,3,Mcp,imah3F_c,1,3,Mbp,imah3F_b,1,3,Mtp,imah3F_t]]
 
,wh1 -> 0.00000000000*10^(0)
,wh2 -> 0.00000000000*10^(0)
,wh3 -> 0.00000000000*10^(0)
,GG -> 0.00000000000*10^(0)
           };

substitutions={
 Rqcdh3->Sqrt[1+aQCDh3*(149/12+aQCDh3*(68.6482-aQCDh3*212.447))]
,Rqcdh2->Sqrt[1+aQCDh2*(149/12+aQCDh2*(68.6482-aQCDh2*212.447))]
,Rqcdh1->Sqrt[1+aQCDh1*(149/12+aQCDh1*(68.6482-aQCDh1*212.447))]
,imah3F_t->EE/MW*Mt/SW*Zh33/tb/2/Mt/2
,imah3F_b->EE/MW*Mb/SW/s2b*(cb*Td3i*Zh32-Td3r*Zh33-sb*Td3i*Zh31+sb^2*Zh33)
 /Mb/2
,imah3F_c->EE/MW*Mc/SW*Zh33/tb/2/Mc/2
,imah2F_t->EE/MW*Mt/SW*Zh23/tb/2/Mt/2
,imah2F_b->EE/MW*Mb/SW/s2b*(cb*Td3i*Zh22-Td3r*Zh23-sb*Td3i*Zh21+sb^2*Zh23)
 /Mb/2
,imah2F_c->EE/MW*Mc/SW*Zh23/tb/2/Mc/2
,imah1F_t->EE/MW*Mt/SW*Zh13/tb/2/Mt/2
,imah1F_b->EE/MW*Mb/SW/s2b*(cb*Td3i*Zh12-Td3r*Zh13-sb*Td3i*Zh11+sb^2*Zh13)
 /Mb/2
,imah1F_c->EE/MW*Mc/SW*Zh13/tb/2/Mc/2
,ah3S_b2->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*Zb21r^2*cb^2*MW^2*Zh31-6*Zb21i^2*cb^
 2*MW^2*Zh31-6*Zb21r^2*cb^2*MW^2*Zh31-4*SW^2*Zb11r^2*cb^2*MW^2*Zh31+4*SW^2*
 Zb21i^2*cb^2*MW^2*Zh31-2*SW^2*Zb21r^2*MW^2*Zh32*s2b+3*Zb21i^2*MW^2*Zh32*
 s2b+3*Zb21r^2*MW^2*Zh32*s2b+2*SW^2*Zb11r^2*MW^2*Zh32*s2b-2*SW^2*Zb21i^2*MW^
 2*Zh32*s2b-12*CW^2*Zb11r*Zb21r*Mb*Zh32*mur+12*CW^2*Zb11r*Zb21i*cb*Mb*Zh33*
 mur+12*CW^2*Zb11r*Zb21i*Mb*Zh32*mui+12*CW^2*Zb11r*Zb21r*cb*Mb*Zh33*mui+12*
 CW^2*Mb^2*Zh31+12*CW^2*Zb11r*Zb21r*Abr*Mb*Zh31+12*CW^2*Zb11r*Zb21i*sb*Abr*
 Mb*Zh33+12*CW^2*Zb11r*Zb21i*Abi*Mb*Zh31-12*CW^2*Zb11r*Zb21r*sb*Abi*Mb*Zh33)
 /12/MSb2^2
,ah3S_b1->1/CW^2*EE/MW/SW/cb*(4*SW^2*Zb21r^2*cb^2*MW^2*Zh31+6*Zb11r^2*cb^2*
 MW^2*Zh31-4*SW^2*Zb11r^2*cb^2*MW^2*Zh31+4*SW^2*Zb21i^2*cb^2*MW^2*Zh31-2*SW^
 2*Zb21r^2*MW^2*Zh32*s2b-3*Zb11r^2*MW^2*Zh32*s2b+2*SW^2*Zb11r^2*MW^2*Zh32*
 s2b-2*SW^2*Zb21i^2*MW^2*Zh32*s2b-12*CW^2*Zb11r*Zb21r*Mb*Zh32*mur+12*CW^2*
 Zb11r*Zb21i*cb*Mb*Zh33*mur+12*CW^2*Zb11r*Zb21r*cb*Mb*Zh33*mui+12*CW^2*
 Zb11r*Zb21i*Mb*Zh32*mui-12*CW^2*Mb^2*Zh31+12*CW^2*Zb11r*Zb21r*Abr*Mb*Zh31+
 12*CW^2*Zb11r*Zb21i*sb*Abr*Mb*Zh33-12*CW^2*Zb11r*Zb21r*sb*Abi*Mb*Zh33+12*
 CW^2*Zb11r*Zb21i*Abi*Mb*Zh31)/12/MSb1^2
,ah3S_t2->1/CW^2*EE/MW/SW/sb*(4*SW^2*Zt21r^2*MW^2*Zh31*s2b-3*Zt21i^2*MW^2*
 Zh31*s2b-3*Zt21r^2*MW^2*Zh31*s2b-4*SW^2*Zt11r^2*MW^2*Zh31*s2b+4*SW^2*Zt21i^
 2*MW^2*Zh31*s2b-8*SW^2*Zt21r^2*sb^2*MW^2*Zh32+6*Zt21i^2*sb^2*MW^2*Zh32+6*
 Zt21r^2*sb^2*MW^2*Zh32+8*SW^2*Zt11r^2*sb^2*MW^2*Zh32-8*SW^2*Zt21i^2*sb^2*
 MW^2*Zh32+12*CW^2*Zt11r*Zt21r*Mt*Zh31*mur-12*CW^2*Zt11r*Zt21i*sb*Mt*Zh33*
 mur-12*CW^2*Zt11r*Zt21i*Mt*Zh31*mui-12*CW^2*Zt11r*Zt21r*sb*Mt*Zh33*mui-12*
 CW^2*Mt^2*Zh32-12*CW^2*Zt11r*Zt21r*Atr*Mt*Zh32-12*CW^2*Zt11r*Zt21i*cb*Atr*
 Mt*Zh33-12*CW^2*Zt11r*Zt21i*Ati*Mt*Zh32+12*CW^2*Zt11r*Zt21r*cb*Ati*Mt*Zh33)
 /12/MSt2^2
,ah3S_t1->(-1/CW^2)*EE/MW/SW/sb*(4*SW^2*Zt21r^2*MW^2*Zh31*s2b+3*Zt11r^2*MW^
 2*Zh31*s2b-4*SW^2*Zt11r^2*MW^2*Zh31*s2b+4*SW^2*Zt21i^2*MW^2*Zh31*s2b-8*SW^
 2*Zt21r^2*sb^2*MW^2*Zh32-6*Zt11r^2*sb^2*MW^2*Zh32+8*SW^2*Zt11r^2*sb^2*MW^2*
 Zh32-8*SW^2*Zt21i^2*sb^2*MW^2*Zh32+12*CW^2*Zt11r*Zt21r*Mt*Zh31*mur-12*CW^2*
 Zt11r*Zt21i*sb*Mt*Zh33*mur-12*CW^2*Zt11r*Zt21r*sb*Mt*Zh33*mui-12*CW^2*
 Zt11r*Zt21i*Mt*Zh31*mui+12*CW^2*Mt^2*Zh32-12*CW^2*Zt11r*Zt21r*Atr*Mt*Zh32-
 12*CW^2*Zt11r*Zt21i*cb*Atr*Mt*Zh33+12*CW^2*Zt11r*Zt21r*cb*Ati*Mt*Zh33-12*
 CW^2*Zt11r*Zt21i*Ati*Mt*Zh32)/12/MSt1^2
,ah3S_sR->1/CW^2*EE/MW/SW/cb*(2*SW^2*cb^2*MW^2*Zh31-SW^2*MW^2*Zh32*s2b-6*CW^
 2*Mq^2*Zh31)/6/MSsR^2
,ah3S_sL->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*cb^2*MW^2*Zh31-6*cb^2*MW^2*Zh31-2*
 SW^2*MW^2*Zh32*s2b+3*MW^2*Zh32*s2b+12*CW^2*Mq^2*Zh31)/12/MSsL^2
,ah3S_dR->1/CW^2*EE/MW/SW/cb*(2*SW^2*cb^2*MW^2*Zh31-SW^2*MW^2*Zh32*s2b-6*CW^
 2*Mq^2*Zh31)/6/MSdR^2
,ah3S_dL->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*cb^2*MW^2*Zh31-6*cb^2*MW^2*Zh31-2*
 SW^2*MW^2*Zh32*s2b+3*MW^2*Zh32*s2b+12*CW^2*Mq^2*Zh31)/12/MSdL^2
,ah3S_cR->(-1/CW^2)*EE/MW/SW/sb*(SW^2*MW^2*Zh31*s2b-2*SW^2*sb^2*MW^2*Zh32+3*
 CW^2*Mc^2*Zh32)/3/MScR^2
,ah3S_cL->1/CW^2*EE/MW/SW/sb*(4*SW^2*MW^2*Zh31*s2b-3*MW^2*Zh31*s2b-8*SW^2*
 sb^2*MW^2*Zh32+6*sb^2*MW^2*Zh32-12*CW^2*Mc^2*Zh32)/12/MScL^2
,ah3S_uR->(-1/CW^2)*EE/MW/SW/sb*(SW^2*MW^2*Zh31*s2b-2*SW^2*sb^2*MW^2*Zh32+3*
 CW^2*Mq^2*Zh32)/3/MSuR^2
,ah3S_uL->1/CW^2*EE/MW/SW/sb*(4*SW^2*MW^2*Zh31*s2b-3*MW^2*Zh31*s2b-8*SW^2*
 sb^2*MW^2*Zh32+6*sb^2*MW^2*Zh32-12*CW^2*Mq^2*Zh32)/12/MSuL^2
,ah2S_b2->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*Zb21r^2*cb^2*MW^2*Zh21-6*Zb21i^2*cb^
 2*MW^2*Zh21-6*Zb21r^2*cb^2*MW^2*Zh21-4*SW^2*Zb11r^2*cb^2*MW^2*Zh21+4*SW^2*
 Zb21i^2*cb^2*MW^2*Zh21-2*SW^2*Zb21r^2*MW^2*Zh22*s2b+3*Zb21i^2*MW^2*Zh22*
 s2b+3*Zb21r^2*MW^2*Zh22*s2b+2*SW^2*Zb11r^2*MW^2*Zh22*s2b-2*SW^2*Zb21i^2*MW^
 2*Zh22*s2b-12*CW^2*Zb11r*Zb21r*Mb*Zh22*mur+12*CW^2*Zb11r*Zb21i*cb*Mb*Zh23*
 mur+12*CW^2*Zb11r*Zb21i*Mb*Zh22*mui+12*CW^2*Zb11r*Zb21r*cb*Mb*Zh23*mui+12*
 CW^2*Mb^2*Zh21+12*CW^2*Zb11r*Zb21r*Abr*Mb*Zh21+12*CW^2*Zb11r*Zb21i*sb*Abr*
 Mb*Zh23+12*CW^2*Zb11r*Zb21i*Abi*Mb*Zh21-12*CW^2*Zb11r*Zb21r*sb*Abi*Mb*Zh23)
 /12/MSb2^2
,ah2S_b1->1/CW^2*EE/MW/SW/cb*(4*SW^2*Zb21r^2*cb^2*MW^2*Zh21+6*Zb11r^2*cb^2*
 MW^2*Zh21-4*SW^2*Zb11r^2*cb^2*MW^2*Zh21+4*SW^2*Zb21i^2*cb^2*MW^2*Zh21-2*SW^
 2*Zb21r^2*MW^2*Zh22*s2b-3*Zb11r^2*MW^2*Zh22*s2b+2*SW^2*Zb11r^2*MW^2*Zh22*
 s2b-2*SW^2*Zb21i^2*MW^2*Zh22*s2b-12*CW^2*Zb11r*Zb21r*Mb*Zh22*mur+12*CW^2*
 Zb11r*Zb21i*cb*Mb*Zh23*mur+12*CW^2*Zb11r*Zb21r*cb*Mb*Zh23*mui+12*CW^2*
 Zb11r*Zb21i*Mb*Zh22*mui-12*CW^2*Mb^2*Zh21+12*CW^2*Zb11r*Zb21r*Abr*Mb*Zh21+
 12*CW^2*Zb11r*Zb21i*sb*Abr*Mb*Zh23-12*CW^2*Zb11r*Zb21r*sb*Abi*Mb*Zh23+12*
 CW^2*Zb11r*Zb21i*Abi*Mb*Zh21)/12/MSb1^2
,ah2S_t2->1/CW^2*EE/MW/SW/sb*(4*SW^2*Zt21r^2*MW^2*Zh21*s2b-3*Zt21i^2*MW^2*
 Zh21*s2b-3*Zt21r^2*MW^2*Zh21*s2b-4*SW^2*Zt11r^2*MW^2*Zh21*s2b+4*SW^2*Zt21i^
 2*MW^2*Zh21*s2b-8*SW^2*Zt21r^2*sb^2*MW^2*Zh22+6*Zt21i^2*sb^2*MW^2*Zh22+6*
 Zt21r^2*sb^2*MW^2*Zh22+8*SW^2*Zt11r^2*sb^2*MW^2*Zh22-8*SW^2*Zt21i^2*sb^2*
 MW^2*Zh22+12*CW^2*Zt11r*Zt21r*Mt*Zh21*mur-12*CW^2*Zt11r*Zt21i*sb*Mt*Zh23*
 mur-12*CW^2*Zt11r*Zt21i*Mt*Zh21*mui-12*CW^2*Zt11r*Zt21r*sb*Mt*Zh23*mui-12*
 CW^2*Mt^2*Zh22-12*CW^2*Zt11r*Zt21r*Atr*Mt*Zh22-12*CW^2*Zt11r*Zt21i*cb*Atr*
 Mt*Zh23-12*CW^2*Zt11r*Zt21i*Ati*Mt*Zh22+12*CW^2*Zt11r*Zt21r*cb*Ati*Mt*Zh23)
 /12/MSt2^2
,ah2S_t1->(-1/CW^2)*EE/MW/SW/sb*(4*SW^2*Zt21r^2*MW^2*Zh21*s2b+3*Zt11r^2*MW^
 2*Zh21*s2b-4*SW^2*Zt11r^2*MW^2*Zh21*s2b+4*SW^2*Zt21i^2*MW^2*Zh21*s2b-8*SW^
 2*Zt21r^2*sb^2*MW^2*Zh22-6*Zt11r^2*sb^2*MW^2*Zh22+8*SW^2*Zt11r^2*sb^2*MW^2*
 Zh22-8*SW^2*Zt21i^2*sb^2*MW^2*Zh22+12*CW^2*Zt11r*Zt21r*Mt*Zh21*mur-12*CW^2*
 Zt11r*Zt21i*sb*Mt*Zh23*mur-12*CW^2*Zt11r*Zt21r*sb*Mt*Zh23*mui-12*CW^2*
 Zt11r*Zt21i*Mt*Zh21*mui+12*CW^2*Mt^2*Zh22-12*CW^2*Zt11r*Zt21r*Atr*Mt*Zh22-
 12*CW^2*Zt11r*Zt21i*cb*Atr*Mt*Zh23+12*CW^2*Zt11r*Zt21r*cb*Ati*Mt*Zh23-12*
 CW^2*Zt11r*Zt21i*Ati*Mt*Zh22)/12/MSt1^2
,ah2S_sR->1/CW^2*EE/MW/SW/cb*(2*SW^2*cb^2*MW^2*Zh21-SW^2*MW^2*Zh22*s2b-6*CW^
 2*Mq^2*Zh21)/6/MSsR^2
,ah2S_sL->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*cb^2*MW^2*Zh21-6*cb^2*MW^2*Zh21-2*
 SW^2*MW^2*Zh22*s2b+3*MW^2*Zh22*s2b+12*CW^2*Mq^2*Zh21)/12/MSsL^2
,ah2S_dR->1/CW^2*EE/MW/SW/cb*(2*SW^2*cb^2*MW^2*Zh21-SW^2*MW^2*Zh22*s2b-6*CW^
 2*Mq^2*Zh21)/6/MSdR^2
,ah2S_dL->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*cb^2*MW^2*Zh21-6*cb^2*MW^2*Zh21-2*
 SW^2*MW^2*Zh22*s2b+3*MW^2*Zh22*s2b+12*CW^2*Mq^2*Zh21)/12/MSdL^2
,ah2S_cR->(-1/CW^2)*EE/MW/SW/sb*(SW^2*MW^2*Zh21*s2b-2*SW^2*sb^2*MW^2*Zh22+3*
 CW^2*Mc^2*Zh22)/3/MScR^2
,ah2S_cL->1/CW^2*EE/MW/SW/sb*(4*SW^2*MW^2*Zh21*s2b-3*MW^2*Zh21*s2b-8*SW^2*
 sb^2*MW^2*Zh22+6*sb^2*MW^2*Zh22-12*CW^2*Mc^2*Zh22)/12/MScL^2
,ah2S_uR->(-1/CW^2)*EE/MW/SW/sb*(SW^2*MW^2*Zh21*s2b-2*SW^2*sb^2*MW^2*Zh22+3*
 CW^2*Mq^2*Zh22)/3/MSuR^2
,ah2S_uL->1/CW^2*EE/MW/SW/sb*(4*SW^2*MW^2*Zh21*s2b-3*MW^2*Zh21*s2b-8*SW^2*
 sb^2*MW^2*Zh22+6*sb^2*MW^2*Zh22-12*CW^2*Mq^2*Zh22)/12/MSuL^2
,ah1S_b2->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*Zb21r^2*cb^2*MW^2*Zh11-6*Zb21i^2*cb^
 2*MW^2*Zh11-6*Zb21r^2*cb^2*MW^2*Zh11-4*SW^2*Zb11r^2*cb^2*MW^2*Zh11+4*SW^2*
 Zb21i^2*cb^2*MW^2*Zh11-2*SW^2*Zb21r^2*MW^2*Zh12*s2b+3*Zb21i^2*MW^2*Zh12*
 s2b+3*Zb21r^2*MW^2*Zh12*s2b+2*SW^2*Zb11r^2*MW^2*Zh12*s2b-2*SW^2*Zb21i^2*MW^
 2*Zh12*s2b-12*CW^2*Zb11r*Zb21r*Mb*Zh12*mur+12*CW^2*Zb11r*Zb21i*cb*Mb*Zh13*
 mur+12*CW^2*Zb11r*Zb21i*Mb*Zh12*mui+12*CW^2*Zb11r*Zb21r*cb*Mb*Zh13*mui+12*
 CW^2*Mb^2*Zh11+12*CW^2*Zb11r*Zb21r*Abr*Mb*Zh11+12*CW^2*Zb11r*Zb21i*sb*Abr*
 Mb*Zh13+12*CW^2*Zb11r*Zb21i*Abi*Mb*Zh11-12*CW^2*Zb11r*Zb21r*sb*Abi*Mb*Zh13)
 /12/MSb2^2
,ah1S_b1->1/CW^2*EE/MW/SW/cb*(4*SW^2*Zb21r^2*cb^2*MW^2*Zh11+6*Zb11r^2*cb^2*
 MW^2*Zh11-4*SW^2*Zb11r^2*cb^2*MW^2*Zh11+4*SW^2*Zb21i^2*cb^2*MW^2*Zh11-2*SW^
 2*Zb21r^2*MW^2*Zh12*s2b-3*Zb11r^2*MW^2*Zh12*s2b+2*SW^2*Zb11r^2*MW^2*Zh12*
 s2b-2*SW^2*Zb21i^2*MW^2*Zh12*s2b-12*CW^2*Zb11r*Zb21r*Mb*Zh12*mur+12*CW^2*
 Zb11r*Zb21i*cb*Mb*Zh13*mur+12*CW^2*Zb11r*Zb21r*cb*Mb*Zh13*mui+12*CW^2*
 Zb11r*Zb21i*Mb*Zh12*mui-12*CW^2*Mb^2*Zh11+12*CW^2*Zb11r*Zb21r*Abr*Mb*Zh11+
 12*CW^2*Zb11r*Zb21i*sb*Abr*Mb*Zh13-12*CW^2*Zb11r*Zb21r*sb*Abi*Mb*Zh13+12*
 CW^2*Zb11r*Zb21i*Abi*Mb*Zh11)/12/MSb1^2
,ah1S_t2->1/CW^2*EE/MW/SW/sb*(4*SW^2*Zt21r^2*MW^2*Zh11*s2b-3*Zt21i^2*MW^2*
 Zh11*s2b-3*Zt21r^2*MW^2*Zh11*s2b-4*SW^2*Zt11r^2*MW^2*Zh11*s2b+4*SW^2*Zt21i^
 2*MW^2*Zh11*s2b-8*SW^2*Zt21r^2*sb^2*MW^2*Zh12+6*Zt21i^2*sb^2*MW^2*Zh12+6*
 Zt21r^2*sb^2*MW^2*Zh12+8*SW^2*Zt11r^2*sb^2*MW^2*Zh12-8*SW^2*Zt21i^2*sb^2*
 MW^2*Zh12+12*CW^2*Zt11r*Zt21r*Mt*Zh11*mur-12*CW^2*Zt11r*Zt21i*sb*Mt*Zh13*
 mur-12*CW^2*Zt11r*Zt21i*Mt*Zh11*mui-12*CW^2*Zt11r*Zt21r*sb*Mt*Zh13*mui-12*
 CW^2*Mt^2*Zh12-12*CW^2*Zt11r*Zt21r*Atr*Mt*Zh12-12*CW^2*Zt11r*Zt21i*cb*Atr*
 Mt*Zh13-12*CW^2*Zt11r*Zt21i*Ati*Mt*Zh12+12*CW^2*Zt11r*Zt21r*cb*Ati*Mt*Zh13)
 /12/MSt2^2
,ah1S_t1->(-1/CW^2)*EE/MW/SW/sb*(4*SW^2*Zt21r^2*MW^2*Zh11*s2b+3*Zt11r^2*MW^
 2*Zh11*s2b-4*SW^2*Zt11r^2*MW^2*Zh11*s2b+4*SW^2*Zt21i^2*MW^2*Zh11*s2b-8*SW^
 2*Zt21r^2*sb^2*MW^2*Zh12-6*Zt11r^2*sb^2*MW^2*Zh12+8*SW^2*Zt11r^2*sb^2*MW^2*
 Zh12-8*SW^2*Zt21i^2*sb^2*MW^2*Zh12+12*CW^2*Zt11r*Zt21r*Mt*Zh11*mur-12*CW^2*
 Zt11r*Zt21i*sb*Mt*Zh13*mur-12*CW^2*Zt11r*Zt21r*sb*Mt*Zh13*mui-12*CW^2*
 Zt11r*Zt21i*Mt*Zh11*mui+12*CW^2*Mt^2*Zh12-12*CW^2*Zt11r*Zt21r*Atr*Mt*Zh12-
 12*CW^2*Zt11r*Zt21i*cb*Atr*Mt*Zh13+12*CW^2*Zt11r*Zt21r*cb*Ati*Mt*Zh13-12*
 CW^2*Zt11r*Zt21i*Ati*Mt*Zh12)/12/MSt1^2
,ah1S_sR->1/CW^2*EE/MW/SW/cb*(2*SW^2*cb^2*MW^2*Zh11-SW^2*MW^2*Zh12*s2b-6*CW^
 2*Mq^2*Zh11)/6/MSsR^2
,ah1S_sL->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*cb^2*MW^2*Zh11-6*cb^2*MW^2*Zh11-2*
 SW^2*MW^2*Zh12*s2b+3*MW^2*Zh12*s2b+12*CW^2*Mq^2*Zh11)/12/MSsL^2
,ah1S_dR->1/CW^2*EE/MW/SW/cb*(2*SW^2*cb^2*MW^2*Zh11-SW^2*MW^2*Zh12*s2b-6*CW^
 2*Mq^2*Zh11)/6/MSdR^2
,ah1S_dL->(-1/CW^2)*EE/MW/SW/cb*(4*SW^2*cb^2*MW^2*Zh11-6*cb^2*MW^2*Zh11-2*
 SW^2*MW^2*Zh12*s2b+3*MW^2*Zh12*s2b+12*CW^2*Mq^2*Zh11)/12/MSdL^2
,ah1S_cR->(-1/CW^2)*EE/MW/SW/sb*(SW^2*MW^2*Zh11*s2b-2*SW^2*sb^2*MW^2*Zh12+3*
 CW^2*Mc^2*Zh12)/3/MScR^2
,ah1S_cL->1/CW^2*EE/MW/SW/sb*(4*SW^2*MW^2*Zh11*s2b-3*MW^2*Zh11*s2b-8*SW^2*
 sb^2*MW^2*Zh12+6*sb^2*MW^2*Zh12-12*CW^2*Mc^2*Zh12)/12/MScL^2
,ah1S_uR->(-1/CW^2)*EE/MW/SW/sb*(SW^2*MW^2*Zh11*s2b-2*SW^2*sb^2*MW^2*Zh12+3*
 CW^2*Mq^2*Zh12)/3/MSuR^2
,ah1S_uL->1/CW^2*EE/MW/SW/sb*(4*SW^2*MW^2*Zh11*s2b-3*MW^2*Zh11*s2b-8*SW^2*
 sb^2*MW^2*Zh12+6*sb^2*MW^2*Zh12-12*CW^2*Mq^2*Zh12)/12/MSuL^2
,ah3F_t->(-EE)/MW*Mt/SW*Zh32/sb/2/Mt
,ah3F_b->(-EE)/MW*Mb/SW/s2b*(sb*Zh31-sb*Td3r*Zh31+cb*Td3r*Zh32+Td3i*Zh33)/Mb
,ah3F_c->(-EE)/MW*Mc/SW*Zh32/sb/2/Mc
,ah2F_t->(-EE)/MW*Mt/SW*Zh22/sb/2/Mt
,ah2F_b->(-EE)/MW*Mb/SW/s2b*(sb*Zh21-sb*Td3r*Zh21+cb*Td3r*Zh22+Td3i*Zh23)/Mb
,ah2F_c->(-EE)/MW*Mc/SW*Zh22/sb/2/Mc
,ah1F_t->(-EE)/MW*Mt/SW*Zh12/sb/2/Mt
,ah1F_b->(-EE)/MW*Mb/SW/s2b*(sb*Zh11-sb*Td3r*Zh11+cb*Td3r*Zh12+Td3i*Zh13)/Mb
,ah1F_c->(-EE)/MW*Mc/SW*Zh12/sb/2/Mc
,MSsR->Sqrt[Md2^2+MW^2/CW^2*(-1/3*SW^2)*c2b]
,MSdR->Sqrt[Md2^2+MW^2/CW^2*(-1/3*SW^2)*c2b]
,MSsL->Sqrt[Mq2^2-MW^2/CW^2*(1/2-1/3*SW^2)*c2b]
,MSdL->Sqrt[Mq2^2-MW^2/CW^2*(1/2-1/3*SW^2)*c2b]
,MScR->Sqrt[Mu2^2+MW^2/CW^2*(2/3*SW^2)*c2b]
,MSuR->Sqrt[Mu2^2+MW^2/CW^2*(2/3*SW^2)*c2b]
,MScL->Sqrt[Mq2^2+MW^2/CW^2*(1/2-2/3*SW^2)*c2b]
,MSuL->Sqrt[Mq2^2+MW^2/CW^2*(1/2-2/3*SW^2)*c2b]
,Ati->aAt*Sin[fiAt*PI/180]
,Atr->aAt*Cos[fiAt*PI/180]
,Abi->aAb*Sin[fiAb*PI/180]
,Abr->aAb*Cos[fiAb*PI/180]
,mui->aMu*Sin[fiMu*PI/180]
,mur->aMu*Cos[fiMu*PI/180]
,PI->3.14159*1
,c2b->cb^2-sb^2
,s2b->2*sb*cb
,cb->Sqrt[1-sb^2]
,sb->tb/Sqrt[1+tb^2]
,MW->MZ*CW
,CW->Sqrt[1-SW^2]
              };

inParticles = {"G","G"}
outParticles = {"G","G"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =0^2;
p2/: SC[p2,p2] =0^2;
p3/: SC[p3,p3] =0^2;
p2/: SC[p2,p3] = -1*(0^2-0^2-0^2-0^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((9*GG^4)/(4));
numerator =(4*SC[p1,p3]^2-4*SC[p1,p3]*SC[p1,p2]+SC[p1,p2]^2);
denominator =(propDen[-p1-p2,0,0]^2);

addToSum;

(*
  Diagram  2 in subprocess 1
*)
totFactor = ((9*GG^4)/(2));
numerator =(2*SC[p1,p3]^2+SC[p1,p3]*SC[p1,p2]-SC[p1,p2]^2);
denominator =(propDen[-p1-p2,0,0]*propDen[-p1+p4,0,0]);

addToSum;

(*
  Diagram  4 in subprocess 1
*)
totFactor = ((24*GG^2*Rqcdh1^2)/(1));
numerator =(4*SC[p1,p3]^3*imLGGh1^2+SC[p1,p3]^3*LGGh1^2-16*SC[p1,p3]^2*
 SC[p1,p2]*imLGGh1^2-3*SC[p1,p3]^2*SC[p1,p2]*LGGh1^2+12*SC[p1,p3]*SC[p1,p2]^
 2*imLGGh1^2+3*SC[p1,p3]*SC[p1,p2]^2*LGGh1^2-SC[p1,p2]^3*LGGh1^2);
denominator =(propDen[-p1-p2,0,0]*propDen[-p1+p4,Mh1,0]);

addToSum;

(*
  Diagram  6 in subprocess 1
*)
totFactor = ((24*GG^2*Rqcdh2^2)/(1));
numerator =(4*SC[p1,p3]^3*imLGGh2^2+SC[p1,p3]^3*LGGh2^2-16*SC[p1,p3]^2*
 SC[p1,p2]*imLGGh2^2-3*SC[p1,p3]^2*SC[p1,p2]*LGGh2^2+12*SC[p1,p3]*SC[p1,p2]^
 2*imLGGh2^2+3*SC[p1,p3]*SC[p1,p2]^2*LGGh2^2-SC[p1,p2]^3*LGGh2^2);
denominator =(propDen[-p1-p2,0,0]*propDen[-p1+p4,Mh2,0]);

addToSum;

(*
  Diagram  8 in subprocess 1
*)
totFactor = ((24*GG^2*Rqcdh3^2)/(1));
numerator =(4*SC[p1,p3]^3*imLGGh3^2+SC[p1,p3]^3*LGGh3^2-16*SC[p1,p3]^2*
 SC[p1,p2]*imLGGh3^2-3*SC[p1,p3]^2*SC[p1,p2]*LGGh3^2+12*SC[p1,p3]*SC[p1,p2]^
 2*imLGGh3^2+3*SC[p1,p3]*SC[p1,p2]^2*LGGh3^2-SC[p1,p2]^3*LGGh3^2);
denominator =(propDen[-p1-p2,0,0]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  9 in subprocess 1
*)
totFactor = ((9*GG^4)/(2));
numerator =(SC[p1,p3]^2-4*SC[p1,p3]*SC[p1,p2]+4*SC[p1,p2]^2);
denominator =(propDen[-p1+p3,0,0]^2);

addToSum;

(*
  Diagram  10 in subprocess 1
*)
totFactor = ((-9*GG^4)/(4));
numerator =(SC[p1,p3]^2-SC[p1,p3]*SC[p1,p2]-2*SC[p1,p2]^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1+p4,0,0]);

addToSum;

(*
  Diagram  11 in subprocess 1
*)
totFactor = ((-24*GG^2*Rqcdh1^2)/(1));
numerator =(8*SC[p1,p3]^2*SC[p1,p2]*imLGGh1^2-4*SC[p1,p3]*SC[p1,p2]^2*
 imLGGh1^2-4*SC[p1,p2]^3*imLGGh1^2-SC[p1,p2]^3*LGGh1^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1-p2,Mh1,wh1]);

addToSum;

(*
  Diagram  13 in subprocess 1
*)
totFactor = ((24*GG^2*Rqcdh1^2)/(1));
numerator =(SC[p1,p3]^3*LGGh1^2-12*SC[p1,p3]^2*SC[p1,p2]*imLGGh1^2-3*
 SC[p1,p3]^2*SC[p1,p2]*LGGh1^2+16*SC[p1,p3]*SC[p1,p2]^2*imLGGh1^2+3*
 SC[p1,p3]*SC[p1,p2]^2*LGGh1^2-4*SC[p1,p2]^3*imLGGh1^2-SC[p1,p2]^3*LGGh1^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1+p4,Mh1,0]);

addToSum;

(*
  Diagram  14 in subprocess 1
*)
totFactor = ((-24*GG^2*Rqcdh2^2)/(1));
numerator =(8*SC[p1,p3]^2*SC[p1,p2]*imLGGh2^2-4*SC[p1,p3]*SC[p1,p2]^2*
 imLGGh2^2-4*SC[p1,p2]^3*imLGGh2^2-SC[p1,p2]^3*LGGh2^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1-p2,Mh2,wh2]);

addToSum;

(*
  Diagram  16 in subprocess 1
*)
totFactor = ((24*GG^2*Rqcdh2^2)/(1));
numerator =(SC[p1,p3]^3*LGGh2^2-12*SC[p1,p3]^2*SC[p1,p2]*imLGGh2^2-3*
 SC[p1,p3]^2*SC[p1,p2]*LGGh2^2+16*SC[p1,p3]*SC[p1,p2]^2*imLGGh2^2+3*
 SC[p1,p3]*SC[p1,p2]^2*LGGh2^2-4*SC[p1,p2]^3*imLGGh2^2-SC[p1,p2]^3*LGGh2^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1+p4,Mh2,0]);

addToSum;

(*
  Diagram  17 in subprocess 1
*)
totFactor = ((-24*GG^2*Rqcdh3^2)/(1));
numerator =(8*SC[p1,p3]^2*SC[p1,p2]*imLGGh3^2-4*SC[p1,p3]*SC[p1,p2]^2*
 imLGGh3^2-4*SC[p1,p2]^3*imLGGh3^2-SC[p1,p2]^3*LGGh3^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1-p2,Mh3,wh3]);

addToSum;

(*
  Diagram  19 in subprocess 1
*)
totFactor = ((24*GG^2*Rqcdh3^2)/(1));
numerator =(SC[p1,p3]^3*LGGh3^2-12*SC[p1,p3]^2*SC[p1,p2]*imLGGh3^2-3*
 SC[p1,p3]^2*SC[p1,p2]*LGGh3^2+16*SC[p1,p3]*SC[p1,p2]^2*imLGGh3^2+3*
 SC[p1,p3]*SC[p1,p2]^2*LGGh3^2-4*SC[p1,p2]^3*imLGGh3^2-SC[p1,p2]^3*LGGh3^2);
denominator =(propDen[-p1+p3,0,0]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  20 in subprocess 1
*)
totFactor = ((128*Rqcdh1^4)/(1));
numerator =(16*SC[p1,p2]^4*imLGGh1^4+8*SC[p1,p2]^4*imLGGh1^2*LGGh1^2+
 SC[p1,p2]^4*LGGh1^4);
denominator =(propDen[-p1-p2,Mh1,wh1]^2);

addToSum;

(*
  Diagram  21 in subprocess 1
*)
totFactor = ((32*Rqcdh1^4)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh1^4+8*SC[p1,p3]^2*SC[p1,p2]^2*
 imLGGh1^2*LGGh1^2+SC[p1,p3]^2*SC[p1,p2]^2*LGGh1^4-32*SC[p1,p3]*SC[p1,p2]^3*
 imLGGh1^4-16*SC[p1,p3]*SC[p1,p2]^3*imLGGh1^2*LGGh1^2-2*SC[p1,p3]*SC[p1,p2]^
 3*LGGh1^4+16*SC[p1,p2]^4*imLGGh1^4+8*SC[p1,p2]^4*imLGGh1^2*LGGh1^2+
 SC[p1,p2]^4*LGGh1^4);
denominator =(propDen[-p1-p2,Mh1,wh1]*propDen[-p1+p4,Mh1,0]);

addToSum;

(*
  Diagram  22 in subprocess 1
*)
totFactor = ((256*Rqcdh2^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p2]^4*imLGGh2^2*imLGGh1^2+8*SC[p1,p2]^4*imLGGh2*
 imLGGh1*LGGh2*LGGh1+SC[p1,p2]^4*LGGh2^2*LGGh1^2);
denominator =(propDen[-p1-p2,Mh1,wh1]*propDen[-p1-p2,Mh2,wh2]);

addToSum;

(*
  Diagram  23 in subprocess 1
*)
totFactor = ((32*Rqcdh2^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^2*imLGGh1^2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh2^2*LGGh1^2+16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2*imLGGh1*
 LGGh2*LGGh1-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh1^2*LGGh2^2+SC[p1,p3]^2*
 SC[p1,p2]^2*LGGh2^2*LGGh1^2-32*SC[p1,p3]*SC[p1,p2]^3*imLGGh2^2*imLGGh1^2+8*
 SC[p1,p3]*SC[p1,p2]^3*imLGGh2^2*LGGh1^2-32*SC[p1,p3]*SC[p1,p2]^3*imLGGh2*
 imLGGh1*LGGh2*LGGh1+8*SC[p1,p3]*SC[p1,p2]^3*imLGGh1^2*LGGh2^2-2*SC[p1,p3]*
 SC[p1,p2]^3*LGGh2^2*LGGh1^2+16*SC[p1,p2]^4*imLGGh2^2*imLGGh1^2-4*SC[p1,p2]^
 4*imLGGh2^2*LGGh1^2+16*SC[p1,p2]^4*imLGGh2*imLGGh1*LGGh2*LGGh1-4*SC[p1,p2]^
 4*imLGGh1^2*LGGh2^2+SC[p1,p2]^4*LGGh2^2*LGGh1^2);
denominator =(propDen[-p1-p2,Mh1,wh1]*propDen[-p1+p4,Mh2,0]);

addToSum;

(*
  Diagram  24 in subprocess 1
*)
totFactor = ((256*Rqcdh3^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p2]^4*imLGGh3^2*imLGGh1^2+8*SC[p1,p2]^4*imLGGh3*
 imLGGh1*LGGh3*LGGh1+SC[p1,p2]^4*LGGh3^2*LGGh1^2);
denominator =(propDen[-p1-p2,Mh1,wh1]*propDen[-p1-p2,Mh3,wh3]);

addToSum;

(*
  Diagram  25 in subprocess 1
*)
totFactor = ((32*Rqcdh3^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*imLGGh1^2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh3^2*LGGh1^2+16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3*imLGGh1*
 LGGh3*LGGh1-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh1^2*LGGh3^2+SC[p1,p3]^2*
 SC[p1,p2]^2*LGGh3^2*LGGh1^2-32*SC[p1,p3]*SC[p1,p2]^3*imLGGh3^2*imLGGh1^2+8*
 SC[p1,p3]*SC[p1,p2]^3*imLGGh3^2*LGGh1^2-32*SC[p1,p3]*SC[p1,p2]^3*imLGGh3*
 imLGGh1*LGGh3*LGGh1+8*SC[p1,p3]*SC[p1,p2]^3*imLGGh1^2*LGGh3^2-2*SC[p1,p3]*
 SC[p1,p2]^3*LGGh3^2*LGGh1^2+16*SC[p1,p2]^4*imLGGh3^2*imLGGh1^2-4*SC[p1,p2]^
 4*imLGGh3^2*LGGh1^2+16*SC[p1,p2]^4*imLGGh3*imLGGh1*LGGh3*LGGh1-4*SC[p1,p2]^
 4*imLGGh1^2*LGGh3^2+SC[p1,p2]^4*LGGh3^2*LGGh1^2);
denominator =(propDen[-p1-p2,Mh1,wh1]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  26 in subprocess 1
*)
totFactor = ((256*Rqcdh1^4)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh1^4+8*SC[p1,p3]^4*imLGGh1^2*LGGh1^2+
 SC[p1,p3]^4*LGGh1^4);
denominator =(propDen[-p1+p3,Mh1,0]^2);

addToSum;

(*
  Diagram  27 in subprocess 1
*)
totFactor = ((16*Rqcdh1^4)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh1^4+8*SC[p1,p3]^4*imLGGh1^2*LGGh1^2+
 SC[p1,p3]^4*LGGh1^4-32*SC[p1,p3]^3*SC[p1,p2]*imLGGh1^4-16*SC[p1,p3]^3*
 SC[p1,p2]*imLGGh1^2*LGGh1^2-2*SC[p1,p3]^3*SC[p1,p2]*LGGh1^4+16*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh1^4+8*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh1^2*LGGh1^2+
 SC[p1,p3]^2*SC[p1,p2]^2*LGGh1^4);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1+p4,Mh1,0]);

addToSum;

(*
  Diagram  28 in subprocess 1
*)
totFactor = ((32*Rqcdh2^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^2*imLGGh1^2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh2^2*LGGh1^2+16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2*imLGGh1*
 LGGh2*LGGh1-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh1^2*LGGh2^2+SC[p1,p3]^2*
 SC[p1,p2]^2*LGGh2^2*LGGh1^2);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1-p2,Mh2,wh2]);

addToSum;

(*
  Diagram  29 in subprocess 1
*)
totFactor = ((512*Rqcdh2^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh2^2*imLGGh1^2+8*SC[p1,p3]^4*imLGGh2*
 imLGGh1*LGGh2*LGGh1+SC[p1,p3]^4*LGGh2^2*LGGh1^2);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1+p3,Mh2,0]);

addToSum;

(*
  Diagram  30 in subprocess 1
*)
totFactor = ((32*Rqcdh2^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh2^2*imLGGh1^2-4*SC[p1,p3]^4*imLGGh2^2*
 LGGh1^2+16*SC[p1,p3]^4*imLGGh2*imLGGh1*LGGh2*LGGh1-4*SC[p1,p3]^4*imLGGh1^2*
 LGGh2^2+SC[p1,p3]^4*LGGh2^2*LGGh1^2-32*SC[p1,p3]^3*SC[p1,p2]*imLGGh2^2*
 imLGGh1^2+8*SC[p1,p3]^3*SC[p1,p2]*imLGGh2^2*LGGh1^2-32*SC[p1,p3]^3*
 SC[p1,p2]*imLGGh2*imLGGh1*LGGh2*LGGh1+8*SC[p1,p3]^3*SC[p1,p2]*imLGGh1^2*
 LGGh2^2-2*SC[p1,p3]^3*SC[p1,p2]*LGGh2^2*LGGh1^2+16*SC[p1,p3]^2*SC[p1,p2]^2*
 imLGGh2^2*imLGGh1^2-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^2*LGGh1^2+16*
 SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2*imLGGh1*LGGh2*LGGh1-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh1^2*LGGh2^2+SC[p1,p3]^2*SC[p1,p2]^2*LGGh2^2*LGGh1^2);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1+p4,Mh2,0]);

addToSum;

(*
  Diagram  31 in subprocess 1
*)
totFactor = ((32*Rqcdh3^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*imLGGh1^2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh3^2*LGGh1^2+16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3*imLGGh1*
 LGGh3*LGGh1-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh1^2*LGGh3^2+SC[p1,p3]^2*
 SC[p1,p2]^2*LGGh3^2*LGGh1^2);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1-p2,Mh3,wh3]);

addToSum;

(*
  Diagram  32 in subprocess 1
*)
totFactor = ((512*Rqcdh3^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh3^2*imLGGh1^2+8*SC[p1,p3]^4*imLGGh3*
 imLGGh1*LGGh3*LGGh1+SC[p1,p3]^4*LGGh3^2*LGGh1^2);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1+p3,Mh3,0]);

addToSum;

(*
  Diagram  33 in subprocess 1
*)
totFactor = ((32*Rqcdh3^2*Rqcdh1^2)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh3^2*imLGGh1^2-4*SC[p1,p3]^4*imLGGh3^2*
 LGGh1^2+16*SC[p1,p3]^4*imLGGh3*imLGGh1*LGGh3*LGGh1-4*SC[p1,p3]^4*imLGGh1^2*
 LGGh3^2+SC[p1,p3]^4*LGGh3^2*LGGh1^2-32*SC[p1,p3]^3*SC[p1,p2]*imLGGh3^2*
 imLGGh1^2+8*SC[p1,p3]^3*SC[p1,p2]*imLGGh3^2*LGGh1^2-32*SC[p1,p3]^3*
 SC[p1,p2]*imLGGh3*imLGGh1*LGGh3*LGGh1+8*SC[p1,p3]^3*SC[p1,p2]*imLGGh1^2*
 LGGh3^2-2*SC[p1,p3]^3*SC[p1,p2]*LGGh3^2*LGGh1^2+16*SC[p1,p3]^2*SC[p1,p2]^2*
 imLGGh3^2*imLGGh1^2-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*LGGh1^2+16*
 SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3*imLGGh1*LGGh3*LGGh1-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh1^2*LGGh3^2+SC[p1,p3]^2*SC[p1,p2]^2*LGGh3^2*LGGh1^2);
denominator =(propDen[-p1+p3,Mh1,0]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  34 in subprocess 1
*)
totFactor = ((128*Rqcdh2^4)/(1));
numerator =(16*SC[p1,p2]^4*imLGGh2^4+8*SC[p1,p2]^4*imLGGh2^2*LGGh2^2+
 SC[p1,p2]^4*LGGh2^4);
denominator =(propDen[-p1-p2,Mh2,wh2]^2);

addToSum;

(*
  Diagram  35 in subprocess 1
*)
totFactor = ((32*Rqcdh2^4)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^4+8*SC[p1,p3]^2*SC[p1,p2]^2*
 imLGGh2^2*LGGh2^2+SC[p1,p3]^2*SC[p1,p2]^2*LGGh2^4-32*SC[p1,p3]*SC[p1,p2]^3*
 imLGGh2^4-16*SC[p1,p3]*SC[p1,p2]^3*imLGGh2^2*LGGh2^2-2*SC[p1,p3]*SC[p1,p2]^
 3*LGGh2^4+16*SC[p1,p2]^4*imLGGh2^4+8*SC[p1,p2]^4*imLGGh2^2*LGGh2^2+
 SC[p1,p2]^4*LGGh2^4);
denominator =(propDen[-p1-p2,Mh2,wh2]*propDen[-p1+p4,Mh2,0]);

addToSum;

(*
  Diagram  36 in subprocess 1
*)
totFactor = ((256*Rqcdh3^2*Rqcdh2^2)/(1));
numerator =(16*SC[p1,p2]^4*imLGGh3^2*imLGGh2^2+8*SC[p1,p2]^4*imLGGh3*
 imLGGh2*LGGh3*LGGh2+SC[p1,p2]^4*LGGh3^2*LGGh2^2);
denominator =(propDen[-p1-p2,Mh2,wh2]*propDen[-p1-p2,Mh3,wh3]);

addToSum;

(*
  Diagram  37 in subprocess 1
*)
totFactor = ((32*Rqcdh3^2*Rqcdh2^2)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*imLGGh2^2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh3^2*LGGh2^2+16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3*imLGGh2*
 LGGh3*LGGh2-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^2*LGGh3^2+SC[p1,p3]^2*
 SC[p1,p2]^2*LGGh3^2*LGGh2^2-32*SC[p1,p3]*SC[p1,p2]^3*imLGGh3^2*imLGGh2^2+8*
 SC[p1,p3]*SC[p1,p2]^3*imLGGh3^2*LGGh2^2-32*SC[p1,p3]*SC[p1,p2]^3*imLGGh3*
 imLGGh2*LGGh3*LGGh2+8*SC[p1,p3]*SC[p1,p2]^3*imLGGh2^2*LGGh3^2-2*SC[p1,p3]*
 SC[p1,p2]^3*LGGh3^2*LGGh2^2+16*SC[p1,p2]^4*imLGGh3^2*imLGGh2^2-4*SC[p1,p2]^
 4*imLGGh3^2*LGGh2^2+16*SC[p1,p2]^4*imLGGh3*imLGGh2*LGGh3*LGGh2-4*SC[p1,p2]^
 4*imLGGh2^2*LGGh3^2+SC[p1,p2]^4*LGGh3^2*LGGh2^2);
denominator =(propDen[-p1-p2,Mh2,wh2]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  38 in subprocess 1
*)
totFactor = ((256*Rqcdh2^4)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh2^4+8*SC[p1,p3]^4*imLGGh2^2*LGGh2^2+
 SC[p1,p3]^4*LGGh2^4);
denominator =(propDen[-p1+p3,Mh2,0]^2);

addToSum;

(*
  Diagram  39 in subprocess 1
*)
totFactor = ((16*Rqcdh2^4)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh2^4+8*SC[p1,p3]^4*imLGGh2^2*LGGh2^2+
 SC[p1,p3]^4*LGGh2^4-32*SC[p1,p3]^3*SC[p1,p2]*imLGGh2^4-16*SC[p1,p3]^3*
 SC[p1,p2]*imLGGh2^2*LGGh2^2-2*SC[p1,p3]^3*SC[p1,p2]*LGGh2^4+16*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh2^4+8*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^2*LGGh2^2+
 SC[p1,p3]^2*SC[p1,p2]^2*LGGh2^4);
denominator =(propDen[-p1+p3,Mh2,0]*propDen[-p1+p4,Mh2,0]);

addToSum;

(*
  Diagram  40 in subprocess 1
*)
totFactor = ((32*Rqcdh3^2*Rqcdh2^2)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*imLGGh2^2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh3^2*LGGh2^2+16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3*imLGGh2*
 LGGh3*LGGh2-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh2^2*LGGh3^2+SC[p1,p3]^2*
 SC[p1,p2]^2*LGGh3^2*LGGh2^2);
denominator =(propDen[-p1+p3,Mh2,0]*propDen[-p1-p2,Mh3,wh3]);

addToSum;

(*
  Diagram  41 in subprocess 1
*)
totFactor = ((512*Rqcdh3^2*Rqcdh2^2)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh3^2*imLGGh2^2+8*SC[p1,p3]^4*imLGGh3*
 imLGGh2*LGGh3*LGGh2+SC[p1,p3]^4*LGGh3^2*LGGh2^2);
denominator =(propDen[-p1+p3,Mh2,0]*propDen[-p1+p3,Mh3,0]);

addToSum;

(*
  Diagram  42 in subprocess 1
*)
totFactor = ((32*Rqcdh3^2*Rqcdh2^2)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh3^2*imLGGh2^2-4*SC[p1,p3]^4*imLGGh3^2*
 LGGh2^2+16*SC[p1,p3]^4*imLGGh3*imLGGh2*LGGh3*LGGh2-4*SC[p1,p3]^4*imLGGh2^2*
 LGGh3^2+SC[p1,p3]^4*LGGh3^2*LGGh2^2-32*SC[p1,p3]^3*SC[p1,p2]*imLGGh3^2*
 imLGGh2^2+8*SC[p1,p3]^3*SC[p1,p2]*imLGGh3^2*LGGh2^2-32*SC[p1,p3]^3*
 SC[p1,p2]*imLGGh3*imLGGh2*LGGh3*LGGh2+8*SC[p1,p3]^3*SC[p1,p2]*imLGGh2^2*
 LGGh3^2-2*SC[p1,p3]^3*SC[p1,p2]*LGGh3^2*LGGh2^2+16*SC[p1,p3]^2*SC[p1,p2]^2*
 imLGGh3^2*imLGGh2^2-4*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*LGGh2^2+16*
 SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3*imLGGh2*LGGh3*LGGh2-4*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh2^2*LGGh3^2+SC[p1,p3]^2*SC[p1,p2]^2*LGGh3^2*LGGh2^2);
denominator =(propDen[-p1+p3,Mh2,0]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  43 in subprocess 1
*)
totFactor = ((128*Rqcdh3^4)/(1));
numerator =(16*SC[p1,p2]^4*imLGGh3^4+8*SC[p1,p2]^4*imLGGh3^2*LGGh3^2+
 SC[p1,p2]^4*LGGh3^4);
denominator =(propDen[-p1-p2,Mh3,wh3]^2);

addToSum;

(*
  Diagram  44 in subprocess 1
*)
totFactor = ((32*Rqcdh3^4)/(1));
numerator =(16*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^4+8*SC[p1,p3]^2*SC[p1,p2]^2*
 imLGGh3^2*LGGh3^2+SC[p1,p3]^2*SC[p1,p2]^2*LGGh3^4-32*SC[p1,p3]*SC[p1,p2]^3*
 imLGGh3^4-16*SC[p1,p3]*SC[p1,p2]^3*imLGGh3^2*LGGh3^2-2*SC[p1,p3]*SC[p1,p2]^
 3*LGGh3^4+16*SC[p1,p2]^4*imLGGh3^4+8*SC[p1,p2]^4*imLGGh3^2*LGGh3^2+
 SC[p1,p2]^4*LGGh3^4);
denominator =(propDen[-p1-p2,Mh3,wh3]*propDen[-p1+p4,Mh3,0]);

addToSum;

(*
  Diagram  45 in subprocess 1
*)
totFactor = ((256*Rqcdh3^4)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh3^4+8*SC[p1,p3]^4*imLGGh3^2*LGGh3^2+
 SC[p1,p3]^4*LGGh3^4);
denominator =(propDen[-p1+p3,Mh3,0]^2);

addToSum;

(*
  Diagram  46 in subprocess 1
*)
totFactor = ((16*Rqcdh3^4)/(1));
numerator =(16*SC[p1,p3]^4*imLGGh3^4+8*SC[p1,p3]^4*imLGGh3^2*LGGh3^2+
 SC[p1,p3]^4*LGGh3^4-32*SC[p1,p3]^3*SC[p1,p2]*imLGGh3^4-16*SC[p1,p3]^3*
 SC[p1,p2]*imLGGh3^2*LGGh3^2-2*SC[p1,p3]^3*SC[p1,p2]*LGGh3^4+16*SC[p1,p3]^2*
 SC[p1,p2]^2*imLGGh3^4+8*SC[p1,p3]^2*SC[p1,p2]^2*imLGGh3^2*LGGh3^2+
 SC[p1,p3]^2*SC[p1,p2]^2*LGGh3^4);
denominator =(propDen[-p1+p3,Mh3,0]*propDen[-p1+p4,Mh3,0]);

addToSum;

finishSum;
