* This is a sub-model of the RICE model for solving the Nash
* equilibrium.

* Assuming that regions are linked only by climate externalities.

CA.FX(N,T) = 0;
NFA.FX(N,T) = 0;

SCALAR TOL  /1.E-10/;

SET ITER        Iteration index                 /1*10/
    SIT         Sub interation index            /1*8/;

PARAMETERS
   WLEVEL(T)
   UTLEVEL(SIT,N,T)
   UTLEVEL1(SIT,N,T)
   TAXRAT(SIT,N,T)
   NPVCC(SIT,N,T)
   DEV(ITER)
   ERR(ITER,N)
   MIA(N,T)
   EIE(ITER,SIT,N,T)
   MIE(ITER,SIT,N,T);

PARAMETERS
   NPROLOR(SIT,N,T)     Produzione lorda
   NPROLORr(SIT,N,T)    Produzione lorda diviso 1000
   NMPR(SIT,N,T)        minor prod da da
   NMOLD(SIT,N,T)       moltiplicatore fu da
   NCCOST(SIT,N,T)      Control cost
   NCDAM(SIT,N,T)       Climate change damages
   NNPVCOST(SIT,N,T)    Net present value of cost
   NNPVDAM(SIT,N,T)     Net present value of damage;

* Redefine the number of countries in a region.
NN0("ROW") = 10;
NN(N) = NN0(N);

DEV("1") = 10;

LOOP(ITER,

* Calculating the reaction path of the 1st region

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
IF(ORD(ITER)=1,
MIE(ITER,"1",N,T)=0;
ELSE
MIE(ITER,"1",N,T)=MIE(ITER-1,"8",N,T)
);
MIU.FX(J,T)$(ORD(J) NE 1)=MIE(ITER,"1",J,T);
IF(ORD(ITER)=1,
EIE(ITER,"1",N,T)=E.L(N,T);
ELSE
EIE(ITER,"1",N,T)=EIE(ITER-1,"8",N,T)
);
E.FX(J,T)$(ORD(J) NE 1)=EIE(ITER,"1",J,T);

SOLVE  RICE MAXIMIZING UT1 USING NLP;

TAXRAT("1","USA",T)$(KK.M("USA",T) NE 0) = -(EE.M("USA",T)/KK.M("USA",T))*1000;

* Calculating the reaction path of the 2nd region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"2",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 2)=MIE(ITER,"1",J,T);
EIE(ITER,"2",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 2)=EIE(ITER,"1",J,T);

SOLVE RICE MAXIMIZING UT2 USING NLP;

TAXRAT("2","JPN",T)$(KK.M("JPN",T) NE 0) = -(EE.M("JPN",T)/KK.M("JPN",T))*1000;

* Calculating the reaction path of the 3rd region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"3",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 3)=MIE(ITER,"1",J,T);
EIE(ITER,"3",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 3)=EIE(ITER,"1",J,T);

SOLVE RICE MAXIMIZING UT3 USING NLP;
TAXRAT("3","EEC",T)$(KK.M("EEC",T) NE 0) = -(EE.M("EEC",T)/KK.M("EEC",T))*1000;

* Calculating the reaction path of the 4th region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"4",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 4)=MIE(ITER,"1",J,T);
EIE(ITER,"4",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 4)=EIE(ITER,"1",J,T);

SOLVE RICE MAXIMIZING UT4 USING NLP;
TAXRAT("4","CHN",T)$(KK.M("CHN",T) NE 0) = -(EE.M("CHN",T)/KK.M("CHN",T))*1000;

* Calculating the reaction path of the 5th region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"5",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 5)=MIE(ITER,"1",J,T);
EIE(ITER,"5",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 5)=EIE(ITER,"1",J,T);

SOLVE RICE MAXIMIZING UT5 USING NLP;
TAXRAT("5","FSU",T)$(KK.M("FSU",T) NE 0) = -(EE.M("FSU",T)/KK.M("FSU",T))*1000;

* Calculating the reaction path of the 6th region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"6",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 6)=MIE(ITER,"1",J,T);
EIE(ITER,"6",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 6)=EIE(ITER,"1",J,T);

SOLVE RICE MAXIMIZING UT6 USING NLP;

TAXRAT("6","ROW",T)$(KK.M("ROW",T) NE 0) = -(EE.M("ROW",T)/KK.M("ROW",T))*1000;

MIE(ITER,"7",N,T)=MIU.L(N,T);
EIE(ITER,"7",N,T)=E.L(N,T);

* Simultaneous mitigation decisions and reset for the next iteration
MIE(ITER,"8","USA",T)=MIE(ITER,"2","USA",T);
MIE(ITER,"8","JPN",T)=MIE(ITER,"3","JPN",T);
MIE(ITER,"8","EEC",T)=MIE(ITER,"4","EEC",T);
MIE(ITER,"8","CHN",T)=MIE(ITER,"5","CHN",T);
MIE(ITER,"8","FSU",T)=MIE(ITER,"6","FSU",T);
MIE(ITER,"8","ROW",T)=MIE(ITER,"7","ROW",T);
EIE(ITER,"8","USA",T)=EIE(ITER,"2","USA",T);
EIE(ITER,"8","JPN",T)=EIE(ITER,"3","JPN",T);
EIE(ITER,"8","EEC",T)=EIE(ITER,"4","EEC",T);
EIE(ITER,"8","CHN",T)=EIE(ITER,"5","CHN",T);
EIE(ITER,"8","FSU",T)=EIE(ITER,"6","FSU",T);
EIE(ITER,"8","ROW",T)=EIE(ITER,"7","ROW",T);

ERR(ITER,N) = SUM((SIT,T), SQRT(SQR(MIE(ITER,SIT,N,T)-MIE(ITER-1,SIT,N,T))));
DEV(ITER) = MAX(ERR(ITER,"USA"),ERR(ITER,"JPN"),ERR(ITER,"CHN"),
  ERR(ITER,"EEC"),ERR(ITER,"FSU"),ERR(ITER,"ROW"));
);

********************************************************************************
NPROLOR("8",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("8",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("8",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("8",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));

NCCOST("8",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
  *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("8",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("8",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("8",N,S));
NNPVDAM("8",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("8",N,S));
NPVCC("8",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));

DISPLAY ERR;
DISPLAY DEV;

UTLEVEL("8","USA",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("USA",S)
  *LOG(C.L("USA",S)/(L("USA",S)*1000)));
UTLEVEL1("8","USA",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("USA",S)*LB(S,"USA")
  *LOG(C.L("USA",S)/(L("USA",S)*1000)));
UTLEVEL("8","JPN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("JPN",S)
  *LOG(C.L("JPN",S)/(L("JPN",S)*1000)));
UTLEVEL1("8","JPN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("JPN",S)*LB(S,"JPN")
  *LOG(C.L("JPN",S)/(L("JPN",S)*1000)));
UTLEVEL("8","EEC",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("EEC",S)
  *LOG(C.L("EEC",S)/(L("EEC",S)*1000)));
UTLEVEL1("8","EEC",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("EEC",S)*LB(S,"EEC")
  *LOG(C.L("EEC",S)/(L("EEC",S)*1000)));
UTLEVEL("8","CHN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("CHN",S)
  *LOG(C.L("CHN",S)/(L("CHN",S)*1000)));
UTLEVEL1("8","CHN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("CHN",S)*LB(S,"CHN")
  *LOG(C.L("CHN",S)/(L("CHN",S)*1000)));
UTLEVEL("8","FSU",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("FSU",S)
  *LOG(C.L("FSU",S)/(L("FSU",S)*1000)));
UTLEVEL1("8","FSU",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("FSU",S)*LB(S,"FSU")
  *LOG(C.L("FSU",S)/(L("FSU",S)*1000)));
UTLEVEL("8","ROW",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("ROW",S)
  *LOG(C.L("ROW",S)/(L("ROW",S)*1000)));
UTLEVEL1("8","ROW",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("ROW",S)*LB(S,"ROW")
  *LOG(C.L("ROW",S)/(L("ROW",S)*1000)));

WLEVEL(T) = UTLEVEL1("8","USA",T)+UTLEVEL1("8","JPN",T)+UTLEVEL1("8","EEC",T)+
   UTLEVEL1("8","CHN",T)+UTLEVEL1("8","FSU",T)+UTLEVEL1("8","ROW",T);


* Writing a short report of convergence status.

FILE DITER;
PUT DITER;
DITER.PC = 4;
DITER.PW = 255;
PUT // @20, "ITERATION OUTCOMES"/;
LOOP(ITER, PUT / ITER.TL:<2;
LOOP(SIT, PUT / SIT.TL:<2;
LOOP(N, PUT / N.TL:<5;
LOOP(T, PUT @(ORD(T)*9-4), MIE(ITER,SIT,N,T):9:6))));
PUTCLOSE DITER;

$include nashout.gms

