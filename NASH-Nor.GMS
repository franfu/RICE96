* This is a sub-model of the RICE model for solving the Nash
* equilibrium.

* Assuming that regions are linked only by climate externalities.

CA.FX(N,T) = 0;
NFA.FX(N,T) = 0;

SCALAR TOL  /1.E-10/;

SET ITER        Iteration index                 /1*5/
    SIT         Sub interation index            /1*7/;

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

* Calculating the reaction path of the 1st region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"1",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 1)=MIE(ITER,"1",J,T);
EIE(ITER,"1",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 1)=EIE(ITER,"1",J,T);

SOLVE  RICE MAXIMIZING UT1 USING NLP;


TAXRAT("1",N,T)$(KK.M(N,T) NE 0) = -(EE.M(N,T)/KK.M(N,T))*1000;
NPROLOR("1",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("1",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("1",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("1",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));


NCCOST("1",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
  *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("1",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("1",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("1",N,S));
NNPVDAM("1",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("1",N,S));
NPVCC("1",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));
UTLEVEL("1","USA",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("USA",S)
  *LOG(C.L("USA",S)/(L("USA",S)*1000)));
UTLEVEL1("1","USA",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("USA",S)*LB(S,"USA")
  *LOG(C.L("USA",S)/(L("USA",S)*1000)));

* Calculating the reaction path of the 2nd region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"2",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 2)=MIE(ITER,"2",J,T);
EIE(ITER,"2",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 2)=EIE(ITER,"2",J,T);
SOLVE RICE MAXIMIZING UT2 USING NLP;


TAXRAT("2",N,T)$(KK.M(N,T) NE 0) = -(EE.M(N,T)/KK.M(N,T))*1000;
NPROLOR("2",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("2",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("2",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("2",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));



NCCOST("2",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
 *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("2",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("2",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("2",N,S));
NNPVDAM("2",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("2",N,S));
NPVCC("2",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));
UTLEVEL("2","JPN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("JPN",S)
  *LOG(C.L("JPN",S)/(L("JPN",S)*1000)));
UTLEVEL1("2","JPN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("JPN",S)*LB(S,"JPN")
  *LOG(C.L("JPN",S)/(L("JPN",S)*1000)));

* Calculating the reaction path of the 3rd region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"3",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 3)=MIE(ITER,"3",J,T);
EIE(ITER,"3",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 3)=EIE(ITER,"3",J,T);
SOLVE RICE MAXIMIZING UT3 USING NLP;


TAXRAT("3",N,T)$(KK.M(N,T) NE 0) = -(EE.M(N,T)/KK.M(N,T))*1000;
NPROLOR("3",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("3",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("3",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("3",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));



NCCOST("3",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
  *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("3",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("3",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("3",N,S));
NNPVDAM("3",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("3",N,S));
NPVCC("3",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));
UTLEVEL("3","EEC",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("EEC",S)
  *LOG(C.L("EEC",S)/(L("EEC",S)*1000)));
UTLEVEL1("3","EEC",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("EEC",S)*LB(S,"EEC")
  *LOG(C.L("EEC",S)/(L("EEC",S)*1000)));

* Calculating the reaction path of the 4th region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"4",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 4)=MIE(ITER,"4",J,T);
EIE(ITER,"4",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 4)=EIE(ITER,"4",J,T);
SOLVE RICE MAXIMIZING UT4 USING NLP;


TAXRAT("4",N,T)$(KK.M(N,T) NE 0) = -(EE.M(N,T)/KK.M(N,T))*1000;
NPROLOR("4",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("4",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("4",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("4",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));


NCCOST("4",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
  *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("4",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("4",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("4",N,S));
NNPVDAM("4",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("4",N,S));
NPVCC("4",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));
UTLEVEL("4","CHN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("CHN",S)
  *LOG(C.L("CHN",S)/(L("CHN",S)*1000)));
UTLEVEL1("4","CHN",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("CHN",S)*LB(S,"CHN")
  *LOG(C.L("CHN",S)/(L("CHN",S)*1000)));

* Calculating the reaction path of the 5th region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"5",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 5)=MIE(ITER,"5",J,T);
EIE(ITER,"5",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 5)=EIE(ITER,"5",J,T);
SOLVE RICE MAXIMIZING UT5 USING NLP;


TAXRAT("5",N,T)$(KK.M(N,T) NE 0) = -(EE.M(N,T)/KK.M(N,T))*1000;
NPROLOR("5",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("5",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("5",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("5",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));


NCCOST("5",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
  *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("5",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("5",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("5",N,S));
NNPVDAM("5",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("5",N,S));
NPVCC("5",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));
UTLEVEL("5","FSU",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("FSU",S)
  *LOG(C.L("FSU",S)/(L("FSU",S)*1000)));
UTLEVEL1("5","FSU",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("FSU",S)*LB(S,"FSU")
  *LOG(C.L("FSU",S)/(L("FSU",S)*1000)));

* Calculating the reaction path of the 6th region.

E.LO(N,T) = 1;
E.UP(N,T) = 450;
MIU.UP(N,T)=1;
MIU.LO(N,T)=0.001;
MIU.FX(N,"1") = 0;
MIE(ITER,"6",N,T)=MIU.L(N,T);
MIU.FX(J,T)$(ORD(J) NE 6)=MIE(ITER,"6",J,T);
EIE(ITER,"6",N,T)=E.L(N,T);
E.FX(J,T)$(ORD(J) NE 6)=EIE(ITER,"6",J,T);
SOLVE RICE MAXIMIZING UT6 USING NLP;


TAXRAT("6",N,T)$(KK.M(N,T) NE 0) = -(EE.M(N,T)/KK.M(N,T))*1000;
NPROLOR("6",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NPROLORr("6",N,T) = AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA);
NMPR("6",N,T) = (-1)*1*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NMOLD("6",N,T) = ( 1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));


NCCOST("6",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)
  *B3(T)*B1(N)*(MIU.L(N,T)**B2(N));
NCDAM("6",N,T) = 1000*AL(N,T)*K.L(N,T)**GAMA*L(N,T)**(1-GAMA)*(1 -
  1/(1+A1(N)*((TE.L(T)/2.5)**A2(N))));
NNPVCOST("6",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCCOST("6",N,S));
NNPVDAM("6",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*NCDAM("6",N,S));
NPVCC("6",N,T) = SUM(S$(ORD(S) LE ORD(T)), RR(S)*C.L(N,S));
UTLEVEL("6","ROW",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("ROW",S)
  *LOG(C.L("ROW",S)/(L("ROW",S)*1000)));
UTLEVEL1("6","ROW",T)=R*SUM(S$(ORD(S) LE ORD(T)),10*RR(S)*L("ROW",S)*LB(S,"ROW")
  *LOG(C.L("ROW",S)/(L("ROW",S)*1000)));

* Reset for the next iteration.

MIE(ITER,"7",N,T)=MIU.L(N,T);
EIE(ITER,"7",N,T)=E.L(N,T);
MIE(ITER+1,"1",N,T)=MIE(ITER,"7",N,T);
EIE(ITER+1,"1",N,T)=EIE(ITER,"7",N,T);

ERR(ITER,N) = SUM((SIT,T), SQRT(SQR(MIE(ITER,SIT,N,T)-MIE(ITER-1,SIT,N,T))));
DEV(ITER) = MAX(ERR(ITER,"USA"),ERR(ITER,"JPN"),ERR(ITER,"CHN"),
  ERR(ITER,"EEC"),ERR(ITER,"FSU"),ERR(ITER,"ROW"));

WLEVEL(T) = UTLEVEL1("1","USA",T)+UTLEVEL1("2","JPN",T)+UTLEVEL1("3","EEC",T)+
   UTLEVEL1("4","CHN",T)+UTLEVEL1("5","FSU",T)+UTLEVEL1("6","ROW",T);
);

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

$include nashoutor.gms

