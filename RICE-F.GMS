$TITLE  THE RICE MODEL: Version: RCE0896

$ontext
                        THE RICE MODEL

        William D. Nordhaus (Yale), Zili Yang (MIT)

This version simulates strategic coalition behavior among regions.

$offtext

OPTION NLP=minos5;

*fare attenzione all'extended horizon: quello orig. e' 1*25

* !!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!1
*per periodi superiori al 27 e' posta una condizione su kk.m
*in modo che sia possibile calcolare i pesi di Negishi

SETS    T      Extended horizon        /1*25/
        TB(T)  Base period
        TT(T)  Terminal period
        N      Regions of the world    /USA, JPN, EEC, CHN, FSU, ROW/
        MD     Categories of initial economic-emission data
        /Y0
         YPC0
         K0
         L0
         E_F0
         E_L0/
        TD     Categories of trend data  /YPC_A, L_A, L_I, A_A, A_I, SIG_I, SIG_A, DEL_TREE /;

ALIAS   (N, J);
ALIAS   (T, S);
TB(T) = YES$(ORD(T) EQ 1);
TT(T) = YES$(ORD(T) EQ CARD(T));
DISPLAY TB, TT;

* assign one-dimension initial values and coefficients.

SCALARS  BET      Elasticity of marginal untility      /0/
         R        Rate of social time pref per year    /.03/
         DELTAM   Removal rate carbon per decade    /.0833/
         DK       Depreciatin rate of capital per year /.1/
         GAMA     Capital elasticity in output  /.25/
         M0       CO2-equiv concentration 65 b.t.c.  /688.4/
         TL0      Lower stratum temp (c) 65    /.1/
         T0       Atmospheric temp (C) 65  /.5/
         ATRET    Marginal atmospheric rentention R  /.64/
         C1       Coefficient for upper level   /.226/
         LAM      Climate feedback factor  /1.41/
         C3       Coeff trans upper to lower level  /.44/
         C4       Coeff of transfer for lower level /.02/

$ontext
The following table contains the data of the economy and emissions in
the initial period. i.e., 1990.

Y0   --- GNP in 1990 US$ (billions)
YPC0 --- GNP per capita in 1990 US$ (thousands)
K0   --- Capital stock in 1990 US$ (billions)
L0   --- Population (millions)
E_F0 --- GHG emissions from fuel uses (million tons of CO2 equivalent)
E_L0 --- GHG emissions from land use changes
$offtext

TABLE INI_DATA(MD,N)   Initial macroeconomic data
        USA      JPN       EEC       CHN       FSU       ROW
Y0    5464.80  2932.06   6828.04    370.02    855.21   4628.62
YPC0    21.83    23.73     18.63      0.33      2.96      1.49
K0   14262.51  8442.25  18435.71   1025.79   2281.90   9842.22
L0     250.37   123.54    366.50   1133.68    289.32   3102.69
E_F0    13.6      2.92      8.72      6.69     10.66     17.0
E_L0      .1      0.        0.        1.36      0.       17.3 ;

$ontext
The following table contains the data of trends of population
growth, productivity growth and  the change of ghg emission vs.
output ratios.

YPC_A --- Asymptotic level of GNP per capita (1990 = 1) not used
L_A   --- Asymptotic level of population (1990 = 1)
L_I   --- Initial population growth rates (percent per decade)
A_A   --- Asymptotic level of productivity (1990 =  1)
A_I   --- Initial productivity growth rates (percent per decade)
SIG_A --- Asymptotic level of (GHG emission / GNP ratios) (1990 = 1)
SIG_I --- Initial change rates of (GHG emission / GNP) ratios
$offtext

TABLE TREND_DATA(TD,N) Trend data
         USA       JPN       EEC       CHN       FSU       ROW
YPC_A   10.00      10.00    100.00     15.00     30.00     50.00
L_A      1.1781     1.0081    1.1713    1.4819    1.2812    2.3262
L_I      5.51       0.65      4.30     10.73      6.14     18.78
A_I     13.10      23.00     15.00     43.50     20.00     32.00
A_A      4.19       3.41      4.23     65.50     14.7      15.0
SIG_I   16.00      17.00     18.00     17.00     20.00     10.00
SIG_A    0.39       0.602     0.547     0.05      0.075     0.23
DEL_TREE  .1         .1        .1        .1        .1        .1 ;

PARAMETERS    GL0(N)            Initial growth rate coeff of labor
              TAXRATE(N,T)      Carbon tax in period T and region N
              ETREE(N,T)        Emissions from deforestation
              DLAB(N)           Decelaration rate of labor growth rate
              GL(N,T)           Growth rate of labor
              L(N,T)            Level of labor force
              A0(N)             Efficiency parameter
              GA0(N)            Initial growth rate coeff of productivity
              DELA(N)           Decelaration rate of productivity growth rate
              GA(N,T)           Growth rate of productivity
              AL(N,T)           Level of productivity
              SIG0(N)           Initial E\Y ratio coeff
              SSIG(N)
              DSIG(N)
              SIGMA(N,T)        Emission-output ratio
              GSIG(N,T)         E\Y deceleration rate      ;


* Population growth trend.
GL0(N) = LOG(TREND_DATA("L_A", N));
DLAB(N) = - LOG(1 - LOG(1 + TREND_DATA("L_I", N) / 100) / GL0(N));
GL(N,T) = GL0(N)*(1-EXP(-DLAB(N)*(ORD(T)-1)));
L(N,T) = INI_DATA("L0",N)*EXP(GL(N,T)) / 1000;

DISPLAY GL0, DLAB;
DISPLAY GL, L;

*  Productivity growth trend.
A0(N) = INI_DATA("Y0",N)/(INI_DATA("K0",N)**GAMA*INI_DATA("L0",N)**(1-GAMA));
GA0(N) = LOG(TREND_DATA("A_A", N));
DELA(N) = - LOG(1 - LOG(1 + TREND_DATA("A_I", N) / 100) / GA0(N));
GA(N,T) = GA0(N)*(1-EXP(-DELA(N)*(ORD(T)-1)));
AL(N,T) = A0(N)*EXP(GA(N,T));

DISPLAY A0, GA0, DELA;
DISPLAY GA, AL;

*  GHG emissions/ output ratio trend.
SSIG(N) = LOG(TREND_DATA("SIG_A",N));
SIG0(N) = 100 * (INI_DATA("E_F0",N)) / INI_DATA("Y0",N);
DSIG(N) =  LOG(1 - LOG(1 + TREND_DATA("SIG_I", N) / 100) / SSIG(N));
GSIG(N,T) = SSIG(N) * (1-EXP(-DSIG(N)*(ORD(T)-1)));
SIGMA(N,T) = SIG0(N) * EXP(GSIG(N,T));

DISPLAY SIG0;
DISPLAY GSIG, SIGMA;

* Deforestation rates
ETREE(N,T) = INI_DATA("E_L0",N)*(1-TREND_DATA("DEL_TREE",N))**(ORD(T)-1);


* Initial welfare weights.
parameter       lb(t,N);

lb(t,N) = 1;

PARAMETERS    A1(N)     Damage function coeff #1 of doubling CO2
                    /USA   0.01102
                     JPN   0.01174
                     EEC   0.01174
                     CHN   0.01523
                     FSU   0.00857
                     ROW   0.02093/
              A2(N)     Damage function coeff #2 of doubling CO2
                    /USA   3.5
                     JPN   3.5
                     EEC   3.5
                     CHN   3.5
                     FSU   3.5
                     ROW   3.5/
              B1(N)     Intercept of control cost function
                    /USA   0.07
                     JPN   0.05
                     EEC   0.05
                     CHN   0.15
                     FSU   0.15
                     ROW   0.10/
              B2(N)     Exponent of control cost function
                    /USA   2.887
                     JPN   2.887
                     EEC   2.887
                     CHN   2.887
                     FSU   2.887
                     ROW   2.887/;

PARAMETERS    RR(T)           Discount factor
              LB(T,N)         Welfare weights
              FORCOTH(T)      Exogenous force
              B3(T)           Cost function adjustment factor
              NN(N)           Number of countries in region N
              NN0(N)
              DELTA(N)        Delta function of coalition
              DUM(T)          Dummy variable 0 except TLAST;

*B3("1") = 1;
*B3("2") = 2;
*B3("3") = 2;
*B3("4") = 1.4;
*B3("5") = 1.2;
*B3("6") = 1.1;
*B3(T)$(ORD(T) GT 6) = 1;

b3(t) = 1;
NN0(N) = 1;
NN(N) = NN0(N);
RR(T) = (1+R)**(10*(1-ORD(T)));
DELTA(N) = 1;
FORCOTH(T) = .78;
FORCOTH(T)$(ORD(T) LT 12) = -.07+(.85/10)*(ORD(T)-1);

DUM(T) = 1$(ORD(T) EQ CARD(T));

VARIABLES     C(N,T)        Consumption
              Y(N,T)        Output
              YS(T)         World aggregate output
              K(N,T)        Capital stock
              I(N,T)        Investment
              RI(N,T)       Interest rate per annum
              MIU(N,T)      Emission control rates
              FORC(T)       Radiative forcing
              TE(T)         Temperature atmosphere C
              TL(T)         Temperature lower ocean C
              E(N,T)        CO2 emissions
              ES(T)         World aggregate CO2 emissions
              M(T)          CO2 concentrations
              NIR(N,T)      Net invisible receipts
              NFA(N,T)      Net amount surplus
              NFY(N,T)      NFA per Y
              CA(N,T)       Current account definition
              CAS(N,T)      Current account surplus
              NRR(T)        Average global rate of return to capital
              UT1
              UT2
              UT3
              UT4
              UT5
              UT6

              UTILITY       The maximand;

POSITIVE VARIABLES C, Y, K, I, MIU, TE, TL, M, E;

EQUATIONS     KK(N,T)      Capital stock balance
              INVEQ(N,T)   Definition of investment
              KB(N,T)      Terminal condition of K
              YY(N,T)      Production function
              YYS(T)       Aggregate production
              KK0(N,T)     Initial capital stock
              EE(N,T)      Emission process
              EES(T)       Aggregate emissions
              MM(T)        CO2 concentration process
              TTE(T)       Temperature raising process
              MM0(T)       Initial condition of CO2 concentration
              TTE0(T)      Initial condition of Temperature
              TTL(T)       Deep ocean temperature process
              TTL0(T)      Initial condition for deep ocean temperature
              FORCE(T)     Radiative forcing equation
              NRREQ(T)     Rate of return definition
              NFAEQ(N,T)   Definition of net foreign assets
              NFACEQ(N,T)  NFA constraint
              CAN(T)       World current account balance
              CASEQ(N,T)   Current account surplus
              CASCON(N,T)  Current account deficit constraint
              UTT1
              UTT2
              UTT3
              UTT4
              UTT5
              UTT6

              UTILI        Definition of maximand;

YY(N,T)..          Y(N,T) =E= AL(N,T)*K(N,T)**GAMA*L(N,T)**(1-GAMA)
                           *(1-B3(T)*B1(N)*(MIU(N,T)**B2(N)))
                           /(1+(A1(N)/NN(N))*((TE(T)/2.5)**A2(N)));
YYS(T)..           YS(T) =E= SUM(N,Y(N,T));
INVEQ(N,T)..       I(N,T) =E= Y(N,T) - C(N,T);
KK0(N,TB)..        K(N,TB) =E= INI_DATA("K0", N) / 1000;
KB(N,TT)..         R*K(N,TT) =L= I(N,TT);
KK(N,T+1)..        K(N,T+1) =E= 10*I(N,T)+(1-DK)**10*K(N,T);
EE(N,T)..          E(N,T) =E= 10*SIGMA(N,T)*(1-MIU(N,T))*AL(N,T)*
                                K(N,T)**GAMA*L(N,T)**(1-GAMA)+ETREE(N,T);
EES(T)..           ES(T) =E= SUM(N,E(N,T));
FORCE(T)..         FORC(T) =E=4.1*(LOG(M(T)/590)/LOG(2))+FORCOTH(T);
MM0(TB)..          M(TB) =E= M0;
MM(T+1)..          M(T+1) =E= 590+ATRET*SUM(N,E(N,T))+(1-DELTAM)*(M(T)-590);
TTE0(TB)..         TE(TB) =E= T0;
TTE(T+1)..         TE(T+1) =E=TE(T)+C1*(FORC(T)-LAM*TE(T)-C3*(TE(T)-TL(T)));
TTL0(TB)..         TL(TB) =E= TL0;
TTL(T+1)..         TL(T+1) =E= TL(T)+C4*(TE(T)-TL(T));
NRREQ(T)..         NRR(T) =E= SUM(N, GAMA*Y(N,T))/SUM(N,K(N,T))-DK;
NFAEQ(N,T+1)..     NFA(N,T+1) =E= NFA(N,T) + 10*CAS(N,T);
CASEQ(N,T)..       CAS(N,T) =E= NRR(T)*NFA(N,T)+CA(N,T);
CASCON(N,T)..      CAS(N,T)+0.1*AL(N,T)*K(N,T)**GAMA*L(N,T)**(1-GAMA)=G=0;
NFACEQ(N,T)..      NFA(N,T)+0.5*Y(N,T) =G= 0;
CAN(T)..           SUM(N, CA(N,T)) =E= 0;

UTT1..          UT1 =E= SUM(T,RR(T)*L("USA",T)*LOG(C("USA",T)/L("USA",T)));
UTT2..          UT2 =E= SUM(T,RR(T)*L("JPN",T)*LOG(C("JPN",T)/L("JPN",T)));
UTT3..          UT3 =E= SUM(T,RR(T)*L("EEC",T)*LOG(C("EEC",T)/L("EEC",T)));
UTT4..          UT4 =E= SUM(T,RR(T)*L("CHN",T)*LOG(C("CHN",T)/L("CHN",T)));
UTT5..          UT5 =E= SUM(T,RR(T)*L("FSU",T)*LOG(C("FSU",T)/L("FSU",T)));
UTT6..          UT6 =E= SUM(T,RR(T)*L("ROW",T)*LOG(C("ROW",T)/L("ROW",T)));


UTILI..         UTILITY =E= R*SUM((T,N), DELTA(N)*LB(T,N)*RR(T)*L(N,T)*
                        LOG(C(N,T)/L(N,T)));


*Initial conditions

*Y.FX(N,TB) = INI_DATA("Y0",N) / 1000;
*K.FX(N,TB) = INI_DATA("K0",N) / 1000;
*C.FX(N,TB) = INI_DATA("C0",N);
NFA.FX(N,TB) = 0;

*Upper and lower bounds

K.LO(N,T) = 0.9*INI_DATA("K0",N) / 1000;
Y.LO(N,T) = 0.9*INI_DATA("Y0",N) / 1000;

Y.UP("USA",T) = 10*INI_DATA("Y0","USA") / 1000;
Y.UP("JPN",T) = 20*INI_DATA("Y0","JPN") / 1000;
Y.UP("EEC",T) = 10*INI_DATA("Y0","EEC") / 1000;
Y.UP("CHN",T) = 250*INI_DATA("Y0","CHN") / 1000;
Y.UP("FSU",T) = 25*INI_DATA("Y0","FSU") / 1000;
Y.UP("ROW",T) = 150*INI_DATA("Y0","ROW") / 1000;

NFA.FX(N,T) = 0;
CA.FX(N,T) = 0;
I.LO(N,T) = 0.001;
NFA.UP(N,TT) = 0;
NFA.LO(N,TT) = 0;
MIU.FX(N,"1") = 0;
TE.LO(T) = 0.9*T0;
TE.UP(T) = 7;
TL.LO(T) = TL0;
TL.UP(T) = 2.5;
M.LO(T) = 100;
*M.LO(T) = M0;
E.LO(N,T) = 1;
*ES.UP(T) = 80;
*ES.UP(T)$(ORD(T) GT 3) = 0.8*85.50;
*E.UP(N,T)$(ORD(T) GT 3) = 0.8*INI_DATA("E0",N);
C.LO(N,T) = 0.5*INI_DATA("Y0",N) / 1000;


OPTION ITERLIM = 99999;
OPTION LIMROW = 0;
OPTION LIMCOL = 0;
OPTION SOLPRINT = OFF;
OPTION RESLIM = 99999;
MODEL RICE /ALL/;
*RICE.OPTFILE = 1;
*FILE OPT MINOS OPTION FILE /MINOS5.OPT/;
*PUT OPT;
*PUT 'BEGIN'/
*'MAJOR ITERATIONS     10000'/
*'END';


* Welfare weights searching.
$ include weigh.gms

* Market solution.
*$include mark-n.gms

* Cooperative solution.
*$include coop-n.gms

*Sequential solution
*$include nash-nor.gms

* Nash solution.
$include nash-n.gms

* Stackelberg solution.
*$include stack-n.gms
