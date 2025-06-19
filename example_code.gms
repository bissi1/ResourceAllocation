$ONTEXT
Bismark Singh
2025
$OFFTEXT

OPTIONS PROFILE=3, RESLIM=42000, LIMROW=5, LIMCOL=5, 
        LP=CPLEX, MIP=CPLEX, RMIP=CPLEX, SOLPRINT=OFF, 
        DECIMALS=8, OPTCR=0.001, OPTCA=0.001, THREADS=8, INTEGER4=0;

*---------------------------------------------
* Runtime recorder
SCALAR start_time, end_time, run_time_total;

*---------------------------------------------
* Sets
SETS
    T  times     /t1*t16/
    C  counties  /c1*c254/
    W  scenarios /w1*w2/ ;

ALIAS (w, ww);
ALIAS (t, tt);

*---------------------------------------------
* Parameters
PARAMETER
    benefit(c,t,w)     benefit
    pop(c,t,w)         population of sick
    total_pop(c,t,w)   total population ;

SCALAR frac_worry;
frac_worry = 0.000;

TABLE available(t,*)   Availability of antivirals
$ONDELIM
$INCLUDE avail_monthly.csv
$OFFDELIM

PARAMETER avail(t);
avail(t) = available(t, 'Amount');

PARAMETER BigM(t);
BigM(t) = SUM(tt$(ORD(tt) LE ORD(t)), avail(tt));

$INCLUDE pop_data.txt
$INCLUDE ben_data.txt
$INCLUDE clear_data.txt

total_pop(c,t,w) = (1 + frac_worry) * pop(c,t,w);

*---------------------------------------------
* Validation
PARAMETER dummy(c,t,w);
dummy(c,t,w) = 1$((benefit(c,t,w) NE 0) AND (total_pop(c,t,w) EQ 0));
DISPLAY dummy;

SCALAR dummy2;
dummy2 = SUM((c,t,w), dummy(c,t,w));
DISPLAY dummy2;

IF (dummy2 > 0, ABORT "benefits and total_pop don't align", dummy2);

* Round down total population
total_pop(c,t,w) = FLOOR(total_pop(c,t,w));

* Zero benefit if no population
benefit(c,t,w)$(total_pop(c,t,w) EQ 0) = 0;

DISPLAY total_pop, benefit;

* Best possible objective value
SCALAR probability;
probability = 1 / CARD(w);

SCALAR best_obj;
best_obj = SUM((c,t,w)$(benefit(c,t,w) > 0), probability * benefit(c,t,w));
DISPLAY best_obj;

*---------------------------------------------
* Begin Model
VARIABLES OBJ;
POSITIVE VARIABLES R(c,t), Q(c,t,w), F(c,t,w);
BINARY VARIABLES X(c,t,w);

EQUATIONS
    OBJECTIVE_TRUE,
    SHEL_TRUE1(c,t,w),
    SHEL_TRUE2(c,t,w),
    RELEASE_TRUE(t),
    FRAC2_TRUE(c,t,w),
    FRAC3_TRUE(c,t,w);

OBJECTIVE_TRUE..
    OBJ =E= probability * SUM((c,t,w)$(total_pop(c,t,w) > 0), benefit(c,t,w) * f(c,t,w));

SHEL_TRUE1(c,t,w)$(total_pop(c,t,w) > 0)..
    q(c,t,w) =E= q(c,t-1,w)$(ORD(t) GE 2) + r(c,t) - f(c,t,w) * total_pop(c,t,w);

SHEL_TRUE2(c,t,w)$(total_pop(c,t,w) = 0)..
    q(c,t,w) =E= q(c,t-1,w)$(ORD(t) GE 2) + r(c,t);

RELEASE_TRUE(t)..
    SUM((c,tt)$(ORD(tt) LE ORD(t)), r(c,tt)) =L= BigM(t);

FRAC2_TRUE(c,t,w)..
    x(c,t,w) =L= f(c,t,w);

FRAC3_TRUE(c,t,w)..
    q(c,t,w) =L= BigM(t) * x(c,t,w);

*---------------------------------------------
* Stopping conditions
SET t_stop(c,t);
PARAMETER dum_stop(c,t);

dum_stop(c,t) = SUM(w, SUM(tt$(ORD(tt) GE ORD(t)), total_pop(c,tt,w)));
t_stop(c,t) = YES$(dum_stop(c,t) EQ 0);
DISPLAY t_stop;

MODEL Antivirals_True /ALL/;
Antivirals_True.OPTFILE = 1;

* Bound and value fixes
f.UP(c,t,w) = 1;
r.UP(c,t) = BigM(t);
q.UP(c,t,w) = BigM(t);

x.FX(c,t,w)$(total_pop(c,t,w) = 0) = 1;
f.FX(c,t,w)$(total_pop(c,t,w) = 0) = 1;
r.FX(c,t)$t_stop(c,t) = 0;

*---------------------------------------------
* Warm start
r.L(c,t) = avail(t) / CARD(c);

*---------------------------------------------
* Solve
SOLVE Antivirals_True USING MIP MAXIMIZING OBJ;
DISPLAY obj.L, r.L, f.L;

* Output to file
FILE true_release_schedule /true_release_schedule.csv/;
true_release_schedule.PC = 5;
true_release_schedule.ND = 3;
true_release_schedule.PW = 2000;

PUT true_release_schedule;
PUT '' LOOP(t, PUT t.TL) PUT /;
LOOP(c, PUT c.TL LOOP(t, PUT r.L(c,t)) PUT /);
PUT obj.L PUT /;
PUTCLOSE true_release_schedule;