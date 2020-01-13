/*************************************************************************

Reference:
Nevalainen et al. (2020)
Frailty modelling under a selective sampling protocol: an application to 
type 1 diabetes related autoantibodies
Submitted

The data is assumed to contain the following variables:

ID 		An identification variable for identifying subjects
AB 		A variable that identifies the autoantibodies (four levels)
AGE1	Lower limit of the interval-censored time to event
AGE2	Upper limit of the interval-censored time to event
RESPONSE indicates whether the autoantibody was negative or positive
GENOTYPE, SEX and MUNIC are binary (0/1) explanatory variables

The data set can contain several rows per individual per autoantibody
corresponding to the times of measurement

*************************************************************************/

* An example frailty model code with three explanatory variables and a piecewise linear log hazard *;

proc nlmixed data = abs (where=(age2>0 and munic ne .));
  parms a1 = -7.4  b11 = 0.29  b21 = 0.069 b31 = 0.0042
        a2 = -8.0  b12 = 0.19  b22 = 0.055 b32 = 0.0052 
        a3 = -9.01 b13 = 0.09  b23 = 0.098 b33 = 0.0052 
        a4 = -7.6  b14 = 0.04  b24 = 0.075 b34 = 0.0029
		_grisk1 = 0.614 _grisk2 = 0.416 _grisk3 = 0.91  _grisk4 = 0.59
        _sex1   = -0.1  _sex2   = -0.31 _sex3   = -0.21 _sex4 = -0.51
        _munic1 = 0.4  _munic2 = 0.12  _munic3 = 0.12  _munic4 = 0.19;

  *** Linear predictors ***;
  eta1 = _grisk1*genotype + _sex1*sex + _munic1*munic;
  eta2 = _grisk2*genotype + _sex2*sex + _munic2*munic;
  eta3 = _grisk3*genotype + _sex3*sex + _munic3*munic;
  eta4 = _grisk4*genotype + _sex4*sex + _munic4*munic;
  *** The part of the integral which lies in the interval (0,2) ***;
  lim1l = min(age1,2);
  lim1u = min(age2,2);
  if ab = "ICA"  then int1 = exp(a1+b11*lim1u)/b11 - exp(a1+b11*lim1l)/b11;
  if ab = "IAA"  then int1 = exp(a2+b12*lim1u)/b12 - exp(a2+b12*lim1l)/b12;
  if ab = "IA2A" then int1 = exp(a3+b13*lim1u)/b13 - exp(a3+b13*lim1l)/b13;
  if ab = "GADA" then int1 = exp(a4+b14*lim1u)/b14 - exp(a4+b14*lim1l)/b14;
  *** The part of the integral which lies in the interval (2,4) ***;
  lim2l = max(2,min(age1,4));
  lim2u = max(2,min(age2,4));
  if ab = "ICA"  then int2 = exp(a1+b11*lim2u+b21*(lim2u-2))/(b11+b21) - exp(a1+b11*lim2l+b21*(lim2l-2))/(b11+b21);
  if ab = "IAA"  then int2 = exp(a2+b12*lim2u+b22*(lim2u-2))/(b12+b22) - exp(a2+b12*lim2l+b22*(lim2l-2))/(b12+b22);
  if ab = "IA2A" then int2 = exp(a3+b13*lim2u+b23*(lim2u-2))/(b13+b23) - exp(a3+b13*lim2l+b23*(lim2l-2))/(b13+b23);
  if ab = "GADA" then int2 = exp(a4+b14*lim2u+b24*(lim2u-2))/(b14+b24) - exp(a4+b14*lim2l+b24*(lim2l-2))/(b14+b24);
  *** Remaining part of the integral ***;
  lim3l = max(4,age1);
  lim3u = max(4,age2);
  if ab = "ICA" then int3 = exp(a1+b11*lim3u+b21*(lim3u-2)+b31*(lim3u-4))/(b11+b21+b31)
                          - exp(a1+b11*lim3l+b21*(lim3l-2)+b31*(lim3l-4))/(b11+b21+b31);
  if ab = "IAA" then int3 = exp(a2+b12*lim3u+b22*(lim3u-2)+b32*(lim3u-4))/(b12+b22+b32)
                          - exp(a2+b12*lim3l+b22*(lim3l-2)+b32*(lim3l-4))/(b12+b22+b32);
  if ab = "IA2A" then int3= exp(a3+b13*lim3u+b23*(lim3u-2)+b33*(lim3u-4))/(b13+b23+b33)
                          - exp(a3+b13*lim3l+b23*(lim3l-2)+b33*(lim3l-4))/(b13+b23+b33);
  if ab = "GADA" then int3= exp(a4+b14*lim3u+b24*(lim3u-2)+b34*(lim3u-4))/(b14+b24+b34)
                          - exp(a4+b14*lim3l+b24*(lim3l-2)+b34*(lim3l-4))/(b14+b24+b34);

  if ab = "ICA"  then pi = exp(-(int1+int2+int3)*exp(eta1+u1));
  if ab = "IAA"  then pi = exp(-(int1+int2+int3)*exp(eta2+(3-g2-g3)*u1)); * have checked that they are all positive *;
  if ab = "IA2A" then pi = exp(-(int1+int2+int3)*exp(eta3+g2*u1));
  if ab = "GADA" then pi = exp(-(int1+int2+int3)*exp(eta4+g3*u1));

  ll = response*log(1-pi)+(1-response)*log(pi);

  random u1 ~ normal(0,exp(var)) subject = id;
  model response ~ general(ll);

  * Tests for differential associations *;
  contrast "sex heterogeneity" _sex2-_sex1, _sex3-_sex1, _sex4-_sex1;
  contrast "grisk heterogeneity" _grisk2-_grisk1, _grisk3-_grisk1, _grisk4-_grisk1;
  contrast "munic heterogeneity" _munic2-_munic1, _munic3-_munic1, _munic4-_munic1;
run;
