options nocenter nonumber nodate ls = 85;
libname in  "K:\Documents\My Papers\DIPP Antibodies\Data";
libname out "K:\Documents\My Papers\DIPP Antibodies\Output";

%macro ABSIM(size,sizeofrs,times);

data pars;
  set _null_;
run;

%do s = 1 %to &times;

proc datasets;
  delete parseparate parfrailty fullfrailty fullnofrailty;
run; quit;

data vastedata;
  epsilon = 0.1; 
  followup = 10;
  do id = 1 to &size;
    H1 = 0;
    u = rannor(441+&s);
	* covariate *;
	x = ranbin(253,1,0.5);
    do t = 0 to followup by epsilon;
      deltat = epsilon;
	  hazard1 = exp(+0.2*x)*exp(-6  + 0.5*(t+epsilon/2) - 0.7*max(0,t+epsilon/2-2) + u);
	  hazard2 = exp(-0.2*x)*exp(-6  + 0.5*(t+epsilon/2) - 0.2*max(0,t+epsilon/2-2) - 0.8*max(0,t+epsilon/2-4) + u);
      ica_pos = ranbin(542+&s,1,hazard1*deltat);
      iaa_pos = ranbin(292+&s,1,hazard2*deltat); 
	  age1 = t;
	  age2 = t+epsilon;
*	  output;
      if t = 0 or ica_pos = 1 or iaa_pos = 1 or age2 > followup - epsilon then output;
    end;
  end;
run;

data iaa ica;
  set vastedata;
  response = ica_pos; ab = "ICA";
  if response ne . then output ica;
  response = iaa_pos; ab = "IAA";
  if response ne . then output iaa;
run;

proc sort data = ica;  by ab id age1;
proc sort data = iaa;  by ab id age1;

data antibodies;
  set ica (where = (response = 1)) iaa (where = (response = 1));
  by ab id age1;
  if first.id;
run;

proc freq data = antibodies;
  table ab;
run;

data antibodies1;
  set antibodies ica (where = (response = 0) in = a) iaa (where = (response = 0) in = b);
  if a or b then age1 = 0;
  if a or b then age2 = followup;
run;

proc sort data = antibodies1 nodupkey;
  by id ab age1 age2 response;
run;

data antibodies2;
  merge antibodies1
        antibodies (where = (ab = "ICA") rename = (age1 = maxtime) keep = id ab age1 age2)
        antibodies (where = (ab = "IAA") rename = (age1 = maxtime) keep = id ab age1 age2);
  by id ab;
  label maxtime = " ";
  if age2 ge maxtime and maxtime ne . and response = 0 then age2 = maxtime;
  drop u hazard1 hazard2 deltat epsilon t;
  format maxtime 6.2;
run;

proc sort data = antibodies2;
  by id ab age1;
run;

* erilliset analyysit kaikille vasta-aineille *;

proc sort data = antibodies2 out = antibodies4;
  by ab id age1;
run;

*** ensin kaksi analyysia täydelliselle aineistolle ***;

proc nlmixed data = antibodies4 (where=(age1<age2 and age1 ne .));
  by ab;
  parms a1 = -5.2 b11 = 0.04 b21 = -0.0003 b31 = 0.0000002 gamma = -0.17;
  *** Covariates assumed constant within each interval ***;
  eta = gamma*x;
  *** The part of the integral which lies in the interval (0,2) ***;
  lim1l = min(age1,2);
  lim1u = min(age2,2);
  int1 = exp(a1+b11*lim1u)/b11 - exp(a1+b11*lim1l)/b11;
  *** The part of the integral which lies in the interval (2,4) ***;
  lim2l = max(2,min(age1,4));
  lim2u = max(2,min(age2,4));
  int2 = exp(a1+b11*lim2u+b21*(lim2u-2))/(b11+b21) - exp(a1+b11*lim2l+b21*(lim2l-2))/(b11+b21);
  *** Remaining part of the integral ***;
  lim3l = max(4,age1);
  lim3u = max(4,age2);
  int3  = exp(a1+b11*lim3u+b21*(lim3u-2)+b31*(lim3u-4))/(b11+b21+b31)
        - exp(a1+b11*lim3l+b21*(lim3l-2)+b31*(lim3l-4))/(b11+b21+b31);
  *** Note that this is not yet the survival function ***;
  pi = exp(-(int1+int2+int3)*exp(eta));
  ll = response*log(1-pi)+(1-response)*log(pi);
  model response ~ general(ll);
  ods output ParameterEstimates = fullnofrailty;
run;

proc sort data = antibodies2;
  by id ab age1;
run;

proc nlmixed data = antibodies2 (where=(age1<age2 and age1 ne .));
  parms a1 = -6.2  b11 = 0.29  b21 = -0.39 b31 = 0.802
        a2 = -4.0  b12 = 0.19  b22 = -0.35 b32 = 0.052 var = 0.1079 gamma1 = 0.11 gamma2 = -0.27;
  midage = (age1+age2)/2;
  *** No covariates ***;
  eta1 = gamma1*x;
  eta2 = gamma2*x;
  *** The part of the integral which lies in the interval (0,2) ***;
  lim1l = min(age1,2);
  lim1u = min(age2,2);
  if ab = "ICA" then int1 = exp(a1+b11*lim1u)/b11 - exp(a1+b11*lim1l)/b11;
  if ab = "IAA" then int1 = exp(a2+b12*lim1u)/b12 - exp(a2+b12*lim1l)/b12;
  *** The part of the integral which lies in the interval (2,4) ***;
  lim2l = max(2,min(age1,4));
  lim2u = max(2,min(age2,4));
  if ab = "ICA" then int2 = exp(a1+b11*lim2u+b21*(lim2u-2))/(b11+b21) - exp(a1+b11*lim2l+b21*(lim2l-2))/(b11+b21);
  if ab = "IAA" then int2 = exp(a2+b12*lim2u+b22*(lim2u-2))/(b12+b22) - exp(a2+b12*lim2l+b22*(lim2l-2))/(b12+b22);
  *** Remaining part of the integral ***;
  lim3l = max(4,age1);
  lim3u = max(4,age2);
  if ab = "ICA" then int3 = exp(a1+b11*lim3u+b21*(lim3u-2)+b31*(lim3u-4))/(b11+b21+b31)
                          - exp(a1+b11*lim3l+b21*(lim3l-2)+b31*(lim3l-4))/(b11+b21+b31);
  if ab = "IAA" then int3 = exp(a2+b12*lim3u+b22*(lim3u-2)+b32*(lim3u-4))/(b12+b22+b32)
                          - exp(a2+b12*lim3l+b22*(lim3l-2)+b32*(lim3l-4))/(b12+b22+b32);
  *** Note that this is not yet the survival function ***;
  if ab = "ICA" then pi = exp(-(int1+int2+int3)*exp(eta1+u1));
  if ab = "IAA" then pi = exp(-(int1+int2+int3)*exp(eta2+u1));
  ll = response*log(1-pi)+(1-response)*log(pi);
  random u1 ~ normal(0,exp(var)) subject = id;
  model response ~ general(ll);
  ods output ParameterEstimates = fullfrailty;
run;

* only obtain IAA data on those who have ICA = 1 or a random sample of individuals who have ICA = 0 *;

data notdeleted1;
  set antibodies2 (where = (response = 1 and ab = "ICA"));
  delindex1 = 0;
  keep id delindex1;
run;

* the data only contains max of two rows of IAA information *;

data notdeleted2;
  set antibodies2 (where = (ab = "IAA" and response = 0));
  if ranuni(91919+&s) < &sizeofrs/100 then delindex2 = 0;
  if delindex2 = 0;
  keep id delindex2;
run;

data antibodies3;
  merge antibodies2 notdeleted1 notdeleted2;
  by id;
  if ab = "ICA" or delindex1 = 0 or delindex2 = 0;
run;

proc nlmixed data = antibodies3 (where=(age1<age2 and age1 ne .));
  parms a1 = -6.2  b11 = 0.29  b21 = -0.39 b31 = 0.802
        a2 = -4.0  b12 = 0.19  b22 = -0.35 b32 = 0.052 var = 0.179
        gamma1 = 0.11 gamma2 = -0.27;
  midage = (age1+age2)/2;
  *** No covariates ***;
  eta1 = gamma1*x;
  eta2 = gamma2*x;
  *** The part of the integral which lies in the interval (0,2) ***;
  lim1l = min(age1,2);
  lim1u = min(age2,2);
  if ab = "ICA" then int1 = exp(a1+b11*lim1u)/b11 - exp(a1+b11*lim1l)/b11;
  if ab = "IAA" then int1 = exp(a2+b12*lim1u)/b12 - exp(a2+b12*lim1l)/b12;
  *** The part of the integral which lies in the interval (2,4) ***;
  lim2l = max(2,min(age1,4));
  lim2u = max(2,min(age2,4));
  if ab = "ICA" then int2 = exp(a1+b11*lim2u+b21*(lim2u-2))/(b11+b21) - exp(a1+b11*lim2l+b21*(lim2l-2))/(b11+b21);
  if ab = "IAA" then int2 = exp(a2+b12*lim2u+b22*(lim2u-2))/(b12+b22) - exp(a2+b12*lim2l+b22*(lim2l-2))/(b12+b22);
  *** Remaining part of the integral ***;
  lim3l = max(4,age1);
  lim3u = max(4,age2);
  if ab = "ICA" then int3 = exp(a1+b11*lim3u+b21*(lim3u-2)+b31*(lim3u-4))/(b11+b21+b31)
                          - exp(a1+b11*lim3l+b21*(lim3l-2)+b31*(lim3l-4))/(b11+b21+b31);
  if ab = "IAA" then int3 = exp(a2+b12*lim3u+b22*(lim3u-2)+b32*(lim3u-4))/(b12+b22+b32)
                          - exp(a2+b12*lim3l+b22*(lim3l-2)+b32*(lim3l-4))/(b12+b22+b32);
  *** Note that this is not yet the survival function ***;
  if ab = "ICA" then pi = exp(-(int1+int2+int3)*exp(eta1+u1));
  if ab = "IAA" then pi = exp(-(int1+int2+int3)*exp(eta2+u1));
  ll = response*log(1-pi)+(1-response)*log(pi);
  random u1 ~ normal(0,exp(var)) subject = id;
  model response ~ general(ll);
  ods output ParameterEstimates = parfrailty;
run;

proc sort data = antibodies3;
  by ab id age1;
run;

proc nlmixed data = antibodies3 (where=(age1<age2 and age1 ne .));
  by ab;
  parms a1 = -7.8 b11 = 0.04 b21 = -0.03 b31 = 0.0000002 gamma = 0.17;
  *** Covariates assumed constant within each interval ***;
  eta = gamma*x;
  *** The part of the integral which lies in the interval (0,2) ***;
  lim1l = min(age1,2);
  lim1u = min(age2,2);
  int1 = exp(a1+b11*lim1u)/b11 - exp(a1+b11*lim1l)/b11;
  *** The part of the integral which lies in the interval (2,4) ***;
  lim2l = max(2,min(age1,4));
  lim2u = max(2,min(age2,4));
  int2 = exp(a1+b11*lim2u+b21*(lim2u-2))/(b11+b21) - exp(a1+b11*lim2l+b21*(lim2l-2))/(b11+b21);
  *** Remaining part of the integral ***;
  lim3l = max(4,age1);
  lim3u = max(4,age2);
  int3  = exp(a1+b11*lim3u+b21*(lim3u-2)+b31*(lim3u-4))/(b11+b21+b31)
        - exp(a1+b11*lim3l+b21*(lim3l-2)+b31*(lim3l-4))/(b11+b21+b31);
  *** Note that this is not yet the survival function ***;
  pi = exp(-(int1+int2+int3)*exp(eta));
  ll = response*log(1-pi)+(1-response)*log(pi);
  model response ~ general(ll);
  ods output ParameterEstimates = parseparate;
run;

data pars;
  set parfrailty (in = a) fullfrailty (in = c) parseparate (in=b) fullnofrailty (in = d) pars;
  if a or c or b or d then simulation = &s;
  if a then model = "wfrailty - deletion    ";
  if c then model = "wfrailty - complete    ";
  if (a or c) and (substr(parameter,3,1) = "1" or parameter = "a1") then ab = "ICA";
  if (a or c) and (substr(parameter,3,1) = "2" or parameter = "a2") then ab = "IAA";
  if (a or c) and parameter = "gamma1" then ab = "ICA";
  if (a or c) and parameter = "gamma2" then ab = "IAA";
  if parameter in ("gamma1","gamma2") then parameter = "gamma";
  if b or d then do;
    if b then model = "nofrailty - deletion";
	if d then model = "nofrailty - complete";
    if parameter ne "gamma" then parameter = substr(parameter,1,2);
    if ab = "ICA" and substr(parameter,1,1) = "b" then parameter = trim(parameter)||"1";
    if ab = "IAA" and substr(parameter,1,1) = "b" then parameter = trim(parameter)||"2";
	if ab = "IAA" and parameter = "a1" then parameter = "a2";
  end;
run;

%end;

proc sort data = pars;
  by ab parameter model;
run;

%let myout = z:\Documents\My Papers\DIPP Antibodies\Output;

data out.sim_par_estimates_&size._&sizeofrs;
  set pars;
run;

ods latex image_dpi = 600 style = journal gpath = "&myout";
ods graphics on / imagename = "absim-IAA-&size.-&sizeofrs.-" imagefmt = jpg reset = index  height = 5cm width = 10cm;

proc sgpanel data = pars (where = (ab = "IAA"));
  panelby parameter / rows = 1 columns = 6 novarname uniscale = row;
  hbox estimate / group = model;
run;

ods graphics off;
ods latex close;

ods latex image_dpi = 600 style = journal gpath = "&myout";
ods graphics on / imagename = "absim-ICA-&size.-&sizeofrs.-" imagefmt = jpg reset = index  height = 5cm width = 10cm;

proc sgpanel data = pars (where = (ab = "ICA"));
  panelby parameter / rows = 1 columns = 5 novarname uniscale = row;
  hbox estimate / group = model;
run;

ods graphics off;
ods latex close;

%mend ABSIM;

***%ABSIM(4000, 20,20);

*%ABSIM(5000, 10,500);
*%ABSIM(5000, 50,500);

* Kokeillaan riippuvalla x:llä *;
%ABSIM(5000, 10,1000);
%ABSIM(5000, 50,1000);


* Tätä ei tosiasiassa tarvita, sillä completet tulee tehtyä osana makroa *;
*%ABSIM(5000,100,3);


