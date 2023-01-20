********************************************************;
*
* M-estimation: a worked example with connections 
* to maximum likelihood estimation 
* XXX et al.
*
* SAS Code prepared by XXX (2023/01/20)
*
********************************************************;

*Data;
data dat_cnts;
input anemia bp ptb n;
datalines;
0 0 0 496
0 0 1 74
0 1 0 113
0 1 1 25
1 0 0 85
1 0 1 15
1 1 0 15
1 1 1 3
;
run;

data dat;
set dat_cnts;
do i=1 to n;
	output;
end;
drop n;
run;

*****************************;
*
* Illustrative example      
*
*****************************;

**************************;
*** MAXIMUM LIKELIHOOD ***;
**************************;

*** Using proc logistic;
proc logistic data=dat outest=mle covout noprint;
model ptb(ref='0')= anemia bp;
run;
proc print data=mle noobs;
var _type_ _name_ intercept anemia bp;
run;

data mle_cov;
set mle;
where _type_="COV";
keep intercept anemia bp;
run;

*Information matrix - Hessian based estimator;
proc iml;
use mle_cov;
	read all var {intercept anemia bp} into cov;
close mle_cov;

infom = inv(cov)/826;
print infom;
quit;

*** Manually;
proc iml;
	use dat;
		read all var {ptb} into ptb;
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

*Save sample size;
n = nrow(ptb);

*Create data matrix;
xmat = j(n,1,1)||anemia||bp;

*Save number of parameters to be estimated;
q = ncol(xmat);

*Function for the log-likelihood;
start logl(beta) global(xmat, ptb);
	p = 1 / (1 + exp(-(xmat * beta`))); 
	logli = ptb#log(p) + (1-ptb)#log(1-p);
	logl = sum(logli);
	return (logl);
 finish logl;

*Optimize to find maximum;
beta = {0,0,0}; *initial values (set to zero);	
opt1 = {1,1}; *options to be used by optimizer, first 1 tells it to find maximum; 
call nlpnra(rc, bhat, "logl", beta, opt1);

*Point estimates;
print bhat;

*Obtain Hessian;
call nlpfdd(func, gradient, hessmat, "logl", bhat); 

*Information matrix;
infomat = -hessmat/n;
print infomat;

*Covariance matrix;
covmat = inv(n*infomat);
print covmat;

quit;


**************************;
***   M-ESTIMATION     ***;
**************************;

proc iml;
	use dat;
		read all var {ptb} into ptb;
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

*Save sample size;
n = nrow(ptb);

*Create data matrix;
xmat = j(n,1,1)||anemia||bp;

*Save number of parameters to be estimated;
q = ncol(xmat);

*Function for the score;
start ef(beta) global(xmat, ptb, p);
	p = 1 / (1 + exp(-(xmat * beta`))); 
	score = xmat` * (ptb - p);
	return (score);
finish ef;

*Root-finding;
beta = {0,0,0}; *initial values (set to zero);	
opt1 = q||1; *options to be used by root finding function; 
call nlplm(rc, bhat, "ef", beta, opt1);
print bhat;

*Bread;
opt2 = q||.||.; *options to be used for calculating Hessian;
call nlpfdd(func, hessmat, na, "ef", bhat, opt2); *using nlpfdd to calculate Hessian at beta hat;
print hessmat;

bread = - hessmat / n; *also the information matrix (Hessian estimator);
print bread;

*Meat;
ef1 = (ptb - p) # xmat[,1]; *note that p was output from our score function ef function above (at bhat);
ef2 = (ptb - p) # xmat[,2]; 
ef3 = (ptb - p) # xmat[,3];
efhat = ef1 || ef2 || ef3; 
meat = efhat` * efhat / n; *also the information matrix (outerproduct estimator);
print meat;

*Sandwich;
covmat = ( inv(bread) * meat * inv(bread)` ) / n; 
print covmat;
quit;


*************************************;
***   STANDARDIZATION EXAMPLE     ***;
*************************************;

proc iml;
	use dat;
		read all var {ptb} into ptb;
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

*Save sample size;
n = nrow(ptb);

*Create data matrix;
xmat = j(n,1,1)||anemia||bp;

*Save number of parameters to be estimated;
q = ncol(xmat) + 3;

*Function for the score;
start ef(theta) global(n, xmat, ptb, bp, p, out_r);
	beta = theta[1:3]`; 
	mu = theta[4:5]`;
	delta = theta[6];

	p = 1 / (1 + exp(-(xmat * beta`))); 
	score = xmat` * (ptb - p); 

	r1 = 1 / (1 + exp(-(beta[1] + beta[2]*1 + beta[3]*bp))) - mu[1]; 
	r0 = 1 / (1 + exp(-(beta[1] + beta[2]*0 + beta[3]*bp))) - mu[2];
	rd = j(n,1,(mu[1] - mu[2]) - delta[1]);

	out_r = r1||r0||rd; *output to use in sandwich later;

	toreturn = score//sum(r1)//sum(r0)//sum(rd); 
	return (toreturn);
finish ef;

*Root-finding;
theta = {-2,0,0,.1,.1,0}; *initial values (set to zero);	
opt1 = q||1; *options to be used by root finding function; 
call nlplm(rc, that, "ef", theta, opt1);
print that;

*Bread;
opt2 = q||.||.; *options to be used for calculating Hessian;
call nlpfdd(func, hessmat, na, "ef", that, opt2); *using nlpfdd to calculate Hessian at beta hat;
print hessmat;

bread = - hessmat / n; *also the information matrix (Hessian estimator);
print bread;

*Meat;
ef1 = (ptb - p) # xmat[,1]; *note that p was output from our score function ef function above (at bhat);
ef2 = (ptb - p) # xmat[,2]; 
ef3 = (ptb - p) # xmat[,3];
efhat = ef1 || ef2 || ef3 || out_r; 
meat = efhat` * efhat / n; *also the information matrix (outerproduct estimator);
print meat;

*Sandwich;
covmat = ( inv(bread) * meat * inv(bread)` ) / n; 
print covmat;
quit;

