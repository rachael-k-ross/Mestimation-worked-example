/*******************************************************************************************************************
M-estimation: introduction and applied examples

Rachael Ross(2023/06/08)
*******************************************************************************************************************/

/***********************************************
Example 1: Logistic Regression */
/***********************************************

/***********************************************
Loading Data */

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


/***********************************************
Rregression by MLE */
proc logistic data=dat;
model ptb(ref='0')= anemia bp;
ods output parameterestimates=ests_mle;
run;

data ests_mle;
set ests_mle(keep=variable estimate stderr);
lcl = estimate - 1.96*stderr;
ucl = estimate + 1.96*stderr;
run;

proc print data=ests_mle noobs;
run;

/***********************************************
M-estimator */
PROC IML;                            /*All steps are completed in PROC IML*/
	*Read data;
	use dat;								/*Open input data from above*/
		read all var {ptb} into ptb;		/*Read in each column as its own vector*/
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

	n = nrow(ptb);                        /*Save number of observations in data */

	/***********************************************
	Defining estimating equation */

	q = 3;								/*Save number parameters to be estimated*/

	START efunc(beta) global(ptb, anemia, bp);							/*Start to define estimating function */ 
		p = 1 / (1 + exp(-(beta[1] + beta[2]*anemia + beta[3]*bp)));	/*Predicted probability of the outcome */ 
		ef_1 = ptb - p;                                                 /*Estimating function for beta_0 (intercept) */ 
		ef_2 = (ptb - p)#anemia;                                        /*Estimating function for beta_1 */ 
		ef_3 = (ptb - p)#bp;                                            /*Estimating function for beta_2 */ 
		ef_mat = ef_1||ef_2||ef_3;
		RETURN(ef_mat);                         						/*Return n by q matrix for estimating functions*/
	FINISH efunc;                       								/*End definition of estimating equation*/

	START eequat(beta);					 								/*Start to define estimating equation (single argument)*/ 
		ef = efunc(beta);
		RETURN(ef[+,]);                  								/*Return column sums, 1 by q vector)*/
	FINISH eequat;                       								/*End definition of estimating equation*/

	/***********************************************
	Root-finding */
	initial = {-2,0,0};                 * Initial parameter values;
			/*Starting values need to be provided. A good starting value is within the plausible range 
			and not close to the bounds. For example, if the parameter is a risk then a starting value of 
			0.5 would be a good choice. For regression, one can generally provide starting values of 0. 
			To increase computational efficiency of M-estimation, subsets of estimating equations can be 
			solved separately and then used as the starting values for the overall estimating equations. 
			For example, in Example 2, we can obtain the point estimates for propensity score model parameters 
			using built-in functions/procedures for logistic regression use those as starting values.*/

	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   beta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	print beta_hat;

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                beta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	print bread;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(beta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	print filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	print sandwich;

	/***********************************************
	Formatting results for output */
	variable = {"Intercept","anemia","bp"};  
	est = beta_hat`;                    
	se = sqrt(vecdiag(sandwich));       /*Extract corresponding SE for each parameter*/
	lcl = est - 1.96*se; 				/*Calculated lcl*/
	ucl = est + 1.96*se;				/*Calculate ucl*/

	PRINT variable est se lcl ucl;     	/*Print information to the Results Viewer*/

	CREATE ests_mest VAR {variable est se lcl ucl};   /*Create an output data set called `out`*/
		APPEND;                         		  	  /*... that includes the parameter estimates, variance, and SE*/
	CLOSE ests_mest;                          		  /*Close the output*/
	QUIT;                                   
RUN;


/***********************************************
Results for logistic regression */

title1 "Logistic regression: MLE";
proc print data=ests_mle noobs; run;

title1 "Logistic regression: M-estimation";
proc print data=ests_mest noobs; run;

title1;


/***********************************************
Example 2: Standardization by IPW */
/***********************************************

/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use dat;								
		read all var {ptb} into ptb;		
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

	n = nrow(ptb);                        

	/***********************************************
	Defining estimating equation */

	q = 5;								

	START efunc(theta) global(ptb, anemia, bp, n);							
		pscore = 1 / (1 + exp(-(theta[1] + theta[2]*bp)));	 /*Predicted propensity score */ 
		ef_1 = anemia - pscore;								 /*Estimating function for alpha_0 */ 
		ef_2 = (anemia - pscore)#bp;						 /*Estimating function for alpha_1 */

		wt = anemia/pscore + (1-anemia)/(1-pscore);   		 /*IPW */
		ef_r0 = (1-anemia)#wt#ptb - theta[3]; 				 /*Estimating function for mu_0 */
		ef_r1 = anemia#wt#ptb - theta[4];					 /*Estimating function for mu_1 */


		ef_rd = j(n,1,(theta[4] - theta[3]) - theta[5]);	 /*Estimating function for delta */

		ef_mat = ef_1||ef_2||ef_r0||ef_r1||ef_rd;
		RETURN(ef_mat);                         						
	FINISH efunc;                       								

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {-2,0,.1,.1,0};           * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"Intercept","bp","risk0","risk1","rd"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_2 VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_2;                          		  
	QUIT;                                   
RUN;

/***********************************************
Example 3: Standardization by g-computation */
/***********************************************

/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use dat;								
		read all var {ptb} into ptb;		
		read all var {anemia} into anemia;
		read all var {bp} into bp;
	close dat;

	n = nrow(ptb);                        

	/***********************************************
	Defining estimating equation */

	q = 6;								

	START efunc(theta) global(ptb, anemia, bp, n);							
		p = 1 / (1 + exp(-(theta[1] + theta[2]*anemia + theta[3]*bp)));			/*Predicted probability of outcome */
		ef_1 = ptb - p;															/*Estimating function for beta_0 */
		ef_2 = (ptb - p)#anemia;												/*Estimating function for beta_1 */
		ef_3 = (ptb - p)#bp;													/*Estimating function for beta_2 */

		ef_r0 = 1 / (1 + exp(-(theta[1] + theta[2]*0 + theta[3]*bp))) - theta[4];	/*Estimating function for mu_0 */
		ef_r1 = 1 / (1 + exp(-(theta[1] + theta[2]*1 + theta[3]*bp))) - theta[5];	/*Estimating function for mu_1 */
		
		ef_rd = j(n,1,(theta[5] - theta[4]) - theta[6]);						/*Estimating function for delta */

		ef_mat = ef_1||ef_2||ef_3||ef_r0||ef_r1||ef_rd;
		RETURN(ef_mat);                         						
	FINISH efunc;                       								

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {-2,0,0,.1,.1,0};         * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"Intercept","anemia","bp","risk0","risk1","rd"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_3 VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_3;                          		  
	QUIT;                                   
RUN;

/***********************************************
Example 4: Data fusion */
/***********************************************

/***********************************************
Loading Data */

data datfusion_cnts;
input r y w n;
datalines;
1 1 1 0 
1 1 0 0
1 0 1 680
1 0 0 270
0 1 1 204
0 1 0 38
0 0 1 18
0 0 0 71
;
run;

data datfusion;
set datfusion_cnts;
do i=1 to n;
	output;
end;
drop n;
run;


/***********************************************
M-estimator */

PROC IML;                            
	*Read data;
	use datfusion;								
		read all var {r} into r;		
		read all var {y} into y;
		read all var {w} into w;
	close datfusion;

	n = nrow(r);                        

	/***********************************************
	Defining estimating equation */

	q = 4;								

	START efunc(theta) global(r, y, w, n);							
		ef_1 = r#(w - theta[1]);							/*Estimating function for misclassified prevalence */
		ef_2 = (1 - r)#y#(w - theta[2]);					/*Estimating function for sensitivity */
		ef_3 = (1 - r)#(1 - y)#((1 - w) - theta[3]);		/*Estimating function for specificity */
		ef_4 = j(n,1,theta[4]*(theta[2] - (1 - theta[3])) - (theta[1] - (1 - theta[3])));	/*Estimating function for prevalence (Rogan Gladen equation) */
		ef_mat = ef_1||ef_2||ef_3||ef_4;
		RETURN(ef_mat);                         						
	FINISH efunc;       

	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {.7,1,1,.7};       * Initial parameter values;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"pr_w","se","sp","pr_y"};  
	est = theta_hat`;                    
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est se lcl ucl;     	

	CREATE ests_4 VAR {variable est se lcl ucl};   
		APPEND;                         		  	  
	CLOSE ests_4;                          		  
	QUIT;                                   
RUN;
