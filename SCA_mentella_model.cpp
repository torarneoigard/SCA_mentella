#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // *** data reading
  DATA_INTEGER(minYear);      // first year of the assessment (population model)
  DATA_INTEGER(minAge);       // first age of the assessment (population model)
  DATA_INTEGER(CatchNrow);    // number of catch-at-age observations
  DATA_INTEGER(SurveyNrow);   // number of survey observations
  DATA_INTEGER(NageBlocks);   // number of age blocks in ageblocks file
  DATA_ARRAY(TotalCatches);   // Total catch in Tonnes
  DATA_ARRAY(CatchAtAge);     // Catch-at-age
  DATA_ARRAY(SurveyIndex);    // Survey indices
  DATA_ARRAY(WeightAtAge);    // year * age matrix of mean individual weight (in kg)
  DATA_ARRAY(MaturityAtAge);  // year * age matrix of mean individual maturity (0 to 1)
  DATA_ARRAY(AgeBlocks);      // predefined age blocks in ageblocks file
  DATA_INTEGER(nYears);       // number of years in the population model
  DATA_INTEGER(nAges);        // number of age in the population model
  DATA_INTEGER(nSurveys);     // number of surveys
  DATA_INTEGER(nSurveysProp); // number of surveys with proportion data
  DATA_VECTOR(SurveyTime);    // timing of the surveys (0 = beginning of the year to 1 = end of the year)
  DATA_VECTOR(SurveyTimeProp);// timing of the surveys with proportion data
  DATA_VECTOR(lowerAgeBoundary); // Lower limits for the a0 parameter of the selectivity curve (for each survey)
  DATA_VECTOR(upperAgeBoundary); // Upper limits for the a0 parameter for the selectivity curve (for each survey)
  
  DATA_INTEGER(REswitch);     // 0 fixed effect, 1 random effect
  // end of reading data

  // *** parameters declaration
  // dimensions and initial values of parameters are set in R
  PARAMETER_VECTOR(logNY1); // numbers of fish in year 1 (1992)
  PARAMETER_VECTOR(logNA1fe); // numbers of fish at age 1 (2y old), fixed effects model
  PARAMETER_VECTOR(DemlogFY); // Log of the demersal fleet fishing mortality
  PARAMETER_VECTOR(PellogFY); // Log of the pelagic fleet fishing mortality
  PARAMETER(pDema50); // probit of the 'a50' parameter for demersal fleet selectivity 
  PARAMETER(Demlogw); // scale parameter for demersal fleet selectivity 
  PARAMETER(pPela50); // probit of the 'a50' parameter for pelagic fleet selectivity
  PARAMETER(Pellogw); // scale parameter for pelagic fleet selectivity
  PARAMETER(DemlogVarLogC); // log of variance of log Catches-at-age for the demersal fleet
  PARAMETER(PellogVarLogC); // log of variance of log Catches-at-age for the pelagic fleet
  PARAMETER_VECTOR(logVarLogI); // log of variances of log Numbers-at-age for the surveys
  PARAMETER_VECTOR(logQSurvey); // log of the scaling coefficients for the surveys
  PARAMETER(Demsplus); // used to calculate demersal fleet selectivity for the +group, not estimated - obsolete?
  PARAMETER(Pelsplus); // used to calculate pelagic fleet selectivity for the +group, not estimated - obsolete?
  PARAMETER_VECTOR(pa0); // probit of the coefficients for the survey selectivity-at-age
  PARAMETER_VECTOR(pa0Prop); // probit of the coefficients for the survey selectivity-at-age for proportion data
  PARAMETER_VECTOR(logb1); // coefficients for the survey selectivity-at-age
  PARAMETER_VECTOR(logb2); // coefficients for the survey selectivity-at-age
  PARAMETER_VECTOR(logb1Prop);
  PARAMETER_VECTOR(logb2Prop);
  PARAMETER(logM2); // log of natural mortality
  PARAMETER(palogNA1);       // RE starts here, probit of the log of initial number (age 1 in year 1)
  PARAMETER(logSigmalogNA1); // variance of the AR process noise
  PARAMETER_VECTOR(ulogNA1); // noise itself
  // end of parameters declaration

  // *** preliminary calculations, indices and log-transformed data
  // probit transformations for bounded a50 coefficients (fleet selectivity)
  Type Dema50=Type(6.0)+(exp(pDema50)/(Type(1.0)+exp(pDema50)))*Type(13.0); // bounded between 6 and 19
  Type Pela50=Type(6.0)+(exp(pPela50)/(Type(1.0)+exp(pPela50)))*Type(13.0); // bounded between 6 and 19
  // probit transformations for bounded a0 coefficients (survey selectivity)
  vector<Type> a0(nSurveys);
  for (int i=0; i<(nSurveys); i++){
    a0(i) = lowerAgeBoundary(i)+(exp(pa0(i))/(Type(1.0)+exp(pa0(i))))*(upperAgeBoundary(i)); // bounded between lower and upper age in the data
  }

  // probit transformations for bounded a0 coefficients (survey selectivity at age for proportions)
  vector<Type> a0Prop(nSurveysProp);
  for (int i=0; i<(nSurveysProp); i++){
    a0Prop(i) = lowerAgeBoundary(i)+(exp(pa0Prop(i))/(Type(1.0)+exp(pa0Prop(i))))*(upperAgeBoundary(i)); // bounded between lower and upper age in the data
  }
  
  // probit transformations for bounded autoregressive recruitment model parameter
  Type alogNA1=Type(2)/(Type(1) + exp(-Type(2)*palogNA1)) - Type(1); // bounded between -1 and 1
  
  Rcout << "Dema50: " << Dema50 << "\n"; // prints out untransformed paramter values on screen
  Rcout << "Pela50: " << Pela50 << "\n";
  Rcout << "a0: " << a0 << "\n";
  
  //RE for NY1  // <- test code to set age distribution in year 1 as an AR(1) process. Not in use for the moment
  //Type SigmalogNY1 = exp(logSigmalogNY1);
  //vector<Type> tmplogNY1(nAges+1);     //Dummy vector initialised at age 1
  //vector<Type> logNY1(nAges);
  //tmplogNY1(0)=initlogNY1;
  //logNY1(0) = initlogNY1A1;
  //for (int i=1; i<nAges; i++){
    //tmplogNY1(i) = alogNY1*tmplogNY1(i-1)+SigmalogNY1*ulogNY1(i-1);
  //  logNY1(i) = alogNY1*logNY1(i-1)+SigmalogNY1*ulogNY1(i-1);
    //logNY1(i-1)=tmplogNY1(i);
  //   }

  // Random effect: creating the AR(1) process (for recruits)
  Type SigmalogNA1 = exp(logSigmalogNA1);
  vector<Type> tmplogNA1(nYears);
  vector<Type> logNA1re(nYears-1);
  tmplogNA1(0)=logNY1(0);
  for (int i=1; i<(nYears); i++){
    tmplogNA1(i) = alogNA1*tmplogNA1(i-1)+SigmalogNA1*ulogNA1(i-1);
    logNA1re(i-1)=tmplogNA1(i);
  }
  
  // recoding catch at age data
  vector <int> CatchYear(CatchNrow);
  vector <int> CatchAge(CatchNrow);
  vector<Type> TotalCatch(CatchNrow);
  vector<Type> DemCatch(CatchNrow);
  vector<Type> PelCatch(CatchNrow);
  vector<Type> TotallogCatch(CatchNrow);
  vector<Type> DemlogCatch(CatchNrow);
  vector<Type> PellogCatch(CatchNrow);
  for (int i=0; i<CatchNrow; i++){
    CatchYear(i) = CppAD::Integer(CatchAtAge(i,0));
    CatchAge(i) = CppAD::Integer(CatchAtAge(i,1));
    TotalCatch(i) = CatchAtAge(i,2);
    PelCatch(i) = CatchAtAge(i,3);
    DemCatch(i) = CatchAtAge(i,4);
    TotallogCatch(i)=log(TotalCatch(i));            // store Total Catch in log
    DemlogCatch(i)=log(DemCatch(i));	            // store Demersal Catch in log
    PellogCatch(i)=log(PelCatch(i)+Type(1));	    // store Pelagic Catch in log (+1 to avoid log(0)=-Inf) - zero catches should be set to NA instead
  }
  
  // recoding survey data
  vector<int> SurveyYear(SurveyNrow);
  vector<int> SurveyAge(SurveyNrow);
  vector<int> Survey(SurveyNrow);
  vector<Type> Index(SurveyNrow);
  vector<Type> logIndex(SurveyNrow);
  for (int i=0; i<SurveyNrow; ++i){
    SurveyYear(i) = CppAD::Integer(SurveyIndex(i,0));
    SurveyAge(i) = CppAD::Integer(SurveyIndex(i,1));
    Survey(i) = CppAD::Integer(SurveyIndex(i,2))-int(1);
    Index(i) = SurveyIndex(i,3);
    logIndex(i)=log(Index(i));		  	    // store survey index in log
  }

  // recoding catch in tonnes data
  vector<Type> CatchInTonnesTotal(nYears);
  vector<Type> logCatchInTonnesTotal(nYears);
  for (int y=0; y<nYears; ++y){
    CatchInTonnesTotal(y) = TotalCatches(y,1);
    logCatchInTonnesTotal(y)=log(CatchInTonnesTotal(y));		  	  	        // store catch in tonnes in log
  }  

  // *** Computation of FA's (fleet selectivities-at-age)
  // Demersal fleet
  Type x;
  Type tanhx;
  vector <Type> DemFA(nAges);
  vector <Type> logitDemFA(nAges); // can be removed?
  DemFA.setZero(); // fill in with zeroes
  for(int a=0;a<(6-minAge);++a){              // loop on ages before 6y, // can be removed?
    logitDemFA(a)=-20;// can be removed?
  }// can be removed?
  for(int a=(6-minAge); a<(nAges-1); ++a){    // loop on ages from age 6 to one before max age
    x=(Type(a+minAge)-Dema50)/exp(Demlogw);
    tanhx=(exp(x)-exp(-x))/(exp(x)+exp(-x));  // need to see if we can find a proper function in cpp for the tanh
    DemFA(a)=Type(0.5)*(Type(1.0)+tanhx);	      // fleet selectivity for ages 6+, demersal
    logitDemFA(a)=log((DemFA(a))/(1-DemFA(a)));// can be removed?
    }
  DemFA(nAges-1)=exp(Demsplus)/(Type(1.0)+exp(Demsplus));// fleet selectivity for the plus group, demersal
  logitDemFA(nAges-1)=log(DemFA(nAges-1)/(1-DemFA(nAges-1)));// can be removed?
  // *** Pelagic fleet
  vector <Type> PelFA(nAges);
  vector <Type> logitPelFA(nAges);// can be removed?
  PelFA.setZero(); // fill in with zeroes
  for(int a=0;a<(7-minAge);++a){              // loop on ages before 7y, // can be removed?
    logitPelFA(a)=-20;// can be removed?
  }// can be removed?
  for(int a=(7-minAge); a<(nAges-1); ++a){    // loop on ages from age 7 to one before max age
    x=(Type(a+minAge)-Pela50)/exp(Pellogw);
    Type tanhx;
    tanhx=(exp(x)-exp(-x))/(exp(x)+exp(-x));
    PelFA(a)=Type(0.5)*(Type(1.0)+tanhx);	      // fleet selectivity for ages 9+, pelagic
    logitPelFA(a)=log((PelFA(a))/(1-PelFA(a)));// can be removed?
    }
  PelFA(nAges-1)=exp(Pelsplus)/(Type(1.0)+exp(Pelsplus));// fleet selectivity for the plus group, pelagic
  logitPelFA(nAges-1)=log(PelFA(nAges-1)/(1-PelFA(nAges-1)));// can be removed?

  // *** Computation of FY's (fleet fishing mortality components by year)
  // Demersal fleet
  vector <Type> DemFY(nYears);
  for (int y=0; y<nYears; ++y){ // loop on years
    DemFY(y)=exp(DemlogFY(y));              // fishing mortality for the demersal fleet
  }
  // Pelagic fleet
  vector <Type> PelFY(nYears);
  PelFY.setZero();
  for (int y=(2006-minYear); y<nYears; ++y){             // loop on years 2006 to last year. (move the start to 2002 after comparison with ADMB results of AFWG2014)
    PelFY(y)=exp(PellogFY(y));                         // fishing mortality for the pelagic fleet
  }
  // *** Separable F's (outer_prod)  // check if this can be done with a propoer outer.prod function
  array<Type> DemF(nYears,nAges);
  array<Type> PelF(nYears,nAges);
  array<Type> F(nYears,nAges);
  for (int y=0; y<nYears ; ++y){
    for (int a=0; a<nAges ; ++a){
      DemF(y,a)=DemFY(y)*DemFA(a);
      PelF(y,a)=PelFY(y)*PelFA(a);
    }
  }
  F=DemF+PelF;						      // Matrix of mortality for both fleets combined

  // *** survey selectivities at age (maybe this needs a little documentation at some point)
  vector<Type> a1(nSurveys);  
  vector<Type> a2(nSurveys);
  vector<Type> b1(nSurveys);
  vector<Type> b2(nSurveys);
  vector<Type> c1(nSurveys);
  vector<Type> c2(nSurveys);
  
  for (int i=0; i<(nSurveys); i++){
    b1(i) = exp(logb1(i));
    b2(i) = exp(logb2(i));
    a1(i) = -b1(i)/(Type(2.0)*a0(i));
    a2(i) = -b2(i)/(Type(2.0)*a0(i));
    c1(i) = b1(i)*b1(i)/(Type(4.0)*a1(i));
    c2(i) = b2(i)*b2(i)/(Type(4.0)*a2(i));
  }

  array <Type> SA(nSurveys,nAges);
  array <Type> logSA(nSurveys,nAges);
  
  Type Ta;
  for (int i=0; i<(nSurveys); i++){
    for(int a=0; a<nAges; ++a){          // loop on ages 
      Ta=Type(a)+Type(2);
      SA(i,a) = ((exp(a1(i)*Ta*Ta+b1(i)*Ta+c1(i))*(Ta<a0(i)))+(exp(a2(i)*Ta*Ta+b2(i)*Ta+c2(i))*(Ta>=a0(i))))*Type(0.999)+Type(0.001);// Survey selectivities
      logSA(i,a) = log(SA(i,a));
    }
  }

  // *** survey selectivity at age for proportion data
  vector<Type> a1P(nSurveysProp);  
  vector<Type> a2P(nSurveysProp);
  vector<Type> b1P(nSurveysProp);
  vector<Type> b2P(nSurveysProp);
  vector<Type> c1P(nSurveysProp);
  vector<Type> c2P(nSurveysProp);
  
  for (int i=0; i<(nSurveysProp); i++){
    b1P(i) = exp(logb1Prop(i));
    b2P(i) = exp(logb2Prop(i));
    a1P(i) = -b1P(i)/(Type(2.0)*a0Prop(i));
    a2P(i) = -b2P(i)/(Type(2.0)*a0Prop(i));
    c1P(i) = b1P(i)*b1P(i)/(Type(4.0)*a1P(i));
    c2P(i) = b2P(i)*b2P(i)/(Type(4.0)*a2P(i));
  }
  
  array <Type> SAProp(nSurveysProp,nAges);  
  for (int i=0; i<(nSurveysProp); i++){
    for(int a=0; a<nAges; ++a){          // loop on ages 
      Ta=Type(a)+Type(2);
      SAProp(i,a) = ((exp(a1P(i)*Ta*Ta+b1P(i)*Ta+c1P(i))*(Ta<a0Prop(i)))+(exp(a2P(i)*Ta*Ta+b2P(i)*Ta+c2P(i))*(Ta>=a0Prop(i))))*Type(0.999)+Type(0.001);// Survey selectivities
    }
  }

  // *** Natural mortality
  Type M2=exp(logM2);      		                      

  // *** Fill in the N matrix
  // initial conditions (year 1 and age 1)
  array<Type> logN(nYears,nAges); // log-numbers in the population matrix
  array<Type> logTriN(nYears+1,nAges+nYears); // extended population matrix (1 extra year and as many extra ages as the number of years of observations)
  array<Type> TriN(nSurveysProp,nYears+1,nAges+nYears);
  TriN.setZero();
  
  for(int a=0; a<nAges; ++a){
    logN(0,a)=logNY1(a);                    // fill in first row of logN with logNY1
    logTriN(0,a)=logNY1(a);                 // fill in first row of logN with logNY1
    //TriN(0,a) = exp(logNY1(a));
    }
 
  for(int y=1; y<nYears; ++y){  // loop on years, i.e. rows (start on second year)
    if(REswitch < 1){ // fixed effects on recruits
      logN(y,0)=logNA1fe(y-1);// fill in first column of logN with logNA1
      logTriN(y,0)=logNA1fe(y-1);// fill in first column of logN with logNA1
      //TriN(y,0)=exp(logNA1fe(y-1));
    }
    
    if(REswitch > 0){ // random effects on recruits
      logN(y,0)=logNA1re(y-1);// fill in first line of logN with logNA1
      logTriN(y,0)=logNA1re(y-1);// fill in first line of logN with logNA1
      //TriN(y,0)=exp(logNA1re(y-1));
    }
    
    for(int a=1; a<(nAges-1); ++a){    // loop on ages (start at second age column and end at one before last)
      logN(y,a)=logN(y-1,a-1)-F(y-1,a-1)-M2;// fill in logN for age a and year y
      logTriN(y,a) = logN(y,a);// fill in logN for age a and year y
      //TriN(y,0)=exp(logNA1re(y-1));
    }
    
    // plus group in the logN matrix
    logN(y,(nAges-1))=log(exp(logN(y-1,(nAges-1))-F(y-1,(nAges-1))-M2)+exp(logN(y-1,(nAges-2))-F(y-1,(nAges-2))-M2)); // +group
    
    // filling the extended population matrix for older age groups, including the +group
    for(int a = (nAges-1); a <(nAges+y); ++a){
      logTriN(y,a) = logTriN(y-1,a-1)-F(y-1,nAges-1)-M2;
      //TriN(y,a) = exp(logTriN(y-1,a-1)-F(y-1,nAges-1)-M2);
    }
  }

  
  // filling the extended population matrix for the last year
  for(int a = 1; a <(nAges+nYears); ++a){
      logTriN(nYears,a) = logTriN(nYears-1,a-1)-F(nYears-1,nAges-1)-M2;
      //TriN(nYears,a) = exp(logTriN(nYears-1,a-1)-F(nYears-1,nAges-1)-M2);
    }

  
  // create a matrix based on the triangular matrix used for the proportions
  array <Type> scalingFactor(nSurveysProp,nYears);
  scalingFactor.setZero(); 
  for(int s=0; s<nSurveysProp; ++s){
    for(int y=0; y<(nYears+1); ++y){;
      for(int a=0; a<(nAges+y); ++a){
	//TriN(y,a) = exp(logTriN(y,a));
	if(a<=(nAges-1)) TriN(s,y,a) = SAProp(s,a)*exp(logTriN(y,a))*exp(-SurveyTimeProp(s)*(M2+F(y,a)));
	if(a>(nAges-1)) TriN(s,y,a) = SAProp(s,nAges-1)*exp(logTriN(y,a))*exp(-SurveyTimeProp(s)*(M2+F(y,nAges-1)));
	scalingFactor(s,y) = scalingFactor(s,y) + TriN(s,y,a);
      }
    }
  }

  // create the proportions matrix
  array<Type> IndexProp(nSurveysProp,nYears+1,nAges+nYears);
  IndexProp.setZero();
  for(int s=0; s<nSurveysProp; ++s){
    for(int y=0; y<(nYears+1); ++y){;
      for(int a=0; a<(nAges+y); ++a){
	IndexProp(s,y,a) = TriN(s,y,a)/scalingFactor(s,y);
      }
    }
  }

  
  // create truncated proportion matrix. Truncated according to the predefined age blocks
  array<Type> IndexPropTruncMat(nSurveysProp,nYears+1,NageBlocks);
  IndexPropTruncMat.setZero();
  for(int s=0; s<nSurveysProp; ++s){
    for (int y=0; y<(nYears+1); ++y){
      //Find max age
      //maxAgeMat <- y - 1992 + 19;
      for (int ageBl = 0; ageBl<(NageBlocks); ++ageBl){
	for(int age = CppAD::Integer(AgeBlocks(ageBl,1)); age < CppAD::Integer(AgeBlocks(ageBl,2)); ++age){
	  IndexPropTruncMat(s,y,ageBl) = IndexPropTruncMat(s,y,ageBl) + IndexProp(s,y,age-minAge);
	}
      }
    }
  }
  
  // *** initialise nll
  Type nll=0;

  // *** predict/estimate Catches and compute nll component
  vector <Type> DemPredlogC(CatchNrow);
  vector <Type> PelPredlogC(CatchNrow);
  array <Type> predCmatrix(nYears,nAges);
  predCmatrix.setZero();                            // fill in with zeroes
  Type DemVarlogC=exp(DemlogVarLogC); 	            // variance of logCatch, demersal
  Type PelVarlogC=exp(PellogVarLogC);	            // variance of logCatch, pelagic
  for(int i=0;i<CatchNrow; ++i){	            // loop on the Catch data
    int sy=CatchYear(i)-minYear;
    int sa=CatchAge(i)-minAge;
    
    if(DemCatch(i)>0){	    		            // Check if catch data exists (i.e. >0)
      DemPredlogC(i)=log(DemF(sy,sa))-log(F(sy,sa)+M2)+log(Type(1.0)-exp(-F(sy,sa)-M2))+logN(sy,sa); // compute predicted logCatch for all years and all ages, demersal
      nll+=-dnorm(DemlogCatch(i),DemPredlogC(i),DemVarlogC,true); // negative log likelihood for catches at age, demersal
      predCmatrix(sy,sa)=exp(DemPredlogC(i)); // add demersal catch to predicted total catch matrix
    }
    if(PelCatch(i)>0){
      PelPredlogC(i)=log(PelF(sy,sa))-log(F(sy,sa)+M2)+log(Type(1.0)-exp(-F(sy,sa)-M2))+logN(sy,sa);  // compute predicted logCatch for all years and all ages, pelagic
      nll+=-dnorm(PellogCatch(i),PelPredlogC(i),sqrt(PelVarlogC),true); // negative log likelihood for catches at age, pelagic
      predCmatrix(sy,sa)+=exp(PelPredlogC(i)); // add pelagic catch to predicted total catch matrix
    }
  } 
  
  Type nll1 = nll; // likelihood component for the catch-at-age
  
  // *** predict survey indices and compute nll component
  vector<Type> VarLogI(nSurveys);			      // variance of logSurvey
  for(int i=0; i<nSurveys; ++i){
    VarLogI(i)=exp(logVarLogI(i));
  }
  vector<Type> predlogI(SurveyNrow);
  for(int i=0; i<SurveyNrow; ++i){ 		      // loop on survey data
    int sy=SurveyYear(i)-minYear; // index of survey year
    int sa=SurveyAge(i)-minAge; // index of survey age
      predlogI(i)=logQSurvey(Survey(i))+log(SA(Survey(i),sa))-(F(sy,sa)+M2)*SurveyTime(Survey(i))+logN(sy,sa); //estimate logSurvey index for year y and age a
      nll+=-dnorm(logIndex(i),predlogI(i),sqrt(VarLogI(Survey(i))),true); //increase total nll estimate by the contribution of survey at year y and age a    
  }
   Type nll2 = nll - nll1; // likelihood component for the survey indices

  // *** predict total catches in tonnes and compute nll component
  array<Type> CatchMatrix(nYears,nAges); // matrix of catches in kg per year/age
  vector<Type> PredTotalCatches(nYears); // vector of catches in tonnes per year
  PredTotalCatches.setZero(); // initialise vector
  for(int y=0; y<nYears; ++y){
    for(int a=0; a<nAges; ++a){
      CatchMatrix(y,a)=predCmatrix(y,a)*WeightAtAge(y,a);
      PredTotalCatches(y)+=CatchMatrix(y,a)/Type(1000); // Total Catch in tonnes
    }
    nll+=-dnorm(logCatchInTonnesTotal(y),log(PredTotalCatches(y)),sqrt(Type(0.001)),true);//increase total nll estimate by the contribution of catches in year y. variance of log-catches in assumed to be 0.001.
  }
  Rcout << "PredTotalCatches " << PredTotalCatches << "\n";
  
  Type nll3 = nll - nll1 - nll2; // likelihood component for the total catches in tonnes
  

/*<- test code to set age distribution in year 1 as an AR(1) process. Not in use for the moment
//Likelihood contribution for the RE logNY1
 for(int i=0;i<(nAges-1);i++)
  {
    nll += -dnorm(ulogNY1(i),Type(0),Type(1),true);
  }
*/

 //Likelihood contribution for the RE logNA1
   for(int i=0;i<(nYears-1);i++)
   {
     nll += -dnorm(ulogNA1(i),Type(0),Type(1),true);
   }
  // nll calculation end here

   
  // rest of the code is only for additional outputs
  // Calculating SSB
  array<Type> SSBmatrix(nYears,nAges);
  vector<Type> logSSB(nYears);
  vector<Type> SSB(nYears);
  SSB.setZero();
  for(int y=0; y<nYears; ++y){
    for(int a=0; a<nAges; ++a){
      SSBmatrix(y,a)=exp(logN(y,a))*MaturityAtAge(y,a)*WeightAtAge(y,a);
      SSB(y)+=SSBmatrix(y,a)/Type(1000); // SSB in tonnes
    }
    logSSB(y)=log(SSB(y));
  }

  // *** addtional outputs
  vector<Type> NY1(nAges); // age distribution in the first year
  NY1=exp(logNY1);
  vector<Type> NA1(nYears-1); // recruitment (number at age 2y) across years
  vector<Type> logNA1(nYears-1); // log-recruitment (number at age 2y) across years
  if(REswitch < 1){
    logNA1=logNA1fe;
  }
  if(REswitch > 0){
    logNA1=logNA1re;
  }
  NA1=exp(logNA1);
  
  vector<Type> RecAge6(nYears); // recruitment at age 6
  for (int y=0; y<nYears; ++y){
    RecAge6(y)=exp(logN(y,6));
  }

  // *** reporting
  ADREPORT(logSSB);              // report logSSB into sdreport
  //ADREPORT(SSB);              // report SSB into sdreport
  ADREPORT(PredTotalCatches); // report total catches in tonnes
  //ADREPORT(logNY1);
  //ADREPORT(NY1);
  ADREPORT(logNA1);
  //ADREPORT(NA1);
  ADREPORT(RecAge6);
  //ADREPORT(PelFY);
  //ADREPORT(DemFY);
  //ADREPORT(PelFA);
  //ADREPORT(DemFA);
  ADREPORT(logitDemFA);
  ADREPORT(logitPelFA);
  ADREPORT(SA);
  ADREPORT(SAProp);
  ADREPORT(logSA);
  ADREPORT(M2);
  ADREPORT(nll1);
  ADREPORT(nll2);
  ADREPORT(nll3);
  //ADREPORT(nll);
  //ADREPORT(logNY1);
  //if(REswitch > 0){
  //  ADREPORT(logNA1re);
  //}
  //if(REswitch > 0){
  //  ADREPORT(alogNA1);
  //}
  ADREPORT(logTriN);
  // returning negative log-likelihood
  return nll;

}
