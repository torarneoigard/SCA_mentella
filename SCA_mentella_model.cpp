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
  DATA_ARRAY(TotalCatches);   // Total catch in Tonnes
  DATA_ARRAY(CatchAtAge);     // Catch-at-age
  DATA_ARRAY(SurveyIndex);    // Survey indices
  DATA_ARRAY(WeightAtAge);    // year * age matrix of mean individual weight (in kg)
  DATA_ARRAY(MaturityAtAge);  // year * age matrix of mean individual maturity (0 to 1)
  DATA_INTEGER(nYears);       // number of years in the population model
  DATA_INTEGER(nAges);        // number of age in the population model
  DATA_INTEGER(nSurveys);     // number of surveys
  DATA_VECTOR(SurveyTime);    // timing of the surveys (0 = beginning of the year to 1 = end of the year)
  DATA_VECTOR(lowerAgeBoundary); // Lower limits for the a0 parameter of the selectivity curve (for each survey)
  DATA_VECTOR(upperAgeBoundary); // Upper limits for the a0 parameter for the selectivity curve (for each survey)
  
  DATA_INTEGER(REswitch);     // 0 fixed effect, 1 random effect
  // end of reading data

  // *** parameters declaration
  // dimensions and initial values of parameters are set in R
  PARAMETER_VECTOR(logNY1);
  PARAMETER_VECTOR(logNA1);
  PARAMETER_VECTOR(DemlogFY);
  PARAMETER_VECTOR(PellogFY);
  PARAMETER(pDema50);
  PARAMETER(Demlogw);
  PARAMETER(pPela50);
  PARAMETER(Pellogw);
  PARAMETER(DemlogVarLogC);
  PARAMETER(PellogVarLogC);
  PARAMETER_VECTOR(logVarLogI);
  PARAMETER_VECTOR(logQSurvey);
  //PARAMETER(logQSurvey1);
  //PARAMETER(logQSurvey2);
  //PARAMETER(logQSurvey3);
  PARAMETER(Demsplus);
  PARAMETER(Pelsplus);
  //PARAMETER(pa0Winter);
  //PARAMETER(logb1Winter);
  //PARAMETER(logb2Winter);
  //PARAMETER(pa0Eco);
  //PARAMETER(logb1Eco);
  //PARAMETER(logb2Eco);
  //PARAMETER(pa0Russian);
  //PARAMETER(logb1Russian);
  //PARAMETER(logb2Russian);
  PARAMETER_VECTOR(pa0);
  PARAMETER_VECTOR(logb1);
  PARAMETER_VECTOR(logb2);
  PARAMETER(logM2);
  //PARAMETER(palogNY1);     // RE starts here
  //PARAMETER(logSigmalogNY1);
  //PARAMETER_VECTOR(ulogNY1);
  PARAMETER(palogNA1);
  PARAMETER(logSigmalogNA1);
  PARAMETER_VECTOR(ulogNA1);
  //PARAMETER(initlogNY1A1);
  // end of parameters declaration

  // *** preliminary calculations, indices and log-transformed data
  // probit transformations for bounded a50 coefficients (selectivity)
  Type Dema50=Type(6.0)+(exp(pDema50)/(Type(1.0)+exp(pDema50)))*Type(13.0); // bounded between 6 and 19
  Type Pela50=Type(6.0)+(exp(pPela50)/(Type(1.0)+exp(pPela50)))*Type(13.0); // bounded between 6 and 19
  //Type a0Winter=Type(2.0)+(exp(pa0Winter)/(Type(1.0)+exp(pa0Winter)))*Type(13.0); // bounded between 2 and 15
  //Type a0Eco=Type(2.0)+(exp(pa0Eco)/(Type(1.0)+exp(pa0Eco)))*Type(13.0); // bounded between 2 and 15
  //Type a0Russian=Type(2.0)+(exp(pa0Russian)/(Type(1.0)+exp(pa0Russian)))*Type(9.0); // bounded between 2 and 11
  //Type alogNY1=Type(2)/(Type(1) + exp(-Type(2)*palogNY1)) - Type(1); // bounded between -1 and 1
  Type alogNA1=Type(2)/(Type(1) + exp(-Type(2)*palogNA1)) - Type(1); // bounded between -1 and 1
  
  vector<Type> a0(nSurveys);
  for (int i=0; i<(nSurveys); i++){
    a0(i) = lowerAgeBoundary(i)+(exp(pa0(i))/(Type(1.0)+exp(pa0(i))))*(upperAgeBoundary(i));
  }
  
  Rcout << "Dema50: " << Dema50 << "\n";
  Rcout << "Pela50: " << Pela50 << "\n";
  //Rcout << "a0Winter: " << a0Winter << "\n";
  //Rcout << "a0Eco: " << a0Eco << "\n";
  //Rcout << "a0Russian: " << a0Russian << "\n";  
  Rcout << "a0: " << a0 << "\n";
  
  //RE for NY1
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
  
 //if(REswitch > 0){   
    Type SigmalogNA1 = exp(logSigmalogNA1);
    vector<Type> tmplogNA1(nYears);
    vector<Type> logNA1re(nYears-1);
    tmplogNA1(0)=logNY1(0);
    for (int i=1; i<(nYears); i++){
      tmplogNA1(i) = alogNA1*tmplogNA1(i-1)+SigmalogNA1*ulogNA1(i-1);
      logNA1re(i-1)=tmplogNA1(i);
    }
 //} 
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
    PellogCatch(i)=log(PelCatch(i)+Type(0.0001));	    // store Pelagic Catch in log (+0.0001 to avoid log(0)=-Inf)
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
  //Int nAgesinSurveys = CppAD::Integer(max(SurveyAge)-min(SurveyAge)+1);

  // recoding catch in tonnes data
  vector<Type> CatchInTonnesTotal(nYears);
  vector<Type> logCatchInTonnesTotal(nYears);
  for (int y=0; y<nYears; ++y){
    CatchInTonnesTotal(y) = TotalCatches(y,1);
    logCatchInTonnesTotal(y)=log(CatchInTonnesTotal(y));		  	  	        // store catch in tonnes in log
  }  

  // *** Computation of FA's
  // Demersal fleet
  Type x;
  Type tanhx;
  vector <Type> DemFA(nAges);
  vector <Type> logitDemFA(nAges);
  DemFA.setZero(); // fill in with zeroes
  for(int a=0;a<(6-minAge);++a){              // loop on ages before 6y
    logitDemFA(a)=-20;
  }
  for(int a=(6-minAge); a<(nAges-1); ++a){    // loop on ages from age 6 to one before max age
    x=(Type(a+minAge)-Dema50)/exp(Demlogw);
    tanhx=(exp(x)-exp(-x))/(exp(x)+exp(-x));
    DemFA(a)=Type(0.5)*(Type(1.0)+tanhx);	      // fleet selectivity for ages 6+, demersal
    logitDemFA(a)=log((DemFA(a))/(1-DemFA(a)));
    }
  DemFA(nAges-1)=exp(Demsplus)/(Type(1.0)+exp(Demsplus));// fleet selectivity for the plus group, demersal
  logitDemFA(nAges-1)=log(DemFA(nAges-1)/(1-DemFA(nAges-1)));
  // *** Pelagic fleet
  vector <Type> PelFA(nAges);
  vector <Type> logitPelFA(nAges);
  PelFA.setZero(); // fill in with zeroes
  for(int a=0;a<(9-minAge);++a){              // loop on ages before 9y
    logitPelFA(a)=-20;
  }
  for(int a=(9-minAge); a<(nAges-1); ++a){    // loop on ages from age 9 to one before max age
    x=(Type(a+minAge)-Pela50)/exp(Pellogw);
    Type tanhx;
    tanhx=(exp(x)-exp(-x))/(exp(x)+exp(-x));
    PelFA(a)=Type(0.5)*(Type(1.0)+tanhx);	      // fleet selectivity for ages 9+, pelagic
    logitPelFA(a)=log((PelFA(a))/(1-PelFA(a)));
    }
  PelFA(nAges-1)=exp(Pelsplus)/(Type(1.0)+exp(Pelsplus));// fleet selectivity for the plus group, pelagic
  logitPelFA(nAges-1)=log(PelFA(nAges-1)/(1-PelFA(nAges-1)));

  // *** Computation of FY's
  vector <Type> DemFY(nYears);
  for (int y=0; y<nYears; ++y){ // loop on years
    DemFY(y)=exp(DemlogFY(y));              // fishing mortality for the demersal fleet
  }
  //Rcout << "DemFY" << DemFY << "\n";
  vector <Type> PelFY(nYears);
  PelFY.setZero();
  for (int y=(2006-minYear); y<nYears; ++y){             // loop on years 2006 to last year. (move the start to 2002 after comparison with ADMB results of AFWG2014)
    PelFY(y)=exp(PellogFY(y));                         // fishing mortality for the pelagic fleet
  }
  //Rcout << "PelFY" << PelFY << "\n";
  // *** Separable F's (outer_prod)
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

  // *** survey selectivities at age
  //Type b1Winter=exp(logb1Winter);
  //Type b1Eco=exp(logb1Eco);
  //Type b1Russian=exp(logb1Russian);
  //Type b2Winter=exp(logb2Winter);
  //Type b2Eco=exp(logb2Eco);
  //Type b2Russian=exp(logb2Russian);
  //Type a1Winter=-b1Winter/(Type(2)*a0Winter);
  //Type a1Eco=-b1Eco/(Type(2)*a0Eco);
  //Type a1Russian=-b1Russian/(Type(2)*a0Russian);
  //Type a2Winter=-b2Winter/(Type(2)*a0Winter);
  //Type a2Eco=-b2Eco/(Type(2)*a0Eco);
  //Type a2Russian=-b2Russian/(Type(2)*a0Russian);
  //Type c1Winter=b1Winter*b1Winter/(Type(4)*a1Winter);
  //Type c1Eco=b1Eco*b1Eco/(Type(4)*a1Eco);
  //Type c1Russian=b1Russian*b1Russian/(Type(4)*a1Russian);
  //Type c2Winter=b2Winter*b2Winter/(Type(4)*a2Winter);
  //Type c2Eco=b2Eco*b2Eco/(Type(4)*a2Eco);
  //Type c2Russian=b2Russian*b2Russian/(Type(4)*a2Russian);
  
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
  
  //Rcout << "logb2Winter" << logb2Winter << "\n";
  Rcout << "logb2" << logb2 << "\n";
  //vector <Type> SAWinter(nAges);
  //vector <Type> SAEco(nAges);
  //vector <Type> SARussian(nAges);
  //vector <Type> logSAWinter(nAges);
  //vector <Type> logSAEco(nAges);
  //vector <Type> logSARussian(nAges);
  
  array <Type> SA(nSurveys,nAges);
  array <Type> logSA(nSurveys,nAges);
  
  Type Ta;
  for (int i=0; i<(nSurveys); i++){
    for(int a=0; a<nAges; ++a){          // loop on ages 
      Ta=Type(a)+Type(2);
      //SAWinter(a)=((exp(a1Winter*Ta*Ta+b1Winter*Ta+c1Winter)*(Ta<a0Winter))+(exp(a2Winter*Ta*Ta+b2Winter*Ta+c2Winter)*(Ta>=a0Winter)))*Type(0.999)+Type(0.001);// Winter survey selectivity
      //logSAWinter(a)=log(SAWinter(a));
      //SAEco(a)=((exp(a1Eco*Ta*Ta+b1Eco*Ta+c1Eco)*(Ta<a0Eco))+(exp(a2Eco*Ta*Ta+b2Eco*Ta+c2Eco)*(Ta>=a0Eco)))*Type(0.999)+Type(0.001);// Ecosystem survey selectivity
      //logSAEco(a)=log(SAEco(a));
      //SARussian(a)=((exp(a1Russian*Ta*Ta+b1Russian*Ta+c1Russian)*(Ta<a0Russian))+(exp(a2Russian*Ta*Ta+b2Russian*Ta+c2Russian)*(Ta>=a0Russian)))*Type(0.999)+Type(0.001);// Russian survey selectivity
      //logSARussian(a)=log(SARussian(a));
      SA(i,a) = ((exp(a1(i)*Ta*Ta+b1(i)*Ta+c1(i))*(Ta<a0(i)))+(exp(a2(i)*Ta*Ta+b2(i)*Ta+c2(i))*(Ta>=a0(i))))*Type(0.999)+Type(0.001);// Winter survey selectivity
      logSA(i,a) = log(SA(i,a));
    }
  }

  // *** Natural mortality
  Type M2=exp(logM2);      		                      

  // *** Fill in the N matrix
  // initial conditions (year 1 and age 1)
  array<Type> logN(nYears,nAges);
  array<Type> logTriN(nYears+1,nAges+nYears);
  //logTriN = Type(-100000000);
  for(int a=0; a<nAges; ++a){
    logN(0,a)=logNY1(a);                    // fill in first column of logN with logNY1
    logTriN(0,a)=logNY1(a);
    }
 
  for(int y=1; y<nYears; ++y){  // loop on years (start on second year)
    if(REswitch < 1){
      logN(y,0)=logNA1(y-1);// fill in first line of logN with logNA1
      logTriN(y,0)=logNA1(y-1);
    }
    if(REswitch > 0){
      logN(y,0)=logNA1re(y-1);// fill in first line of logN with logNA1
      logTriN(y,0)=logNA1re(y-1);
    }
    
    for(int a=1; a<(nAges-1); ++a){    // loop on ages (start at second age column and end at one before last)
      logN(y,a)=logN(y-1,a-1)-F(y-1,a-1)-M2;// fill in logN for age a and year y
      logTriN(y,a) = logN(y,a);
    }
    logN(y,(nAges-1))=log(exp(logN(y-1,(nAges-1))-F(y-1,(nAges-1))-M2)+exp(logN(y-1,(nAges-2))-F(y-1,(nAges-2))-M2)); // +group
    
    
    for(int a = (nAges-1); a <(nAges+y); ++a){
      logTriN(y,a) = logTriN(y-1,a-1)-F(y-1,nAges-1)-M2;
    }
  }
  
  //The last year
  for(int a = 1; a <(nAges+nYears); ++a){
      logTriN(nYears,a) = logTriN(nYears-1,a-1)-F(nYears-1,nAges-1)-M2;
    }
  
  //vector<Type> replogTriN = logTriN(nYears,);
  
  //Rcout << "N" << exp(logN) << "\n";


  // *** initialise nll
  Type nll=0;

  // *** predict Catches and compute nll component
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
  Rcout << "DemVarlogC " << DemVarlogC << "\n";

  Rcout << "PelVarlogC " << PelVarlogC << "\n";

  Rcout << "nnl1  " << nll << "\n";
  //Rcout << "Survey: " << Survey << "\n";
  Type nll1 = nll;
  
  // *** predict survey indices and compute nll component
  //vector <Type> logQSurvey(3);
  //logQSurvey(0)=logQSurvey1;
  //logQSurvey(1)=logQSurvey2;
  //logQSurvey(2)=logQSurvey3;
  //Type Sel=0;                                          // survey selectivity
  vector<Type> VarLogI(nSurveys);			      // variance of logSurvey
  for(int i=0; i<nSurveys; ++i){
    VarLogI(i)=exp(logVarLogI(i));
  }
  vector<Type> predlogI(SurveyNrow);
  for(int i=0; i<SurveyNrow; ++i){ 		      // loop on survey data
    int sy=SurveyYear(i)-minYear;
    int sa=SurveyAge(i)-minAge;
      /*if(Survey(i)<1){	      	   	     	              // check if Winter survey 
            Sel=SAWinter(sa);                        // set selectivity for Winter survey
          }
      if((Survey(i)>0)&(Survey(i)<2)){	      	   	     	              // check if Ecosystem survey 
            Sel=SAEco(sa);                           // set selectivity for Ecosystem survey
          }
      if(Survey(i)>1){                                       // check if Russian survey
          Sel=SARussian(sa);                         // set selectivity for Russian Survey
          }
       */
      //predlogI(i)=logQSurvey(Survey(i))+log(Sel)-(F(sy,sa)+M2)*SurveyTime(Survey(i))+logN(sy,sa); //estimate logSurvey index for year y and age a
      predlogI(i)=logQSurvey(Survey(i))+log(SA(Survey(i),sa))-(F(sy,sa)+M2)*SurveyTime(Survey(i))+logN(sy,sa); //estimate logSurvey index for year y and age a
      
      nll+=-dnorm(logIndex(i),predlogI(i),sqrt(VarLogI(Survey(i))),true); //increase total nll estimate by the contribution of survey at year y and age a    
  }

   Type nll2 = nll - nll1;
   Rcout << "nnl2  " << nll2 << "\n";

  // *** predict total catches and compute nll component
  array<Type> CatchMatrix(nYears,nAges);
  vector<Type> PredTotalCatches(nYears);
  PredTotalCatches.setZero();
  for(int y=0; y<nYears; ++y){
    for(int a=0; a<nAges; ++a){
      CatchMatrix(y,a)=predCmatrix(y,a)*WeightAtAge(y,a);
      PredTotalCatches(y)+=CatchMatrix(y,a)/Type(1000); // Total Catch in tonnes
    }
    nll+=-dnorm(logCatchInTonnesTotal(y),log(PredTotalCatches(y)),sqrt(Type(0.001)),true);//increase total nll estimate by the contribution of catches in year y. variance of log-catches in assumed to be 0.001.
  }
  Rcout << "PredTotalCatches " << PredTotalCatches << "\n";
  
  Type nll3 = nll - nll1 - nll2;
  Rcout << "nnl3  " << nll3 << "\n";
/*
//Likelihood contribution for the RE logNY1
 for(int i=0;i<(nAges-1);i++)
  {
    nll += -dnorm(ulogNY1(i),Type(0),Type(1),true);
  }
*/

 //Likelihood contribution for the RE logNA1
// if(REswitch > 0){
   for(int i=0;i<(nYears-1);i++)
   {
     nll += -dnorm(ulogNA1(i),Type(0),Type(1),true);
   }
// }

   
  // *** Calculating SSB
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
  vector<Type> NY1(nAges); // first age
  NY1=exp(logNY1);
  vector<Type> NA1(nYears-1); // first year
  if(REswitch < 1){
    NA1=exp(logNA1);
  }
  if(REswitch > 0){
    NA1=exp(logNA1re);
  }
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
  //ADREPORT(logNA1);
  //ADREPORT(NA1);
  ADREPORT(RecAge6);
  //ADREPORT(PelFY);
  //ADREPORT(DemFY);
  //ADREPORT(PelFA);
  //ADREPORT(DemFA);
  ADREPORT(logitDemFA);
  ADREPORT(logitPelFA);
  ADREPORT(SA);
  ADREPORT(logSA);
  //NEED TO FIX ADREPORT
  //ADREPORT(SAWinter);
  //ADREPORT(SAEco);
  //ADREPORT(SARussian);
  
  
  //ADREPORT(logSAWinter);
  //ADREPORT(logSAEco);
  //ADREPORT(logSARussian);
  ADREPORT(M2);
  ADREPORT(nll1);
  ADREPORT(nll2);
  ADREPORT(nll3);
  //ADREPORT(logNY1);
  if(REswitch > 0){
    ADREPORT(logNA1re);
  }
  //ADREPORT(alogNY1);
  if(REswitch > 0){
    ADREPORT(alogNA1);
  }
   
  //ADREPORT(nll);
  ADREPORT(logTriN);
  // returning negative log-likelihood
  return nll;

}
