// multilikeihood model for malaria disaggregation
// Given point parasite surveys and country/admin level cases data, estimate pixel level malaria prevalence.
// 
// Author: Tim Lucas, Ewan Cameron
// Date: 2018-05-24
// For comparison paper

// Requirements
// To be compiled/optimised with TMB - https://cran.r-project.org/web/packages/TMB/index.html
// Uses a mesh built with INLA - http://www.r-inla.org/
// hopefully small enough that tmbstan can give good posterior

// Table of contents
// Data: Spatial field mesh and matrices, polygon data, point data
// Read in parameter starting values and set priors: regression slopes, field parameters, time series parameters.
// Measure size of objects
// Likelihood from parameters and priors
// Calculate fields
// Likelihood from data: point data, polygon data.




#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

using namespace R_inla;
using namespace density;
using namespace Eigen;

// ------------------------------------------------------------------------ //
// Spatial field data
// ------------------------------------------------------------------------ //

// The A matrices are for projecting the mesh to a point for the pixel and point data respectively.
// spde is the spde model itself.
DATA_SPARSE_MATRIX(Apixel);
DATA_SPARSE_MATRIX(Apoint);
DATA_STRUCT(spde, spde_t);

// Indicate which datasets to use.
// If neither, just priors. So don't do that.
DATA_IVECTOR(use_points);
DATA_IVECTOR(use_polygons);

// ------------------------------------------------------------------------ //
// Polygon level data
// ------------------------------------------------------------------------ //

// Pixel data. 
// All in long format so that if a pixel is in multiple polygons or multiple years, it will be represented by multiple rows.
// Environmental/other covariates data matrix
DATA_MATRIX(x);

// Population fraction of a polygon in a pixel. i.e. pixel i makes contains 0.01 of the population of polygon j
DATA_VECTOR(xpop);

// two col matrix with start end indices for each shape case. 
DATA_IARRAY(startendindex); 

// Shape data. Cases and region id.
DATA_VECTOR(polygon_cases);



// ------------------------------------------------------------------------ //
// Point level data
// ------------------------------------------------------------------------ //

// Point data. Covariate matrix, number of cases and denominators.
DATA_MATRIX(pointx);
DATA_VECTOR(pointcases);
DATA_VECTOR(pointtested);



// ------------------------------------------------------------------------ //
// Parameters
// ------------------------------------------------------------------------ //

// regression slopes
// (log of) empirical mean incidence to guide intercept
PARAMETER(intercept); // intercept
PARAMETER_VECTOR(slope); 



DATA_SCALAR(priormean_intercept); // = -4.0; 
DATA_SCALAR(priorsd_intercept);// = 2.0
DATA_SCALAR(priormean_slope); // = 0.0;
DATA_SCALAR(priorsd_slope); // = 1.0;

// spde hyperparameters
// tau defines strength of random field. 
// kappa defines distance within which points in field affect each other. 
PARAMETER(log_tau);
PARAMETER(log_kappa);

// Priors on spde hyperparameters
//   kappa -- i.e. exp(priormean_log_kappa) -- set as approximately the width of the region being studied.
//   This implies prior belief in a fairly flat field.
//   tau -- exp(priormean_log_tau) -- set to close to zero. Betas on regression coefficients have priors of 0 so this is reasonable.
//Type priormean_log_kappa = -3;
//Type priorsd_log_kappa   = 0.5;
//Type priormean_log_tau   = -0.50;
//Type priorsd_log_tau     = 2.0;
DATA_SCALAR(priormean_log_kappa);
DATA_SCALAR(priorsd_log_kappa);
DATA_SCALAR(priormean_log_tau);
DATA_SCALAR(priorsd_log_tau);

// Convert hyperparameters to natural scale
Type tau = exp(log_tau);
Type kappa = exp(log_kappa);


// Space-time random effect parameters
// matrix logit_pr_offset [nrows = n_mesh, col=n_years].
PARAMETER_VECTOR(nodemean);

// Prevalence to incidence conversion parameters 
DATA_VECTOR(prev_inc_par); // length: 3 


// get number of data points to loop over
// y (cases) length
int n = polygon_cases.size(); 
// Number of pixels
int pixn = x.rows();
// Number of point surveys
int pointn = pointcases.size();


// ------------------------------------------------------------------------ //
// Likelihood from priors
// ------------------------------------------------------------------------ //



// Initialise negative log likelihood then calc nll of the gaussian markov random field and AR time series
Type nll = 0.0;


// Likelihood of slope parameters given priors
nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
for(int s = 0; s < slope.size(); s++){
  nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
}



// Likelihood of hyperparameters for field
nll -= dnorm(log_kappa, priormean_log_kappa, priorsd_log_kappa, true);
nll -= dnorm(log_tau, priormean_log_tau, priorsd_log_tau, true);


// Build spde matrix
SparseMatrix<Type> Q = Q_spde(spde, kappa);

// Likelihood of the random field.
nll += SCALE(GMRF(Q), 1.0/tau)(nodemean);

Type nll1 = nll;







// ------------------------------------------------------------------------ //
// Likelihood from data
// ------------------------------------------------------------------------ //


if(use_points[1] == 1){
  // Calculate field for point data
  vector<Type> logit_prevalence_point_field;
  logit_prevalence_point_field = Apoint * nodemean;
  
  // Point data likelihood
  vector<Type> point_linear_pred(pointn);
  point_linear_pred = intercept + pointx*slope +
    logit_prevalence_point_field.array();
  point_linear_pred = invlogit(point_linear_pred);
  
  // Hopefully vectorised dbinom.
  nll -= sum(dnbinom(pointcases, pointtested, point_linear_pred, true));
}






// Polygon level likelihood
// For each i in n = cases.size()
//   Find the pixels with correct OBJECTID
//   Sum rate from each pixel
//   Calculate likelihood of binomial.


if(use_polygons[1] == 1){
  
  
  
  // Calculate field for pixel data
  vector<Type> logit_prevalence_pixel_field(pixn);
  logit_prevalence_pixel_field = Apixel * nodemean;
  
  vector<Type> pixel_linear_pred(pixn);
  pixel_linear_pred = intercept + x*slope +
    logit_prevalence_pixel_field.array();
  
  REPORT(pixel_linear_pred);
  
  // recalculate startendindices to be in the form start, n
  startendindex.col(1) = startendindex.col(1) - startendindex.col(0) + 1;
  
  
  vector<Type> inshape_prev;
  vector<Type> inshape_incidencerate;
  vector<Type> inshape_incidence;
  vector<Type> inshape_pop;
  Type shapeincidence = 0.0;
  Type shapepop = 0.0;
  
  vector<Type> reportinc(n);
  vector<Type> reportnll(n);
  vector<Type> reportinccount(n);
  vector<Type> reportpop(n);
  //vector<Type> reportpopprev(n);
  vector<Type> reportprev(n);
  vector<Type> min_inshape_prev(n);
  
  //For each shape use startendindex to find sum of pixel incidence rates
  for (int s = 0; s < n; s++) {
    // Sum pixel risks (raster + field 
  
    // Create logit prevalence
    inshape_prev = pixel_linear_pred.segment(startendindex(s, 0), startendindex(s, 1)).array();
    inshape_prev = invlogit(inshape_prev);
    reportprev[s] = sum(inshape_prev) / inshape_prev.size();
    min_inshape_prev[s] = min(inshape_prev);
      
      
    // Push through ewans prevalence to incidence rate model
    inshape_incidencerate = inshape_prev * prev_inc_par[0] +
                            inshape_prev.pow(2) * prev_inc_par[1] +
                            inshape_prev.pow(3) * prev_inc_par[2];
    
    // Calculate pixel incidence and then polyogn incidence
    inshape_incidence = (inshape_incidencerate * xpop.segment(startendindex(s, 0), startendindex(s, 1)).array());
    shapeincidence = sum(inshape_incidence);
  
    // extract pixel pop and then sum to polygon population
    inshape_pop = xpop.segment(startendindex(s, 0), startendindex(s, 1));
    shapepop = sum(inshape_pop);
  
    reportinccount[s] = shapeincidence;
    reportinc[s] = 1000 * shapeincidence / shapepop;
    reportpop[s] = shapepop;
  
    
   
      nll -= dpois(polygon_cases[s], shapeincidence, true); 
      reportnll[s] = -dpois(polygon_cases[s], shapeincidence, true); 
  }
  REPORT(reportinc);
  REPORT(reportinccount);
  REPORT(reportpop);
  REPORT(reportprev);
  REPORT(polygon_cases);
  REPORT(min_inshape_prev);
  REPORT(reportnll);
}

Type nll3 = nll;
REPORT(nll1);
REPORT(nll3);


return nll;
}
