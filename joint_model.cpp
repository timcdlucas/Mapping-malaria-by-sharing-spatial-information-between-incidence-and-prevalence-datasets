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
DATA_SCALAR(use_points);
DATA_SCALAR(use_polygons);

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

// Pixel i (containing PR data) is in polygon j
// Length of PR data and same order as pointcases
// Values correspond to elements of polygon data in the same order as polygon_cases
DATA_IVECTOR(pointtopolygonmap);


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



// polygon iid effect
PARAMETER_VECTOR(iideffect);
PARAMETER(iideffect_log_tau);
Type iideffect_tau = exp(iideffect_log_tau);
Type iideffect_sd = 1 / sqrt(iideffect_tau);

Type iideffect_mean = 0.0;

// Priors on iid random effect for polygons
DATA_SCALAR(prior_iideffect_sd_max);
DATA_SCALAR(prior_iideffect_sd_prob);

// point iid effect
PARAMETER_VECTOR(iideffect_pr);
PARAMETER(iideffect_pr_log_tau);
Type iideffect_pr_tau = exp(iideffect_pr_log_tau);
Type iideffect_pr_sd = 1 / sqrt(iideffect_pr_tau);

Type iideffect_pr_mean = 0.0;

// Priors on iid random effect for points
DATA_SCALAR(prior_iideffect_pr_sd_max);
DATA_SCALAR(prior_iideffect_pr_sd_prob);


// spde hyperparameters
// sigma defines strength of random field. 
// rho defines range (corr = 0.1 after rho space). 
//   Might need to switch to log.
PARAMETER(log_sigma);
PARAMETER(log_rho);
Type sigma = exp(log_sigma);
Type rho = exp(log_rho);


// Priors on spde hyperparameters
//   kappa -- i.e. exp(priormean_log_kappa) -- set as approximately the width of the region being studied.
//   This implies prior belief in a fairly flat field.
//   tau -- exp(priormean_log_tau) -- set to close to zero. Betas on regression coefficients have priors of 0 so this is reasonable.
//Type priormean_log_kappa = -3;
//Type priorsd_log_kappa   = 0.5;
//Type priormean_log_tau   = -0.50;
//Type priorsd_log_tau     = 2.0;
DATA_SCALAR(prior_rho_min);
DATA_SCALAR(prior_rho_prob);
DATA_SCALAR(prior_sigma_max);
DATA_SCALAR(prior_sigma_prob);

// Convert hyperparameters to natural scale
// todo
Type kappa = sqrt(8) / rho;
Type nu = 1;
// nu = 1
// gamma(nu + d / 2) = gamma(2) = 1
// gamma(nu) = 1
// So gamma(nu + d / 2)(4pi)^{d/2) / gamma(nu) is just 4 * pi
Type tau = sigma * pow(kappa, nu) * sqrt(4 * M_PI);


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


// Likelihood of hyperparameter of polygon iid random effect.
//nll -= dgamma(iideffect_sd, prior_iideffect_sd_shape, prior_iideffect_sd_scale, true);
Type lambda = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
Type pcdensityiid = lambda / 2 * pow(iideffect_tau, -3/2) * exp( - lambda * pow(iideffect_tau, -1/2));
nll -= log(pcdensityiid);

// Likelihood of random effect for polygons
for(int p = 0; p < iideffect.size(); p++) {
  nll -= dnorm(iideffect[p], iideffect_mean, iideffect_sd, true);
}

// Likelihood of hyperparameter of point iid random effect.
Type lambda_pr = -log(prior_iideffect_pr_sd_prob) / prior_iideffect_pr_sd_max;
Type pcdensityiid_pr = lambda_pr / 2 * pow(iideffect_pr_tau, -3/2) * exp( - lambda_pr * pow(iideffect_pr_tau, -1/2));
nll -= log(pcdensityiid_pr);

// Likelihood of random effect for points
for(int p = 0; p < iideffect_pr.size(); p++) {
  nll -= dnorm(iideffect_pr[p], iideffect_pr_mean, iideffect_pr_sd, true);
}

// Likelihood of hyperparameters for field

Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;

Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;

Type pcdensity = lambdatilde1 * lambdatilde2 * pow(rho, -2) * exp(-lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma);

nll -= log(pcdensity);


// Build spde matrix
SparseMatrix<Type> Q = Q_spde(spde, kappa);

// Likelihood of the random field.
nll += SCALE(GMRF(Q), 1.0/tau)(nodemean);

Type nll1 = nll;


Type polygon_weight = 1.0;
Type point_weight = 1.0;

if(use_polygons != 1){
  polygon_weight = 0;
}
if(use_points != 1){
  point_weight = 0;
}



// ------------------------------------------------------------------------ //
// Likelihood from data
// ------------------------------------------------------------------------ //


// Calculate field for point data
vector<Type> logit_prevalence_point_field;
logit_prevalence_point_field = Apoint * nodemean;

// Point data likelihood
vector<Type> point_linear_pred(pointn);
point_linear_pred = intercept + pointx*slope +
  logit_prevalence_point_field.array();

// Hopefully vectorised dbinom.
vector<Type> reportnllpoint(pointn);

for(int q = 0; q < pointn; q++){
  point_linear_pred[q] = point_linear_pred[q] + iideffect[pointtopolygonmap[q] - 1] + iideffect_pr[q]; // Check index of array
  point_linear_pred[q] = invlogit(point_linear_pred[q]);
  nll -= point_weight * dbinom(pointcases[q], pointtested[q], point_linear_pred[q], true);
  reportnllpoint[q] = -point_weight * dbinom(pointcases[q], pointtested[q], point_linear_pred[q], true);
}
REPORT(reportnllpoint);
Type nll2 = nll;






// Polygon level likelihood
// For each i in n = cases.size()
//   Find the pixels with correct OBJECTID
//   Sum rate from each pixel
//   Calculate likelihood of binomial.



  

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
  inshape_prev = pixel_linear_pred.segment(startendindex(s, 0), startendindex(s, 1)).array() + iideffect[s]; // Check about elementwise addition
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

  
 
  nll -= polygon_weight * dpois(polygon_cases[s], shapeincidence, true); 
  reportnll[s] = - polygon_weight * dpois(polygon_cases[s], shapeincidence, true); 
}
REPORT(reportinc);
REPORT(reportinccount);
REPORT(reportpop);
REPORT(reportprev);
REPORT(polygon_cases);
REPORT(min_inshape_prev);
REPORT(reportnll);


Type nll3 = nll;
REPORT(nll1);
REPORT(nll2);
REPORT(nll3);


return nll;
}
