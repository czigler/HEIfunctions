#include <iostream>
#include <string>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#include "covmodel.h"

extern"C" {

SEXP mvCov(SEXP coordsD_r, SEXP n_r, SEXP m_r, SEXP Psi_r, SEXP V_r, SEXP theta_r, SEXP covModel_r){
  
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int mm = m*m;
  int nn = n*n;
  int nm = n*m;
  int nmnm = nm*nm;

  int i, j, k, l, info;
 
  //BLAS and LAPACK vars
  char const *lower = "L";
  char const *upper = "U";
  char const *ntran = "N";
  char const *ytran = "T";
  char const *rside = "R";
  char const *lside = "L";
  double one = 1.0;
  double negOne = -1.0;
  double zero = 0.0;
  int incOne = 1;

  double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
  double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
  double *A = (double *) R_alloc(mm, sizeof(double));
  
  SEXP C_r;
  PROTECT(C_r = allocMatrix(REALSXP, nm, nm)); 
  
  int nPramPtr = 1;
  string covModel = CHAR(STRING_ELT(covModel_r,0));
  
  void (covmodel::*cov1ParamPtr)(double, double &, double &) = NULL; 
  void (covmodel::*cov2ParamPtr)(double, double, double &, double&) = NULL;
  
  if(covModel == "exponential"){
    cov1ParamPtr = &covmodel::exponential;
  }else if(covModel == "spherical"){
    cov1ParamPtr = &covmodel::spherical;
  }else if(covModel == "gaussian"){
    cov1ParamPtr = &covmodel::gaussian;
  }else if(covModel == "matern"){
    cov2ParamPtr = &covmodel::matern;
    nPramPtr = 2;
  }else{
    error("c++ error: cov.model is not correctly specified");
  }
  
  //my covmodel object for calling cov function
  covmodel *covModelObj = new covmodel;
  
  //Get A
  F77_NAME(dcopy)(&mm, REAL(V_r), &incOne, A, &incOne);
  F77_NAME(dpotrf)(lower, &m, A, &m, &info); 
  if(info != 0){error("mvCov: Cholesky failed (1)\n");}

  //clear upper tri
  for(i = 0; i < m-1; i++){
    for(j = i+1; j < m; j++){
      A[j*m+i] = 0.0;
    }
  }
  
  //make C
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      
      zeros(tmp_mm, mm);
      
      for(k = 0; k < m; k++){
	if(nPramPtr == 1)
	  (covModelObj->*cov1ParamPtr)(REAL(theta_r)[k], tmp_mm[k*m+k], REAL(coordsD_r)[j*n+i]);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(REAL(theta_r)[k], REAL(theta_r)[m+k], tmp_mm[k*m+k], REAL(coordsD_r)[j*n+i]);
      }
      
      F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
      
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  REAL(C_r)[((j*m+l)*nm)+(i*m+k)] = tmp_mm[l*m+k];
	  tmp_mm[l*m+k] = 0.0; //zero out
	}
      }
    }
  }
     
  for(i = 0; i < n; i++){
    for(k = 0; k < m; k++){
      for(l = 0; l < m; l++){
	REAL(C_r)[(i*m+l)*nm+(i*m+k)] += REAL(Psi_r)[l*m+k];
      }
    }
  }
  
  UNPROTECT(1);

  return(C_r);
}

}
