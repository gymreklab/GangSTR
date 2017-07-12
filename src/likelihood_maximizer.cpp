/*
Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>
and Nima Mousavi (mousavi@ucsd.edu)

This file is part of GangSTR.

GangSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GangSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GangSTR.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <nlopt.hpp>
// #include <nlopt.h>

#include <gsl/gsl_multimin.h>
#include "src/likelihood_maximizer.h"
#include "src/mathops.h"
#include <iostream>

using namespace std;
//////////
int count = 0;
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
  ++count;
    if (grad) {
        cerr<< "No grad!"<<endl;
        return 0.0;
    }
    // return 2 * sqrt(x[1]);
    return pow(x[0], 2) - 0.5 * pow(pow(x[0] - 1, 2) - 3, 2);
}
typedef struct {
    double a, b;
} my_constraint_data;
double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }
////
struct nlopt_data{
  LikelihoodMaximizer* lm_ptr;
  int read_len, motif_len, ref_count;
  nlopt_data(int Read_Len, int Motif_Len, int Ref_Count, LikelihoodMaximizer* LM_OBJ) : 
    read_len(Read_Len), motif_len(Motif_Len), ref_count(Ref_Count), lm_ptr(LM_OBJ) {
      // cout<<endl<<"constr LM\t"<<LM_OBJ->GetEnclosingDataSize();
      // cout<<endl<<"constr lm\t"<<lm_ptr->GetEnclosingDataSize();
    }
};

void dummy_test(LikelihoodMaximizer* LM_OBJ){
    cout<<endl<<"Dummyyy\t"<<LM_OBJ->GetEnclosingDataSize()<<endl;
}

double nloptNegLikelihood(unsigned n, const double *x, double *grad, void *data)
{
  ++count;
  if (grad) {
        cerr<< "No grad!"<<endl;
        return 0.0;
  }
  nlopt_data *d = (nlopt_data *) data;
  int read_len  = d -> read_len;
  int motif_len = d -> motif_len;
  int ref_count = d -> ref_count; 
  LikelihoodMaximizer* lm_ptr = d -> lm_ptr;
  // cout<<lm_ptr->GetSpanningDataSize();
  double A = x[0], B = x[1];
  double gt_ll;
  if(!lm_ptr->GetGenotypeNegLogLikelihood(A, B, read_len, motif_len, ref_count, &gt_ll))
    return -1.0;
  else{
    return gt_ll;
  }
}
////

double gslNegLikelihood(const gsl_vector *v, void *params)
{
  double A, B;
  LikelihoodMaximizer *p = (LikelihoodMaximizer *)params;
  double gt_ll;
  A = gsl_vector_get(v, 0);
  B = gsl_vector_get(v, 1);
 
  if(!p[0].GetGenotypeNegLogLikelihood(A, B, p[0].options->read_len, p[0].options->motif_len, p[0].options->ref_count, &gt_ll))
    return -1.0;
  else
    return gt_ll;
}
////////////

LikelihoodMaximizer::LikelihoodMaximizer(const Options& _options) {
  options = &_options;
  enclosing_class_.SetOptions(*options);
  frr_class_.SetOptions(*options);
  spanning_class_.SetOptions(*options);
}

// // Copy constructor:
// LikelihoodMaximizer::LikelihoodMaximizer(const LikelihoodMaximizer& lm_obj){
//   options = lm_obj.options;
//   enclosing_class_.SetOptions(*options);
//   frr_class_.SetOptions(*options);
//   spanning_class_.SetOptions(*options);
// }

void LikelihoodMaximizer::Reset() {
  enclosing_class_.Reset();
  frr_class_.Reset();
  spanning_class_.Reset();
}

void LikelihoodMaximizer::AddEnclosingData(const int32_t& data) {
  enclosing_class_.AddData(data);
}
void LikelihoodMaximizer::AddSpanningData(const int32_t& data) {
  spanning_class_.AddData(data);
}
void LikelihoodMaximizer::AddFRRData(const int32_t& data) {
  frr_class_.AddData(data);
}
std::size_t LikelihoodMaximizer::GetEnclosingDataSize() {
  return enclosing_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetSpanningDataSize() {
  return spanning_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetFRRDataSize() {
  return frr_class_.GetDataSize();
}

bool LikelihoodMaximizer::GetGenotypeNegLogLikelihood(const int32_t& allele1,
						      const int32_t& allele2,
						      const int32_t& read_len,
						      const int32_t& motif_len,
						      const int32_t& ref_count,
						      double* gt_ll) {
  double frr_ll, span_ll, encl_ll = 0.0;
  frr_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, &frr_ll);
  spanning_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, &span_ll);
  enclosing_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, &encl_ll);
  // cout<<endl<<frr_ll<<"\t"<<span_ll<<"\t"<<encl_ll;
  *gt_ll = -1*(options->frr_weight*frr_ll +
	       options->spanning_weight*span_ll +
	       options->enclosing_weight*encl_ll);
}



bool LikelihoodMaximizer::OptimizeLikelihood(const int32_t& read_len, const int32_t& motif_len,
					     const int32_t& ref_count,
					     int32_t* allele1, int32_t* allele2) {
  // //////////// NLOPT C
  // double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
  // nlopt_opt opt;
  // opt = nlopt_create(NLOPT_LN_COBYLA, 2); /* algorithm and dimensionality */
  // nlopt_set_lower_bounds(opt, lb);
  // nlopt_set_min_objective(opt, myfunc, NULL);

  // my_constraint_data data[2] = { {2,0}, {-1,1} };
  // nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
  // nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);

  // nlopt_set_xtol_rel(opt, 1e-4);

  // double x[2] = { 1.234, 5.678 };  /* some initial guess */
  // double minf; /* the minimum objective value, upon return */

  // if (nlopt_optimize(opt, x, &minf) < 0) {
  //     printf("nlopt failed!\n");
  // }
  // else {
  //  printf("found minimum after %d evaluations\n", count);
  //     printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
  // }

  // nlopt_destroy(opt);

  // //////////// NLOPT C

  // //////////// NLOPT C++ n = 1

  nlopt::opt opt(nlopt::LN_COBYLA, 2);

  std::vector<double> lb(2);
  lb[0] = 33;
  lb[1] = 33;
  opt.set_lower_bounds(lb);

  std::vector<double> ub(2);
  ub[0] = 200;
  ub[1] = 200;
  opt.set_upper_bounds(ub);

  nlopt_data data[1] = nlopt_data(read_len, motif_len, ref_count, this);
  

  opt.set_min_objective(nloptNegLikelihood, data);    // Change to max for maximization

  opt.set_xtol_rel(1e-4);

  std::vector<double> xx(2);
  xx[0] = 40;
  xx[1] = 50;
  double minf;
  nlopt::result result = opt.optimize(xx, minf);

  cout<<result<<"\t"<<xx[0]<<","<<xx[1]<<"\t"<<minf<<endl;
  
  // dummy_test(this);
  // cout<<data[0].lm_ptr->GetSpanningDataSize()<<endl;
  // cout<<data[0].read_len<<endl;
  const double yy[2] = {50, 70};
  cout<<endl<<"(50, 70): "<<nloptNegLikelihood(2, yy, NULL, data)<<endl;
  double ans = 0.111;
  this->frr_class_.GetClassLogLikelihood(50, 70, 100, 3, 13, &ans);
  cout<<endl<<this->frr_class_.GetDataSize()<<endl;
  cout<<endl<<"(50, 70) No wrapper: "<<ans<<endl;
  // const double uu[2] = {52, 78};
  // cout<<endl<<"(52, 78): "<<nloptNegLikelihood(2, uu, NULL, data)<<endl;
  // const double ii[2] = {40, 80};
  // cout<<endl<<"(40, 80): "<<nloptNegLikelihood(2, ii, NULL, data)<<endl;
  // const double oo[2] = {50, 78};
  // cout<<endl<<"(50, 78): "<<nloptNegLikelihood(2, oo, NULL, data)<<endl;
  // //////////// NLOPT C++ n = 1


  // //////////// NLOPT C++

  // nlopt::opt opt(nlopt::LN_COBYLA, 2);

  // std::vector<double> lb(2);
  // lb[0] = -HUGE_VAL; lb[1] = 0;
  // opt.set_lower_bounds(lb);

  // opt.set_min_objective(myfunc, NULL);    // Change to max for maximization

  // my_constraint_data data[2] = { {2,0}, {-1,1} };
  // opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
  // opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);

  // opt.set_xtol_rel(1e-4);

  // std::vector<double> xx(2);
  // xx[0] = 1.234; xx[1] = 5.678;
  // double minf;
  // nlopt::result result = opt.optimize(xx, minf);

  // cout<<result<<"\t"<<minf<<endl;
  // cout<<"Hi\n";
  // //////////// NLOPT C++

  // //////////// GSL Optimizer
  // size_t minim_dim = 2;
  // const gsl_multimin_fminimizer_type * minim_type = gsl_multimin_fminimizer_nmsimplex2;
  // gsl_multimin_fminimizer * minim_handle = gsl_multimin_fminimizer_alloc (minim_type, minim_dim);

  // gsl_vector *ss, *x;
  // gsl_multimin_function neg_log_like;

  // size_t iter = 0;
  // int status;
  // double size;
  // //LikelihoodMaximizer* par = new LikelihoodMaximizer(options); // Params = [read_len, motif_len, ref_count]
  // /* Starting point */
  // x = gsl_vector_alloc (2);
  // gsl_vector_set (x, 0, 10.0);
  // gsl_vector_set (x, 1, 30.0);

  // /* Set initial step sizes to 1 */
  // ss = gsl_vector_alloc (2);
  // gsl_vector_set_all (ss, 4.0);

  // /* Initialize method and iterate */
  // neg_log_like.n = 2;
  // neg_log_like.f = &(dummy_func);
  // neg_log_like.params = this;

  // minim_handle = gsl_multimin_fminimizer_alloc (minim_type, 2);
  // gsl_multimin_fminimizer_set (minim_handle, &neg_log_like, x, ss);

  // do
  //   {
  //     iter++;
  //     status = gsl_multimin_fminimizer_iterate(minim_handle);
      
  //     if (status) 
  //       break;

  //     size = gsl_multimin_fminimizer_size (minim_handle);
  //     status = gsl_multimin_test_size (size, 1e-2);

  //     if (status == GSL_SUCCESS)
  //       {
  //         printf ("converged to minimum at\n");
  //       }

  //     printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
  //             iter,
  //             gsl_vector_get (minim_handle->x, 0), 
  //             gsl_vector_get (minim_handle->x, 1), 
  //             minim_handle->fval, size);
  //   }
  // while (status == GSL_CONTINUE && iter < 100);
  
  // gsl_vector_free(x);
  // gsl_vector_free(ss);
  // gsl_multimin_fminimizer_free (minim_handle);

  // cout<<gsl_vector_get (minim_handle->x, 0)<<"\t"<<gsl_vector_get (minim_handle->x, 1);
  // *allele1 = 2;
  // //////////// GSL Optimizer
  return false;
}
LikelihoodMaximizer::~LikelihoodMaximizer() {}
