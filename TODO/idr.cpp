// **************************************************************************
//
//    PARALUTION   www.paralution.com
//
//    Copyright (C) 2015  PARALUTION Labs UG (haftungsbeschränkt) & Co. KG
//                        Am Hasensprung 6, 76571 Gaggenau
//                        Handelsregister: Amtsgericht Mannheim, HRA 706051
//                        Vertreten durch:
//                        PARALUTION Labs Verwaltungs UG (haftungsbeschränkt)
//                        Am Hasensprung 6, 76571 Gaggenau
//                        Handelsregister: Amtsgericht Mannheim, HRB 721277
//                        Geschäftsführer: Dimitar Lukarski, Nico Trost
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// **************************************************************************



// PARALUTION version 1.1.0 


#include "../../utils/def.hpp"
#include "idr.hpp"
#include "../iter_ctrl.hpp"

#include "../../base/local_matrix.hpp"
#include "../../base/matrix_formats_ind.hpp"
#include "../../base/local_stencil.hpp"
#include "../../base/local_vector.hpp"

#include "../../utils/log.hpp"
#include "../../utils/allocate_free.hpp"
#include "../../utils/math_functions.hpp"

#include <math.h>
#include <time.h>

namespace paralution {

template <class OperatorType, class VectorType, typename ValueType>
IDR<OperatorType, VectorType, ValueType>::IDR() {

  LOG_DEBUG(this, "IDR::IDR()",
            "default constructor");

  this->s_ = 4;
  this->kappa_ = ValueType(0.7f);

  this->fhost_ = NULL;
  this->Mhost_ = NULL;

  this->g_ = NULL;
  this->u_ = NULL;
  this->P_ = NULL;

}

template <class OperatorType, class VectorType, typename ValueType>
IDR<OperatorType, VectorType, ValueType>::~IDR() {

  LOG_DEBUG(this, "IDR::~IDR()",
            "destructor");

  this->Clear();

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::Print(void) const {

  if (this->precond_ == NULL) { 

    LOG_INFO("IDR(" << this->s_ << ") solver");

  } else {

    LOG_INFO("IDR(" << this->s_ << ") solver, with preconditioner:");
    this->precond_->Print();

  }
  
}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::PrintStart_(void) const {

  if (this->precond_ == NULL) { 

    LOG_INFO("IDR(" << this->s_ << ") (non-precond) linear solver starts");

  } else {

    LOG_INFO("PIDR(" << this->s_ << ") solver starts, with preconditioner:");
    this->precond_->Print();

  }

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::PrintEnd_(void) const {

  if (this->precond_ == NULL) { 

    LOG_INFO("IDR(" << this->s_ << ") (non-precond) ends");

  } else {

    LOG_INFO("PIDR(" << this->s_ << ") ends");

  }

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::Build(void) {

  LOG_DEBUG(this, "IDR::Build()",
            this->build_ <<
            " #*# begin");

  if (this->build_ == true)
    this->Clear();

  assert(this->build_ == false);
  this->build_ = true;

  assert(this->op_ != NULL);
  assert(this->op_->get_nrow() == this->op_->get_ncol());
  assert(this->op_->get_nrow() > 0);

  if (this->s_ > this->op_->get_nrow())
    this->s_ = this->op_->get_nrow();

  if (this->precond_ != NULL) {

    this->precond_->SetOperator(*this->op_);

    this->precond_->Build();

    this->z_.CloneBackend(*this->op_);
    this->z_.Allocate("z", this->op_->get_nrow());

  }

  allocate_host(this->s_, &this->fhost_);
  allocate_host(this->s_*this->s_, &this->Mhost_);

  this->r_.CloneBackend(*this->op_);
  this->r_.Allocate("r", this->op_->get_nrow());

  this->v_.CloneBackend(*this->op_);
  this->v_.Allocate("v", this->op_->get_nrow());

  this->t_.CloneBackend(*this->op_);
  this->t_.Allocate("t", this->op_->get_nrow());

  this->g_ = new VectorType*[this->s_];
  this->u_ = new VectorType*[this->s_];
  this->P_ = new VectorType*[this->s_];

  for (int i=0; i<this->s_; ++i) {

    this->g_[i] = new VectorType;
    this->g_[i]->CloneBackend(*this->op_);
    this->g_[i]->Allocate("g", this->op_->get_nrow());

    this->u_[i] = new VectorType;
    this->u_[i]->CloneBackend(*this->op_);
    this->u_[i]->Allocate("u", this->op_->get_nrow());

    this->P_[i] = new VectorType;
    this->P_[i]->CloneBackend(*this->op_);
    this->P_[i]->Allocate("P", this->op_->get_nrow());
    this->P_[i]->SetRandom(-2.0, 2.0, int((i+1)*time(NULL)));

  }

  // Build ONB out of P using modified Gram-Schmidt algorithm
  for (int k=0; k<this->s_; ++k) {
    for (int j=0; j<k; ++j)
      this->P_[k]->AddScale(*this->P_[j] , ValueType(-1.0)*this->P_[j]->Dot(*this->P_[k]));
    this->P_[k]->Scale(ValueType(1.0) / this->P_[k]->Norm());
  }

  LOG_DEBUG(this, "IDR::Build()",
            this->build_ <<
            " #*# end");

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::Clear(void) {

  LOG_DEBUG(this, "IDR::Clear()",
            this->build_);

  if (this->build_ == true) {

    this->r_ .Clear();
    this->v_ .Clear();
    this->t_ .Clear();

    for (int i=0; i<this->s_; ++i) {

      delete this->u_[i];
      delete this->g_[i];
      delete this->P_[i];

    }

    delete [] this->u_;
    delete [] this->g_;
    delete [] this->P_;

    free_host(&this->fhost_);
    free_host(&this->Mhost_);

    if (this->precond_ != NULL) {

      this->precond_->Clear();
      this->precond_   = NULL;

      this->z_.Clear();

    }

    this->iter_ctrl_.Clear();

    this->build_ = false;

  }

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::SetShadowSpace(const int s) {

  LOG_DEBUG(this, "IDR:SetShadowSpace()",
            s);

  assert(this->build_ == false);
  assert(s > 0);
  assert(this->op_ != NULL);
  assert(s <= this->op_->get_nrow());

  this->s_ = s;

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::MoveToHostLocalData_(void) {

  LOG_DEBUG(this, "IDR::MoveToHostLocalData_()",
            this->build_);  

  if (this->build_ == true) {

    this->r_.MoveToHost();
    this->v_.MoveToHost();
    this->t_.MoveToHost();

    for (int i=0; i<this->s_; ++i) {

      this->u_[i]->MoveToHost();
      this->g_[i]->MoveToHost();
      this->P_[i]->MoveToHost();

    }

    if (this->precond_ != NULL) {

      this->z_.MoveToHost();

    }

  }

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::MoveToAcceleratorLocalData_(void) {

  LOG_DEBUG(this, "IDR::MoveToAcceleratorLocalData_()",
            this->build_);

  if (this->build_ == true) {

    this->r_.MoveToAccelerator();
    this->v_.MoveToAccelerator();
    this->t_.MoveToAccelerator();

    for (int i=0; i<this->s_; ++i) {

      this->u_[i]->MoveToAccelerator();
      this->g_[i]->MoveToAccelerator();
      this->P_[i]->MoveToAccelerator();

    }

    if (this->precond_ != NULL) {

      this->z_.MoveToAccelerator();

    }

  }

}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::SolveNonPrecond_(const VectorType &rhs,
                                                                      VectorType *x) {

  LOG_DEBUG(this, "IDR::SolveNonPrecond_()",
            " #*# begin");

  assert(x != NULL);
  assert(x != &rhs);
  assert(this->op_  != NULL);
  assert(this->precond_  == NULL);
  assert(this->build_ == true);

  const OperatorType *op = this->op_;

  VectorType *r = &this->r_;
  VectorType *v = &this->v_;
  VectorType *t = &this->t_;
  VectorType **g = this->g_;
  VectorType **u = this->u_;
  VectorType **P = this->P_;

  const int s = this->s_;
  const ValueType kappa = this->kappa_;

  ValueType *fhost = this->fhost_;
  ValueType *Mhost = this->Mhost_;

  ValueType alpha;
  ValueType beta;
  ValueType omega = ValueType(1.0);
  ValueType rho;

  // initial residual r = b - Ax
  op->Apply(*x, r);
  r->ScaleAdd(ValueType(-1.0), rhs);

  // use for |b-Ax0|
  ValueType res_norm = this->Norm(*r);

  if (this->iter_ctrl_.InitResidual(paralution_abs(res_norm)) == false) {

      LOG_DEBUG(this, "::SolveNonPrecond_()",
            " #*# end");

      return;
  }

  // g = u = v = 0, M = I
  for (int i=0; i<s*s; ++i)
    Mhost[i] = ValueType(0.0);

  for (int i=0; i<s; ++i) {
    g[i]->Zeros();
    u[i]->Zeros();
    Mhost[DENSE_IND(i,i,s,s)] = ValueType(1.0);
  }

  // f = P^T * r
  for (int i=0; i<s; ++i)
    fhost[i] = P[i]->Dot(*r);

  for (int k=0; k<s; ++k) {

    LocalVector<ValueType> f;
    ValueType *fsub = fhost+k;
    f.SetDataPtr(&fsub, "f", s-k);

    LocalVector<ValueType> c;
    c.Allocate("c", s-k);

    ValueType *subMat = NULL;
    allocate_host((s-k)*(s-k), &subMat);
    for (int i=k; i<s; ++i)
      for (int j=k; j<s; ++j)
        subMat[DENSE_IND((i-k),(j-k),(s-k),(s-k))] = Mhost[DENSE_IND(i,j,s,s)];

    LocalMatrix<ValueType> M;
    M.ConvertToDENSE();
    M.SetDataPtrDENSE(&subMat, "M", s-k, s-k);
    subMat = NULL;

    // Solve Mc = f
    M.QRDecompose();
    M.QRSolve(f, &c);

    f.LeaveDataPtr(&fsub);

    ValueType *chost = NULL;
    c.LeaveDataPtr(&chost);

    v->Zeros();

    // v = r - sum(c_i g_i)
    for (int i=k; i<s; ++i)
      v->AddScale(*g[i], chost[i-k]);

    v->ScaleAdd(ValueType(-1.0), *r);

    // u_k = sum(c_i u_i)
    for (int i=k; i<s; ++i) {
      if (k == i)
        u[k]->Scale(chost[i-k]);
      else
        u[k]->AddScale(*u[i], chost[i-k]);
    }

    free_host(&chost);

    // u_k = u_k + v
    u[k]->ScaleAdd(ValueType(1.0), *v);

    // g_k = A u_k
    op->Apply(*u[k], g[k]);

    for (int i=0; i<k; ++i) {

      // alpha = P^T_i * g_k / M_ii
      alpha = P[i]->Dot(*g[k]) / Mhost[DENSE_IND(i,i,s,s)];

      // g_k = g_k - alpha * g_i
      g[k]->AddScale(*g[i], ValueType(-1.0)*alpha);

      // u_k = u_k - alpha * u_i
      u[k]->AddScale(*u[i], ValueType(-1.0)*alpha);

    }

    for (int i=k; i<s; ++i) {

      // M_i_k = P^T_i * g_k
      Mhost[DENSE_IND(i,k,s,s)] = P[i]->Dot(*g[k]);

      if (Mhost[DENSE_IND(k,k,s,s)] == ValueType(0.0)) {
        LOG_INFO("IDR(s) break down ; M(k,k) == 0.0");
        FATAL_ERROR(__FILE__, __LINE__);
      }

    }

    // beta = f_k / M_k_k
    beta = fhost[k] / Mhost[DENSE_IND(k,k,s,s)];

    // r = r - beta * g_k
    r->AddScale(*g[k], ValueType(-1.0)*beta);

    // x = x + beta * u_k
    x->AddScale(*u[k], beta);

    // Residual norm
    res_norm = this->Norm(*r);

    // Check inner loop for convergence
    if (this->iter_ctrl_.CheckResidual(paralution_abs(res_norm), this->index_))
      break;

    // f_i = f_i - beta * M_i_k
    for (int i=k+1; i<s; ++i)
      fhost[i] -= beta * Mhost[DENSE_IND(i,k,s,s)];

  }

  // t = Ar
  op->Apply(*r, t);

  // omega = (t,r) / ||t||^2
  ValueType rt = t->Dot(*r);
  ValueType nt = t->Norm();

  rt /= nt;

  rho = paralution_abs(rt / res_norm);
  omega = rt / nt;

  if (rho < kappa)
    omega = omega*kappa/rho;

  if (omega == ValueType(0.0)) {
    LOG_INFO("IDR(s) break down ; omega == 0.0");
    FATAL_ERROR(__FILE__, __LINE__);
  }

  // x = x + omega * r
  x->AddScale(*r, omega);

  // r = r - omega * t
  r->AddScale(*t, ValueType(-1.0)*omega);


  // Residual norm to check outer loop convergence
  res_norm = this->Norm(*r);

  while (!this->iter_ctrl_.CheckResidual(paralution_abs(res_norm), this->index_)) {

    // f = P^T * r
    for (int i=0; i<s; ++i)
      fhost[i] = P[i]->Dot(*r);

    for (int k=0; k<s; ++k) {

      LocalVector<ValueType> f;
      ValueType *fsub = fhost+k;
      f.SetDataPtr(&fsub, "f", s-k);

      LocalVector<ValueType> c;
      c.Allocate("c", s-k);

      ValueType *subMat = NULL;
      allocate_host((s-k)*(s-k), &subMat);
      for (int i=k; i<s; ++i)
        for (int j=k; j<s; ++j)
          subMat[DENSE_IND((i-k),(j-k),(s-k),(s-k))] = Mhost[DENSE_IND(i,j,s,s)];

      LocalMatrix<ValueType> M;
      M.ConvertToDENSE();
      M.SetDataPtrDENSE(&subMat, "M", s-k, s-k);
      subMat = NULL;

      // Solve Mc = f
      M.QRDecompose();
      M.QRSolve(f, &c);

      f.LeaveDataPtr(&fsub);

      ValueType *chost = NULL;
      c.LeaveDataPtr(&chost);

      v->Zeros();

      // v = r - sum(c_i g_i)
      for (int i=k; i<s; ++i)
        v->AddScale(*g[i], chost[i-k]);

      v->ScaleAdd(ValueType(-1.0), *r);

      // v = omega * v
      v->Scale(omega);

      // u_k = sum(c_i u_i)
      for (int i=k; i<s; ++i) {
        if (k == i)
          u[k]->Scale(chost[i-k]);
        else
          u[k]->AddScale(*u[i], chost[i-k]);
      }

      free_host(&chost);

      // u_k = u_k + v
      u[k]->ScaleAdd(ValueType(1.0), *v);

      // g_k = A u_k
      op->Apply(*u[k], g[k]);

      for (int i=0; i<k; ++i) {

        // alpha = P^T_i * g_k / M_ii
        alpha = P[i]->Dot(*g[k]) / Mhost[DENSE_IND(i,i,s,s)];

        // g_k = g_k - alpha * g_i
        g[k]->AddScale(*g[i], ValueType(-1.0)*alpha);

        // u_k = u_k - alpha * u_i
        u[k]->AddScale(*u[i], ValueType(-1.0)*alpha);

      }

      for (int i=k; i<s; ++i) {

        // M_i_k = P^T_i * g_k
        Mhost[DENSE_IND(i,k,s,s)] = P[i]->Dot(*g[k]);

        if (Mhost[DENSE_IND(k,k,s,s)] == ValueType(0.0)) {
          LOG_INFO("IDR(s) break down ; M(k,k) == 0.0");
          FATAL_ERROR(__FILE__, __LINE__);
        }

      }

      // beta = f_k / M_k_k
      beta = fhost[k] / Mhost[DENSE_IND(k,k,s,s)];

      // r = r - beta * g_k
      r->AddScale(*g[k], ValueType(-1.0)*beta);

      // x = x + beta * u_k
      x->AddScale(*u[k], beta);

      // Residual norm
      res_norm = this->Norm(*r);

      // Check inner loop for convergence
      if (this->iter_ctrl_.CheckResidual(paralution_abs(res_norm), this->index_))
        break;

      // f_i = f_i - beta * M_i_k
      for (int i=k+1; i<s; ++i)
        fhost[i] -= beta * Mhost[DENSE_IND(i,k,s,s)];

    }

    // t = Ar
    op->Apply(*r, t);

    // omega = (t,r) / ||t||^2
    ValueType rt = t->Dot(*r);
    ValueType nt = t->Norm();

    rt /= nt;

    rho = paralution_abs(rt / res_norm);
    omega = rt / nt;

    if (rho < kappa)
      omega = omega*kappa/rho;

    if (omega == ValueType(0.0)) {
      LOG_INFO("IDR(s) break down ; omega == 0.0");
      FATAL_ERROR(__FILE__, __LINE__);
    }

    // x = x + omega * r
    x->AddScale(*r, omega);

    // r = r - omega * t
    r->AddScale(*t, ValueType(-1.0)*omega);


    // Residual norm to check outer loop convergence
    res_norm = this->Norm(*r);

  }

  LOG_DEBUG(this, "IDR::SolveNonPrecond_()",
            " #*# end");


}

template <class OperatorType, class VectorType, typename ValueType>
void IDR<OperatorType, VectorType, ValueType>::SolvePrecond_(const VectorType &rhs,
                                                                   VectorType *x) {

  LOG_DEBUG(this, "IDR::SolvePrecond_()",
            " #*# begin");

  assert(x != NULL);
  assert(x != &rhs);
  assert(this->op_  != NULL);
  assert(this->precond_ != NULL);
  assert(this->build_ == true);

  const OperatorType *op = this->op_;

  VectorType *r = &this->r_;
  VectorType *v = &this->v_;
  VectorType *z = &this->z_;
  VectorType *t = &this->t_;
  VectorType **g = this->g_;
  VectorType **u = this->u_;
  VectorType **P = this->P_;

  const int s = this->s_;
  const ValueType kappa = this->kappa_;

  ValueType *fhost = this->fhost_;
  ValueType *Mhost = this->Mhost_;

  ValueType alpha;
  ValueType beta;
  ValueType omega = ValueType(1.0);
  ValueType rho;

  // initial residual r = b - Ax
  op->Apply(*x, r);
  r->ScaleAdd(ValueType(-1.0), rhs);

  // use for |b-Ax0|
  ValueType res_norm = this->Norm(*r);

  if (this->iter_ctrl_.InitResidual(paralution_abs(res_norm)) == false) {

      LOG_DEBUG(this, "::SolvePrecond_()",
            " #*# end");

      return;
  }

  // g = u = v = 0, M = I
  for (int i=0; i<s*s; ++i)
    Mhost[i] = ValueType(0.0);

  for (int i=0; i<s; ++i) {
    g[i]->Zeros();
    u[i]->Zeros();
    Mhost[DENSE_IND(i,i,s,s)] = ValueType(1.0);
  }

  // f = P^T * r
  for (int i=0; i<s; ++i)
    fhost[i] = P[i]->Dot(*r);

  for (int k=0; k<s; ++k) {

    LocalVector<ValueType> f;
    ValueType *fsub = fhost+k;
    f.SetDataPtr(&fsub, "f", s-k);

    LocalVector<ValueType> c;
    c.Allocate("c", s-k);

    ValueType *subMat = NULL;
    allocate_host((s-k)*(s-k), &subMat);
    for (int i=k; i<s; ++i)
      for (int j=k; j<s; ++j)
        subMat[DENSE_IND((i-k),(j-k),(s-k),(s-k))] = Mhost[DENSE_IND(i,j,s,s)];

    LocalMatrix<ValueType> M;
    M.ConvertToDENSE();
    M.SetDataPtrDENSE(&subMat, "M", s-k, s-k);
    subMat = NULL;

    // Solve Mc = f
    M.QRDecompose();
    M.QRSolve(f, &c);

    f.LeaveDataPtr(&fsub);

    ValueType *chost = NULL;
    c.LeaveDataPtr(&chost);

    v->Zeros();

    // v = r - sum(c_i g_i)
    for (int i=k; i<s; ++i)
      v->AddScale(*g[i], chost[i-k]);

    v->ScaleAdd(ValueType(-1.0), *r);

    // Preconditioning Mz = v
    this->precond_->SolveZeroSol(*v, z);

    // u_k = sum(c_i u_i)
    for (int i=k; i<s; ++i) {
      if (k == i)
        u[k]->Scale(chost[i-k]);
      else
        u[k]->AddScale(*u[i], chost[i-k]);
    }

    free_host(&chost);

    // u_k = u_k + z
    u[k]->ScaleAdd(ValueType(1.0), *z);

    // g_k = A u_k
    op->Apply(*u[k], g[k]);

    for (int i=0; i<k; ++i) {

      // alpha = P^T_i * g_k / M_ii
      alpha = P[i]->Dot(*g[k]) / Mhost[DENSE_IND(i,i,s,s)];

      // g_k = g_k - alpha * g_i
      g[k]->AddScale(*g[i], ValueType(-1.0)*alpha);

      // u_k = u_k - alpha * u_i
      u[k]->AddScale(*u[i], ValueType(-1.0)*alpha);

    }

    for (int i=k; i<s; ++i) {

      // M_i_k = P^T_i * g_k
      Mhost[DENSE_IND(i,k,s,s)] = P[i]->Dot(*g[k]);

      if (Mhost[DENSE_IND(k,k,s,s)] == ValueType(0.0)) {
        LOG_INFO("IDR(s) break down ; M(k,k) == 0.0");
        FATAL_ERROR(__FILE__, __LINE__);
      }

    }

    // beta = f_k / M_k_k
    beta = fhost[k] / Mhost[DENSE_IND(k,k,s,s)];

    // r = r - beta * g_k
    r->AddScale(*g[k], ValueType(-1.0)*beta);

    // x = x + beta * u_k
    x->AddScale(*u[k], beta);

    // Residual norm
    res_norm = this->Norm(*r);

    // Check inner loop for convergence
    if (this->iter_ctrl_.CheckResidual(paralution_abs(res_norm), this->index_))
      break;

    // f_i = f_i - beta * M_i_k
    for (int i=k+1; i<s; ++i)
      fhost[i] -= beta * Mhost[DENSE_IND(i,k,s,s)];

  }

  // Mv = r
  this->precond_->SolveZeroSol(*r, v);

  // t = Av
  op->Apply(*v, t);

  // omega = (t,r) / ||t||^2
  ValueType rt = t->Dot(*r);
  ValueType nt = t->Norm();

  rt /= nt;

  rho = paralution_abs(rt / res_norm);
  omega = rt / nt;

  if (rho < kappa)
    omega = omega*kappa/rho;

  if (omega == ValueType(0.0)) {
    LOG_INFO("IDR(s) break down ; omega == 0.0");
    FATAL_ERROR(__FILE__, __LINE__);
  }

  // r = r - omega * t
  r->AddScale(*t, ValueType(-1.0)*omega);

  // x = x + omega * v
  x->AddScale(*v, omega);

  // Residual norm to check outer loop convergence
  res_norm = this->Norm(*r);

  while (!this->iter_ctrl_.CheckResidual(paralution_abs(res_norm), this->index_)) {

    // f = P^T * r
    for (int i=0; i<s; ++i)
      fhost[i] = P[i]->Dot(*r);

    for (int k=0; k<s; ++k) {

      LocalVector<ValueType> f;
      ValueType *fsub = fhost+k;
      f.SetDataPtr(&fsub, "f", s-k);

      LocalVector<ValueType> c;
      c.Allocate("c", s-k);

      ValueType *subMat = NULL;
      allocate_host((s-k)*(s-k), &subMat);
      for (int i=k; i<s; ++i)
        for (int j=k; j<s; ++j)
          subMat[DENSE_IND((i-k),(j-k),(s-k),(s-k))] = Mhost[DENSE_IND(i,j,s,s)];

      LocalMatrix<ValueType> M;
      M.ConvertToDENSE();
      M.SetDataPtrDENSE(&subMat, "M", s-k, s-k);
      subMat = NULL;

      // Solve Mc = f
      M.QRDecompose();
      M.QRSolve(f, &c);

      f.LeaveDataPtr(&fsub);

      ValueType *chost = NULL;
      c.LeaveDataPtr(&chost);

      v->Zeros();

      // v = r - sum(c_i g_i)
      for (int i=k; i<s; ++i)
        v->AddScale(*g[i], chost[i-k]);

      v->ScaleAdd(ValueType(-1.0), *r);

      // Preconditioning Mz = v
      this->precond_->SolveZeroSol(*v, z);

      // z = omega * z
      z->Scale(omega);

      // u_k = sum(c_i u_i)
      for (int i=k; i<s; ++i) {
        if (k == i)
          u[k]->Scale(chost[i-k]);
        else
          u[k]->AddScale(*u[i], chost[i-k]);
      }

      free_host(&chost);

      // u_k = u_k + z
      u[k]->ScaleAdd(ValueType(1.0), *z);

      // g_k = A u_k
      op->Apply(*u[k], g[k]);

      for (int i=0; i<k; ++i) {

        // alpha = P^T_i * g_k / M_ii
        alpha = P[i]->Dot(*g[k]) / Mhost[DENSE_IND(i,i,s,s)];

        // g_k = g_k - alpha * g_i
        g[k]->AddScale(*g[i], ValueType(-1.0)*alpha);

        // u_k = u_k - alpha * u_i
        u[k]->AddScale(*u[i], ValueType(-1.0)*alpha);

      }

      for (int i=k; i<s; ++i) {

        // M_i_k = P^T_i * g_k
        Mhost[DENSE_IND(i,k,s,s)] = P[i]->Dot(*g[k]);

        if (Mhost[DENSE_IND(k,k,s,s)] == ValueType(0.0)) {
          LOG_INFO("IDR(s) break down ; M(k,k) == 0.0");
          FATAL_ERROR(__FILE__, __LINE__);
        }

      }

      // beta = f_k / M_k_k
      beta = fhost[k] / Mhost[DENSE_IND(k,k,s,s)];

      // r = r - beta * g_k
      r->AddScale(*g[k], ValueType(-1.0)*beta);

      // x = x + beta * u_k
      x->AddScale(*u[k], beta);

      // Residual norm
      res_norm = this->Norm(*r);

      // Check inner loop for convergence
      if (this->iter_ctrl_.CheckResidual(paralution_abs(res_norm), this->index_))
        break;

      // f_i = f_i - beta * M_i_k
      for (int i=k+1; i<s; ++i)
        fhost[i] -= beta * Mhost[DENSE_IND(i,k,s,s)];

    }

    // Mv = r
    this->precond_->SolveZeroSol(*r, v);

    // t = Av
    op->Apply(*v, t);

    // omega = (t,r) / ||t||^2
    ValueType rt = t->Dot(*r);
    ValueType nt = t->Norm();

    rt /= nt;

    rho = paralution_abs(rt / res_norm);
    omega = rt / nt;

    if (rho < kappa)
      omega = omega*kappa/rho;

    if (omega == ValueType(0.0)) {
      LOG_INFO("IDR(s) break down ; omega == 0.0");
      FATAL_ERROR(__FILE__, __LINE__);
    }

    // r = r - omega * t
    r->AddScale(*t, ValueType(-1.0)*omega);

    // x = x + omega * v
    x->AddScale(*v, omega);

    // Residual norm to check outer loop convergence
    res_norm = this->Norm(*r);

  }

  LOG_DEBUG(this, "::SolvePrecond_()",
            " #*# end");

}


template class IDR< LocalMatrix<double>, LocalVector<double>, double >;
template class IDR< LocalMatrix<float>,  LocalVector<float>, float >;
#ifdef SUPPORT_COMPLEX
template class IDR< LocalMatrix<std::complex<double> >, LocalVector<std::complex<double> >, std::complex<double> >;
template class IDR< LocalMatrix<std::complex<float> >,  LocalVector<std::complex<float> >,  std::complex<float> >;
#endif

template class IDR< LocalStencil<double>, LocalVector<double>, double >;
template class IDR< LocalStencil<float>,  LocalVector<float>, float >;
#ifdef SUPPORT_COMPLEX
template class IDR< LocalStencil<std::complex<double> >, LocalVector<std::complex<double> >, std::complex<double> >;
template class IDR< LocalStencil<std::complex<float> >,  LocalVector<std::complex<float> >,  std::complex<float> >;
#endif

}
