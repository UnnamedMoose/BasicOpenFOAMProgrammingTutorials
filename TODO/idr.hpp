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


#ifndef PARALUTION_KRYLOV_IDR_HPP_
#define PARALUTION_KRYLOV_IDR_HPP_

#include "../solver.hpp"

#include <vector>

namespace paralution {

/// IDR(s) - Induced Dimension Reduction method, taken from "An Elegant IDR(s)
/// Variant that Efficiently Exploits Biorthogonality Properties" by Martin B.
/// van Gijzen and Peter Sonneveld, Delft University of Technology
template <class OperatorType, class VectorType, typename ValueType>
class IDR : public IterativeLinearSolver<OperatorType, VectorType, ValueType> {

public:

  IDR();
  virtual ~IDR();

  virtual void Print(void) const;

  virtual void Build(void);
  virtual void Clear(void);

  /// Set the size of the Shadow Space
  virtual void SetShadowSpace(const int s);

protected:

  virtual void SolveNonPrecond_(const VectorType &rhs,
                                VectorType *x);
  virtual void SolvePrecond_(const VectorType &rhs,
                             VectorType *x);

  virtual void PrintStart_(void) const;
  virtual void PrintEnd_(void) const;

  virtual void MoveToHostLocalData_(void);
  virtual void MoveToAcceleratorLocalData_(void);

private:

  int s_;

  ValueType kappa_;

  ValueType *fhost_, *Mhost_;

  VectorType r_, v_, z_, t_;
  VectorType **g_, **u_, **P_;

};


}

#endif // PARALUTION_KRYLOV_IDR_HPP_
