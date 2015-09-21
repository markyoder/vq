// Copyright (c) 2012-2014 Eric M. Heien, Michael K. Sachs, John B. Rundle
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "GreensFunctions.h"

#ifndef _GREENS_FUNC_CALC_H_
#define _GREENS_FUNC_CALC_H_

/*!
 Initializes the simulation greens values using the user specified Greens function and parameters.
 */
class GreensInit : public SimPlugin {
    public:
        virtual std::string name(void) const {
            return "Greens function calculation";
        }
        virtual void initDesc(const SimFramework *_sim) const {};

        virtual void dryRun(SimFramework *_sim);
        virtual void init(SimFramework *_sim);
        //
        // yoder: also get some greens stats:
        void getGreensStats(Simulation *sim, double &shear_min, double &shear_max, double &shear_mean, double &normal_min, double &normal_max, double &normal_mean);
        void getGreensDiagStats(Simulation *sim, double &shear_diag_min, double &shear_diag_max, double &shear_diag_mean, double &normal_diag_min, double &normal_diag_max, double &normal_diag_mean, double &shear_offdiag_min, double &shear_offdiag_max, double &shear_offdiag_mean, double &normal_offdiag_min, double &normal_offdiag_max, double &normal_offdiag_mean);

};

#endif
