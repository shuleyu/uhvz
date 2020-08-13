#ifndef ASU_CQ
#define ASU_CQ

#include<cmath>
#include<vector>

#include<EvenSampledSignal.hpp>
#include<PNormErr.hpp>

/*************************************************
 * This C++ template returns the "Compare Quality"
 * of two input signals.
 *
 * Because it's very specific, it's not included
 * in the CPP library.
 *
 * input(s):
 * const vector<T> &p  ----  Input array.
 *
 * output(s):
 * double cq  ----  comparison result.
 *
 * Shule Yu
 * Feb 12 2020
 *
 * Key words: comparison quality
*************************************************/

std::vector<double> CalculateCQ(const EvenSampledSignal &dataTrace, const EvenSampledSignal &modelTrace, const double &compareLength){

    // Quality 1. cross-correlation.

    double cc1 = dataTrace.CrossCorrelation(0, compareLength, modelTrace, 0, compareLength, 1, std::make_pair(0,0)).second;

    // take amplitude into consideration.

    double ed = sqrt(dataTrace.SumArea(0, compareLength, 2));
    double es = sqrt(modelTrace.SumArea(0, compareLength, 2));

    double cc2 = cc1 * std::min(ed, es) / std::max(ed, es);

    // convert [ -1(worst) ~ 1(best) ] to [ 0(worst) ~ 1(best) ]
    double cc = (1 + cc2) / 2;


    // Quality 2. Norm-2 difference.

    // |x-y|^2 / |y|^2
    // currently using this:
    double nn2 = PNormErr(modelTrace.GetAmp(0, compareLength), dataTrace.GetAmp(0, compareLength), 2);

    // should be this?
    double nn2x = PNormErr(dataTrace.GetAmp(0, compareLength), modelTrace.GetAmp(0, compareLength), 2);


    // [ 0(best) ~ inf(worst) ] to [ 1(best) ~ 0(worst) ].
    double nd = 1.0 / (1.0 + nn2);
    double ndx = 1.0 / (1.0 + nn2x);

    return std::vector<double> {cc, nd, ndx, cc1, cc2, nn2, nn2x};
}

#endif
