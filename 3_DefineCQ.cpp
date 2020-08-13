#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<string>
#include<unistd.h>
#include<iomanip>
#include<set>

#include<EvenSampledSignal.hpp>
#include<MariaDB.hpp>
#include<GMTPlotSignal.hpp>
#include<ShellExec.hpp>
#include<Float2String.hpp>

#include "CalculateCQ.hpp"

using namespace std;

/********************************************************
 *
 * This code demostrate our definition of the comparison
 * quality tackles all possible scenarios.
 *
********************************************************/


// Inputs. ---------------------------------

const double signalLen = 40, compareLen = 10, dt = 0.025, sigma = 1.2;
const double plotWidth = 5, plotHeight = 3;
const vector<double> amps{1, 0.75, 0.5, 0.33333, 0, -0.33333, -0.5, -0.75, -1};
const vector<double> times{-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6};
const vector<double> stretches{0.5, 0.6666, 0.8, 1, 1.25, 1.5};
const vector<double> dcs{-0.1, 0.0, 0.1};

// -----------------------------------------


int main(){


    // Create plot data.

    EvenSampledSignal gaussianSignalCenteredOnPeak;

    gaussianSignalCenteredOnPeak = EvenSampledSignal(GaussianSignal((size_t)floor(signalLen / dt), dt, sigma), dt, 0);
    gaussianSignalCenteredOnPeak.FindPeakAround(signalLen / 2.0);
    gaussianSignalCenteredOnPeak.NormalizeToPeak();
    gaussianSignalCenteredOnPeak.ShiftTimeReferenceToPeak();


    EvenSampledSignal dataSignal, modelSignal;
    vector<EvenSampledSignal> data, model;

    dataSignal = gaussianSignalCenteredOnPeak;
    dataSignal.ShiftTime(compareLen / 2.0);

    // ---- Amplitude varying.

//     for (auto amp: amps) {
//         data.push_back(dataSignal);
//         model.push_back(dataSignal * amp);
//     }

    // ---- Time shift.
//     for (auto time: times) {
//         data.push_back(dataSignal);
//         model.push_back(dataSignal);
//         model.back().ShiftTime(time);
//     }

    // ---- Stretch / shrink.

//     for (auto stretch: stretches) {
//         data.push_back(dataSignal);
//         model.push_back(dataSignal.Stretch(stretch));
//     }

    // ---- DC shift

    for (auto dc: dcs) {
        data.push_back(dataSignal);
        model.push_back(dataSignal + dc);
    }

    // Plot
    double YSIZE = plotHeight * data.size() + 1, XSIZE = 2 + plotWidth;
    string outfile = GMT::BeginEasyPlot(XSIZE, YSIZE, ShellExec("pwd", true) + "/" + string(__FILE__));
    GMT::MoveReferencePoint(outfile, "-Xf1i -Yf" + to_string(YSIZE - 0.5) + "i");
    vector<GMT::Text> texts;
    double cq;
    vector<double> compareResult;


    for (size_t i = 0; i < data.size(); ++i) {

        GMT::MoveReferencePoint(outfile, "-Y-" + to_string(plotHeight) + "i");

        dataSignal = data[i];
        modelSignal = model[i];


        GMT::psxy(outfile, vector<double> {0, 0, compareLen, compareLen}, vector<double> {-1, 1, 1, -1}, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -Glightblue -O -K");
        GMT::psxy(outfile, vector<double> {0, 20}, vector<double> {0, 0}, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -W1p,gray,- -O -K");
        GMT::psxy(outfile, modelSignal, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -W2p,orange -O -K");
        GMT::psxy(outfile, dataSignal, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -W1p,black -O -K");
        GMT::psbasemap(outfile, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -Bxa5f1 -Bx+l\"Time (sec.)\" -Bya0.5f0.1 -By+l\"Amplitude\" -BWS -O -K");

        if ( i > 0 ) {
            GMT::MoveReferencePoint(outfile, "-Y" + to_string(i * plotHeight) + "i");
            GMT::psxy(outfile, modelSignal, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -W2p,orange -O -K");
            GMT::MoveReferencePoint(outfile, "-Y-" + to_string(i * plotHeight) + "i");
        }


        compareResult = CalculateCQ(dataSignal, modelSignal, compareLen);
        cq = compareResult[0] * compareResult[1];

        texts.clear();
        texts.push_back(GMT::Text(compareLen + 0.5, 1.0, "CQ: " + Float2String(cq, 3), 12 , "LT"));
        texts.push_back(GMT::Text(compareLen + 0.5, 0.75, "CC: " + Float2String(compareResult[3], 3) + ", CC': " + Float2String(compareResult[4], 3), 12 , "LT"));
        texts.push_back(GMT::Text(compareLen + 0.5, 0.5, "CC\": " + Float2String(compareResult[0], 3), 12 , "LT"));
        texts.push_back(GMT::Text(compareLen + 0.5, 0.25, "N2: " + Float2String(compareResult[5], 3) + ", N2\": " + Float2String(compareResult[1], 3), 12 , "LT"));
        GMT::pstext(outfile, texts, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight*0.8) + "i -R0/15/-1.1/1.1 -N -O -K");

    }

    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile,__FILE__);

    return 0;
}
