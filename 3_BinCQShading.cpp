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

using namespace std;

// Inputs. -----------------------------

const string modelingTable = "gen2CA_D.ModelingResult_Subtract_cutoff70";
const string binTable = "gen2CA_D.Bins";
const string propertyTable="gen2CA_D.Properties";

const double plotHeight = 4.901875, plotWidth = 35;
const vector<double> rhoDimensions{-5, -2.5, 0, 2.5, 5, 7.5, 10};
const double xinc = 1, yinc = 1;

const vector<size_t> plotTheseBins {1,4,9,13,14,15,16,17,18,22,23,41};
//const vector<size_t> plotTheseBins {};

// ------------------------------------

int main(){


    // For each bin, get the location and bin number.
    auto binInfo = MariaDB::Select("bin from "+binTable);


    // Get modelname - property map.
    map<string, vector<double>> modelNameToProperties;

    auto modelInfo = MariaDB::Select("modelname, thickness, dvs, drho from "+propertyTable);
    for (size_t i=0; i<modelInfo.NRow(); ++i) {
        modelNameToProperties[ modelInfo.GetString("modelname")[i] ] = {modelInfo.GetDouble("thickness")[i], modelInfo.GetDouble("dvs")[i], modelInfo.GetDouble("drho")[i]};
    }

    // Plot
    double YSIZE = plotHeight * binInfo.NRow() + 1;
    if (!plotTheseBins.empty()) {
        YSIZE = plotHeight * plotTheseBins.size() + 1;
    }
    double XSIZE = 1 + plotWidth;
    string outfile = GMT::BeginEasyPlot(XSIZE, YSIZE, ShellExec("pwd",true) + "/" + string(__FILE__));
    GMT::MoveReferencePoint(outfile, "-Xf0.5i -Yf" + to_string(YSIZE-0.5) + "i");

    ofstream fpout("tmp.cpt");
    fpout <<
"0   white   1   black\n"
"B   white\n"
"F   black\n"
"N   white";
    fpout.close();

    // Plot bin result.
    for (size_t i=0; i<binInfo.NRow(); ++i) {

        if (!plotTheseBins.empty()) {

            if (find(plotTheseBins.begin(), plotTheseBins.end(), i + 1) == plotTheseBins.end()) {
                continue;
            }
        }


        GMT::MoveReferencePoint(outfile, "-Y-" + to_string(plotHeight) + "i");

        int binN = binInfo.GetInt("bin")[i];

        // ModelName is in the form: "ModelType_2015xxx"
        auto modelingResult = MariaDB::Select("modelName, cq from " + modelingTable + " where bin=" + to_string(binN));

        // Find each density category data grid.
        vector<vector<vector<double>>> plotData(rhoDimensions.size());
        for (size_t j=0; j<modelingResult.NRow(); ++j) {

            auto properties = modelNameToProperties [ modelingResult.GetString("modelName")[j] ];
            size_t k = distance(rhoDimensions.begin(), lower_bound(rhoDimensions.begin(), rhoDimensions.end(), properties[2]-0.01));
            plotData[k].push_back({properties[1], properties[0], modelingResult.GetDouble("cq")[j]});
        }


        // plot image.
        for (size_t j=0; j<rhoDimensions.size(); ++j) {

            GMT::MoveReferencePoint(outfile, "-Xf" + to_string(0.5 + plotWidth / rhoDimensions.size() * j ) + "i");

            GMT::grdimage(outfile, plotData[j], xinc, yinc, "-JX" + to_string(plotWidth / rhoDimensions.size() * 0.8 ) + "i/" + to_string(plotHeight*0.8) + "i -R-30.5/20.5/0.5/50.5 -Ctmp.cpt -O -K");

            GMT::psbasemap(outfile, "-J -R -Bxa5 -Bya5 -BWSne -O -K");
        }

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-1, 0.8, "Bin "+to_string(binN),12,"LT"));
        GMT::MoveReferencePoint(outfile, "-Xf0.5i");
        GMT::pstext(outfile,texts,"-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight)+"i -R-1/1/-1/1 -N -O -K");

    }

    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile,__FILE__);

    return 0;
}
