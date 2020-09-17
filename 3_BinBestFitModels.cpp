#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<unistd.h>
#include<iomanip>
#include<set>

#include<EvenSampledSignal.hpp>
#include<MariaDB.hpp>
#include<GMTPlotSignal.hpp>
#include<ShellExec.hpp>

using namespace std;

// Inputs. -----------------------------

const size_t topN= 10;

const double dCQThreshold=0.45;

const string modelingTable = "gen2CA_D.ModelingResult_Subtract";
const string binTable = "gen2CA_D.Bins";
const string propertyTable = "gen2CA_D.Properties";

const double plotHeight = 4.901875, plotWidth = 5;

// must be a full set.
const vector<double> rhoDimensions{-5, -2.5, 0, 2.5, 5, 7.5, 10};
const vector<string> penColor{"darkred", "cyan", "black", "lightblue", "darkblue", "lightgreen", "darkgreen"};

// const vector<size_t> plotTheseBins {};
const vector<size_t> plotTheseBins {28,1,4,9,13,14,15,16,18,22,23,38,40,41};


// ------------------------------------

int main() {

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
        auto modelingResult = MariaDB::Select("modelName, cq from " + modelingTable + " where bin=" + to_string(binN) + " and modelName not like \"Lamella_%\" order by cq desc");

        // Find PREM cq.
        double premCQ = 0;
        for (size_t j = 0; j < modelingResult.NRow(); ++j) {
            if (modelingResult.GetString("modelName")[j] == "PREM_201500000000") {
                premCQ = modelingResult.GetDouble("cq")[j];
                break;
            }
        }

        // Find each density category data grid.
        vector<vector<vector<double>>> goodData(rhoDimensions.size()), badData = goodData;
        vector<size_t> cnt(rhoDimensions.size(), 0);
        for (size_t j = 0; j < modelingResult.NRow(); ++j) {

            auto properties = modelNameToProperties [ modelingResult.GetString("modelName")[j] ];

            size_t k = distance(rhoDimensions.begin(), lower_bound(rhoDimensions.begin(), rhoDimensions.end(), properties[2]-0.01));
            if (k == rhoDimensions.size()) {
                continue;
            }
            double cq =  modelingResult.GetDouble("cq")[j];

            if (cnt[k] >= topN) {
                continue;
            }

            ++cnt[k];

            if (cq < premCQ + dCQThreshold) {
                badData[k].push_back({properties[1], properties[0], 0.05});
            }
            else {
                goodData[k].push_back({properties[1], properties[0], cq/6.0});
            }

        }



        // plot image.

        GMT::MoveReferencePoint(outfile, "-Xf0.5i");

        for (size_t k = 0; k < rhoDimensions.size(); ++k) {
            GMT::psxy(outfile, badData[k],  "-JX" + to_string(plotWidth * 0.8 ) + "i/" + to_string(plotHeight * 0.8) + "i -R-30.5/20.5/0.5/50.5 -Sx -W1p," + penColor[k] + " -O -K");
        }
        for (size_t k = 0; k < rhoDimensions.size(); ++k) {

            GMT::psxy(outfile, goodData[k], "-JX" + to_string(plotWidth * 0.8 ) + "i/" + to_string(plotHeight * 0.8) + "i -R-30.5/20.5/0.5/50.5 -Sc -W1p," + penColor[k] + " -O -K");
        }

        GMT::psbasemap(outfile, "-J -R -Bxa5 -Bya5 -BWSne -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-1, 0.8, "Bin "+to_string(binN),12,"LT"));
        GMT::MoveReferencePoint(outfile, "-Xf0.5i");
        GMT::pstext(outfile, texts, "-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight)+"i -R-1/1/-1/1 -N -O -K");

    }

    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile,__FILE__);

    return 0;
}

//     // Plot best fit model properties for the significant bins.
//     vector<GMT::Text> texts;
//     for (size_t i=0;i<Data.size();++i)
//         if (Data[i].dcq>dCQThreshold)
//             texts.push_back(GMT::Text{Data[i].dvs,Data[i].h,to_string(Data[i].binN),12,"CM"});
//
//
//
//     GMT::MoveReferencePoint(outfile,"-Y-6i");
//     GMT::psbasemap(outfile,"-Jx0.14i -R-35/25/0/35 -Bxa5g5f1 -Bya5g5f1 -Bx+l\"dVs (%)\" -By+l\"Thickness (km)\" -BWS -O -K");
//     GMT::pstext(outfile,texts,"-J -R -O -K");
//
//     for (size_t i=0;i<Data.size();++i) {
//         if (Data[i].dcq<=dCQThreshold) continue;
//         vector<vector<double>> best5Data;
//         for (auto item:Data[i].best5Properties)
//             best5Data.push_back({item.first,item.second,double(Data[i].binN),0.05});
//         GMT::psxy(outfile,best5Data,"-J -R -Ctmp.cpt -Sc -O -K");
//     }
//
