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

const string modelingTable="gen2CA_D.ModelingResult_Subtract_cutoff70";
const string binTable="gen2CA_D.Bins";
const double plotHeight = 4.901875 / 2.0, plotWidth = 4, binSize=0.05;
//const vector<size_t> plotTheseBins {1,4,9,13,14,15,16,17,18,22,23,41};
const vector<size_t> plotTheseBins {};

// ------------------------------------

size_t getIndex(double val){
    return (size_t)floor(val/binSize);
}

int main(){


    // For each bin, get the location and bin number.
    auto binInfo=MariaDB::Select("bin from "+binTable);


    // Plot
    double YSIZE=plotHeight * binInfo.NRow() + 1, XSIZE= 1 + plotWidth;
    string outfile=GMT::BeginEasyPlot(XSIZE,YSIZE,ShellExec("pwd",true)+"/"+string(__FILE__));
    GMT::MoveReferencePoint(outfile,"-Xf0.5i -Yf"+to_string(YSIZE-0.5)+"i");
        

    // Plot bin result.

    for (size_t i=0; i<binInfo.NRow(); ++i) {

        if (!plotTheseBins.empty()) {

            if (find(plotTheseBins.begin(), plotTheseBins.end(), i + 1) == plotTheseBins.end()) {
                continue;
            }
        }

        GMT::MoveReferencePoint(outfile,"-Y-"+to_string(plotHeight)+"i");

        int binN=binInfo.GetInt("bin")[i];

        // ModelName is in the form: "ModelType_2015xxx"
        auto res=MariaDB::Select("modelName, cq from " + modelingTable + " where bin=" + to_string(binN));

        // Find each category number
        double premCQ=-1;
        vector<double> ulvzVal(getIndex(1.0),0), uhvzVal=ulvzVal;

        for (size_t j=0; j<res.NRow(); ++j) {
            double x=res.GetDouble("cq")[j];
            size_t index=getIndex(x);
            if (res.GetString("modelName")[j]=="PREM_201500000000"){
                premCQ=x;
            }
            else if (res.GetString("modelName")[j].find("ULVZ")!=string::npos){
                ++ulvzVal[index];
            }
            else {
                ++uhvzVal[index];
            }
        }

        // plot histogram.
        double x1=0, x2=binSize;
        for (size_t j=0; j<ulvzVal.size(); ++j) {
            ulvzVal[j]=ulvzVal[j]/res.NRow()*100;
            uhvzVal[j]=uhvzVal[j]/res.NRow()*100+ulvzVal[j];

            GMT::psxy(outfile, vector<double> {x1, x1, x2, x2}, vector<double> {0, ulvzVal[j], ulvzVal[j], 0}, "-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight*0.8)+"i -R0/1/0/50 -W0.247p,black -Gred -O -K");
            GMT::psxy(outfile, vector<double> {x1, x1, x2, x2}, vector<double> {ulvzVal[j], uhvzVal[j], uhvzVal[j], ulvzVal[j]}, "-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight*0.8)+"i -R0/1/0/50 -W0.247p,black -Gblue -O -K");
            x1=x2;
            x2+=binSize;
        }

        GMT::psxy(outfile, vector<double> {premCQ, premCQ}, vector<double> {0, 100}, "-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight*0.8)+"i -R0/1/0/50 -W2p,orange -O -K");
        GMT::psbasemap(outfile, "-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight*0.8)+"i -R0/1/0/50 -Bx -By  -Bws -O -K");
        //GMT::psbasemap(outfile, "-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight*0.8)+"i -R0/1/0/50 -Bxa0.5f0.1 -Bx+l\"cq\" -Bya10f10 -By+l\"\%\" -BWS -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-1, 0.8, "Bin "+to_string(binN),12,"LT"));
        GMT::pstext(outfile,texts,"-JX"+to_string(plotWidth)+"i/"+to_string(plotHeight)+"i -R-1/1/-1/1 -N -O -K");

    }

    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile,__FILE__);

    return 0;
}
