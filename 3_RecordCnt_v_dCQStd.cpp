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

#include<MariaDB.hpp>
#include<GMT.hpp>
#include<AvrStd.hpp>
#include<CreateGrid.hpp>
#include<Interpolate.hpp>

using namespace std;

// Inputs. -----------------------------

const double maxThickness = 50;
const double dRhoMin = -10, dRhoMax = 20;

const double plotHeight = 6, plotWidth = 5;

const string modelingTable = "gen2CA_D.ModelingResult_Subtract";
const string binTable = "gen2CA_D.Bins";
const string propertyTable = "gen2CA_D.Properties";

// ------------------------------------


struct BinResult{

    int binN;

    // Best fit model.

    vector<double> h, dvs, drho, dcq, cnt;

    double h_bar, dvs_bar, drho_bar, h_std, dvs_std, drho_std, premCQ, cnt_bar;

    BinResult (int num){

        binN=num;
    }
};

int main(){

    // Get bin info.
    auto binInfo = MariaDB::Select("bin from " + binTable);

    vector<BinResult> Data;

    for (size_t i = 0; i < binInfo.NRow(); ++i) {

        Data.push_back(BinResult(binInfo.GetInt("bin")[i]));
    }

    // Get model properties.

    auto modelInfo = MariaDB::Select("modelName, thickness, dvs, drho from " + propertyTable + " where modelName like \"U%\" order by modelName");

    size_t bestFitCnt = (size_t)round(modelInfo.NRow() / 100);


    // For each bin, get the best fit model.
    for (size_t i = 0; i < Data.size(); ++i) {

        auto res = MariaDB::Select("modelName, cq, stackTraceCnt from " + modelingTable + " where bin = " + to_string(Data[i].binN) + " and (modelName like \"U%\" or modelName like \"PREM%\") order by cq desc");

        // Find prem result.
        for (size_t j = 0; j < res.NRow(); ++j) {
            if (res.GetString("modelName")[j] == "PREM_201500000000"){

                Data[i].premCQ = res.GetDouble("cq")[j];
                break;
            }
        }

        // Find top 1% best fit models.
        for (size_t j=0; j<res.NRow(); ++j) {

            if (Data[i].h.size() >= bestFitCnt) {

                break;
            }

            string modelName = res.GetString("modelName")[j];

            if (modelName.find("PREM") != string::npos) {

                continue;
            }


            // get model property.
            auto it = lower_bound(modelInfo.GetString("modelName").begin(), modelInfo.GetString("modelName").end(), modelName);

            if (it == modelInfo.GetString("modelName").end() || *it != res.GetString("modelName")[j]) {

                throw runtime_error("Can't find model property for " + res.GetString("modelName")[j] + " for bin " + to_string(Data[i].binN));
            }
            else {

                size_t index = distance(modelInfo.GetString("modelName").begin(), it);

                if ( dRhoMin <= modelInfo.GetDouble("drho")[index] &&
                     modelInfo.GetDouble("drho")[index] <= dRhoMax && 
                     modelInfo.GetDouble("thickness")[index] <= maxThickness ){


                    Data[i].h.push_back(modelInfo.GetDouble("thickness")[index]);
                    Data[i].dvs.push_back(modelInfo.GetDouble("dvs")[index]);
                    Data[i].drho.push_back(modelInfo.GetDouble("drho")[index]);

                    Data[i].dcq.push_back(res.GetDouble("cq")[j] - Data[i].premCQ);
                    Data[i].cnt.push_back(res.GetInt("stackTraceCnt")[j]);
                }
            }
        }
    }

    // Calculate averages.

    for (auto &item: Data) {

        auto res = AvrStd(item.h);
        item.h_bar = res.first;
        item.h_std= res.second;

        res = AvrStd(item.dvs);
        item.dvs_bar = res.first;
        item.dvs_std= res.second;

        res = AvrStd(item.drho);
        item.drho_bar = res.first;
        item.drho_std= res.second;

        res = AvrStd(item.cnt);
        item.cnt_bar = res.first;

        // sort according to dcq.

        auto sortIndex = SortWithIndex(item.dcq.begin(), item.dcq.end(), std::greater<double> ());

        ReorderUseIndex(item.h.begin(), item.h.end(), sortIndex);
        ReorderUseIndex(item.dvs.begin(), item.dvs.end(), sortIndex);
        ReorderUseIndex(item.drho.begin(), item.drho.end(), sortIndex);
        ReorderUseIndex(item.cnt.begin(), item.cnt.end(), sortIndex);
    }


    // Plot.

    double XSIZE = 2 + plotWidth, YSIZE = 1 + plotHeight * 3;

    string outfile = GMT::BeginEasyPlot(XSIZE, YSIZE);

    // Plot height std vs cnt.

    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf" + to_string(YSIZE - plotHeight) + "i");

    GMT::psbasemap(outfile, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight * 0.7) + "i -R0/2000/0/20 -Bxa500g500 -Bya10g10 -Bx+l\"Averaged stacked trace count.\" -By+l\"Thickness Standard Deviation\" -BWSne -O -K");

    vector<GMT::Text> texts;

    for (size_t i = 0; i < Data.size(); ++i) {

        string color = (Data[i].dvs[0] > 0 ? "blue" : "red");

        GMT::psxy(outfile, vector<double> {Data[i].cnt_bar}, vector<double> {Data[i].h_std}, "-J -R -Sc0.1i -Glight" + color + " -O -K");

        texts.push_back(GMT::Text(Data[i].cnt_bar, Data[i].h_std, to_string(Data[i].binN), 9, "LB"));
    }

    GMT::pstext(outfile, texts, "-J -R -N -O -K");



    // Plot dvs std vs cnt.

    GMT::MoveReferencePoint(outfile,"-Y" + to_string(-plotHeight) + "i");

    GMT::psbasemap(outfile, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight * 0.7) + "i -R0/2000/0/20 -Bxa500g500 -Bya10g10 -Bx+l\"Averaged stacked trace count.\" -By+l\"dVs Standard Deviation\" -BWSne -O -K");

    texts.clear();

    for (size_t i = 0; i < Data.size(); ++i) {

        string color = (Data[i].dvs[0] > 0 ? "blue" : "red");

        GMT::psxy(outfile, vector<double> {Data[i].cnt_bar}, vector<double> {Data[i].dvs_std}, "-J -R -Sc0.1i -Glight" + color + " -O -K");

        texts.push_back(GMT::Text(Data[i].cnt_bar, Data[i].dvs_std, to_string(Data[i].binN), 9, "LB"));
    }

    GMT::pstext(outfile, texts, "-J -R -N -O -K");




    // Plot dvs std vs cnt.

    GMT::MoveReferencePoint(outfile,"-Y" + to_string(-plotHeight) + "i");

    GMT::psbasemap(outfile, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight * 0.7) + "i -R0/2000/0/10 -Bxa500g500 -Bya10g10 -Bx+l\"Averaged stacked trace count.\" -By+l\"Density Standard Deviation\" -BWSne -O -K");

    texts.clear();

    for (size_t i = 0; i < Data.size(); ++i) {

        string color = (Data[i].dvs[0] > 0 ? "blue" : "red");

        GMT::psxy(outfile, vector<double> {Data[i].cnt_bar}, vector<double> {Data[i].drho_std}, "-J -R -Sc0.1i -Glight" + color + " -O -K");

        texts.push_back(GMT::Text(Data[i].cnt_bar, Data[i].drho_std, to_string(Data[i].binN), 9, "LB"));
    }

    GMT::pstext(outfile, texts, "-J -R -N -O -K");


    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile, __FILE__);

    return 0;
}
