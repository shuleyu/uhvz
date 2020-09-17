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

const string modelingTable = "gen2CA_D.ModelingResult_Subtract_cutoff70";
const string binTable = "gen2CA_D.Bins";
const string propertyTable = "gen2CA_D.Properties";

// ------------------------------------


struct BinResult{

    int binN;

    // Best fit model.

    vector<double> h, dvs, drho, dcq;

    double h_bar, dvs_bar, drho_bar, h_std, dvs_std, drho_std, premCQ;

    BinResult (int num){
        binN=num;
    }
};

int main(){

    // For each bin, get the location and bin number.
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

        // ModelName is in the form: "ModelType_2015xxx"
        auto res = MariaDB::Select("modelName, cq from " + modelingTable + " where bin = " + to_string(Data[i].binN) + " and (modelName like \"U%\" or modelName like \"PREM%\") order by cq desc");

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
                }
            }
        }
    }

    // Calculate averages.

    for (auto &item: Data) {

        auto res = AvrStd(item.h);
        item.h_bar = res.first;
        item.h_std = res.second;

        res = AvrStd(item.dvs);
        item.dvs_bar = res.first;
        item.dvs_std = res.second;

        res = AvrStd(item.drho);
        item.drho_bar = res.first;
        item.drho_std = res.second;

        // sort according to dcq.

        auto sortIndex = SortWithIndex(item.dcq.begin(), item.dcq.end(), std::greater<double> ());

        ReorderUseIndex(item.h.begin(), item.h.end(), sortIndex);
        ReorderUseIndex(item.dvs.begin(), item.dvs.end(), sortIndex);
        ReorderUseIndex(item.drho.begin(), item.drho.end(), sortIndex);
    }

    // read in bk's data.

    vector<double> diamond_drho, diamond_dvs, diamond_vol;
    double rho, vs, vol;
    ifstream fpin("BK.txt");

    while (fpin >> rho >> vs >> vol){
        diamond_drho.push_back(rho);
        diamond_dvs.push_back(vs);
        diamond_vol.push_back(vol);
    }

    vector<double> diamondVol_xx_drho = CreateGrid(-7.5, 0, 100, 0);
    vector<double> diamondVol_xx_dvs = CreateGrid(0, 20, 100, 0);

    auto diamondVol_yy_drho = Interpolate(diamond_drho, diamond_vol, diamondVol_xx_drho, true);
    auto diamondVol_yy_dvs = Interpolate(diamond_dvs, diamond_vol, diamondVol_xx_dvs, true);


    fpin.close();


    // Plot.

    double XSIZE = 2 + plotWidth, YSIZE = 1 + plotHeight * (Data.size() + 1);

    string outfile = GMT::BeginEasyPlot(XSIZE, YSIZE);


    for (size_t i = 0; i < Data.size(); ++i) {

        GMT::MoveReferencePoint(outfile,"-Xf1i -Yf" + to_string(1 + (Data.size() - i) * plotHeight) + "i");

        GMT::psbasemap(outfile, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight * 0.7) + "i -R-7.5/12.5/-30.5/20.5 -Bxa2.5g2.5 -Bya5g5 -Bx+ldrho(\%) -By+ldvs(%) -BWSne -O -K");

        GMT::psxy(outfile, vector<double>{0, 0}, vector<double>{-30.5, 20.5}, "-J -R -W0.5p,black -O -K");
        GMT::psxy(outfile, vector<double>{-7.5, 12.5}, vector<double>{0, 0}, "-J -R -W0.5p,black -O -K");

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(-10, 25, "Bin " + to_string(Data[i].binN), 12, "LB"));
        GMT::pstext(outfile, texts, "-J -R -N -O -K");

        GMT::psxy(outfile, diamond_drho, diamond_dvs, "-J -R -W2p,lightgreen -O -K");



        for (size_t j = 0; j < Data[i].h.size(); ++j) {

            if (Data[i].dcq[j] >= 0.45) {

                GMT::psxy(outfile, vector<double> {Data[i].drho[j]}, vector<double> {Data[i].dvs[j]}, "-J -R -Sc0.1i -Gblack -O -K");
            }
            else if (0.35 < Data[i].dcq[j] && Data[i].dcq[j] < 0.45) {

                GMT::psxy(outfile, vector<double> {Data[i].drho[j]}, vector<double> {Data[i].dvs[j]}, "-J -R -Sc0.05i -W0.5p,gray -O -K");
            }
            else if (0.25 < Data[i].dcq[j] && Data[i].dcq[j] < 0.35) {

                GMT::psxy(outfile, vector<double> {Data[i].drho[j]}, vector<double> {Data[i].dvs[j]}, "-J -R -Sc0.05i -W1p,black,. -O -K");
            }
            else if (0.15 < Data[i].dcq[j] && Data[i].dcq[j] < 0.25) {

                GMT::psxy(outfile, vector<double> {Data[i].drho[j]}, vector<double> {Data[i].dvs[j]}, "-J -R -Sc0.05i -W0.5p,red -O -K");
            }

        }





        GMT::psxy(outfile, vector<double> {Data[i].drho_bar - Data[i].drho_std, Data[i].drho_bar + Data[i].drho_std}, vector<double> {Data[i].dvs_bar, Data[i].dvs_bar}, "-J -R -W1p,blue -O -K");

        GMT::psxy(outfile, vector<double> {Data[i].drho_bar, Data[i].drho_bar}, vector<double> {Data[i].dvs_bar - Data[i].dvs_std, Data[i].dvs_bar + Data[i].dvs_std}, "-J -R -W1p,blue -O -K");

        if (Data[i].dvs_bar > 0) {

            GMT::psbasemap(outfile, "-JX-" + to_string(plotWidth / 20 * 7.5) + "i/1i -R" + to_string(diamondVol_yy_drho.back()) + "/" + to_string(diamondVol_yy_drho[0]) + "/-1/1 -Bxa5 -Bx+lvol(\%) -BS -Y0.5i -O -K");
        }
        else {

            GMT::psbasemap(outfile, "-JX0.01i/" + to_string(plotHeight * 0.7 / 51 * 20.5) + "i -R-1/1/" + to_string(diamondVol_yy_dvs[0]) + "/" + to_string(diamondVol_yy_dvs.back()) + " -Bya5 -By+lvol(\%) -BE -X3.75i -Y" + to_string(plotHeight * 0.7 / 51 * 30.5) + "i -O -K");
        }
    }

    // All bins (error bar only).


    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf1i");

    GMT::psbasemap(outfile, "-JX" + to_string(plotWidth) + "i/" + to_string(plotHeight * 0.7) + "i -R-7.5/12.5/-30.5/20.5 -Bxa2.5g2.5 -Bya5g5 -Bx+ldrho(\%) -By+ldvs(%) -BWSne -O -K");

    GMT::psxy(outfile, vector<double>{0, 0}, vector<double>{-30.5, 20.5}, "-J -R -W0.5p,black -O -K");
    GMT::psxy(outfile, vector<double>{-7.5, 12.5}, vector<double>{0, 0}, "-J -R -W0.5p,black -O -K");

    vector<GMT::Text> texts;
    texts.push_back(GMT::Text(-10, 25, "All Bins", 12, "LB"));
    GMT::pstext(outfile, texts, "-J -R -N -O -K");

    GMT::psxy(outfile, diamond_drho, diamond_dvs, "-J -R -W2p,lightgreen -O -K");

//     for (size_t i = 0; i < Data.size(); ++i) {
// 
//         if (Data[i].binN == 1) {
//             Data[i].h[0] = 50;
//             Data[i].dvs[0] = 2;
//             Data[i].drho[0] = -2.5;
//         }
//         if (Data[i].binN == 4) {
//             Data[i].h[0] = 39;
//             Data[i].dvs[0] = 12;
//             Data[i].drho[0] = -2.5;
//         }
//         else if (Data[i].binN == 9) {
//             Data[i].h[0] = 41;
//             Data[i].dvs[0] = 11;
//             Data[i].drho[0] = -5;
//         }
//         else if (Data[i].binN == 13) {
//             Data[i].h[0] = 37;
//             Data[i].dvs[0] = 4;
//             Data[i].drho[0] = -5;
//         }
//         else if (Data[i].binN == 14) {
//             Data[i].h[0] = 47;
//             Data[i].dvs[0] = 5;
//             Data[i].drho[0] = 7.5;
//         }
//         else if (Data[i].binN == 15) {
//             Data[i].h[0] = 49;
//             Data[i].dvs[0] = 1;
//             Data[i].drho[0] = -5;
//         }
//         else if (Data[i].binN == 16) {
//             Data[i].h[0] = 50;
//             Data[i].dvs[0] = 2;
//             Data[i].drho[0] = -5;
//         }
//         else if (Data[i].binN == 17) {
//             Data[i].h[0] = 33;
//             Data[i].dvs[0] = 9;
//             Data[i].drho[0] = 7.5;
//             Data[i].dcq[0] = 0.46;
//         }
//         else if (Data[i].binN == 18) {
//             Data[i].h[0] = 38;
//             Data[i].dvs[0] = 4;
//             Data[i].drho[0] = -5;
//         }
//         else if (Data[i].binN == 22) {
//             Data[i].h[0] = 50;
//             Data[i].dvs[0] = 18;
//             Data[i].drho[0] = 5;
//         }
//         else if (Data[i].binN == 23) {
//             Data[i].h[0] = 50;
//             Data[i].dvs[0] = 6;
//             Data[i].drho[0] = 5;
//         }
//         else if (Data[i].binN == 28) {
//             Data[i].h[0] = 6;
//             Data[i].dvs[0] = -16;
//             Data[i].drho[0] = 10;
//         }
//         else if (Data[i].binN == 41) {
//             Data[i].h[0] = 46;
//             Data[i].dvs[0] = 2;
//             Data[i].drho[0] = 2.5;
//         }
// 
//         else if (Data[i].binN == 19 || Data[i].binN == 38 || Data[i].binN == 40) {
// 
//             Data[i].dcq[0] = 0.44;
//         }
// cout << Data[i].binN << " - " << Data[i].dcq[0] << endl;
// 
//     }

    auto cmp = [](const BinResult &d1, const BinResult &d2){
        if (d1.dcq == d2.dcq) {
            return d1.binN < d2.binN;
        }
        return d1.dcq > d2.dcq;
    };
    auto cmp2 = [](const BinResult &d1, const BinResult &d2){
        return d1.binN < d2.binN;
    };

    sort(Data.begin(), Data.end(), cmp2);


    for (size_t i = 0; i < Data.size(); ++i) {

        string color = "red";
        if (Data[i].dvs[0] > 0) {
            color = "blue";
        }

        if (Data[i].dcq[0] > 0.45) {

cout << "Group 1 - Bin: " << Data[i].binN << " " << Data[i].dcq[0] << " " << Data[i].h[0] << " " << Data[i].h_std << " " << Data[i].dvs[0] << " " << Data[i].dvs_std << " " << Data[i].drho[0] << " " << Data[i].drho_std << endl;

            GMT::psxy(outfile, vector<double> {Data[i].drho[0] - Data[i].drho_std, Data[i].drho[0] + Data[i].drho_std}, vector<double> {Data[i].dvs[0], Data[i].dvs[0]}, "-J -R -W2p," + color + " -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0], Data[i].drho[0]}, vector<double> {Data[i].dvs[0] - Data[i].dvs_std, Data[i].dvs[0] + Data[i].dvs_std}, "-J -R -W2p," + color + " -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0]}, vector<double> {Data[i].dvs[0]}, "-J -R -Sc0.10i -G" + color + " -O -K");

        }
        else if (0.40 < Data[i].dcq[0] && Data[i].dcq[0] <= 0.45) {

cout << "Group 2 - Bin: " << Data[i].binN << " " << Data[i].dcq[0] << " " << Data[i].h[0] << " " << Data[i].h_std << " " << Data[i].dvs[0] << " " << Data[i].dvs_std << " " << Data[i].drho[0] << " " << Data[i].drho_std << endl;

            GMT::psxy(outfile, vector<double> {Data[i].drho[0] - Data[i].drho_std, Data[i].drho[0] + Data[i].drho_std}, vector<double> {Data[i].dvs[0], Data[i].dvs[0]}, "-J -R -W1.234p,light" + color + " -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0], Data[i].drho[0]}, vector<double> {Data[i].dvs[0] - Data[i].dvs_std, Data[i].dvs[0] + Data[i].dvs_std}, "-J -R -W1.234p,light" + color + " -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0]}, vector<double> {Data[i].dvs[0]}, "-J -R -Sc0.05i -Glight" + color + " -O -K");

        }
        else if (0.35 < Data[i].dcq[0] && Data[i].dcq[0] <= 0.40) {

cout << "Group 3 - Bin: " << Data[i].binN << " " << Data[i].dcq[0] << " " << Data[i].h[0] << " " << Data[i].h_std << " " << Data[i].dvs[0] << " " << Data[i].dvs_std << " " << Data[i].drho[0] << " " << Data[i].drho_std << endl;

            GMT::psxy(outfile, vector<double> {Data[i].drho[0] - Data[i].drho_std, Data[i].drho[0] + Data[i].drho_std}, vector<double> {Data[i].dvs[0], Data[i].dvs[0]}, "-J -R -W1.123p,light" + color + ",- -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0], Data[i].drho[0]}, vector<double> {Data[i].dvs[0] - Data[i].dvs_std, Data[i].dvs[0] + Data[i].dvs_std}, "-J -R -W1.123p,light" + color + ",- -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0]}, vector<double> {Data[i].dvs[0]}, "-J -R -Sc0.1i -W1.123p,black -O -K");
        }
        else if (0.30 < Data[i].dcq[0] && Data[i].dcq[0] <= 0.35) {

cout << "Group 4 - Bin: " << Data[i].binN << " " << Data[i].dcq[0] << " " << Data[i].h[0] << " " << Data[i].h_std << " " << Data[i].dvs[0] << " " << Data[i].dvs_std << " " << Data[i].drho[0] << " " << Data[i].drho_std << endl;

            GMT::psxy(outfile, vector<double> {Data[i].drho[0] - Data[i].drho_std, Data[i].drho[0] + Data[i].drho_std}, vector<double> {Data[i].dvs[0], Data[i].dvs[0]}, "-J -R -W0.5p,light" + color + ",. -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0], Data[i].drho[0]}, vector<double> {Data[i].dvs[0] - Data[i].dvs_std, Data[i].dvs[0] + Data[i].dvs_std}, "-J -R -W0.5p,light" + color + ",. -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0]}, vector<double> {Data[i].dvs[0]}, "-J -R -Sc0.05i -W1p,black -O -K");
        }

        else {

cout << "Group 5 - Bin: " << Data[i].binN << " " << Data[i].dcq[0] << " " << Data[i].h[0] << " " << Data[i].h_std << " " << Data[i].dvs[0] << " " << Data[i].dvs_std << " " << Data[i].drho[0] << " " << Data[i].drho_std << endl;

            GMT::psxy(outfile, vector<double> {Data[i].drho[0] - Data[i].drho_std, Data[i].drho[0] + Data[i].drho_std}, vector<double> {Data[i].dvs[0], Data[i].dvs[0]}, "-J -R -W0.234p,light" + color + ",. -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0], Data[i].drho[0]}, vector<double> {Data[i].dvs[0] - Data[i].dvs_std, Data[i].dvs[0] + Data[i].dvs_std}, "-J -R -W0.234p,light" + color + ",. -O -K");

            GMT::psxy(outfile, vector<double> {Data[i].drho[0]}, vector<double> {Data[i].dvs[0]}, "-J -R -Sc0.05i -W0.234p,black -O -K");

        }



    }

    for (size_t i = 0; i < Data.size(); ++i) {

        vector<GMT::Text> texts;
        texts.push_back(GMT::Text(Data[i].drho[0], Data[i].dvs[0], to_string(Data[i].binN), 6, "CM"));
        GMT::pstext(outfile, texts, "-J -R -N -O -K");
    }

    GMT::psbasemap(outfile, "-JX-" + to_string(plotWidth / 20 * 7.5) + "i/1i -R" + to_string(diamondVol_yy_drho.back()) + "/" + to_string(diamondVol_yy_drho[0]) + "/-1/1 -Bxa5 -Bx+lvol(\%) -BS -Y0.5i -O -K");


    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile, __FILE__);

    return 0;
}
