#include<iostream>
#include<vector>
#include<algorithm>
#include<set>

#include<ShellExecVec.hpp>
#include<ShellExec.hpp>
#include<EvenSampledSignal.hpp>
#include<SACSignals.hpp>
#include<GMTPlotSignal.hpp>
#include<MariaDB.hpp>
#include<Float2String.hpp>
#include<GetHomeDir.hpp>

/**********************************************************************************
 *
 * Run this code on t041.DATA -- Use subtraction instead of deconvolution.
 *
 * This will generate S_ESW-subtraced waveforms of S and ScS center on peaks.
 *
 * Time window is -100 ~ 100 second.
 *
 *********************************************************************************/

using namespace std;

// Inputs. ------------------------------------

const bool makePlots=true;

const string homeDir=GetHomeDir();

const int beginIndex=1, endIndex=1;
const string infoTable="gen2CA_D.Master_a14";
const double filterCornerLow=0.033, filterCornerHigh=0.3, dt=0.025;
const double cutSourceT1=-100, cutSourceT2=100;
const double cutResultT1=-100, cutResultT2=100;

// Outputs. ------------------------------------

const string dirPrefix=homeDir+"/PROJ/t041.CA_D/Subtract";
const string outputDB="gen2CA_D", outputTable="Subtract";

// --------------------------------------------

int main(int argc, char **argv){

    auto eqNames=MariaDB::Select("eq from "+infoTable+" group by eq order by eq");
    string outfile="";
    vector<vector<string>> sqlData(6);

    for (size_t Index=0; Index<endIndex-beginIndex+1; ++Index) {

        const string eqName=eqNames.GetString("eq")[beginIndex+Index-1];

        cout << "Processing data (subtraction): " << eqName << " ..." << endl;

        /****************************
         *
         * 1. Read in waveform.
         *
        ****************************/

        auto dataInfo=MariaDB::Select("pairname as pn, concat(dirPrefix,'/',File) as file, Peak_S, Peak_ScS, stnm from "+infoTable+" where eq="+eqName);

        SACSignals Data(dataInfo.GetString("file"));

        Data.Interpolate(dt);
        Data.RemoveTrend();
        Data.HannTaper(20);
        Data.Butterworth(filterCornerLow,filterCornerHigh);



        // find S peak and shift time reference to it.
        Data.ShiftTime(Data.GetTravelTimes("S"));
        Data.FindPeakAround(dataInfo.GetDouble("Peak_S"),2);
        Data.ShiftTimeReferenceToPeak();
        Data.FlipPeakUp();
        Data.NormalizeToPeak();



        /********************************************
         *
         * 2. Make S ESW stack.
         *
        ********************************************/

        SACSignals sESWData=Data;
        sESWData.CheckAndCutToWindow(cutSourceT1-10,cutSourceT2+10);

        auto sESW=sESWData.XCorrStack(0,-15,15,2).second.first;
        sESW.FindPeakAround(0,10);
        sESW.ShiftTimeReferenceToPeak();
        sESW.CheckAndCutToWindow(cutSourceT1,cutSourceT2);
        sESW.FlipPeakUp();
        sESW.NormalizeToPeak();
        sESW.HannTaper(20);


        // Optional: make S ESW again, this time, stretch/shrink each S to match S ESW, then stack.


        sESWData.StretchToFit(sESW,-13,13,-0.3,0.3,0.25,true); // stretch, then compare waveform Amp_WinDiff.

        sESW=sESWData.XCorrStack(0,-15,15,2).second.first;
        sESW.FindPeakAround(0,10);
        sESW.ShiftTimeReferenceToPeak();
        sESW.CheckAndCutToWindow(cutSourceT1,cutSourceT2);
        sESW.FlipPeakUp();
        sESW.NormalizeToPeak();
        sESW.HannTaper(20);


        /************************************************
         *
         * 3. Modifiy S ESW to match S waveforms.
         *    Strip modified S ESW from each S waveform.
         *
        ************************************************/

        SACSignals beforeSStrip=Data;

        // Properly modify the S_ESW to look like S.
        vector<EvenSampledSignal> modifiedToFitS;
        for (size_t i=0; i<dataInfo.NRow(); ++i) {
            modifiedToFitS.push_back(sESW.StretchToFit(beforeSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true)); // stretch, then compare waveform Amp_WinDiff.
            //modifiedToFitS.push_back(sESW.StretchToFit(beforeSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true,1)); // stretch, then compare waveform Amp_Diff.
            //modifiedToFitS.push_back(sESW.StretchToFitHalfWidth(beforeSStrip.GetData()[i])); // stretch to fit the half-width.
        }
        auto SXCTimeShift=beforeSStrip.CrossCorrelation(-15,15,modifiedToFitS,-15,15).first;
        auto afterSStrip=beforeSStrip;
        afterSStrip.StripSignal(modifiedToFitS,SXCTimeShift);

        beforeSStrip.CheckAndCutToWindow(cutResultT1,cutResultT2);
        afterSStrip.CheckAndCutToWindow(cutResultT1,cutResultT2);



        /************************************************
         *
         * 4. Modifiy S ESW to match ScS waveforms.
         *    Strip modified S ESW from each ScS waveform.
         *
        ************************************************/

        SACSignals beforeScSStrip=Data;


        beforeScSStrip.ShiftTime(beforeScSStrip.GetTravelTimes("ScS"));
        beforeScSStrip.FindPeakAround(dataInfo.GetDouble("Peak_ScS"),2);
        beforeScSStrip.ShiftTimeReferenceToPeak();
        beforeScSStrip.FlipPeakUp();
        beforeScSStrip.NormalizeToPeak();


        // Properly modify the S ESW to look like ScS.
        vector<EvenSampledSignal> modifiedToFitScS;
        for (size_t i=0; i<dataInfo.NRow(); ++i) {
            modifiedToFitScS.push_back(sESW.StretchToFit(beforeScSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true)); // stretch, then compare waveform Amp_WinDiff.
            //modifiedToFitScS.push_back(sESW.StretchToFit(beforeScSStrip.GetData()[i],-13,13,-0.3,0.3,0.25,true,1)); // stretch, then compare waveform Amp_Diff.
            //modifiedToFitScS.push_back(sESW.StretchToFitHalfWidth(beforeScSStrip.GetData()[i])); // stretch to fit the half-width.
        }
        auto ScSXCTimeShift=beforeScSStrip.CrossCorrelation(-10,10,modifiedToFitScS,-10,10).first;
        auto afterScSStrip=beforeScSStrip;
        afterScSStrip.StripSignal(modifiedToFitScS,ScSXCTimeShift);

        beforeScSStrip.CheckAndCutToWindow(cutResultT1,cutResultT2);
        afterScSStrip.CheckAndCutToWindow(cutResultT1,cutResultT2);



        // Check peak finding error.


        for (size_t i=0;i<dataInfo.NRow();++i) {

            size_t index=beforeSStrip.GetData()[i].GetPeak();
            if (beforeSStrip.GetData()[i].GetAmp()[index-1]>beforeSStrip.GetData()[i].GetAmp()[index] ||
                beforeSStrip.GetData()[i].GetAmp()[index+1]>beforeSStrip.GetData()[i].GetAmp()[index] )
                cout << "S Peak finding error for: " << beforeSStrip.GetData()[i].GetFileName() << endl;

            index=beforeScSStrip.GetData()[i].GetPeak();
            if (beforeScSStrip.GetData()[i].GetAmp()[index-1]>beforeScSStrip.GetData()[i].GetAmp()[index] ||
                beforeScSStrip.GetData()[i].GetAmp()[index+1]>beforeScSStrip.GetData()[i].GetAmp()[index] )
                cout << "S Peak finding error for: " << beforeScSStrip.GetData()[i].GetFileName() << endl;

            sqlData[0].push_back(dataInfo.GetString("pn")[i]);
            sqlData[1].push_back(to_string(-beforeSStrip.GetTravelTimes("S",{i})[0]));
            sqlData[2].push_back(to_string(-beforeScSStrip.GetTravelTimes("ScS",{i})[0]));
        }

        /******************
         *
         * 5. Output.
         *
        ******************/

        // Output Stripped waveforms.
        ShellExec("mkdir -p "+dirPrefix+"/"+eqName);
        for (size_t i=0;i<dataInfo.NRow();++i) {
            sqlData[3].push_back(eqName+"/"+dataInfo.GetString("stnm")[i]+".SStripped");
            sqlData[4].push_back(eqName+"/"+dataInfo.GetString("stnm")[i]+".ScSStripped");
            sqlData[5].push_back(dirPrefix);
            afterSStrip.GetData()[i].OutputToFile(dirPrefix+"/"+sqlData[3].back());
            afterScSStrip.GetData()[i].OutputToFile(dirPrefix+"/"+sqlData[4].back());
        }

        // Plot.
        if (!makePlots) continue;

        for (size_t i=0; i<dataInfo.NRow(); ++i) {
            string stnm=dataInfo.GetString("stnm")[i];
            double dist=afterSStrip.GetMData()[i].gcarc;

            if (i%17==0) { // A New page.
                if (outfile.empty()) outfile=GMT::BeginEasyPlot(69.5,40);
                else GMT::NewPage(outfile);
                GMT::MoveReferencePoint(outfile,"-Xf1i -Yf37.2i");
            }
            else GMT::MoveReferencePoint(outfile,"-Y-2.3i");

            // Plot to verify sESW.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1.05 -Bxa10 -Bya0.5 -BWSne -O -K -Xf1i");
            GMT::psxy(outfile,vector<double> {0,0},vector<double> {-2,2},"-J -R -W0.5p,red,- -O -K");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,darkyellow -O -K");

            vector<GMT::Text> texts{GMT::Text(-30,0.9,to_string(beginIndex+Index)+". "+eqName,12,"LT")};
            GMT::pstext(outfile,texts,"-J -R -N -O -K");


            // Plot to verify S Strip.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1.05 -Bxa10 -Bya0.5 -BWSne -O -K -Xf14.5i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,vector<double> {beforeSStrip.GetTravelTimes("S",{i})[0]},vector<double> {0},"-J -R -Sy0.05i -W1p,red -O -K");
            // GMT::psxy(outfile,sESWData,i,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,beforeSStrip,i,"-J -R -W1p,black -O -K");
            modifiedToFitS[i].ShiftTime(SXCTimeShift[i]);
            GMT::psxy(outfile,modifiedToFitS[i],"-J -R -W1p,cyan -O -K");
            GMT::psxy(outfile,afterSStrip,i,"-J -R -W1p,green -O -K");
            GMT::psxy(outfile,vector<double> {0},vector<double> {1},"-J -R -Sc0.05i -Gblue -W0p -O -K");

            texts.clear();
            texts.push_back(GMT::Text(-30,0.9,stnm+"("+Float2String(dist,2)+")",12,"LT"));
            GMT::pstext(outfile,texts,"-J -R -N -O -K");


            // Plot S Strip result.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1.05 -Bxa10 -Bya0.5 -BWSne -O -K -Xf28i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,afterSStrip,i,"-J -R -W1p,black -O -K");


            // Plot to verify ScS Strip.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1.05 -Bxa10 -Bya0.5 -BWSne -O -K -Xf41.5i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,vector<double> {beforeScSStrip.GetTravelTimes("ScS",{i})[0]},vector<double> {0},"-J -R -Sy0.05i -W1p,red -O -K");
            GMT::psxy(outfile,beforeScSStrip,i,"-J -R -W1p,black -O -K");
            modifiedToFitScS[i].ShiftTime(ScSXCTimeShift[i]);
            GMT::psxy(outfile,modifiedToFitScS[i],"-J -R -W1p,cyan -O -K");
            GMT::psxy(outfile,afterScSStrip,i,"-J -R -W1p,green -O -K");
            GMT::psxy(outfile,vector<double> {0},vector<double> {1},"-J -R -Sc0.05i -Gblue -W0p -O -K");


            // Plot ScS Strip result.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1.05 -Bxa10 -Bya0.5 -BWSne -O -K -Xf55i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,afterScSStrip,i,"-J -R -W1p,black -O -K");
        }
    }

//     MariaDB::Query("create database if not exists "+outputDB);
//     MariaDB::Query("drop table if exists "+outputDB+"."+outputTable);
//     MariaDB::Query("create table "+outputDB+"."+outputTable+" (PairName varchar(30) not null unique primary key, Peak_S double, Peak_ScS double, SStripped varchar(200), ScSStripped varchar(200), dirPrefix varchar(200))");
//     MariaDB::LoadData(outputDB,outputTable,vector<string> {"pairname","Peak_S","Peak_ScS","SStripped", "ScSStripped", "dirPrefix"}, sqlData);

    if (!makePlots) return 0;

    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile,__FILE__);

    return 0;
}
