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

/*
 * Run this code on t039.Cx. - different source depth synthetics.
 *
 * The source is S stack -> modified according to ScS Stack (for individual events).
 *
 * This will generate a folder containing all cherry picked synthetic traces, deconed and frsed.
 * You can then swap in the synthetics as data to see how a 1D structure show up in the 1D modeling result.
*/

using namespace std;

// Inputs. ------------------------------------

const bool makePlots=false;

const string homeDir=GetHomeDir();
const string synDataDir=homeDir+"/PROJ/t039.reflRecoverStructure/t039.C0";
const string outputDIR=homeDir+"/PROJ/t041.reflRecoverStructure/C0/WaterFRS";
const string outputDIR2=homeDir+"/PROJ/t041.reflRecoverStructure/C0/WaterDecon";
const int beginIndex=1,endIndex=58;

const string infoTable="gen2CA_D.Master_a14";
const double filterCornerLow=0.033,filterCornerHigh=0.3;
const double cutDeconSourceT1=-60,cutDeconSourceT2=60;
const double cutDeconSignalT1=-120,cutDeconSiganlT2=120;
const double cutDeconResultT1=-100,cutDeconResultT2=100;
const double waterLevel=0.1,sigma=1.27398,deconFilterCornerLow=0.03,dt=0.025;

// --------------------------------------------

int main(int argc, char **argv){

    auto eqNames=MariaDB::Select("eq from "+infoTable+" group by eq order by eq");
    string outfile="";

    for (size_t Index=0; Index<endIndex-beginIndex+1; ++Index) {

        const string modelName=to_string(201500000000+beginIndex+Index);
        const string modelFolder=synDataDir+"/"+modelName;
        const string eqName=eqNames.GetString("eq")[beginIndex+Index-1];

        cout << "Processing synthetics: " << modelName << " -> " << eqName << " ..." << endl;

        /*************************************
        * 1. Read in waveform.               *
        *************************************/

        SACSignals Data(ShellExecVec("ls "+modelFolder+"/*.THT.sac"));

        Data.SortByGcarc();
        Data.Interpolate(dt);
        Data.RemoveTrend();
        Data.HannTaper(20);
        Data.Butterworth(filterCornerLow,filterCornerHigh);


        // find S peak and shift time reference to it.
        Data.FindPeakAround(Data.GetTravelTimes("S"),10);
        Data.ShiftTimeReferenceToPeak();
        Data.FlipPeakUp();
        Data.NormalizeToPeak();


        // Cherry pick: to find a set resemble the data.
        // Happens before ESW, decon source etc ..  or after them?
//         auto dataInfo=MariaDB::Select("stnm,gcarc from "+infoTable+" where eq="+eqName);
// 
//         const auto &dist=Data.GetDistances();
//         vector<size_t> selectTheseRecords;
//         for (auto item:dataInfo.GetDouble("gcarc")) {
//             auto it=lower_bound(dist.begin(),dist.end(),item);
//             if (it==dist.end()) it=prev(it);
//             if (it!=dist.begin() && fabs(*prev(it)-item)<fabs(*it-item)) it=prev(it);
//             selectTheseRecords.push_back(distance(dist.begin(),it));
//         }
//         Data=SACSignals(Data,selectTheseRecords);


        /********************************************
        * 2. Make decon source:                     *
        * stack S and ScS                           *
        * stretch/shrink S esw to look like ScS esw *
        ********************************************/

        SACSignals partialData;

        // stack s.

        partialData=Data;

        partialData.CheckAndCutToWindow(cutDeconSourceT1-10,cutDeconSourceT2+10);
        auto sESW=partialData.XCorrStack(0,-15,15,2).second.first;

        sESW.FindPeakAround(0,10);
        sESW.ShiftTimeReferenceToPeak();
        sESW.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
        sESW.FlipPeakUp();
        sESW.NormalizeToPeak();
        sESW.HannTaper(10);

        // stack scs

        partialData=Data;
        partialData.FindPeakAround(partialData.GetTravelTimes("ScS"),10);
        partialData.ShiftTimeReferenceToPeak();
        partialData.FlipPeakUp();
        partialData.NormalizeToPeak();


        partialData.CheckAndCutToWindow(cutDeconSourceT1-10,cutDeconSourceT2+10);
        auto scsESW=partialData.XCorrStack(0,-15,15,2).second.first;

        scsESW.FindPeakAround(0,10);
        scsESW.ShiftTimeReferenceToPeak();
        scsESW.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
        scsESW.FlipPeakUp();
        scsESW.NormalizeToPeak();
        scsESW.HannTaper(10);


        // Properly modify the S ESW to look like ScS ESW: make the decon source.
        EvenSampledSignal deconSource=sESW.StretchToFit(scsESW,-13,13,-0.1,0.2,0.25);

        deconSource.FindPeakAround(0,4);
        deconSource.ShiftTimeReferenceToPeak();
        deconSource.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
        deconSource.NormalizeToPeak();
        deconSource.FlipPeakUp();
        deconSource.HannTaper(10);


        // Cherry pick: to find a set resemble the data.
        // Happens here ... or before ESW, decon source etc.
        auto dataInfo=MariaDB::Select("stnm,gcarc from "+infoTable+" where eq="+eqName);

        const auto &dist=Data.GetDistances();
        vector<size_t> selectTheseRecords;
        for (auto item:dataInfo.GetDouble("gcarc")) {
            auto it=lower_bound(dist.begin(),dist.end(),item);
            if (it==dist.end()) it=prev(it);
            if (it!=dist.begin() && fabs(*prev(it)-item)<fabs(*it-item)) it=prev(it);
            selectTheseRecords.push_back(distance(dist.begin(),it));
        }
        Data=SACSignals(Data,selectTheseRecords);


        /*************
        * 3. Strip S *
        *************/

        // keep one trace to verify the S strip process.
        SACSignals beforeStrip, afterStrip;
        if (makePlots) beforeStrip=Data;

        // Strip S from waveform (choice 1: not the modified sESW)
//         vector<EvenSampledSignal> modifiedForStrip;
//         for (size_t i=0; i<Data.Size(); ++i) modifiedForStrip.push_back(sESW);


        // Strip S from waveform (choice 2: first stretch to best fit S)
        vector<EvenSampledSignal> modifiedForStrip;
        for (size_t i=0; i<Data.Size(); ++i) {
            //modifiedForStrip.push_back(sESW.StretchToFit(Data.GetData()[i],-13,13,-0.1,0.2,0.25));
            modifiedForStrip.push_back(sESW.StretchToFitHalfWidth(Data.GetData()[i]));
        }


        auto XCTimeShift=Data.CrossCorrelation(-15,15,modifiedForStrip,-15,15).first;
        Data.StripSignal(modifiedForStrip,XCTimeShift);



        if (makePlots) afterStrip=Data;


        /************
        * 4. Decon  *
        ************/


        // Find ScS from the stripped waveform.
        Data.FindPeakAround(Data.GetTravelTimes("ScS"),10);
        Data.ShiftTimeReferenceToPeak();
        Data.FlipPeakUp();
        Data.NormalizeToPeak();
        Data.CheckAndCutToWindow(cutDeconSignalT1,cutDeconSiganlT2);
        Data.HannTaper(20);


        // keep traces to verify the decon process.
        SACSignals beforeDecon;
        if (makePlots) beforeDecon=Data;


        // Decon.
        Data.WaterLevelDecon(deconSource,waterLevel);
        Data.GaussianBlur(sigma);
        Data.Butterworth(deconFilterCornerLow,10000);

        Data.FindPeakAround(0,5);
        Data.ShiftTimeReferenceToPeak();
        Data.FlipPeakUp();
        Data.NormalizeToPeak();
        Data.CheckAndCutToWindow(cutDeconResultT1,cutDeconResultT2);


        // Output Deconed waveforms.
        ShellExec("mkdir -p "+outputDIR2+"/"+eqName);
        for (size_t i=0;i<dataInfo.NRow();++i) {
            string outfileName=outputDIR2+"/"+eqName+"/"+dataInfo.GetString("stnm")[i]+".trace";
            ofstream fpout(outfileName);
            fpout << Data.GetData()[i];
            fpout.close();
        }


        /*********
        * 5. FRS *
        *********/

        SACSignals FRS=Data;

        FRS.FlipReverseSum(0);
        FRS.CheckAndCutToWindow(0,15-dt*0.8);

        // Measure the peak time and amplitude.
        FRS.FindPeakAround(7.5,7.5);

        // Output FRS waveforms.
        ShellExec("mkdir -p "+outputDIR);
        for (size_t i=0;i<dataInfo.NRow();++i) {
            string outfileName=outputDIR+"/"+eqName+"_"+dataInfo.GetString("stnm")[i]+".frs";
            ofstream fpout(outfileName);
            fpout << FRS.GetData()[i];
            fpout.close();
        }


        // Plot.
        if (!makePlots) continue;

        for (size_t i=0; i<dataInfo.NRow(); ++i) {
            string stnm=dataInfo.GetString("stnm")[i];
            double dist=dataInfo.GetDouble("gcarc")[i];
            double originalDist=Data.GetMData()[i].gcarc;

            if (i%17==0) { // A New page.
                if (outfile.empty()) outfile=GMT::BeginEasyPlot(69.5,40);
                else GMT::NewPage(outfile);
                GMT::MoveReferencePoint(outfile,"-Xf1i -Yf37.2i");
            }
            else GMT::MoveReferencePoint(outfile,"-Y-2.3i");

            // Plot to verify sESW.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf1i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,darkyellow -O -K");
            vector<GMT::Text> texts{GMT::Text(-30,0.9,eqName,12,"LT")};
            GMT::pstext(outfile,texts,"-J -R -N -O -K");

            // Plot to verify scsESW.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf14.5i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,scsESW,"-J -R -W1p,black -O -K");

            // Plot to verify the decon source.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf28i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,scsESW,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,yellow -O -K");
            GMT::psxy(outfile,deconSource,"-J -R -W1p,red -O -K");

            // Plot to verify the stripping.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf41.5i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,beforeStrip,i,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,modifiedForStrip[i],"-J -R -W1p,cyan -O -K");
            GMT::psxy(outfile,afterStrip,i,"-J -R -W1p,green -O -K");
            texts.clear();
            texts.push_back(GMT::Text(-30,0.9,stnm+"("+Float2String(dist,2)+" -> "+Float2String(originalDist,2)+")",12,"LT"));
            GMT::pstext(outfile,texts,"-J -R -N -O -K");

            // Plot to verify deconed trace.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf55i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,beforeDecon,i,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,deconSource,"-J -R -W1p,red -O -K");
            GMT::psxy(outfile,Data,i,"-J -R -W1p,purple -O -K");
        }
    }

    if (!makePlots) return 0;

    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile,__FILE__);

    return 0;
}
