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

/*
 * Run this code on t052.Cx.
 * This will generate a folder containing all cherry picked synthetic traces, deconed and frsed.
 * You can then swap in the synthetics as data to see how a 2D structure show up in the 1D modeling result.
*/

using namespace std;

// Inputs. ------------------------------------

const bool makePlots=false;

const string sourceDataDir="/home/shule/PROJ/t052.PREMX.PREM/201500000000"; // Use PREMX @500km S or ScS stack as the source?
const string synDataDir="/home/shule/PROJ/t052.shaxiRecoverStructure/t052.C7.ULVZ";
const string outputDIR="/home/shule/PROJ/t041.SHAXI_C7/WaterFRS";
const string outputDIR2="/home/shule/PROJ/t041.SHAXI_C7/WaterDecon";
const int beginIndex=1,endIndex=58;

const string infoTable="gen2CA_D.Master_a14";
const double filterCornerLow=0.033,filterCornerHigh=0.3;
const double cutDeconSourceT1=-60,cutDeconSourceT2=60;
const double cutDeconSignalT1=-120,cutDeconSiganlT2=120;
const double cutDeconResultT1=-100,cutDeconResultT2=100;
// const double waterLevel=0.1,sigma=1.69864,deconFilterCornerLow=0.03,dt=0.025; // width=4.
const double waterLevel=0.1,sigma=1.35891,deconFilterCornerLow=0.03,dt=0.025; // width=3.2.
// const double waterLevel=0.1,sigma=1.31645,deconFilterCornerLow=0.03,dt=0.025; // width=3.1.
//const double waterLevel=0.1,sigma=1.27398,deconFilterCornerLow=0.03,dt=0.025; // width=3.

// --------------------------------------------

int main(int argc, char **argv){

    /*****************************************
    * Make decon source.
    * Read in PREMX.
    * Stack to make sESW, scsESW.
    * Modify sESW to look like scsESW.
    * Use the modified sESW as decon source?
    * Or use the scsESW (50~80 deg) as the decond source?
    *****************************************/

    // read in prem data.
    SACSignals premData(ShellExecVec("ls "+sourceDataDir+"/*.THT.sac"));

    premData.SortByGcarc();
    premData.Interpolate(dt);
    premData.RemoveTrend();
    premData.HannTaper(20);
    premData.Butterworth(filterCornerLow,filterCornerHigh);

    // find S peak and shift time reference to the peak.
    premData.FindPeakAround(premData.GetTravelTimes("S"),10);
    premData.ShiftTimeReferenceToPeak();
    premData.FlipPeakUp();
    premData.NormalizeToPeak();


    // make sESW.
    SACSignals partialData;

    partialData=premData;
    //partialData.CheckDist(40,50);// It is observed in SHAXI: ScS is uneven, the most approximate S waveform (to ScS) distance is between 40~50.

    partialData.CheckAndCutToWindow(cutDeconSourceT1-10,cutDeconSourceT2+10);
    auto sESW=partialData.XCorrStack(0,-15,15,2).second.first;

    sESW.FindPeakAround(0,10);
    sESW.ShiftTimeReferenceToPeak();
    sESW.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
    sESW.FlipPeakUp();
    sESW.NormalizeToPeak();
    sESW.HannTaper(10);


    // make scsESW.

    partialData=premData;
    partialData.CheckDist(50,80); // sS is coming in from 40 to 50.


    // strip S and sS to make scsESW, if use scsESW as the decon source.

    auto XCTimeShift=partialData.CrossCorrelation(-15,15,sESW,-15,15).first;
    partialData.StripSignal(sESW,XCTimeShift);

    partialData.FindPeakAround(partialData.GetTravelTimes("sS"),10);
    partialData.ShiftTimeReferenceToPeak();
    partialData.FlipPeakUp();
    partialData.NormalizeToPeak();

    XCTimeShift=partialData.CrossCorrelation(-15,15,sESW,-15,15).first;
    partialData.StripSignal(sESW,XCTimeShift);


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

    //EvenSampledSignal deconSource=sESW.StretchToFit(scsESW,-13,13,-0.1,0.2,0.25); // use modified sESW as decon source.
    EvenSampledSignal deconSource=scsESW; // use scs ESW as decon source.

    deconSource.FindPeakAround(0,4);
    deconSource.ShiftTimeReferenceToPeak();
    deconSource.CheckAndCutToWindow(cutDeconSourceT1,cutDeconSourceT2);
    deconSource.NormalizeToPeak();
    deconSource.FlipPeakUp();
    deconSource.HannTaper(10);



    /*****************************************************************************
    * For each synthetic model...
    * 1. Read in synthetic ScS waveform. Preprocess it (Band-pass, cherry-pick)
    * 2. Strip S.
    * 3. Decon.
    * 4. FRS.
    * 5. UpdateTables.
    *****************************************************************************/

    auto eqNames=MariaDB::Select("eq from "+infoTable+" group by eq order by eq");

    vector<string> outfiles;
    size_t maskCnt=0;

    for (size_t Index=0; Index<endIndex-beginIndex+1; ++Index) {

        const string modelName=to_string(201500000000+beginIndex+Index);
        const string modelFolder=synDataDir+"/"+modelName;
        const string eqName=eqNames.GetString("eq")[beginIndex+Index-1];

        cout << "Processing synthetics: " << modelName << " -> " << eqName << " ..." << endl;

        /*************************************
        * 1. Read in synthetic ScS waveform. *
        *    Select the right distance set.  *
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
        auto dataInfo=MariaDB::Select("stnm,gcarc from "+infoTable+" where eq="+eqName);

        const auto &dist=Data.GetDistances();
        set<size_t> selectTheseRecords;
        for (auto item:dataInfo.GetDouble("gcarc")) {
            auto it=lower_bound(dist.begin(),dist.end(),item);
            if (it==dist.end()) it=prev(it);
            if (it!=dist.begin() && fabs(*prev(it)-item)<fabs(*it-item)) it=prev(it);
            selectTheseRecords.insert(distance(dist.begin(),it));
        }
        Data=SACSignals(Data,selectTheseRecords);
        Data.SortByGcarc();

        /*************
        * 2. Strip S *
        *************/

        // keep one trace to verify the S strip process.
        SACSignals beforeStrip, afterStrip;
        if (makePlots) beforeStrip=Data;


        // Strip S from waveform (notice: not the modified sESW)
        auto XCTimeShift=Data.CrossCorrelation(-15,15,sESW,-15,15).first;
        Data.StripSignal(sESW,XCTimeShift);

        // optional, strip sS from SHAXI because sS on SHAXI is huge.
//         Data.FindPeakAround(Data.GetTravelTimes("sS"),10);
//         Data.ShiftTimeReferenceToPeak();
//         XCTimeShift=Data.CrossCorrelation(-15,15,sESW,-15,15).first;
//         Data.StripSignal(sESW,XCTimeShift);
//         Data.FindPeakAround(Data.GetTravelTimes("S"),10);
//         Data.ShiftTimeReferenceToPeak();


        if (makePlots) afterStrip=Data;


        /***********
        * 3. Decon *
        ***********/

        // Find ScS from the stripped waveform.
        Data.FindPeakAround(Data.GetTravelTimes("ScS"),10);
        Data.ShiftTimeReferenceToPeak();
        Data.FlipPeakUp();
        Data.NormalizeToPeak();
        Data.CheckAndCutToWindow(cutDeconSignalT1,cutDeconSiganlT2);


        // keep one trace to verify the decon process.
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



        // Output Deconed waveforms. Cherry-pick happens here.
        Data.SortByGcarc();
        ShellExec("mkdir -p "+outputDIR2+"/"+eqName);
        for (size_t i=0;i<dataInfo.NRow();++i) {
            string outfileName=outputDIR2+"/"+eqName+"/"+dataInfo.GetString("stnm")[i]+".trace";
            double dist=dataInfo.GetDouble("gcarc")[i];
            size_t j=Data.FindByGcarc(dist)[0];

            // Mask SHAXI sS, SS and other traffic phases. (Their side lobes have large amplitude!)
            // Or strip sS before decon?
            if ( fabs(Data.GetTravelTimes("ScS")[j]-Data.GetTravelTimes("sS")[j])<30 ||
                 fabs(Data.GetTravelTimes("ScS")[j]-Data.GetTravelTimes("S")[j])<30  ||
                 fabs(Data.GetTravelTimes("ScS")[j]-Data.GetTravelTimes("SS")[j])<30 ){
                ++maskCnt;
                Data.Mask(-100000,100000,j);
            }
                
            ofstream fpout(outfileName);
            fpout << Data.GetData()[j];
            fpout.close();
        }


        /*********
        * 4. FRS *
        *********/

        SACSignals FRS=Data;

        FRS.FlipReverseSum(0);
        FRS.CheckAndCutToWindow(0,15-dt*0.8);

        // Measure the peak time and amplitude.
        FRS.FindPeakAround(7.5,7.5);
        FRS.SortByGcarc();

        // Output FRS waveforms.

        ShellExec("mkdir -p "+outputDIR);
        for (size_t i=0;i<dataInfo.NRow();++i) {
            string outfileName=outputDIR+"/"+eqName+"_"+dataInfo.GetString("stnm")[i]+".frs";
            double dist=dataInfo.GetDouble("gcarc")[i];
            ofstream fpout(outfileName);
            fpout << FRS.GetData()[FRS.FindByGcarc(dist)[0]];
            fpout.close();
        }

        // Plot.
        if (!makePlots) continue;

        string outfile;

        for (size_t i=0; i<dataInfo.NRow(); ++i) {
            string stnm=dataInfo.GetString("stnm")[i];
            double dist=dataInfo.GetDouble("gcarc")[i];

            size_t j=Data.FindByGcarc(dist)[0];
            double originalDist=Data.GetMData()[j].gcarc;

            if (i%17==0) { // A New page.
                outfile=GMT::BeginEasyPlot(69.5,40);
                outfiles.push_back(outfile);
                GMT::MoveReferencePoint(outfile,"-Xf1i -Yf37.2i");
            }
            else {
                GMT::MoveReferencePoint(outfile,"-Y-2.3i");
            }

            // Plot to verify sESW.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf1i");
            GMT::psxy(outfile,vector<double> {-100,100},vector<double> {0,0},"-J -R -W1p,gray,- -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,yellow -O -K");

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
            GMT::psxy(outfile,beforeStrip,j,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,sESW,"-J -R -W1p,yellow -O -K");
            GMT::psxy(outfile,afterStrip,j,"-J -R -W1p,green -O -K");
            vector<GMT::Text> texts{GMT::Text(-30,0.9,Float2String(originalDist,2)+" -> "+stnm+"("+Float2String(dist,2)+")",12,"LT")};
            GMT::pstext(outfile,texts,"-J -R -N -O -K");

            // Plot to verify deconed trace.
            GMT::psbasemap(outfile,"-JX13i/2i -R-100/100/-1/1 -Bxa10 -Bya0.5 -BWSne -O -K -Xf55i");
            GMT::psxy(outfile,beforeDecon,j,"-J -R -W1p,black -O -K");
            GMT::psxy(outfile,deconSource,"-J -R -W1p,red -O -K");
            GMT::psxy(outfile,Data,j,"-J -R -W1p,purple -O -K");
        }

    }
    cout << "In total, " << maskCnt << " traces are mask due to traffic phases (sS, SS, S)" << endl;

    if (!makePlots) return 0;

    for (auto file: outfiles) {
        GMT::SealPlot(file);
        ShellExec("cat "+file+" >> tmp.ps");
        remove(file.c_str());
    }
    string pdffile=__FILE__;
    ShellExec("ps2pdf tmp.ps "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
    remove("tmp.ps");

    return 0;
}
