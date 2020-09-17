#include<iostream>
#include<string>
#include<map>

#include<MariaDB.hpp>
#include<GMT.hpp>
#include<ShellExec.hpp>

using namespace std;

// Inputs . --------------

const string infoTable1="gen2CA_D.Master_a14", infoTable2="gen2CA_D.Bins", infoTable3="gen2CA_D.CenterDists";
const string phase1="ScS", phase2="sS"; // Will calcualte dT = phase2 - phase1

// -----------------------

int main(){


    auto recordDT=MariaDB::Select(phase2 + "-" + phase1 +" as dT ,pairname as pn from "+infoTable1); 

    map<string,double> mapDT;
    for (size_t i=0; i<recordDT.NRow(); ++i)
        mapDT[recordDT.GetString("pn")[i]]=recordDT.GetDouble("dT")[i];


    // For each bin, get the pairname, then get the 

    vector<vector<double>> plotData;
    auto binInfo=MariaDB::Select("bin, lon, lat from "+infoTable2); 

    for (size_t i=0; i<binInfo.NRow(); ++i) {
        string binN=to_string(binInfo.GetInt("bin")[i]);
        auto dataInfo=MariaDB::Select("pairname as pn from "+infoTable3+" where dist_"+binN+" >-0.5"); 

        double val=0;
        for (size_t j=0; j<dataInfo.NRow(); ++j) {
            double x=mapDT[dataInfo.GetString("pn")[j]];
            if(x<30) {
                val=1;
                break;
            }
        }
        plotData.push_back({binInfo.GetDouble("lon")[i], binInfo.GetDouble("lat")[i], val, 0.2});
    }


    double XSIZE=10,scale=(XSIZE-2)/35;
    string outfile=GMT::BeginEasyPlot(XSIZE, XSIZE);
    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf1.5i");
    GMT::psbasemap(outfile,"-Jx"+to_string(scale)+"id/"+to_string(scale)+"id -R-99/-64/-5/27 -Bxa10f5 -Bya10f5 -Bx+lLongitude -By+lLatitude -BWSne -O -K");
    GMT::pscoast(outfile,"-J -R -W0.5p,black -O -K");
    GMT::makecpt("-Cgray -T0.0/1.0 -I > tmp.cpt");
    GMT::psxy(outfile,plotData,"-J -R -Sc -Ctmp.cpt -W0p -O -K");
    GMT::psscale(outfile,"-Ctmp.cpt -D"+to_string(scale*35/2.0)+"i/-0.5i/2i/0.1ih -O -K -Bxa1 -Bx+lHas("+phase2+"-"+phase1+"<30sec)");

    // Make PDF.
    GMT::SealPlot(outfile);
    string pdffile=__FILE__;
    ShellExec("ps2pdf "+outfile+" "+pdffile.substr(0,pdffile.find_last_of("."))+".pdf");
    remove("tmp.cpt");
    remove(outfile.c_str());


    return 0;
}
