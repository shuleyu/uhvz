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

#include<SACSignals.hpp>
#include<MariaDB.hpp>
#include<GMT.hpp>
#include<AvrStd.hpp>
#include<WayPoint.hpp>
#include<ShellExec.hpp>

using namespace std;

int main(){


    // data.
    auto HitLocation=MariaDB::Select("hitlo, hitla from gen2CA_D.Master_a14 where shift_gcarc <= 70");

    auto BinInfo=MariaDB::Select("bin, lon, lat, radius from gen2CA_D.Bins");


    // plot.
    string outfile = GMT::BeginEasyPlot(12, 12);

    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf1i");
    GMT::psbasemap(outfile,"-Jx0.15id -R-101/-64/-5/27 -Bxa5 -Bya5 -Bx+lLongitude -By+lLatitude -BWSne -O -K");
    GMT::pscoast(outfile,"-J -R -A10000 -W0.213p,black -O -K");

    // Plot hit points.
    GMT::psxy(outfile, HitLocation.GetDouble("hitlo"), HitLocation.GetDouble("hitla"), "-J -R -W0 -Sc0.01i -Ggray -O -K");

    // Plot bins.
    int ret = 0;
    vector<GMT::Text> texts;
    for (size_t i=0;i<BinInfo.NRow();++i) {

        // more data.
        auto Cnt = MariaDB::Select("count(*) as cnt from gen2CA_D.CenterDists as A join gen2CA_D.Master_a14 as B on A.pairname = B.pairname where A.dist_" + to_string(BinInfo.GetInt("bin")[i]) + " > -0.5 and B.shift_gcarc <=70");

        // Plot a circle.
        vector<double> c_lon,c_lat;
        for (double az=0;az<=360;az+=0.1) {
            c_lon.push_back(0);
            c_lat.push_back(0);
            tie(c_lon.back(),c_lat.back())=WayPoint(BinInfo.GetDouble("lon")[i],BinInfo.GetDouble("lat")[i],az,BinInfo.GetDouble("radius")[i]);
        }
        GMT::psxy(outfile,c_lon, c_lat, "-J -R -W0.12p,black -O -K");

        // Add texts.
        texts.push_back(GMT::Text(BinInfo.GetDouble("lon")[i],BinInfo.GetDouble("lat")[i],to_string(Cnt.GetInt("cnt")[0]),8,"CM"));

        ret += Cnt.GetInt("cnt")[0];
    }

    cout << ret / BinInfo.NRow() << endl;

    // Plot number of traces in each bin.
    GMT::pstext(outfile,texts,"-J -R -N -O -K");

    GMT::ps2pdf(outfile, __FILE__);

    return 0;

}
