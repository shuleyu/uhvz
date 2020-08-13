#include <iostream>
#include <fstream>
#include <vector>

#include <GMT.hpp>

using namespace std;

int main(){

    vector<double> drho, dvs;
    double rho, vs;
    string s;
    ifstream fpin("Dan.txt");

    vector<GMT::Text> texts;
    while (fpin >> s >> rho >> vs){

        texts.push_back(GMT::Text(rho + 0.3, vs + 1, s, 8, "LB"));
        drho.push_back(rho);
        dvs.push_back(vs);
    }
    fpin.close();

    string outfile = GMT::BeginEasyPlot(10, 10);

    GMT::MoveReferencePoint(outfile,"-Xf1i -Yf1i");

    GMT::psbasemap(outfile, "-JX5.4i/4i -R-10/12.5/-30.5/20.5 -Bxa2.5g2.5 -Bya5g5 -Bx+ldrho(\%) -By+ldvs(%) -BWSne -O -K");

    GMT::psxy(outfile, vector<double>{0, 0}, vector<double>{-30.5, 20.5}, "-J -R -W0.5p,black -O -K");

    GMT::psxy(outfile, vector<double>{-10, 12.5}, vector<double>{0, 0}, "-J -R -W0.5p,black -O -K");

    GMT::psxy(outfile, drho, dvs, "-J -R -Sc0.1i -Gblack -O -K");

    GMT::pstext(outfile, texts, "-J -R -N -O -K");

    // Make PDF.
    GMT::SealPlot(outfile);
    GMT::ps2pdf(outfile, __FILE__);

    return 0;
}
