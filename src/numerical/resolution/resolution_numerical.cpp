#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"


//#define NDEBUG


using namespace std;


// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("iter_conflict");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

double bound(double const val)
{
    return(val < 0 ? 0 : (val > 1.0 ? 1.0 : val));
}
double ent(double p)
{
    if (p == 0)
    {
        return(0);
    }
    if (p==1.0)
    {
        return(0);
    }

    return(- p * (log2(p)+(1-p) * log2(1-p)));
}
double siginfo(double p1, double p2)
{
    return(1 - (((p1 + p2)/2) * ent(p1/ ( p1 + p2)) + (1 - ((p1 + p2)/2)) * ent((1-p1)/ ((1-p1)+(1-p2)))));
}

int main(int argc, char **argv)
{
    double s1, s1tplus1;
    double s2, s2tplus1;
    double qNS, qNStplus1;
    double qS, qStplus1;

    double p1 = atof(argv[1]);
    double n = atof(argv[2]);
    double d = atof(argv[3]);
    double sigma12 = atof(argv[4]);
    double sigma21 = atof(argv[5]);

    // eigenvector related thiings
    double ev, v1, v2, u1;
    double dWds1, dWds2, dWdqS, dWdqNS;

    DataFile << "s1;s2;qS;qNS;p1;d;n;c1;c2;sigma12;sigma21;ev;u1;v1;v2;siginfo;" << endl;

    // first offspring strategies
    for (double c2 = 0; c2 < 1.02; c2+=0.02)
    {
        for (double c1 = 0; c1 < 1.02; c1+=0.02)
        {
            s1 = 0.1;
            s2 = 0.9;
            qS = 0.9;
            qNS = 0.1;

            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                // first calculate eigenvectors et al 
                EIGENVECS

                // then evaluate selection gradients mother
                SELGRADS

                // update values
                s1tplus1 = s1 + 0.01 * dWds1;
                s2tplus1 = s2 + 0.01 * dWds2;
                qNStplus1 = qNS + 0.01 * dWdqNS;
                qStplus1 = qS + 0.01 * dWdqS;

                s1tplus1 = bound(s1tplus1);
                s2tplus1 = bound(s2tplus1);
                qStplus1 = bound(qStplus1);
                qNStplus1 = bound(qNStplus1);

                if (
                        (
                        fabs(s2tplus1 - s2) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08 &&
                        fabs(qNStplus1 - qNS) < 1e-08 &&
                        fabs(qStplus1 - qS) < 1e-08
                       )
                        ||
                        iter == 10000000 - 1
                )
                {
                    s2 = s2tplus1;
                    s1 = s1tplus1;
                    qNS = qNStplus1;
                    qS = qStplus1;
                    DataFile <<  s1 << ";" << s2 << ";" << qS << ";" << qNS << ";" << p1 << ";" << d << ";" << n << ";" << c1 << ";" << c2 << ";" << sigma12 << ";" << sigma21 << ";" << ev << ";" << u1 << ";" << v1 << ";" << v2 << ";" << siginfo(s1,s2) << ";" << endl;
                    break;
                }
                
                s2 = s2tplus1;
                s1 = s1tplus1;
                qNS = qNStplus1;
                qS = qStplus1;
            }
        }
    }
}
