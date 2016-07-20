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

int main(int argc, char **argv)
{
    double s1, s1tplus1;
    double s2, s2tplus1;

    double p1 = atof(argv[1]);
    double n = atof(argv[2]);
    double d = atof(argv[3]);
    double sigma12 = atof(argv[4]);
    double sigma21 = atof(argv[5]);

    // eigenvector related thiings
    double ev, v1, v2, u1;
    double dWds1mother, dWds2mother, dWds1offspring, dWds2offspring;

    DataFile << "type;s1;s2;p1;d;n;c1;c2;sigma12;sigma21;ev;u1;v1;v2;" << endl;

    // first offspring strategies
    for (double c2 = 0; c2 <= 1.0; c2+=0.02)
    {
        for (double c1 = 0; c1 <= 1.0; c1+=0.02)
        {
            s1 = 0.5;
            s2 = 0.5;

            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                // first calculate eigenvectors et al 
                EIGENVECS

                // then evaluate selection gradients mother
                SELGRADS_MOM

                // update values
                s1tplus1 = s1 + 0.01 * dWds1mother;
                s2tplus1 = s2 + 0.01 * dWds2mother;

                if (s1tplus1 < 0)
                {
                    s1tplus1 = 0;
                } else if (s1tplus1 > 1.0)
                {
                    s1tplus1 = 1.0;
                }
                
                if (s2tplus1 < 0)
                {
                    s2tplus1 = 0;
                } else if (s2tplus1 > 1.0)
                {
                    s2tplus1 = 1.0;
                }

                if (
                        (fabs(s2tplus1 - s2) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08)
                        ||
                        iter == 10000000 - 1
                )
                {
                    s2 = s2tplus1;
                    s1 = s1tplus1;
                    DataFile << "mother;" << s1 << ";" << s2 << ";" << p1 << ";" << d << ";" << n << ";" << c1 << ";" << c2 << ";" << sigma12 << ";" << sigma21 << ";" << ev << ";" << u1 << ";" << v1 << ";" << v2 << ";" << endl;
                    break;
                }
                
                s1 = s1tplus1;
                s2 = s2tplus1;
            }

            s1 = 0.5;
            s2 = 0.5;

            
            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                // first calculate eigenvectors et al 
                EIGENVECS

                // then evaluate selection gradients offspring
                SELGRADS_OFF

                // update values
                s1tplus1 = s1 + 0.01 * dWds1offspring;
                s2tplus1 = s2 + 0.01 * dWds2offspring;

                if (s1tplus1 < 0)
                {
                    s1tplus1 = 0;
                } else if (s1tplus1 > 1.0)
                {
                    s1tplus1 = 1.0;
                }
                
                if (s2tplus1 < 0)
                {
                    s2tplus1 = 0;
                } else if (s2tplus1 > 1.0)
                {
                    s2tplus1 = 1.0;
                }

                if (
                        (fabs(s2tplus1 - s2) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08)
                        ||
                        iter == 10000000 - 1
                )
                {
                    s2 = s2tplus1;
                    s1 = s1tplus1;
                    DataFile << "offspring;" << s1 << ";" << s2 << ";" << p1 << ";" << d << ";" << n << ";" << c1 << ";" << c2 << ";" << sigma12 << ";" << sigma21 << ";" << ev << ";" << u1 << ";" << v1 << ";" << v2 << ";" << endl;
                    break;
                }
                
                s2 = s2tplus1;
                s1 = s1tplus1;
            }

        }
    }
}
