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
    double q1, q1tplus1;
    double q2, q2tplus1;

    double p1 = atof(argv[1]);
    double p2 = 1-p1;
    double c = atof(argv[2]);
    double n = atof(argv[3]);

    DataFile << "q1;q2;p1;d;n;c;" << endl;

        for (double d = 0; d <= 1.0; d+=0.02)
        {
            q1 = 0.5;
            q2 = 0.5;

            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                q1tplus1 = 0.01 * 
                    (
                     c*(-(((-1 + d)*p(1))/
                (1 + d*(-1 + p(1) + p(2)) + 
                  c*(-1 + q(1) + d*
                      (1 + p(1)*(-1 + q(1)) - q(1) + p(2)*(-1 + q(2))))))
               - (Power(-1 + d,2)*p(1)*(p(1) + p(2))*
                (1 + c*(-1 + q(1))))/
              (n*Power(1 + d*(-1 + p(1) + p(2)) + 
                  c*(-1 + q(1) + d*
                      (1 + p(1)*(-1 + q(1)) - q(1) + p(2)*(-1 + q(2)))),
                 2)*(1 - (Power(-1 + d,2)*(-1 + n)*p(1)*
                     Power(1 + c*(-1 + q(1)),2))/
                   (n*Power(1 + d*(-1 + p(1) + p(2)) + 
                       c*(-1 + q(1) + 
                          d*(1 + p(1)*(-1 + q(1)) - q(1) + 
                             p(2)*(-1 + q(2)))),2)) - 
                  (Power(-1 + d,2)*(-1 + n)*p(2)*Power(-1 + c*q(2),2))/
                   (n*Power(-1 + c*q(2) + 
                       d*(p(1)*(-1 + c*q(1)) + (-1 + p(2))*(-1 + c*q(2)))
                       ,2)))) + d*(p(1)/
                 (1 + d*(-1 + p(1) + p(2)) + 
                   c*(-1 + q(1) + d*
                       (1 + p(1)*(-1 + q(1)) - q(1) + p(2)*(-1 + q(2)))))
                  + p(2)/
                 (-1 + c*q(2) + d*(p(1)*(-1 + c*q(1)) + 
                      (-1 + p(2))*(-1 + c*q(2)))))) 
                      );

                q2tplus1 = 0.01 * (
        c*(-(((-1 + d)*p(2))/
                (-1 + c*q(2) + d*(p(1)*(-1 + c*q(1)) + 
                     (-1 + p(2))*(-1 + c*q(2))))) + 
             (Power(-1 + d,2)*p(2)*(p(1) + p(2))*(1 - c*q(2)))/
              (n*Power(-1 + c*q(2) + 
                  d*(p(1)*(-1 + c*q(1)) + (-1 + p(2))*(-1 + c*q(2))),2)*
                (1 - (Power(-1 + d,2)*(-1 + n)*p(1)*
                     Power(1 + c*(-1 + q(1)),2))/
                   (n*Power(1 + d*(-1 + p(1) + p(2)) + 
                       c*(-1 + q(1) + 
                          d*(1 + p(1)*(-1 + q(1)) - q(1) + 
                             p(2)*(-1 + q(2)))),2)) - 
                  (Power(-1 + d,2)*(-1 + n)*p(2)*Power(-1 + c*q(2),2))/
                   (n*Power(-1 + c*q(2) + 
                       d*(p(1)*(-1 + c*q(1)) + (-1 + p(2))*(-1 + c*q(2)))
                       ,2)))) + d*(p(1)/
                 (1 + d*(-1 + p(1) + p(2)) + 
                   c*(-1 + q(1) + d*
                       (1 + p(1)*(-1 + q(1)) - q(1) + p(2)*(-1 + q(2)))))
                  + p(2)/
                 (-1 + c*q(2) + d*(p(1)*(-1 + c*q(1)) + 
                      (-1 + p(2))*(-1 + c*q(2))))))

                        );

                if (fabs(q2tplus1 - q2) < 1e-08 &&
                        fabs(q1tplus1 - q1) < 1e-08
                )
                {
                    q2 = q2tplus1;
                    q1 = q1tplus1;
                    DataFile << q1 << ";" << q2 << ";" << p1 << ";" << d << ";" << n << ";" << c << ";" << endl;
                    break;
                }
                
                q2 = q2tplus1;
                q1 = q1tplus1;
            }
        }
}
