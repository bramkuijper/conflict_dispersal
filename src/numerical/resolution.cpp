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

void bound01(double &val)
{
    if (val < 0.0001)
    {
        val = 0.0001;
    } else if (val > 0.9999)
    {
        val = 0.9999;
    }
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
    double qS, qStplus1;
    double qNS, qNStplus1;

    double p1 = atof(argv[1]);
    double p2 = 1-p1;
    double n = atof(argv[2]);
    double d = atof(argv[3]);

    DataFile << "type;q1;q2;qS;qNS;s1;s2;p1;d;n;c1;c2;siginfo;" << endl;

    // first offspring strategies
    for (double c2 = 0; c2 <= 1.0; c2+=0.02)
    {
        for (double c1 = 0; c1 <= 1.0; c1+=0.02)
        {
            qNS = 0.1;
            qS = 0.9;
            s1 = 0.9;
            s2 = 0.1;

            for (size_t iter = 0; iter < 10000000; ++iter)
            {
                s1tplus1 = s1 + 0.01 *
(
-((pow(1 - d,2)*(-qNS + qS + (qNS - qS)*(1 - c1))*p1*(p1 + p2)*
        (qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)))/
      (n*pow((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
          d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
             p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2))),2)*
        (1 - (pow(-1 + d,2)*(-1 + n)*p2*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2,2))/
           (n*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2 + 
               d*(p1*(1 - qS*c2*s1) + qNS*c2*(p1*(-1 + s1) + (-1 + p2)*(-1 + s2)) - 
                  (-1 + p2)*(-1 + qS*c2*s2)),2)) - 
          (pow(-1 + d,2)*(-1 + n)*p1*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1),2))/
           (n*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1) + 
               d*(1 - p1 - p2 + c1*(-1 + p1 + p2 + qS*s1 - qS*p1*s1 + 
                     qNS*(1 + p1*(-1 + s1) - s1 + p2*(-1 + s2)) - qS*p2*s2)),2))))) + 
   ((1 + (p1 + p2)/(n*(1 - (pow(-1 + d,2)*(-1 + n)*p2*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2,2))/
              (n*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2 + 
                  d*(p1*(1 - qS*c2*s1) + qNS*c2*(p1*(-1 + s1) + (-1 + p2)*(-1 + s2)) - 
                     (-1 + p2)*(-1 + qS*c2*s2)),2)) - 
             (pow(-1 + d,2)*(-1 + n)*p1*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1),2))/
              (n*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1) + 
                  d*(1 - p1 - p2 + c1*(-1 + p1 + p2 + qS*s1 - qS*p1*s1 + 
                        qNS*(1 + p1*(-1 + s1) - s1 + p2*(-1 + s2)) - qS*p2*s2)),2)))))*
      (((1 - d)*(-qNS + qS + (qNS - qS)*(1 - c1))*p1)/
         ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
           d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
              p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
        d*(((-qNS + qS + (qNS - qS)*(1 - c1))*p1)/
            ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
              d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
                 p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
           ((qNS - qS + (-qNS + qS)*(1 - c2))*p2)/
            ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
              d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
                 p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))))))/2.

);

                s2tplus1 = s2  + 0.01 *(
-((pow(1 - d,2)*(qNS - qS + (-qNS + qS)*(1 - c2))*p2*(p1 + p2)*
        (1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))/
      (n*pow((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
          d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
             p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2))),2)*
        (1 - (pow(-1 + d,2)*(-1 + n)*p2*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2,2))/
           (n*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2 + 
               d*(p1*(1 - qS*c2*s1) + qNS*c2*(p1*(-1 + s1) + (-1 + p2)*(-1 + s2)) - 
                  (-1 + p2)*(-1 + qS*c2*s2)),2)) - 
          (pow(-1 + d,2)*(-1 + n)*p1*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1),2))/
           (n*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1) + 
               d*(1 - p1 - p2 + c1*(-1 + p1 + p2 + qS*s1 - qS*p1*s1 + 
                     qNS*(1 + p1*(-1 + s1) - s1 + p2*(-1 + s2)) - qS*p2*s2)),2))))) + 
   ((1 + (p1 + p2)/(n*(1 - (pow(-1 + d,2)*(-1 + n)*p2*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2,2))/
              (n*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2 + 
                  d*(p1*(1 - qS*c2*s1) + qNS*c2*(p1*(-1 + s1) + (-1 + p2)*(-1 + s2)) - 
                     (-1 + p2)*(-1 + qS*c2*s2)),2)) - 
             (pow(-1 + d,2)*(-1 + n)*p1*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1),2))/
              (n*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1) + 
                  d*(1 - p1 - p2 + c1*(-1 + p1 + p2 + qS*s1 - qS*p1*s1 + 
                        qNS*(1 + p1*(-1 + s1) - s1 + p2*(-1 + s2)) - qS*p2*s2)),2)))))*
      (((1 - d)*(qNS - qS + (-qNS + qS)*(1 - c2))*p2)/
         ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
           d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
              p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))) + 
        d*(((-qNS + qS + (qNS - qS)*(1 - c1))*p1)/
            ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
              d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
                 p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
           ((qNS - qS + (-qNS + qS)*(1 - c2))*p2)/
            ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
              d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
                 p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))))))/2.
                        );

                qStplus1 = qS + 0.01 * (

((1 - d)*p1*(s1 - (1 - c1)*s1))/
    ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
      d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
         p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
   ((1 - d)*p2*(-s2 + (1 - c2)*s2))/
    ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
      d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
         p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))) + 
   ((p1 + p2)*(-((pow(1 - d,2)*p1*(s1 - (1 - c1)*s1)*
             (qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)))/
           pow((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
             d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
                p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2))),2)) - 
        (pow(1 - d,2)*p2*(-s2 + (1 - c2)*s2)*
           (1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))/
         pow((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
           d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
              p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2))),2)))/
    (n*(1 - (pow(-1 + d,2)*(-1 + n)*p2*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2,2))/
         (n*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2 + 
             d*(p1*(1 - qS*c2*s1) + qNS*c2*(p1*(-1 + s1) + (-1 + p2)*(-1 + s2)) - 
                (-1 + p2)*(-1 + qS*c2*s2)),2)) - 
        (pow(-1 + d,2)*(-1 + n)*p1*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1),2))/
         (n*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1) + 
             d*(1 - p1 - p2 + c1*(-1 + p1 + p2 + qS*s1 - qS*p1*s1 + 
                   qNS*(1 + p1*(-1 + s1) - s1 + p2*(-1 + s2)) - qS*p2*s2)),2)))) + 
   d*((p1*(s1 - (1 - c1)*s1))/
       ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
         d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
            p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
      (p2*(-s1 + (1 - c2)*s1))/
       ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
         d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
            p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2))))) + 
   d*((p1*(s2 - (1 - c1)*s2))/
       ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
         d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
            p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
      (p2*(-s2 + (1 - c2)*s2))/
       ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
         d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
            p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))))

                );

                qNStplus1 = qNS + 0.01 * (
((1 - d)*p1*(1 + (1 - c1)*(-1 + s1) - s1))/
    ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
      d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
         p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
   ((1 - d)*p2*(-1 + (1 - c2)*(1 - s2) + s2))/
    ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
      d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
         p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))) + 
   ((p1 + p2)*(-((pow(1 - d,2)*p1*(1 + (1 - c1)*(-1 + s1) - s1)*
             (qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)))/
           pow((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
             d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
                p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2))),2)) - 
        (pow(1 - d,2)*p2*(-1 + (1 - c2)*(1 - s2) + s2)*
           (1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))/
         pow((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
           d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
              p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2))),2)))/
    (n*(1 - (pow(-1 + d,2)*(-1 + n)*p2*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2,2))/
         (n*pow(1 + qNS*c2*(-1 + s2) - qS*c2*s2 + 
             d*(p1*(1 - qS*c2*s1) + qNS*c2*(p1*(-1 + s1) + (-1 + p2)*(-1 + s2)) - 
                (-1 + p2)*(-1 + qS*c2*s2)),2)) - 
        (pow(-1 + d,2)*(-1 + n)*p1*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1),2))/
         (n*pow(-1 + c1*(1 + qNS*(-1 + s1) - qS*s1) + 
             d*(1 - p1 - p2 + c1*(-1 + p1 + p2 + qS*s1 - qS*p1*s1 + 
                   qNS*(1 + p1*(-1 + s1) - s1 + p2*(-1 + s2)) - qS*p2*s2)),2)))) + 
   d*((p1*(1 + (1 - c1)*(-1 + s1) - s1))/
       ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
         d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
            p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
      (p2*(-1 + (1 - c2)*(1 - s1) + s1))/
       ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
         d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
            p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2))))) + 
   d*((p1*(1 + (1 - c1)*(-1 + s2) - s2))/
       ((1 - d)*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
         d*(p1*(qNS*(1 - s1) + qS*s1 + (1 - c1)*(1 - qNS*(1 - s1) - qS*s1)) + 
            p2*(qNS*(1 - s2) + qS*s2 + (1 - c1)*(1 - qNS*(1 - s2) - qS*s2)))) + 
      (p2*(-1 + (1 - c2)*(1 - s2) + s2))/
       ((1 - d)*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)) + 
         d*(p1*(1 - qNS*(1 - s1) - qS*s1 + (1 - c2)*(qNS*(1 - s1) + qS*s1)) + 
            p2*(1 - qNS*(1 - s2) - qS*s2 + (1 - c2)*(qNS*(1 - s2) + qS*s2)))))
                );

                bound01(qNStplus1);
                bound01(qStplus1);
                bound01(s1tplus1);
                bound01(s2tplus1);


                if (
                        (fabs(qNStplus1 - qNS) < 1e-08 &&
                        fabs(qStplus1 - qS) < 1e-08 &&
                        fabs(s1tplus1 - s1) < 1e-08 &&
                        fabs(s2tplus1 - s2) < 1e-08
                        )
                        ||
                        iter == 10000000 - 1
                )
                {
                    qS = qStplus1;
                    qNS = qNStplus1;
                    s1 = s1tplus1;
                    s2 = s2tplus1;
                    DataFile << "offspring;" << s1 * qS + (1-s1) * qNS << ";" << s2 * qS + (1-s2) * qNS << ";" << qS << ";" << qNS << ";" << s1 << ";" << s2 << ";" << p1 << ";" << d << ";" << n << ";" << c1 << ";" << c2 << ";" << siginfo(s1, s2) << ";" <<  endl;
                    break;
                }
                
                qS = qStplus1;
                qNS = qNStplus1;
                s1 = s1tplus1;
                s2 = s2tplus1;
            }
        }
    }
}
