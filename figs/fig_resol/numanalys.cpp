#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include "bramauxiliary.h"

using namespace std;

string filename("iter_predict_adapt");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

double entropy(double const val)
{
    return(val == 0 ? 0 : val == 1.0 ? 0 : -(val * log2(val) + (1-val) * log2(1-val)));
}

double bound(double const val)
{
    return(val < 0.001 ? 0.001 : val > 0.999 ? 0.999 : val);
}

int main(int argc, char **argv) 
{
    double dNS,dS,s1,s2, dwmatds1, dwmatds2, dwdS, dwdNS;
    double dNS_t1,dS_t1,s1_t1,s2_t1;

    double nvals[1] = { 5 };

    cout << "time;p1;n;c1;c2;s1;s2;dS;dNS;d1;d2;siginfo;" << endl;

    for (double p1 = 0.5; p1 <= .5; p1 += 0.15)
    {
        for (int n_i = 0; n_i < 1; ++n_i)
        {
            double n = nvals[n_i];

            for (double c1 = 0.0; c1 < 1.01; c1 += 0.01)
            {
                for (double c2 = 0.0; c2 < 1.01; c2 += 0.01)
                {
                    dNS = 0.9;
                    dS = 0.1;
                    s1 = 0.1;
                    s2 = 0.9;

                    if (c1 < c2)
                    {
                        dNS = 0.9;
                        dS = 0.1;
                        s1 = 0.5;
                        s2 = 0.5;
                    }

                    int time = 0;
            
                    for (time = 0; time < 100000; ++time)
                    {
                        dwmatds1 = 
                           -(((dNS - dS)*(-((1 - p1)/n) - p1/n)*(1 - dNS*(1 - s1) - dS*s1))/
      (pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
          (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)*
        (-1 + ((-1 + n)*p1*pow(1 - dNS*(1 - s1) - dS*s1,2))/
           (n*pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
               (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) + 
          ((-1 + n)*(1 - p1)*pow(1 - dNS*(1 - s2) - dS*s2,2))/
           (n*pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
               (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2))))) + 
   ((1 + (-((1 - p1)/n) - p1/n)/
         (-1 + ((-1 + n)*p1*pow(1 - dNS*(1 - s1) - dS*s1,2))/
            (n*pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
                (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) + 
           ((-1 + n)*(1 - p1)*pow(1 - dNS*(1 - s2) - dS*s2,2))/
            (n*pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
                (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2))))*
      ((dNS - dS)/(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
        (-dNS + dS)*(1 - c1)*(p1/
            (1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
              (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
           (1 - p1)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
              (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)))))/2.; 
                            
                        dwmatds2 = 

-(((dNS - dS)*(-((1 - p1)/n) - p1/n)*(1 - dNS*(1 - s2) - dS*s2))/
      (pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
          (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)*
        (-1 + ((-1 + n)*p1*pow(1 - dNS*(1 - s1) - dS*s1,2))/
           (n*pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
               (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) + 
          ((-1 + n)*(1 - p1)*pow(1 - dNS*(1 - s2) - dS*s2,2))/
           (n*pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
               (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2))))) + 
   ((1 + (-((1 - p1)/n) - p1/n)/
         (-1 + ((-1 + n)*p1*pow(1 - dNS*(1 - s1) - dS*s1,2))/
            (n*pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
                (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) + 
           ((-1 + n)*(1 - p1)*pow(1 - dNS*(1 - s2) - dS*s2,2))/
            (n*pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
                (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2))))*
      ((dNS - dS)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
        (-dNS + dS)*(1 - c2)*(p1/
            (1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
              (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
           (1 - p1)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
              (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)))))/2.;

                        dwdS = 
-(s1/(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
        (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2))) - 
   s2/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
      (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
   ((-((1 - p1)/n) - p1/n)*((s1*(1 - dNS*(1 - s1) - dS*s1))/
         pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2) + 
        (s2*(1 - dNS*(1 - s2) - dS*s2))/
         pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)))/
    (-1 + ((-1 + n)*p1*pow(1 - dNS*(1 - s1) - dS*s1,2))/
       (n*pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) + 
      ((-1 + n)*(1 - p1)*pow(1 - dNS*(1 - s2) - dS*s2,2))/
       (n*pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2))) + 
   (1 - c1)*s1*(p1/
       (1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
      (1 - p1)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2))) + 
   (1 - c2)*s2*(p1/
       (1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
      (1 - p1)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)));


                        dwdNS = 
   (-1 + s1)/(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
      (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
   (-1 + s2)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
      (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
   ((-((1 - p1)/n) - p1/n)*(-(((-1 + s1)*(1 - dNS*(1 - s1) - dS*s1))/
           pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
             (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) - 
        ((-1 + s2)*(1 - dNS*(1 - s2) - dS*s2))/
         pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)))/
    (-1 + ((-1 + n)*p1*pow(1 - dNS*(1 - s1) - dS*s1,2))/
       (n*pow(1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2)) + 
      ((-1 + n)*(1 - p1)*pow(1 - dNS*(1 - s2) - dS*s2,2))/
       (n*pow(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
           (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2),2))) + 
   (1 - c1)*(1 - s1)*(p1/
       (1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
      (1 - p1)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2))) + 
   (1 - c2)*(1 - s2)*(p1/
       (1 - dNS*(1 - s1) - dS*s1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)) + 
      (1 - p1)/(1 + (1 - c1)*p1*(dNS*(1 - s1) + dS*s1) - dNS*(1 - s2) - dS*s2 + 
         (1 - c2)*(1 - p1)*(dNS*(1 - s2) + dS*s2)));
    
    
                        s1_t1 = bound(s1 + .2*dwmatds1);
                        s2_t1 = bound(s2 + .2*dwmatds2);
                        dS_t1 = bound(dS + .2*dwdS);
                        dNS_t1 = bound(dNS + .2*dwdNS);

                        if (fabs(s1_t1 - s1) < 1e-05 &&
                            fabs(s2_t1 - s2) < 1e-05 &&
                            fabs(dS_t1 - dS) < 1e-05 &&
                            fabs(dNS_t1 - dNS) < 1e-05) 
                        {
                            s1 = s1_t1;
                            s2 = s2_t1;
                            dS = dS_t1;
                            dNS = dNS_t1;
                            break;
                        }

                        s1 = s1_t1;
                        s2 = s2_t1;
                        dS = dS_t1;
                        dNS = dNS_t1;
                    }

                    double d1 = s1 * dS + (1-s1) * dNS;
                    double d2 = s2 * dS + (1-s2) * dNS;
                    double signinfo = 1 - (((s1 + s2)/2.0) * entropy(s1/(s1+s2)) + (1-((s1+s2)/2.0))*entropy((1-s1)/((1-s1)+(1-s2))));

                    if (c1 == c2)
                    {
                        signinfo = 0.0;
                    }
                    cout << time << ";" << p1 << ";" << n << ";" << c1 << ";" << c2 << ";" << s1 << ";" << s2 << ";" << dS << ";" << dNS << ";" << d1 <<  ";" << d2 << ";" << signinfo << ";" << endl;
                }
            }
        }
    }
}
