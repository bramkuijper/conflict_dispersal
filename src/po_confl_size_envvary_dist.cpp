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
#include <gsl/gsl_histogram.h>
#include "bramauxiliary.h"


//#define NDEBUG

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const int Npatches = 2000; 
const int Npp = 4; // individuals per patch
const int numgen = 50000;
const int Clutch = 50;

// mutation rates
double mu = 0;
double sdmu = 0;

double init_s0 = 0;
double init_s1 = 0;
double init_dS = 0;
double init_dNS = 0;

const int nbins = 200;

// dispersal costs 
double c[2] = { 0, 0 };

// frequency of patch 1
double p = 0;

// tally of dispersers
int Ndisp = 0;

// runtime for stats
time_t total_time; 

int generation = 0;

int seed = -1;

// skip the number of generations in output
// to prevent output files from becoming too large
int skip = 10;
int dist_skip = 100;

// haploid individual
struct Individual
{
    // signal dependent dispersal rates
    double d[2];

    // environment-dependent signalling probs
    double s[2];
};

struct Patch
{
    Individual locals[Npp]; // all the local breeders

    // philopatric offspring
    Individual phils[2 * Npp * Clutch];     
    // (note that dispersing offspring
    // ends up in global pool, see below)

    // total number of kids 
    int Nkids; 

    // variable that allows for correct
    // rounding of the number of immigrants
    // per patch (see below)
    int immigrant_bonus;

    // local environmental state
    bool envt;
};

// generate the complete population
Patch MetaPop[Npatches];
Individual Dispersers[Npatches * Npp * Clutch * 2];

// give the outputfile a unique name
// by using the create_filename function (see 
// "bramauxiliary.h")
string filename("sim_conflict");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

ofstream DistFile((filename_new + "dist").c_str());  // output file 
// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]); // mutation prob
    sdmu = atof(argv[2]); // mutational distribution
    c[0] = atof(argv[3]); 
    c[1] = atof(argv[4]); // survival slope
    p = atof(argv[5]); // dispersal at t=0
    init_s0 = atof(argv[6]);
    init_s1 = atof(argv[7]);
    init_dS = atof(argv[8]);
    init_dNS = atof(argv[9]);
}

void init_pop()
{
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set the seed to the random number generator
    // stupidly enough, this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // go through all patches
    for (int i = 0; i < Npatches; ++i)
    {
        // initialize the population with an environment at the
        // desired frequency
        MetaPop[i].envt = gsl_rng_uniform(r) < p;

        for (int j = 0; j < Npp; ++j)
        {
            MetaPop[i].locals[j].d[0] = init_dS;
            MetaPop[i].locals[j].d[1] = init_dNS;
            MetaPop[i].locals[j].s[0] = init_dS;
            MetaPop[i].locals[j].s[1] = init_dNS;
        }
    }
}

// mutate an allele
double Mut(double h)
{
    h+=gsl_rng_uniform(r) < mu ? gsl_ran_gaussian(r, sdmu) : 0;
    h = h < 0 ? 0 : h > 1.0 ? 1.0 : h;
    return(h);
}

// allocate a kid and give it genes
void Create_Kid(Individual &mother, Individual &father, Individual &Kid)
{
    for (int i = 0; i < 2; ++i)
    {
        Kid.d[i] = gsl_rng_uniform(r) < 0.5 ? Mut(mother.d[i]) : Mut(father.d[i]);
        Kid.s[i] = gsl_rng_uniform(r) < 0.5 ? Mut(mother.s[i]) : Mut(father.s[i]);
    }
}

// mate and create kids across all patches...
void make_juveniles()
{
    Ndisp = 0;

    // variable to hold the mother's signal
    double maternal_signal = 0;
    // and the ensuing offspring's dispersal probability
    double pdispersal = 0;

    // generate offspring on each patch
    for (int i = 0; i < Npatches; ++i)
    {
        // reset local offspring counter
        MetaPop[i].Nkids = 0;

        // variable that takes into account
        // any beyond-the-minimum number of immigrants
        // (see later)
        MetaPop[i].immigrant_bonus = 0;

        for (int j = 0; j < Npp; ++j)
        {
            // find random father within patch
            int father = gsl_rng_uniform_int(r, Npp);

            // express the mother's signal
            maternal_signal = MetaPop[i].locals[j].s[MetaPop[i].envt];

            // create kids and reduce parental resources
            for (int k = 0; k < Clutch; ++k)
            {
                Individual Kid; 

                Create_Kid(MetaPop[i].locals[j],MetaPop[i].locals[father],Kid);

                // effective offspring dispersal probability is given by
                // P(mother gives signal) * P(dispersal upon receiving a signal)
                // or
                // P(mother does not give a signal) * P(dispersal when signal is absent)
                pdispersal = gsl_rng_uniform(r) < maternal_signal ? Kid.d[0] : Kid.d[1];

                if (gsl_rng_uniform(r) < pdispersal)
                {
                    if (gsl_rng_uniform(r) > c[MetaPop[i].envt])
                    {
                        Dispersers[Ndisp++] = Kid;
                    }
                }
                else
                {
                    MetaPop[i].phils[MetaPop[i].Nkids++] = Kid;
                }
            } //clutch
        }//Npp
    }//Npatches

    assert(Ndisp < Npatches * Npp * Clutch * 2);
}

// replacement of adults with juveniles
void replace_adults()
{
    // okay, we have Ndisp/Npatches dispersers per patch
    int dispersers_per_patch = floor((double) Ndisp / Npatches);
    //cout << "Ndisp: " << Ndisp << ", Npatches: " << Npatches << " disp per patch: " << dispersers_per_patch << endl;

    // however, we need to correctly round this rational number
    // to the actual number of immigrants for 
    // a given patch. To this end, 
    // we just randomly distribute the dispersers that remain after the
    // previous rounding over the different patches
    for (int i = 0; i < Ndisp - Npatches * dispersers_per_patch; ++i)
    {
        // randomly picked patch receives additional immigrant
        MetaPop[gsl_rng_uniform_int(r,Npatches)].immigrant_bonus++;
    }

    
    // now replace local breeders on each patch
    for (int i = 0; i < Npatches; ++i)
    {
        if (dispersers_per_patch + MetaPop[i].immigrant_bonus > 0)
        {
            assert(Ndisp > 0);
        }

        int arriving_immigrants = dispersers_per_patch + MetaPop[i].immigrant_bonus;

        for (int j = 0; j < Npp; ++j)
        {
            assert((double) arriving_immigrants / (arriving_immigrants + MetaPop[i].Nkids) >= 0.0 && (double) arriving_immigrants / (arriving_immigrants + MetaPop[i].Nkids) <= 1.0);
            if (gsl_rng_uniform(r) < (double) arriving_immigrants / (arriving_immigrants + MetaPop[i].Nkids))
            {
                int rand_disp = gsl_rng_uniform_int(r,Ndisp);
                MetaPop[i].locals[j] = Dispersers[rand_disp];
                Dispersers[rand_disp] = Dispersers[Ndisp-1];
                --Ndisp;
                --arriving_immigrants;
            }
            else
            {
                int rand_phil = gsl_rng_uniform_int(r,MetaPop[i].Nkids);
                MetaPop[i].locals[j] = MetaPop[i].phils[rand_phil];
                MetaPop[i].phils[rand_phil] = MetaPop[i].phils[MetaPop[i].Nkids-1];
                --MetaPop[i].Nkids;
            }
        }

        // remove all remaining immigrants from the global dispersal pool
        for (int j = 0; j < arriving_immigrants; ++j)
        {
            int rand_disp = gsl_rng_uniform_int(r,Ndisp);
            Dispersers[rand_disp] = Dispersers[Ndisp-1];
            --Ndisp;
        }

        // randomly allocate environments for the following round
        // to prevent any local correlations building up
        MetaPop[i].envt = gsl_rng_uniform(r) < p;
    }

    assert(Ndisp==0);
}

void write_data_headers()
{
    DataFile << "generation;s0;s1;dsignal;dnosignal;d0;d1;vars0;vars1;vardsignal;vardnosignal;" << endl;
}

void write_data()
{
    double means0 = 0;
    double means1 = 0;
    double meand0 = 0;
    double meand1 = 0;
    double sss0 = 0;
    double sss1 = 0;
    double ssd0 = 0;
    double ssd1 = 0;

    for (int i = 0; i < Npatches; ++i)
    {
        for (int j = 0; j < Npp; ++j)
        {
            meand0 += MetaPop[i].locals[j].d[0];
            meand1 += MetaPop[i].locals[j].d[1];
            means0 += MetaPop[i].locals[j].s[0];
            means1 += MetaPop[i].locals[j].s[1];
            ssd0 += MetaPop[i].locals[j].d[0] * MetaPop[i].locals[j].d[0];
            ssd1 += MetaPop[i].locals[j].d[1] * MetaPop[i].locals[j].d[1];
            sss0 += MetaPop[i].locals[j].s[0] * MetaPop[i].locals[j].s[0];
            sss1 += MetaPop[i].locals[j].s[1] * MetaPop[i].locals[j].s[1];
        }
    }
    
    meand0 /= (Npatches * Npp);
    meand1 /= (Npatches * Npp);
    means0 /= (Npatches * Npp);
    means1 /= (Npatches * Npp);
        
    DataFile << generation << ";" << means0 << ";" << means1 << ";" 
                                << meand0 << ";" << meand1 << ";"
                                << means0 * meand0 + (1-means0) * meand1 << ";" << means1 * meand0 + (1-means1) * meand1 << ";"
                                << sss0 / (Npatches * Npp)  - means0 * means0 << ";" 
                                << sss1 / (Npatches * Npp)  - means1 * means1 << ";" 
                                << ssd0 / (Npatches * Npp)  - meand0 * meand0 << ";" 
                                << ssd1 / (Npatches * Npp)  - meand1 * meand1 << ";" << endl;
}

void write_dist_headers()
{
    DistFile << "generation;bin;dS;dNS;s0;s1;" << endl;    
}

// write the distribution 
void write_dist()
{
    gsl_histogram * histogram_d0 = gsl_histogram_alloc(nbins);
    gsl_histogram * histogram_d1 = gsl_histogram_alloc(nbins);
    gsl_histogram * histogram_s0 = gsl_histogram_alloc(nbins);
    gsl_histogram * histogram_s1 = gsl_histogram_alloc(nbins);

    gsl_histogram_set_ranges_uniform(histogram_d0, 0, 1.0001);
    gsl_histogram_set_ranges_uniform(histogram_d1, 0, 1.0001);
    gsl_histogram_set_ranges_uniform(histogram_s0, 0, 1.0001);
    gsl_histogram_set_ranges_uniform(histogram_s1, 0, 1.0001);

    for (int i = 0; i < Npatches; ++i)
    {
        for (int j = 0; j < Npp; ++j)
        {
            gsl_histogram_increment(histogram_d0,MetaPop[i].locals[j].d[0]);
            gsl_histogram_increment(histogram_d1,MetaPop[i].locals[j].d[1]);
            gsl_histogram_increment(histogram_s0,MetaPop[i].locals[j].s[0]);
            gsl_histogram_increment(histogram_s1,MetaPop[i].locals[j].s[1]);
        }
    }

    for (size_t i = 0; i < gsl_histogram_bins(histogram_d0); ++i)
    {
        DistFile << generation << ";" << (double) i / gsl_histogram_bins(histogram_d0) << ";" << gsl_histogram_get(histogram_d0, i) << ";"
            << gsl_histogram_get(histogram_d1, i) << ";"
            << gsl_histogram_get(histogram_s0, i) << ";"
            << gsl_histogram_get(histogram_s1, i) << ";" << endl;
    }

    DistFile << generation << ";1.0;0;0;0;0;" << endl;
    
    gsl_histogram_free(histogram_d0);
    gsl_histogram_free(histogram_d1);
    gsl_histogram_free(histogram_s0);
    gsl_histogram_free(histogram_s1);
}

void write_parameters()
{
    DataFile << "patch;" << Npatches << endl
                << "npp;" << Npp << endl
                << "numgen;" << numgen << endl
                << "mu;" << mu << endl
                << "sdmu;" << sdmu << endl
                << "c0;" << c[0] << endl
                << "c1;" << c[1] << endl
                << "p;" << p << endl
                << "runtime;" << total_time << endl;
}


int main(int argc, char * argv[])
{
    init_arguments(argc,argv);
    init_pop();

    write_data_headers();
    write_dist_headers();

    for (generation = 0; generation < numgen; ++generation)
    {
        make_juveniles();

        if (generation % skip == 0)
        {
            write_data();
        }

        if (generation % dist_skip == 0)
        {
            write_dist();
        }

        replace_adults();
    }

    write_data();
    write_dist();
    write_parameters();
}
