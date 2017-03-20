#ifndef APPGRM_H
#define APPGRM_H

#include "main.h"

class AppGRM
{
public:
    int run(int argc, char *argv[]);

    static void calc_grm_IBS(const Genotype &gt, vector< vector<double> > &x);

    static void calc_grm_EIGENSTRAT(const Genotype &gt, vector< vector<double> > &x);

    static void calc_grm_VanRaden1(const Genotype &gt, vector< vector<double> > &x);

    static void calc_grm_VanRaden2(const Genotype &gt, vector< vector<double> > &x);

private:
    int perform();

    void load_genotype();

private:
    struct Params
    {
        string vcf;
        string out;
        string method;
    };

    Params m_par;
    Genotype m_gt;
    SquareData m_grm;
};

#endif // APPGRM_H
