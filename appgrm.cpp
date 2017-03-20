#include <cmath>
#include <memory>
#include <limits>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include "appgrm.h"
#include "cmdline.h"
#include "vcfio.h"

namespace {

template<typename T>
int count_shared_allele(T a, T b, T c, T d)
{
    if (a == b)
        return static_cast<int>(a == c) + (a == d);

    if (c == d)
        return static_cast<int>(c == a) + (c == b);

    return static_cast<int>(a == c) + (a == d) + (b == c) + (b == d);
}

std::pair<size_t,size_t> match_kernel_haplo(size_t n, const allele_t *x, const allele_t *y)
{
    size_t a = 0, b = 0;
    for (size_t i = 0; i < n; ++i) {
        if (x[i] && y[i]) {
            ++a;
            if (x[i] == y[i])
                ++b;
        }
    }
    return std::make_pair(a, b);
}

std::pair<size_t,size_t> match_kernel_diplo(size_t n, const allele_t *x, const allele_t *y)
{
    size_t a = 0, b = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t j1 = i*2, j2 = i*2+1;
        if (x[j1] && x[j2] && y[j1] && y[j2]) {
            a += 2;
            b += count_shared_allele(x[j1], x[j2], y[j1], y[j2]);
        }
    }
    return std::make_pair(a, b);
}


void recode_012_haploid(const vector<allele_t> &g, vector<double> &x)
{
    size_t n = g.size();
    static const allele_t ref = 1;

    x.assign(n, std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < n; ++i) {
        if ( g[i] )
            x[i] = g[i] == ref ? 2.0 : 0.0;
    }
}

void recode_012_diploid(const vector<allele_t> &g, vector<double> &x)
{
    size_t n = g.size() / 2;
    static const allele_t ref = 1;

    x.assign(n, std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < n; ++i) {
        auto a = g[i*2], b = g[i*2+1];
        if (a && b)
            x[i] = static_cast<int>(a == ref) + static_cast<int>(b == ref);
    }
}

bool is_biallelic(const Genotype &gt)
{
    for (auto &v : gt.allele)
        if (v.size() > 2)
            return false;
    return true;
}

} // namespace

int AppGRM::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("grm [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--out", "output file prefix", "appgrm.out");
    cmd->add("--method", "IBS/EIGENSTRAT/VanRaden1/VanRaden2", "IBS");

    if (argc < 2) {
        cmd->help();
        return 1;
    }

    cmd->parse(argc, argv);

    m_par.vcf = cmd->get("--vcf");
    m_par.out = cmd->get("--out");
    m_par.method = cmd->get("--method");

    cmd.reset();

    std::transform(m_par.method.begin(), m_par.method.end(), m_par.method.begin(), ::toupper);

    int info = perform();

    return info;
}

int AppGRM::perform()
{
    load_genotype();

    if ( m_gt.loc.empty() || m_gt.ind.size() < 2 ) {
        std::cerr << "ERROR: not enough observations\n";
        return 1;
    }

    if (m_par.method != "IBS") {
        if ( ! is_biallelic(m_gt) )
            std::cerr << "WARNING: method " << m_par.method
                      << " is not compatible with multi-allelic marker\n";
    }

    if (m_par.method == "IBS") {
        calc_grm_IBS(m_gt, m_grm.dat);
    }
    else if (m_par.method == "EIGENSTRAT") {
        calc_grm_EIGENSTRAT(m_gt, m_grm.dat);
    }
    else if (m_par.method == "VANRADEN1") {
        calc_grm_VanRaden1(m_gt, m_grm.dat);
    }
    else if (m_par.method == "VANRADEN2") {
        calc_grm_VanRaden2(m_gt, m_grm.dat);
    }
    else {
        std::cerr << "ERROR: invalid method: " << m_par.method << "\n";
        return 1;
    }

    m_grm.ind = m_gt.ind;

    std::ofstream ofs(m_par.out);

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << m_par.out << "\n";
        return 1;
    }

    auto n = m_grm.ind.size();

    for (size_t i = 0; i < n; ++i) {
        ofs << m_grm.ind[i];
        for (size_t j = 0; j < n; ++j)
            ofs << "\t" << m_grm.dat[i][j];
        ofs << "\n";
    }

    return 0;
}

void AppGRM::load_genotype()
{
    if ( m_par.vcf.empty() )
        return;

    std::cerr << "INFO: reading genotype file...\n";

    int info = read_vcf(m_par.vcf, m_gt);

    if (info != 0) {
        m_gt.loc.clear();
        m_gt.ind.clear();
        m_gt.dat.clear();
    }

    std::cerr << "INFO: " << m_gt.ind.size() << " individuals and "
              << m_gt.loc.size() << " loci were observed\n";
}

void AppGRM::calc_grm_IBS(const Genotype &gt, vector< vector<double> > &x)
{
    static const size_t nb = 1000;

    size_t m = 0;
    size_t n = gt.ind.size();
    vector< vector<size_t> > z(n, vector<size_t>(n, 0));

    x.assign(n, vector<double>(n, 0.0));

    if (gt.ploidy == 1) {
        vector< vector<allele_t> > dat(n, vector<allele_t>(nb));

        for (auto &v : gt.dat) {
            if (m < nb) {
                for (size_t i = 0; i < n; ++i)
                    dat[i][m] = v[i];
                ++m;
                continue;
            }
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    if (j <= i)
                        continue;
                    auto p = match_kernel_haplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
            for (size_t i = 0; i < n; ++i)
                dat[i][0] = v[i];
            m = 1;
        }

        if (m > 0) {
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    auto p = match_kernel_haplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
        }
    }
    else {
        vector< vector<allele_t> > dat(n, vector<allele_t>(nb*2));

        for (auto &v : gt.dat) {
            if (m < nb) {
                for (size_t i = 0; i < n; ++i) {
                    dat[i][m*2] = v[i*2];
                    dat[i][m*2+1] = v[i*2+1];
                }
                ++m;
                continue;
            }
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    if (j <= i)
                        continue;
                    auto p = match_kernel_diplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
            for (size_t i = 0; i < n; ++i) {
                dat[i][0] = v[i*2];
                dat[i][1] = v[i*2+1];
            }
            m = 1;
        }

        if (m > 0) {
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i+1; j < n; ++j) {
                    auto p = match_kernel_diplo(m, dat[i].data(), dat[j].data());
                    z[i][j] += p.first;
                    z[j][i] += p.second;
                }
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        x[i][i] = 1.0;
        for (size_t j = i+1; j < n; ++j) {
            if (z[i][j] != 0) {
                auto a = static_cast<double>(z[j][i]) / z[i][j];
                x[i][j] = x[j][i] = a;
            }
        }
    }
}

// Price A.L. et al., Nat Genet, 2006, 38(8): 904-9

void AppGRM::calc_grm_EIGENSTRAT(const Genotype &gt, vector< vector<double> > &x)
{
    size_t n = gt.ind.size();
    std::vector<double> z;

    x.assign(n, vector<double>(n, 0.0));

    for (auto &v : gt.dat) {
        if (gt.ploidy == 1)
            recode_012_haploid(v, z);
        else
            recode_012_diploid(v, z);

        size_t nz = 0;
        double sz = 0.0;
        for (auto e : z) {
            if (e == e) {
                ++nz;
                sz += e;
            }
        }

        auto p = (sz + 1) / (2 * nz + 2);
        auto mu = sz / nz;
        auto sd = std::sqrt( p * (1 - p) );

        for (auto &e : z)
            e = e != e ? 0.0 : (e - mu) / sd;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j)
                x[i][j] += z[i] * z[j];
        }
    }

    double d = 0.0;
    for (size_t i = 0; i < n; ++i)
        d += x[i][i];
    d /= (n - 1);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j)
            x[i][j] = x[j][i] = x[i][j] / d;
    }
}

// VanRaden P.M., J Dairy Sci, 2008, 91(11): 4414-23

void AppGRM::calc_grm_VanRaden1(const Genotype &gt, vector< vector<double> > &x)
{
    size_t n = gt.ind.size();
    std::vector<double> z;

    x.assign(n, vector<double>(n, 0.0));

    double scal = 0.0;
    for (auto &v : gt.dat) {
        if (gt.ploidy == 1)
            recode_012_haploid(v, z);
        else
            recode_012_diploid(v, z);

        size_t nz = 0;
        double sz = 0.0;
        for (auto e : z) {
            if (e == e) {
                ++nz;
                sz += e;
            }
        }

        auto p = 0.5 * sz / nz;
        auto mu = 2 * p;
        scal += 2 * p * (1 - p);

        for (auto &e : z)
            e = e != e ? 0.0 : e - mu;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j)
                x[i][j] += z[i] * z[j];
        }
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j)
            x[i][j] = x[j][i] = x[i][j] / scal;
    }
}

// VanRaden P.M., J Dairy Sci, 2008, 91(11): 4414-23

void AppGRM::calc_grm_VanRaden2(const Genotype &gt, vector< vector<double> > &x)
{
    size_t n = gt.ind.size();
    std::vector<double> z;

    x.assign(n, vector<double>(n, 0.0));

    for (auto &v : gt.dat) {
        if (gt.ploidy == 1)
            recode_012_haploid(v, z);
        else
            recode_012_diploid(v, z);

        size_t nz = 0;
        double sz = 0.0;
        for (auto e : z) {
            if (e == e) {
                ++nz;
                sz += e;
            }
        }

        auto p = 0.5 * sz / nz;
        auto mu = 2 * p;
        auto sd = std::sqrt( 2 * p * (1 - p) );

        for (auto &e : z)
            e = e != e ? 0.0 : (e - mu) / sd;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j)
                x[i][j] += z[i] * z[j];
        }
    }

    size_t m = gt.loc.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j)
            x[i][j] = x[j][i] = x[i][j] / m;
    }
}
