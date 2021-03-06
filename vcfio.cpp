#include <ctime>
#include <limits>
#include <fstream>
#include <iostream>
#include "vcfio.h"
#include "strsplit.h"
#include "util.h"

namespace {

const int NA = std::numeric_limits<int>::min();

// ret < 0: error; ret > 0 : ploidy
int parse_gt(const Token &t, int &a, int &b)
{
    auto beg = t.data();
    auto end = std::find(t.data(), t.data() + t.size(), ':');
    size_t len = end - beg;

    if (len == 0)
        return -1;

    if (len == 1 && beg[0] == '.')
        return 1;

    if (len == 3 && beg[0] == '.' && beg[2] == '.' && (beg[1] == '/' || beg[1] == '|'))
        return 2;

    if (std::count(beg, end, '/') + std::count(beg, end, '|') > 1) {
        std::cerr << "ERROR: unsuppored polyploidy genotype: " << string(beg,end) << "\n";
        return -1;
    }

    size_t pos = 0;
    while (pos < len && beg[pos] != '/' && beg[pos] != '|')
        ++pos;
    if (pos == 0 || pos == len - 1)
        return -1;

    bool ok = false;

    if (pos == len) {
        a = number<int>(string(beg,end), &ok);
        if ( ! ok || a < 0 )
            return -1;
        return 1;
    }

    string gt(beg, beg + pos);
    if (gt != ".") {
        a = number<int>(gt, &ok);
        if ( ! ok || a < 0 )
            return -1;
    }

    gt.assign(beg + pos + 1, end);
    if (gt != ".") {
        b = number<int>(gt, &ok);
        if ( ! ok || b < 0 )
            return -1;
    }

    return 2;
}

string datetime_now()
{
    char buf[16];
    auto t = std::time(nullptr);
    std::strftime(buf, sizeof(buf), "%Y%m%d%H%M%S", std::localtime(&t));
    return string(buf);
}

} // namespace

int read_vcf(const string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;
    string ver;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        if (ln == 1 && line.compare(0,13,"##fileformat=") == 0) {
            ver = line.substr(13);
            if (ver != "VCFv4.0" && ver != "VCFv4.1" && ver != "VCFv4.2") {
                std::cerr << "ERROR: unsupported VCF format version: " << ver << "\n";
                return 1;
            }
        }

        if (line.compare(0,7,"##INFO=") == 0)
            continue;

        if (line.compare(0,9,"##FILTER=") == 0)
            continue;

        if (line.compare(0,9,"##FORMAT=") == 0)
            continue;

        if (line.compare(0,9,"##contig=") == 0)
            continue;

        if (line.compare(0,2,"##") == 0)
            continue;

        if (line.compare(0,1,"#") == 0) {
            vector<string> vs;
            strsplit(delim, line.begin(), line.end(), vs);

            if (vs.size() != 8 && vs.size() < 10) {
                std::cerr << "ERROR: incorrect number of columns in header line: " << vs.size() << "\n";
                return 1;
            }

            if (vs[0] != "#CHROM" || vs[1] != "POS" || vs[2] != "ID" || vs[3] != "REF" ||
                vs[4] != "ALT" || vs[5] != "QUAL" || vs[6] != "FILTER" || vs[7] != "INFO") {
                std::cerr << "ERROR: incorrect column names in header line: "
                          << vs[0] << "\t" << vs[1] << "\t" << vs[2] << "\t" << vs[3] << "\t"
                          << vs[4] << "\t" << vs[5] << "\t" << vs[6] << "\t" << vs[7] << "\n";
                return 1;
            }

            if (vs.size() > 9) {
                if (vs[8] != "FORMAT") {
                    std::cerr << "ERROR: FORMAT is required at 9th field in header line: " << vs[8] << "\n";
                    return 1;
                }
                gt.ind.assign(vs.begin() + 9, vs.end());
            }

            break;
        }

        std::cerr << "ERROR: #CHROME header line is required\n";
        return 1;
    }

    if ( ver.empty() ) {
        std::cerr << "ERROR: ##fileformat field is required and must be the first line\n";
        return 1;
    }

    int ploidy = -1;
    size_t ncols = gt.ind.empty() ? 8 : gt.ind.size() + 9;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<Token> vt;
        strsplit(delim, line.begin(), line.end(), vt);

        if (vt.size() != ncols) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vt.size() << " != "
                      << ncols << "): " << filename << "\n";
            return 1;
        }

        gt.chr.push_back(string(vt[0]));
        gt.pos.push_back(number<int>(string(vt[1])));

        if (vt[2] == Token(".",1))
            gt.loc.push_back(string(vt[0]) + "_" + string(vt[1]));
        else
            gt.loc.push_back(string(vt[2]));

        if ( gt.ind.empty() )
            continue;

        auto format = string(vt[8]);
        if (format.compare(0,2,"GT") != 0) {
            std::cerr << "ERROR: GT is required and must be the first sub-field of FORMAT: " << format << "\n";
            return 1;
        }

        auto ref = string(vt[3]);
        auto alt = string(vt[4]);
        std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
        std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);

        vector<string> as;
        as.push_back(ref);
        strsplit([](char c) { return c == ','; }, alt.begin(), alt.end(), as);
        gt.allele.push_back(as);
        int nas = as.size();

        vector<allele_t> v;
        for (size_t i = 9; i < ncols; ++i) {
            int a = NA, b = NA;
            int info = parse_gt(vt[i], a, b);
            if (info < 0 || a >= nas || b >= nas) {
                std::cerr << "ERROR: invalid genotype data: " << string(vt[i]) << "\n";
                return 1;
            }

            if (ploidy == -1)
                ploidy = info;

            if (info != ploidy) {
                std::cerr << "ERROR: inconsistent ploidy (" << ploidy << ") genotype: " << string(vt[i]) << "\n";
                return 1;
            }

            v.push_back(a == NA ? 0 : a + 1);
            if (ploidy == 2)
                v.push_back(b == NA ? 0 : b + 1);
        }

        gt.dat.push_back(v);
    }

    gt.dist.assign(gt.loc.size(), 0.0);
    gt.ploidy = ploidy;

    return 0;
}

int write_vcf(const Genotype &gt, const string &filename, bool force_diploid)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = gt.loc.size(), n = gt.ind.size();
    bool haploid = gt.ploidy == 1;

    ofs << "##fileformat=VCFv4.2\n";
    ofs << "##filedate=" << datetime_now() << "\n";

    ofs << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if ( ! gt.ind.empty() ) {
        ofs << "\tFORMAT";
        for (size_t i = 0; i < n; ++i)
            ofs << "\t" << gt.ind[i];
    }
    ofs << "\n";

    int p = std::numeric_limits<allele_t>::max() + 1;
    vector<string> gs(p);
    gs[0] = ".";
    for (int i = 1; i < p; ++i)
        gs[i] = std::to_string(i-1);

    string line;

    for (size_t j = 0; j < m; ++j) {
        line.clear();

        line.append(gt.chr[j]).append("\t");
        line.append(std::to_string(gt.pos[j])).append("\t");
        line.append(gt.loc[j]).append("\t");

        if ( ! gt.allele[j].empty() ) {
            line.append(gt.allele[j][0]).append("\t");
            auto na = gt.allele[j].size();
            if (na > 1) {
                for (size_t k = 1; k < na; ++k)
                    line.append(gt.allele[j][k]).append(",");
                line.pop_back();
            }
            else
                line.append(".");
        }
        else
            line.append(".\t.");

        line.append("\t.\t.\t.");

        if ( gt.ind.empty() ) {
            ofs << line << "\n";
            continue;
        }

        line.append("\tGT");

        if ( haploid ) {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i];
                line.append("\t").append(gs[a]);
                if ( force_diploid )
                    line.append("/").append(gs[a]);
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                auto a = gt.dat[j][i*2], b = gt.dat[j][i*2+1];
                line.append("\t").append(gs[a]).append("/").append(gs[b]);
            }
        }

        ofs << line << "\n";
    }

    return 0;
}
