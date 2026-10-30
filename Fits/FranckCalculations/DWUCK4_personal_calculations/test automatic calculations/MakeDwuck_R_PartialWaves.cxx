// =============================================================================
//  MakeDwuckConv.C  --  DWUCK4 convergence test at FIXED Ex:
//                       scan of the integration radius and of the maximum
//                       partial wave.
//
//  Same philosophy as MakeDwuck.C: your working script (g.s. case) is the
//  template and is edited IN PLACE, column by column, so every character that
//  is not explicitly changed stays byte-identical.
//
//  What changes between runs (Ex, OMP, Q, Sn, normalisation are the same in all):
//     - integration radius R (fm) : line 5, cols 1-3  (I3, the "30" of " 30  1  1  3")
//     - maximum partial wave Lmax : line 6, cols 17-24 (F8.0, the "40.0")
//     - tabulated distorted-wave potentials (the four "300." blocks): resized to
//       N = R/DR points (300 = 30 fm / 0.1 fm in your template), see kResizeTables.
//
//  Usage (ROOT): PhysOMP.cpp must be compiled/loaded first (it is reused, not copied)
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuckConv.C+'                              // defaults, runs
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuckConv.C+(4,0.,{20,30,40},{30,40,50})'  // 2*Jf, Ex, R list, Lmax list
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuckConv.C+(4,0.,{30},{20,30,40,50,60})'  // only Lmax scan (R fixed)
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuckConv.C+(4,0.,{},{},false)'            // only write, don't run
//
//  Output:
//     <kOutRoot>/conv_0p3-2_2+_gs/R30_L40/run.sh  -> <NOM>.res, fort.7, fort7_2col.dat, xs_dw4.n, run.log
//     <kOutRoot>/conv_0p3-2_2+_gs/convergence.dat -> summary table (read from each xs_dw4.n)
//  Each (R, Lmax) has its own folder, so the fixed names written by dwuck4 / normdw4
//  never collide.  Calling with run=false after running by hand regenerates only the summary
//  (it also rewrites the run.sh files, which is harmless).
// =============================================================================
#include "PhysOMP.h" // your existing OMP classes (PhysOMP::KoningDelaroche)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dwc
{

// ----------------------------- USER CONFIG -----------------------------------
const std::string kTemplate = "template_dwuck.sh"; // your original g.s. script
const std::string kBase = "Li7_dp";                // NOM = kBase_<tag>_R<R>_L<L>_kSuffix
const std::string kSuffix = "ADKDKD_LEA";

const std::string kOutRoot = "."; // folders are created here

// The ONE final state of this study (default values; can be passed as arguments)
const int kTwoJf = 4;      // 2*Jf of the final state (4 -> 2)
const double kEx = 3.0;    // MeV

// Scan lists (default values; can be passed as arguments).  The largest R and the
// largest Lmax are taken as the reference for the relative differences.
const std::vector<int> kRList = {30, 35, 40, 45, 50, 55, 60}; // fm   (integer: the field is I3)
const std::vector<int> kLList = {40, 45, 50, 55, 60, 65, 70};     // Lmax

// Resize the four tabulated-potential blocks to N = R/DR points in every run.
// (Your template has 300 points = 30 fm / DR=0.1.)  Set to false if DWUCK4 does not
// need the table length to match the integration radius.
const bool kResizeTables = true;

// Transferred neutron (written in the quantum-number line and in line 5)
const int kN = 0;
const int kL = 1;     // p
const int kTwoJ = 3;  // 2*j  (3 -> p3/2)   -> same as the template: " 30  1  1  3" and "0. 1. 3. 1. 50. 0."
const int kTwoS = 1;
const int kTwoJi = 3;
const int kParityI = -1;

// Title card (line 3) relabelling
const bool kRelabelTitle = true;
const std::string kTitleOld = "8Li(gs) 2+";
const std::string kTitleNucl = "8Li";

// Template line numbers (1-based, as in your editor)
const int kLnTransfer = 5; // " 30  1  1  3"
const int kLnRmax = 5;     // integration radius, cols 1-3 (I3)
const int kColRmax = 1;
const int kLnLmax = 6;     // "     0.1     0.0    40.0          0.7457": DR, ?, Lmax, ?, finite range
const int kColLmax = 17;   //   Lmax = 3rd field (F8.0)
const int kLnDR = 6;       //   DR   = 1st field (cols 1-8), only read
const int kLnQ = 257;
const int kLnOMP1 = 258;
const int kLnOMP2 = 259;
const int kLnOMP3 = 260;
const int kLnBound = 261;
const int kLnQN = 264;
const int kLnNorm = 278;

const int kDecEx = 3;
const int kDecOMP = 3;
const char* kNormFmt = "%.5f";

// Kinematics used only to hand the exit-channel energy to your OMP functions
const double kEdLab = 14.532, kMd = 2.0136, kMt = 7.0144, kMp = 1.0073, kMB = 8.0208;

// ------------------------- YOUR PHYSICS (same as MakeDwuck.C) ------------------
struct DwuckOMP
{
    double V, rV, aV;
    double W, rW, aW;
    double Wd, rWd, aWd;
    double Vso, rVso, aVso;
    double Wso, rWso, aWso;
};

const int kResZ = 3, kResA = 8; // 8Li

// Thin subclass, only to read the parameters the existing class already computed.
// Conversion to DWUCK4: volume -> sign change; surface -> x4; spin-orbit -> sign change and x4 on (L.S / 2)
struct KDView : public PhysOMP::KoningDelaroche
{
    KDView(int Z, int A, double Elab) : PhysOMP::KoningDelaroche(Z, A, Elab) {}
    DwuckOMP ToDwuck() const
    {
        DwuckOMP p;
        p.V = -fVr;
        p.rV = frv;
        p.aV = fav;
        p.W = -fWv;
        p.rW = frw;
        p.aW = faw;
        p.Wd = 4. * fWs;
        p.rWd = frs;
        p.aWd = fas;
        p.Vso = -4. * (fVso / 2);
        p.rVso = frvso;
        p.aVso = favso;
        p.Wso = -4. * (fWso / 2);
        p.rWso = frwso;
        p.aWso = fawso;
        return p;
    }
};

DwuckOMP ExitOMP(double /*Ex*/, double Ep)
{
    KDView kd(kResZ, kResA, Ep);
    return kd.ToDwuck();
}

// N = (125.2)^2 / 10^3 * (2Jf+1) / (2Ji+1) / (2j+1)
double Normalization(int twoJf)
{
    return 125.2 * 125.2 / 1000. * (twoJf + 1.) / (kTwoJi + 1.) / (kTwoJ + 1.);
}

// --------------------------------- helpers -----------------------------------
const double kPi = 3.14159265358979323846;

double ExitEnergyEquiv(double ex, double qgs)
{
    double ecmIn = kEdLab * kMt / (kMd + kMt);
    double ecmOut = ecmIn + (qgs - ex);
    return ecmOut * (kMp + kMB) / kMB;
}

typedef std::vector<std::string> VS;
struct Lines
{
    VS v;
    std::string& at(int n) { return v.at(n - 1); } // 1-based access
    const std::string& at(int n) const { return v.at(n - 1); }
};

std::string Trim(const std::string& s)
{
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

VS Tokens(const std::string& s)
{
    std::istringstream is(s);
    VS t;
    std::string w;
    while(is >> w)
        t.push_back(w);
    return t;
}

// Fortran-like Fw.d, right-justified.  dec==0 prints "3." like your file.
std::string Fmt(double v, int dec, size_t w)
{
    if(std::fabs(v) < 0.5 * std::pow(10., -dec))
        v = 0.;
    char buf[64];
    if(dec == 0)
        std::snprintf(buf, sizeof buf, "%.0f.", v);
    else
        std::snprintf(buf, sizeof buf, "%.*f", dec, v);
    std::string s(buf);
    if(s.size() > w)
        throw std::runtime_error("value does not fit in field: " + s);
    return std::string(w - s.size(), ' ') + s;
}

// Overwrite one F8 field of a template line in place (col is 1-based)
void Put(Lines& L, int ln, int col, double val, int dec, size_t w = 8)
{
    std::string& s = L.at(ln);
    size_t c0 = (size_t)(col - 1);
    if(s.size() < c0 + w)
        s.resize(c0 + w, ' ');
    std::string old = s.substr(c0, w);
    if(Fmt(std::atof(old.c_str()), dec, w) != old)
        std::cerr << "WARNING line " << ln << " cols " << col << "-" << col + (int)w - 1 << ": template field '" << old
                  << "' is not in the format being written\n";
    s.replace(c0, w, Fmt(val, dec, w));
}

// Overwrite one integer field (Iw) of a template line in place
void PutInt(Lines& L, int ln, int col, int val, size_t w = 3)
{
    std::string& s = L.at(ln);
    size_t c0 = (size_t)(col - 1);
    if(s.size() < c0 + w)
        s.resize(c0 + w, ' ');
    std::string old = s.substr(c0, w);
    if(std::to_string(std::atoi(old.c_str())) != Trim(old))
        std::cerr << "WARNING line " << ln << " cols " << col << "-" << col + (int)w - 1 << ": template field '" << old
                  << "' is not an integer\n";
    char b[32];
    std::snprintf(b, sizeof b, "%*d", (int)w, val);
    std::string f(b);
    if(f.size() > w)
        throw std::runtime_error("integer does not fit in field: " + f);
    s.replace(c0, w, f);
}

// Ex field: last field of the card (cols 73-80), with explicit sign
void PutEx(Lines& L, int ln, double val)
{
    std::string& s = L.at(ln);
    if(s.size() < 80)
        s.resize(80, ' ');
    std::string old = s.substr(72, 8);
    if(std::fabs(std::atof(old.c_str())) > 0.)
        std::cerr << "WARNING line " << ln << ": template Ex field is '" << old
                  << "', expected 0 (the template should be the g.s.)\n";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%+.*f", kDecEx, val);
    std::string f(buf);
    if(f.size() > 8)
        throw std::runtime_error("Ex does not fit in the field: " + f);
    s.replace(72, 8, std::string(8 - f.size(), ' ') + f);
}

void Require(Lines& L, int ln, std::initializer_list<const char*> needles)
{
    for(const char* n : needles)
        if(L.at(ln).find(n) == std::string::npos)
            throw std::runtime_error("template line " + std::to_string(ln) + " does not contain '" + n +
                                     "' (line numbers shifted?):\n  " + L.at(ln));
}

Lines ReadTemplate(const std::string& fn)
{
    std::ifstream in(fn.c_str());
    if(!in)
        throw std::runtime_error("cannot open template " + fn);
    Lines L;
    std::string s;
    while(std::getline(in, s))
    {
        if(!s.empty() && s[s.size() - 1] == '\r')
            s.erase(s.size() - 1);
        L.v.push_back(s);
    }
    return L;
}

std::string Tag(double ex)
{
    if(std::fabs(ex) < 1e-9)
        return "gs";
    char b[32];
    std::snprintf(b, sizeof b, "Ex%.3f", ex);
    return b;
}

std::string SpinLabel(int twoJ)
{
    return (twoJ % 2 == 0) ? std::to_string(twoJ / 2) : std::to_string(twoJ) + "-2";
}
std::string ShellLabel()
{
    const char* lett = "spdfghik";
    if(kL < 0 || kL > 7)
        throw std::runtime_error("l out of range for the shell label");
    return std::to_string(kN) + lett[kL] + SpinLabel(kTwoJ);
}
std::string StateLabel(int twoJf)
{
    const int pf = kParityI * ((kL % 2 == 0) ? 1 : -1);
    return SpinLabel(twoJf) + (pf > 0 ? "+" : "-");
}

// ------------------- tabulated distorted-wave potentials ----------------------
// Each block in the template is:   "<x> <x>"  card  /  "  NPT.   k."  /  NPT values, 5 per line (5E16.7)
// A block header is recognised as: a line with exactly two tokens that both end in '.',
// whose NEXT line holds an exponent ('E').  (The "+8.  -1.0" card before it does not qualify.)
struct TabBlock
{
    int hdr; // 1-based line of the "NPT.  k." card
    int npt;
};

std::vector<TabBlock> FindTables(const Lines& L)
{
    std::vector<TabBlock> r;
    for(int ln = 1; ln < (int)L.v.size(); ++ln)
    {
        VS t = Tokens(L.at(ln));
        if(t.size() != 2 || t[0].empty() || t[1].empty() || t[0][t[0].size() - 1] != '.' ||
           t[1][t[1].size() - 1] != '.')
            continue;
        if(L.at(ln + 1).find('E') == std::string::npos)
            continue;
        TabBlock b;
        b.hdr = ln;
        b.npt = (int)std::lround(std::atof(t[0].c_str()));
        r.push_back(b);
    }
    return r;
}

// Set every table block to nNew points: zero-pad if longer, cut if shorter.
// Done LAST in Build (it changes the number of lines).  Blocks are processed bottom-up.
int ResizeTables(Lines& L, int nNew)
{
    if(nNew <= 0 || nNew % 5 != 0)
        throw std::runtime_error("N = R/DR = " + std::to_string(nNew) + " must be a positive multiple of 5");
    std::vector<TabBlock> tb = FindTables(L);
    if(tb.empty())
        throw std::runtime_error("no tabulated-potential blocks found in the template");

    std::string zero;
    for(int i = 0; i < 5; ++i)
        zero += "   0.0000000E+00";

    int changed = 0;
    for(size_t k = tb.size(); k-- > 0;)
    {
        const TabBlock& b = tb[k];
        if(b.npt % 5 != 0 || b.hdr + b.npt / 5 > (int)L.v.size())
            throw std::runtime_error("table block at line " + std::to_string(b.hdr) + " is not 5 values per line");
        if(b.npt == nNew)
            continue;
        const int oldLn = b.npt / 5, newLn = nNew / 5;
        if(newLn > oldLn)
            L.v.insert(L.v.begin() + b.hdr + oldLn, newLn - oldLn, zero); // data = 0-based [hdr, hdr+oldLn)
        else
            L.v.erase(L.v.begin() + b.hdr + newLn, L.v.begin() + b.hdr + oldLn);
        Put(L, b.hdr, 1, nNew, 0);
        ++changed;
    }
    return changed;
}

// ------------------------------ main worker ----------------------------------
std::string DirFor(const std::string& baseDir, int rmax, int lmax)
{
    return baseDir + "/R" + std::to_string(rmax) + "_L" + std::to_string(lmax);
}

// Writes <baseDir>/R<rmax>_L<lmax>/run.sh and returns that folder
std::string Build(const Lines& tpl, double ex, int twoJf, int rmax, int lmax, const std::string& baseDir, bool verbose)
{
    Lines L = tpl;

    // ---- sanity checks: are we really editing the lines we think? -------------
    Require(L, 1, {"NOM="});
    Require(L, kLnLmax, {"0.7457"});
    Require(L, kLnQ, {"1.0073", "8.0208"});
    Require(L, kLnOMP1, {"+1."});
    Require(L, kLnOMP2, {"+2."});
    Require(L, kLnOMP3, {"-4."});
    Require(L, kLnBound, {"1.0087", "7.0144"});
    Require(L, kLnQN, {"50."});
    if(Trim(L.at(kLnNorm - 1)) != "xs_dw4.n")
        throw std::runtime_error("line before the normalisation is not 'xs_dw4.n'");
    if(L.at(kLnTransfer).size() < 12)
        throw std::runtime_error("line 5 is too short");

    const double qgs = std::atof(L.at(kLnQ).substr(0, 8).c_str());
    const double dr = std::atof(L.at(kLnDR).substr(0, 8).c_str());
    if(dr <= 0.)
        throw std::runtime_error("DR read from line " + std::to_string(kLnDR) + " is not positive");

    const std::string tag = Tag(ex);
    const std::string nom =
        kBase + "_" + tag + "_R" + std::to_string(rmax) + "_L" + std::to_string(lmax) + "_" + kSuffix;

    // ---- line 1: unique output name ------------------------------------------
    L.at(1) = "NOM=" + nom;

    // ---- line 3: title card ---------------------------------------------------
    if(kRelabelTitle)
    {
        std::string& t = L.at(3);
        size_t pos = t.find(kTitleOld);
        char e[32];
        std::snprintf(e, sizeof e, "%.3f", ex);
        std::string nw =
            kTitleNucl + "(" + (std::fabs(ex) < 1e-9 ? std::string("gs") : std::string(e)) + ") " + StateLabel(twoJf);
        if(pos == std::string::npos)
            std::cerr << "WARNING: '" << kTitleOld << "' not found in line 3, title left as is\n";
        else if(t.size() - kTitleOld.size() + nw.size() > 80)
            std::cerr << "WARNING: new title would exceed 80 columns, title left as is\n";
        else
            t.replace(pos, kTitleOld.size(), nw);
    }

    // ---- line 5: (R, 1, L, 2J) -----------------------------------------------
    PutInt(L, kLnRmax, kColRmax, rmax);
    {
        std::string& s5 = L.at(kLnTransfer);
        char b[8];
        std::snprintf(b, sizeof b, "%3d", kL);
        s5.replace(6, 3, b);
        std::snprintf(b, sizeof b, "%3d", kTwoJ);
        s5.replace(9, 3, b);
    }

    // ---- line 6: Lmax ---------------------------------------------------------
    Put(L, kLnLmax, kColLmax, lmax, 1);

    // ---- lines 257 / 261: Ex ---------------------------------------------------
    PutEx(L, kLnQ, -ex);
    PutEx(L, kLnBound, +ex);

    // ---- lines 258-260: exit-channel OMP --------------------------------------
    DwuckOMP p = ExitOMP(ex, ExitEnergyEquiv(ex, qgs));
    const int d = kDecOMP;
    Put(L, kLnOMP1, 9, p.V, d);
    Put(L, kLnOMP1, 17, p.rV, d);
    Put(L, kLnOMP1, 25, p.aV, d);
    Put(L, kLnOMP1, 41, p.W, d);
    Put(L, kLnOMP1, 49, p.rW, d);
    Put(L, kLnOMP1, 57, p.aW, d);
    Put(L, kLnOMP2, 41, p.Wd, d);
    Put(L, kLnOMP2, 49, p.rWd, d);
    Put(L, kLnOMP2, 57, p.aWd, d);
    Put(L, kLnOMP3, 9, p.Vso, d);
    Put(L, kLnOMP3, 17, p.rVso, d);
    Put(L, kLnOMP3, 25, p.aVso, d);
    Put(L, kLnOMP3, 41, p.Wso, d);
    Put(L, kLnOMP3, 49, p.rWso, d);
    Put(L, kLnOMP3, 57, p.aWso, d);

    // ---- quantum-number line --------------------------------------------------
    Put(L, kLnQN, 1, kN, 0);
    Put(L, kLnQN, 9, kL, 0);
    Put(L, kLnQN, 17, kTwoJ, 0);
    Put(L, kLnQN, 25, kTwoS, 0);

    // ---- normalisation for normdw4 -------------------------------------------
    {
        char b[64];
        std::snprintf(b, sizeof b, kNormFmt, Normalization(twoJf));
        L.at(kLnNorm) = b;
    }

    // ---- show what changed (BEFORE resizing tables: that shifts line numbers) --
    const int nPts = (int)std::lround(rmax / dr);
    std::printf("=== Ex = %.3f MeV, R = %d fm, Lmax = %d (N = %d pts)\n", ex, rmax, lmax, nPts);
    for(size_t i = 0; i < L.v.size(); ++i)
        if(L.v[i] != tpl.v[i] && (verbose || (int)i + 1 == kLnRmax || (int)i + 1 == kLnLmax))
            std::printf("%4zu|%s\n", i + 1, L.v[i].c_str());

    // ---- tabulated potentials: N = R/DR points (LAST: it changes the line count) ----
    if(kResizeTables)
    {
        int n = ResizeTables(L, nPts);
        if(n > 0)
            std::printf("     (%d tabulated-potential blocks resized to %d points)\n", n, nPts);
    }

    // ---- write ------------------------------------------------------------------
    const std::string dir = DirFor(baseDir, rmax, lmax);
    if(std::system(("mkdir -p '" + dir + "'").c_str()) != 0)
        throw std::runtime_error("cannot create " + dir);
    {
        std::ofstream out((dir + "/run.sh").c_str());
        if(!out)
            throw std::runtime_error("cannot write " + dir + "/run.sh");
        for(size_t i = 0; i < L.v.size(); ++i)
            out << L.v[i] << "\n";
    }
    return dir;
}

// ------------------------------- summary --------------------------------------
struct Row
{
    int R, L;
    std::string dir;
    bool ok;
    double pk, thpk, integ; // peak xs, its angle, 2pi * int xs sin(theta) dtheta
};

// Reads "angle xs" pairs from xs_dw4.n (lines that do not parse are skipped)
bool Analyse(Row& r)
{
    std::ifstream in((r.dir + "/xs_dw4.n").c_str());
    if(!in)
        return false;
    std::vector<double> th, xs;
    std::string s;
    while(std::getline(in, s))
    {
        double a, b;
        if(std::sscanf(s.c_str(), "%lf %lf", &a, &b) == 2)
        {
            th.push_back(a);
            xs.push_back(b);
        }
    }
    if(th.size() < 3)
        return false;
    r.pk = xs[0];
    r.thpk = th[0];
    r.integ = 0.;
    for(size_t i = 0; i < th.size(); ++i)
    {
        if(xs[i] > r.pk)
        {
            r.pk = xs[i];
            r.thpk = th[i];
        }
        if(i + 1 < th.size())
        {
            const double a0 = th[i] * kPi / 180., a1 = th[i + 1] * kPi / 180.;
            r.integ += 0.5 * (xs[i] * std::sin(a0) + xs[i + 1] * std::sin(a1)) * (a1 - a0);
        }
    }
    r.integ *= 2. * kPi;
    return true;
}

void Summary(std::vector<Row>& rows, const std::string& baseDir)
{
    int rRef = 0, lRef = 0;
    for(size_t i = 0; i < rows.size(); ++i)
    {
        rRef = std::max(rRef, rows[i].R);
        lRef = std::max(lRef, rows[i].L);
    }
    const Row* ref = nullptr;
    for(size_t i = 0; i < rows.size(); ++i)
    {
        rows[i].ok = Analyse(rows[i]);
        if(rows[i].ok && rows[i].R == rRef && rows[i].L == lRef)
            ref = &rows[i];
    }

    std::ofstream out((baseDir + "/convergence.dat").c_str());
    char b[256];
    std::snprintf(b, sizeof b, "# reference: R=%d fm, Lmax=%d%s", rRef, lRef, ref ? "" : "   (NOT available: no differences)");
    out << b << "\n# R[fm]  Lmax   sigma_max     theta_max   sigma_int    d_max[%]   d_int[%]\n";
    std::cout << "\n" << b << "\n  R[fm]  Lmax   sigma_max   theta_max   sigma_int   d_max[%]  d_int[%]\n";
    for(size_t i = 0; i < rows.size(); ++i)
    {
        const Row& r = rows[i];
        if(!r.ok)
        {
            std::snprintf(b, sizeof b, "%7d %5d   (no xs_dw4.n yet / not readable)", r.R, r.L);
        }
        else
        {
            const double dm = ref ? 100. * (r.pk - ref->pk) / ref->pk : 0.;
            const double di = ref ? 100. * (r.integ - ref->integ) / ref->integ : 0.;
            std::snprintf(b, sizeof b, "%7d %5d  %12.5E %8.1f  %12.5E %9.3f %9.3f", r.R, r.L, r.pk, r.thpk, r.integ,
                          dm, di);
        }
        out << b << "\n";
        std::cout << b << "\n";
    }
    std::cout << "Summary written to " << baseDir << "/convergence.dat" << std::endl;
}

} // namespace dwc

// ROOT entry point (function name = file name).  One Jf, one Ex, lists of R and Lmax.
void MakeDwuck_R_PartialWaves(int twoJf = dwc::kTwoJf, double ex = dwc::kEx, const std::vector<int>& rIn = std::vector<int>(),
                   const std::vector<int>& lIn = std::vector<int>(), bool run = true)
{
    const std::vector<int> rList = rIn.empty() ? dwc::kRList : rIn;
    const std::vector<int> lList = lIn.empty() ? dwc::kLList : lIn;
    try
    {
        dwc::Lines tpl = dwc::ReadTemplate(dwc::kTemplate);
        const std::string baseDir =
            dwc::kOutRoot + "/conv_" + dwc::ShellLabel() + "_" + dwc::StateLabel(twoJf) + "_" + dwc::Tag(ex);
        std::cout << "Folder: " << baseDir << "  (" << rList.size() * lList.size() << " runs)" << std::endl;

        std::vector<dwc::Row> rows;
        bool first = true;
        for(size_t i = 0; i < rList.size(); ++i)
            for(size_t j = 0; j < lList.size(); ++j)
            {
                dwc::Row r;
                r.R = rList[i];
                r.L = lList[j];
                r.ok = false;
                r.pk = r.thpk = r.integ = 0.;
                r.dir = dwc::Build(tpl, ex, twoJf, r.R, r.L, baseDir, first);
                first = false;
                if(run)
                {
                    std::cout << "running " << r.dir << "/run.sh ..." << std::endl;
                    int rc = std::system(("cd '" + r.dir + "' && bash run.sh > run.log 2>&1").c_str());
                    if(rc != 0)
                        std::cerr << "WARNING: run in " << r.dir << " returned " << rc << " (see run.log)\n";
                }
                rows.push_back(r);
            }
        dwc::Summary(rows, baseDir);
    }
    catch(const std::exception& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
}