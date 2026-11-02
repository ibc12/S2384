// =============================================================================
//  MakeDwuck_R_PartialWaves.C -- DWUCK4 convergence test in RMAX and in the number
//  of partial waves (LMAX), for every DRF (step) folder, at FIXED Ex.
//  ALL combinations  step x RMAX x LMAX  are generated; ANY RMAX is accepted for ANY step.
//
//  Layout:
//     S01/   bash_template.sh      <- template already on the grid DRF=0.1
//            R40/L30/run.sh ... R120/L60/run.sh   (+ <NOM>.res, fort.7, xs_dw4.n, run.log)
//            R40/convergence.dat ...
//     S02/   bash_template.sh      <- DRF=0.2
//     S03/   bash_template.sh      <- DRF=0.3
//     convergence_all.dat          <- every (step,RMAX) at its largest LMAX, vs global reference
//
//  WHAT IS EDITED in each copy of the template:
//     - line 1   : NOM=...  (unique name)
//     - line 3   : title card relabelled (Ex, spin)
//     - line 5   : cols 1-3 = LMAX ; cols 7-12 = LTR, 2*JTR
//     - line 6   : cols 17-24 = RMAX                      <-- scanned variable (any value)
//     - tabulated potentials resized to N = ceil(RMAX/DRF) rounded up to a multiple of 5
//     - Ex, exit-channel OMP, bound-state quantum numbers, normalisation
//
//  WHAT IS ONLY CHECKED:
//     - line 6, DRF (cols 1-8) must equal the step of the folder.
//
//  Nothing is skipped: if RMAX/DRF is not an integer, or N > kKmax, a WARNING is printed
//  and the run is still generated (set kEnforceKmax=true to skip N > kKmax instead).
//
//  Usage (ROOT):
//   root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck_R_PartialWaves.C+'                      // defaults
//   root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck_R_PartialWaves.C+(4,3.,{30,40,50},{40,80,120})'
//   root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck_R_PartialWaves.C+(4,3.,{},{},false)'    // only write
//   root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck_R_PartialWaves.C+(4,3.,{},{},true,"S02")' // one step only
// =============================================================================
#include "PhysOMP.h"

#include <algorithm>
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

inline std::vector<int> Range(int a, int b, int s)
{
    std::vector<int> v;
    for(int x = a; x <= b; x += s)
        v.push_back(x);
    return v;
}

// ----------------------------- USER CONFIG -----------------------------------
struct Step
{
    std::string dir; // folder with the template; R<rmax>/L<lmax> sub-folders are created inside
    double step;     // DRF (fm)
};

const std::vector<Step> kSteps = {{"S01", 0.1}, {"S02", 0.2}, {"S03", 0.3}};

// Default scans (can be passed as arguments).  The largest L is the reference inside each (step,R).
const std::vector<int> kRList = Range(40, 120, 10); // 40, 50, ..., 120
const std::vector<int> kLList = Range(30, 60, 5);   // 30, 35, ..., 60

const std::string kTemplate = "bash_template.sh"; // template file name INSIDE each step folder
const std::string kBase = "Li7_dp";               // NOM = kBase_<tag>_<step>_R<R>_L<L>_kSuffix
const std::string kSuffix = "ADKDKD_LEA";

const std::string kOutRoot = "."; // the step folders live here

const int kTwoJf = 4;   // 2*Jf of the final state
const double kEx = 3.0; // MeV

// Resize the four tabulated-potential blocks to N points in every run.
const bool kResizeTables = true;
// Allow CUTTING the tables when RMAX is smaller than the template's range (same grid: r_i = i*DRF).
const bool kAllowCut = true;
// Nominal maximum number of radial points in DWUCK4 (KMAX in FORMF).
const int kKmax = 400;
// false: N > kKmax only prints a WARNING and the run is generated anyway.
// true : such (step,R) pairs are skipped.
const bool kEnforceKmax = false;

// Transferred neutron
const int kN = 0;
const int kL = 1;
const int kTwoJ = 3;
const int kTwoS = 1;
const int kTwoJi = 3;
const int kParityI = -1;

// Title card (line 3) relabelling
const bool kRelabelTitle = true;
const std::string kTitleOld = "8Li(gs) 2+";
const std::string kTitleNucl = "8Li";

// Template line numbers (1-based) BEFORE the potential tables
const int kLnCard3 = 5;
const int kColLmax = 1;
const int kLnCard4 = 6;
const int kColDRF = 1;
const int kColRmax = 17;
const int kDecRmax = 4; // RMAX written as F8.4 (e.g. "120.0000")

// Lines AFTER the tables, counted from the last data line of the last table
const int kOffQ = 2;
const int kOffOMP1 = 3;
const int kOffOMP2 = 4;
const int kOffOMP3 = 5;
const int kOffBound = 6;
const int kOffQN = 9;

const int kDecEx = 3;
const int kDecOMP = 3;
const char* kNormFmt = "%.5f";

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

double Kd()
{
    const double hbarc = 197.327, amu = 931.495;
    const double ecm = kEdLab * kMt / (kMd + kMt);
    const double mu = kMd * kMt / (kMd + kMt);
    return std::sqrt(2. * mu * amu * ecm) / hbarc;
}

// Number of table points for (R,step): ceil(R/step), rounded UP to a multiple of 5 (5E16.7 format).
int GridPoints(int R, double step)
{
    int n = (int)std::ceil(R / step - 1e-6);
    if(n < 5)
        n = 5;
    return ((n + 4) / 5) * 5;
}

// Text describing anything unusual about the (R,step) pair ("" if the grid is clean)
std::string GridWarn(int R, double step)
{
    std::string w;
    const double ratio = R / step;
    const int n = GridPoints(R, step);
    char b[128];
    if(std::fabs(ratio - std::lround(ratio)) > 1e-6)
    {
        std::snprintf(b, sizeof b, "RMAX/DRF = %.4f is not an integer (tables padded to N=%d) ", ratio, n);
        w += b;
    }
    else if(std::lround(ratio) != n)
    {
        std::snprintf(b, sizeof b, "RMAX/DRF = %ld is not a multiple of 5 (tables padded to N=%d) ", std::lround(ratio),
                      n);
        w += b;
    }
    if(n > kKmax)
    {
        std::snprintf(b, sizeof b, "N=%d > KMAX=%d ", n, kKmax);
        w += b;
    }
    return w;
}

bool Skip(int R, double step)
{
    return kEnforceKmax && GridPoints(R, step) > kKmax;
}

typedef std::vector<std::string> VS;
struct Lines
{
    VS v;
    std::string& at(int n) { return v.at(n - 1); }
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

double Field(const Lines& L, int ln, int col, size_t w = 8)
{
    std::string s = L.at(ln);
    const size_t c0 = (size_t)(col - 1);
    if(s.size() < c0 + w)
        s.resize(c0 + w, ' ');
    return std::atof(s.substr(c0, w).c_str());
}

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

// Same as Put but without the "format differs" warning (used for RMAX, which is expected to change)
void PutQuiet(Lines& L, int ln, int col, double val, int dec, size_t w = 8)
{
    std::string& s = L.at(ln);
    size_t c0 = (size_t)(col - 1);
    if(s.size() < c0 + w)
        s.resize(c0 + w, ' ');
    s.replace(c0, w, Fmt(val, dec, w));
}

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
    if(ln < 1 || ln > (int)L.v.size())
        throw std::runtime_error("template line " + std::to_string(ln) + " does not exist (template too short?)");
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
std::string StateLabel(int twoJf)
{
    const int pf = kParityI * ((kL % 2 == 0) ? 1 : -1);
    return SpinLabel(twoJf) + (pf > 0 ? "+" : "-");
}

// ------------------- tabulated distorted-wave potentials ----------------------
struct TabBlock
{
    int hdr;
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

struct Pos
{
    int q, omp1, omp2, omp3, bound, qn, norm;
};

Pos Locate(const Lines& L)
{
    const std::vector<TabBlock> tb = FindTables(L);
    if(tb.empty())
        throw std::runtime_error("no tabulated-potential blocks found in the template");
    const int e = tb.back().hdr + tb.back().npt / 5;
    Pos p;
    p.q = e + kOffQ;
    p.omp1 = e + kOffOMP1;
    p.omp2 = e + kOffOMP2;
    p.omp3 = e + kOffOMP3;
    p.bound = e + kOffBound;
    p.qn = e + kOffQN;
    p.norm = 0;
    for(int ln = p.qn + 1; ln < (int)L.v.size(); ++ln)
        if(Trim(L.at(ln)) == "xs_dw4.n")
        {
            p.norm = ln + 1;
            break;
        }
    if(p.norm == 0)
        throw std::runtime_error("line 'xs_dw4.n' (before the normalisation) not found after the tables");
    if(p.norm > (int)L.v.size())
        throw std::runtime_error("template ends right after 'xs_dw4.n': no normalisation line");
    return p;
}

// Set every table block to nNew points: zero-pad if longer, cut if shorter (if allowed).
int ResizeTables(Lines& L, int nNew)
{
    if(nNew <= 0 || nNew % 5 != 0)
        throw std::runtime_error("N = " + std::to_string(nNew) + " must be a positive multiple of 5");
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
            L.v.insert(L.v.begin() + b.hdr + oldLn, newLn - oldLn, zero);
        else
        {
            if(!kAllowCut)
            {
                for(int i = b.hdr + newLn; i < b.hdr + oldLn; ++i)
                {
                    const VS t = Tokens(L.v[i]);
                    for(size_t q = 0; q < t.size(); ++q)
                        if(std::atof(t[q].c_str()) != 0.)
                            throw std::runtime_error("table at line " + std::to_string(b.hdr) + " has " +
                                                     std::to_string(b.npt) +
                                                     " points with non-zero values beyond point " +
                                                     std::to_string(nNew) + " (kAllowCut=false)");
                }
            }
            L.v.erase(L.v.begin() + b.hdr + newLn, L.v.begin() + b.hdr + oldLn);
        }
        Put(L, b.hdr, 1, nNew, 0);
        ++changed;
    }
    return changed;
}

// ------------------------------ main worker ----------------------------------
std::string StepDir(const Step& s)
{
    return kOutRoot + "/" + s.dir;
}

std::string BaseDir(const Step& s, int R)
{
    return StepDir(s) + "/R" + std::to_string(R);
}

std::string DirFor(const Step& s, int R, int lmax)
{
    return BaseDir(s, R) + "/L" + std::to_string(lmax);
}

// Writes <step>/R<R>/L<lmax>/run.sh and returns that folder.
std::string Build(const Lines& tpl, const Step& s, int R, double ex, int twoJf, int lmax, bool verbose)
{
    Lines L = tpl;
    const Pos p = Locate(L);

    Require(L, 1, {"NOM="});
    Require(L, kLnCard4, {"0.7457"});
    Require(L, p.q, {"1.0073", "8.0208"});
    Require(L, p.omp1, {"+1."});
    Require(L, p.omp2, {"+2."});
    Require(L, p.omp3, {"-4."});
    Require(L, p.bound, {"1.0087", "7.0144"});
    Require(L, p.qn, {"50."});
    if(L.at(kLnCard3).size() < 12)
        throw std::runtime_error("line 5 (card set 3) is too short");

    // the template must be on the DRF of this folder (RMAX is NOT checked: it is written)
    const double drT = Field(L, kLnCard4, kColDRF);
    if(std::fabs(drT - s.step) > 1e-6)
        throw std::runtime_error("template of " + s.dir + " has DRF=" + std::to_string(drT) +
                                 " but this folder needs DRF=" + std::to_string(s.step) + " (re-grid the template)");
    const int nPts = GridPoints(R, s.step);

    const int lWarn = (int)(Kd() * R + 12.5);
    if(lmax + 1 > lWarn)
        std::cerr << "NOTE " << s.dir << " R" << R << " L" << lmax << ": ADWUCK4 should print warning 9905 (LMAX+1 > "
                  << lWarn << " = int(k_d*R+12.5)); k_d=" << Kd() << " fm^-1\n";

    const double qgs = Field(L, p.q, 1);

    const std::string tag = Tag(ex);
    const std::string nom =
        kBase + "_" + tag + "_" + s.dir + "_R" + std::to_string(R) + "_L" + std::to_string(lmax) + "_" + kSuffix;

    // ---- line 1 ----
    L.at(1) = "NOM=" + nom;

    // ---- line 3: title ----
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

    // ---- line 5: LMAX, NLTR, LTR, 2*JTR ----
    PutInt(L, kLnCard3, kColLmax, lmax);
    {
        std::string& s5 = L.at(kLnCard3);
        char b[8];
        std::snprintf(b, sizeof b, "%3d", kL);
        s5.replace(6, 3, b);
        std::snprintf(b, sizeof b, "%3d", kTwoJ);
        s5.replace(9, 3, b);
    }

    // ---- line 6: RMAX (any value) ----
    PutQuiet(L, kLnCard4, kColRmax, (double)R, kDecRmax);

    // ---- Ex fields ----
    PutEx(L, p.q, -ex);
    PutEx(L, p.bound, +ex);

    // ---- exit-channel OMP ----
    DwuckOMP o = ExitOMP(ex, ExitEnergyEquiv(ex, qgs));
    const int d = kDecOMP;
    Put(L, p.omp1, 9, o.V, d);
    Put(L, p.omp1, 17, o.rV, d);
    Put(L, p.omp1, 25, o.aV, d);
    Put(L, p.omp1, 41, o.W, d);
    Put(L, p.omp1, 49, o.rW, d);
    Put(L, p.omp1, 57, o.aW, d);
    Put(L, p.omp2, 41, o.Wd, d);
    Put(L, p.omp2, 49, o.rWd, d);
    Put(L, p.omp2, 57, o.aWd, d);
    Put(L, p.omp3, 9, o.Vso, d);
    Put(L, p.omp3, 17, o.rVso, d);
    Put(L, p.omp3, 25, o.aVso, d);
    Put(L, p.omp3, 41, o.Wso, d);
    Put(L, p.omp3, 49, o.rWso, d);
    Put(L, p.omp3, 57, o.aWso, d);

    // ---- bound-state quantum numbers ----
    Put(L, p.qn, 1, kN, 0);
    Put(L, p.qn, 9, kL, 0);
    Put(L, p.qn, 17, kTwoJ, 0);
    Put(L, p.qn, 25, kTwoS, 0);

    // ---- normalisation ----
    {
        char b[64];
        std::snprintf(b, sizeof b, kNormFmt, Normalization(twoJf));
        L.at(p.norm) = b;
    }

    std::printf("=== %s: Ex = %.3f MeV, RMAX = %d fm, DRF = %.3f fm (N = %d), LMAX = %d\n", s.dir.c_str(), ex, R,
                s.step, nPts, lmax);
    for(size_t i = 0; i < L.v.size(); ++i)
        if(L.v[i] != tpl.v[i] && (verbose || (int)i + 1 == kLnCard3 || (int)i + 1 == kLnCard4))
            std::printf("%4zu|%s\n", i + 1, L.v[i].c_str());

    // ---- tabulated potentials (LAST: it changes the line count) ----
    if(kResizeTables)
    {
        int n = ResizeTables(L, nPts);
        if(n > 0)
            std::printf("     (%d tabulated-potential blocks resized to %d points)\n", n, nPts);
    }

    const std::string dir = DirFor(s, R, lmax);
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
    std::string stepDir;
    double step;
    int R;
    int L;
    std::string dir;
    bool ok;
    double pk, thpk, integ;
};

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

// Per-(step,R) table (reference = largest L).
void Summary(std::vector<Row>& rows, const Step& s, int R)
{
    int lRef = 0;
    for(size_t i = 0; i < rows.size(); ++i)
        lRef = std::max(lRef, rows[i].L);
    const Row* ref = nullptr;
    for(size_t i = 0; i < rows.size(); ++i)
    {
        rows[i].ok = Analyse(rows[i]);
        if(rows[i].ok && rows[i].L == lRef)
            ref = &rows[i];
    }

    std::ofstream out((BaseDir(s, R) + "/convergence.dat").c_str());
    char b[256];
    std::snprintf(b, sizeof b, "# %s R%d: RMAX=%d fm, DRF=%.3f fm; reference Lmax=%d%s", s.dir.c_str(), R, R, s.step,
                  lRef, ref ? "" : "   (NOT available: no differences)");
    out << b << "\n# Lmax   sigma_max     theta_max   sigma_int    d_max[%]   d_int[%]\n";
    std::cout << "\n" << b << "\n  Lmax   sigma_max   theta_max   sigma_int   d_max[%]  d_int[%]\n";
    for(size_t i = 0; i < rows.size(); ++i)
    {
        const Row& r = rows[i];
        if(!r.ok)
            std::snprintf(b, sizeof b, "%6d   (no xs_dw4.n yet / not readable)", r.L);
        else
        {
            const double dm = ref ? 100. * (r.pk - ref->pk) / ref->pk : 0.;
            const double di = ref ? 100. * (r.integ - ref->integ) / ref->integ : 0.;
            std::snprintf(b, sizeof b, "%6d  %12.5E %8.1f  %12.5E %9.3f %9.3f", r.L, r.pk, r.thpk, r.integ, dm, di);
        }
        out << b << "\n";
        std::cout << b << "\n";
    }
    std::cout << "Summary written to " << BaseDir(s, R) << "/convergence.dat" << std::endl;
}

// ALL (step,R,L) rows vs a global reference: largest RMAX, then smallest step, then largest L
void WriteAll(const std::vector<Row>& all)
{
    int iref = -1;
    for(size_t i = 0; i < all.size(); ++i)
    {
        if(!all[i].ok)
            continue;
        if(iref < 0)
        {
            iref = (int)i;
            continue;
        }
        const Row& a = all[i];
        const Row& b = all[(size_t)iref];
        if(a.R > b.R || (a.R == b.R && (a.step < b.step - 1e-9 || (std::fabs(a.step - b.step) < 1e-9 && a.L > b.L))))
            iref = (int)i;
    }

    std::ofstream out((kOutRoot + "/convergence_all.dat").c_str());
    char b[256];
    if(iref >= 0)
        std::snprintf(b, sizeof b, "# all runs; reference: %s R=%d DRF=%.3f Lmax=%d", all[(size_t)iref].stepDir.c_str(),
                      all[(size_t)iref].R, all[(size_t)iref].step, all[(size_t)iref].L);
    else
        std::snprintf(b, sizeof b, "# all runs; reference: (none available)");
    out << b << "\n# step   DRF    RMAX  Lmax   sigma_max     sigma_int    d_max[%]   d_int[%]\n";
    std::cout << "\n" << b << "\n  step   DRF    RMAX  Lmax   sigma_max     sigma_int    d_max[%]  d_int[%]\n";
    for(size_t i = 0; i < all.size(); ++i)
    {
        const Row& r = all[i];
        if(!r.ok)
            std::snprintf(b, sizeof b, "%-6s %5.2f %6d %5d   (no data)", r.stepDir.c_str(), r.step, r.R, r.L);
        else
        {
            const double dm = iref >= 0 ? 100. * (r.pk - all[(size_t)iref].pk) / all[(size_t)iref].pk : 0.;
            const double di = iref >= 0 ? 100. * (r.integ - all[(size_t)iref].integ) / all[(size_t)iref].integ : 0.;
            std::snprintf(b, sizeof b, "%-6s %5.2f %6d %5d  %12.5E %12.5E %9.3f %9.3f", r.stepDir.c_str(), r.step, r.R,
                          r.L, r.pk, r.integ, dm, di);
        }
        out << b << "\n";
        std::cout << b << "\n";
    }
    std::cout << "Summary written to " << kOutRoot << "/convergence_all.dat" << std::endl;
}

} // namespace dwc

// ROOT entry point (function name = file name).
// lIn, rIn : lists of LMAX and RMAX (empty -> defaults kLList, kRList).  All combinations are run.
// only     : name of one step folder (e.g. "S02") to process just that one.
void MakeDwuck_R_PartialWaves(int twoJf = dwc::kTwoJf, double ex = dwc::kEx,
                              const std::vector<int>& lIn = std::vector<int>(),
                              const std::vector<int>& rIn = std::vector<int>(), bool run = true,
                              const std::string& only = "")
{
    const std::vector<int> lList = lIn.empty() ? dwc::kLList : lIn;
    const std::vector<int> rList = rIn.empty() ? dwc::kRList : rIn;
    try
    {
        std::vector<dwc::Row> all;
        int nStep = 0;
        for(size_t is = 0; is < dwc::kSteps.size(); ++is)
        {
            const dwc::Step& s = dwc::kSteps[is];
            if(!only.empty() && only != s.dir)
                continue;
            ++nStep;
            dwc::Lines tpl = dwc::ReadTemplate(dwc::StepDir(s) + "/" + dwc::kTemplate);

            for(size_t ir = 0; ir < rList.size(); ++ir)
            {
                const int R = rList[ir];
                const int n = dwc::GridPoints(R, s.step);
                const std::string w = dwc::GridWarn(R, s.step);
                if(dwc::Skip(R, s.step))
                {
                    std::cout << "SKIP " << s.dir << " (DRF=" << s.step << ") RMAX=" << R << ": " << w
                              << "(kEnforceKmax=true)" << std::endl;
                    continue;
                }
                if(!w.empty())
                    std::cout << "WARNING " << s.dir << " (DRF=" << s.step << ") RMAX=" << R << ": " << w << std::endl;
                std::cout << "Case " << s.dir << " R" << R << ": DRF=" << s.step << " fm, N=" << n << ", "
                          << lList.size() << " runs" << std::endl;

                std::vector<dwc::Row> rows;
                bool first = true;
                for(size_t j = 0; j < lList.size(); ++j)
                {
                    dwc::Row r;
                    r.stepDir = s.dir;
                    r.step = s.step;
                    r.R = R;
                    r.L = lList[j];
                    r.ok = false;
                    r.pk = r.thpk = r.integ = 0.;
                    r.dir = dwc::Build(tpl, s, R, ex, twoJf, r.L, first);
                    first = false;
                    if(run)
                    {
                        std::cout << "running " << r.dir << "/run.sh ..." << std::endl;
                        int rc =
                            std::system(("cd '" + r.dir + "' && rm -f fort.7 && bash run.sh > run.log 2>&1").c_str());
                        if(rc != 0)
                            std::cerr << "WARNING: run in " << r.dir << " returned " << rc << " (see run.log)\n";
                    }
                    rows.push_back(r);
                }
                dwc::Summary(rows, s, R); // fills rows[i].ok and writes convergence.dat
                all.insert(all.end(), rows.begin(), rows.end());
            }
        }
        if(nStep == 0)
            std::cerr << "ERROR: no step folder matches '" << only << "'" << std::endl;
        else if(!all.empty())
            dwc::WriteAll(all);
    }
    catch(const std::exception& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
}