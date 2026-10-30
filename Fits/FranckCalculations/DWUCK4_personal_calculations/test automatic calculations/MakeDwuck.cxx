// =============================================================================
//  MakeDwuck.C  --  generate DWUCK4 run scripts for a list of Ex values
//
//  It reads YOUR working script (the g.s. case) as a template and edits it
//  IN PLACE, column by column (the cards are 80-col, 8-wide fields), so every
//  character that is not explicitly changed stays byte-identical.
//
//  One call = ONE final-state spin Jf, for a list of Ex values.
//
//  Usage (ROOT): PhysOMP.cpp must be compiled/loaded first (it is reused, not copied)
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck.C+'                       // defaults below, runs
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck.C+(4,{0.,0.5})'           // 2*Jf, Ex list, runs
//     root -l -b -q -e '.L PhysOMP.cpp+' 'MakeDwuck.C+(4,{0.,0.5},false)'     // only write, don't run
//
//  Output: one folder per configuration+Jf, and inside it one folder per Ex:
//     <kOutRoot>/0p3-2_2+/gs/        run.sh  -> after running: <NOM>.res, fort.7,
//     <kOutRoot>/0p3-2_2+/Ex0.500/           fort7_2col.dat, xs_dw4.n (same names always), run.log
//  Each Ex has its own working directory, so the fixed names written by dwuck4 /
//  normdw4 never collide and the folders can be run in any order (or in parallel).
// =============================================================================
#include "PhysOMP.h" // your existing OMP classes (PhysOMP::KoningDelaroche)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dw
{

// ----------------------------- USER CONFIG -----------------------------------
const std::string kTemplate = "template_dwuck.sh"; // your original g.s. script
const std::string kBase = "Li7_dp";                // NOM = kBase_<tag>_kSuffix
const std::string kSuffix = "ADKDKD_LEA";

const std::string kOutRoot = "."; // folders are created here

// Final states of this run: ONE Jf, several Ex  (default values; can be passed as arguments)
const int kTwoJf = 4;                                // 2*Jf of the final state (4 -> 2)
const std::vector<double> kExList = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0}; // MeV (example list)

// Transferred neutron (written in the quantum-number line and in line 5)
const int kN = 0; // 1st field of the quantum-number line (as in your file)
const int kL = 0;
const int kTwoJ = 1;     // 2*j  (3 -> j=3/2)   -> line 5 last field + QN line
const int kTwoS = 1;     // 2*s of the neutron
const int kTwoJi = 3;    // 2*Ji of the target (7Li: 3/2)
const int kParityI = -1; // parity of the target (7Li: 3/2-); Pf = Pi * (-1)^l

// Title card (line 3) relabelling, so the .res header says which state it is
const bool kRelabelTitle = true;
const std::string kTitleOld = "8Li(gs) 2+"; // text present in the template
const std::string kTitleNucl = "8Li";

// Q(g.s.) and Sn(g.s.) stay as they are in the template (lines 257 / 261, cols 1-8):
// DWUCK adds Ex to them by itself.  The macro only READS them (sanity check, kinematics).
// What it writes is Ex, in the last field (cols 73-80) of those two lines.

// Template line numbers (1-based, as in your editor)
const int kLnTransfer = 5; // " 30  1  1  3"
const int kLnQ = 257;      // exit-channel card (Q value)
const int kLnOMP1 = 258;   // "+1."  volume
const int kLnOMP2 = 259;   // "+2."  surface
const int kLnOMP3 = 260;   // "-4."  spin-orbit
const int kLnBound = 261;  // bound neutron card (-Sn)
const int kLnQN = 264;     // n, l, 2J, 2s  (the line just before the "9")
const int kLnNorm = 278;   // normalisation read by normdw4

const int kDecEx = 3;          // decimals of the Ex field (cols 73-80) in lines 257 / 261
const int kDecOMP = 3;         // decimals for OMP parameters
const char* kNormFmt = "%.5f"; // format of the normalisation line

// Kinematics used only to hand the exit-channel energy to your OMP functions
const double kEdLab = 14.532, kMd = 2.0136, kMt = 7.0144, kMp = 1.0073, kMB = 8.0208;

// ------------------------- YOUR PHYSICS GOES HERE -----------------------------
// Exit-channel optical potential, in the DWUCK4 sign/units convention of the
// template.  Each group = (depth, radius, diffuseness).
struct DwuckOMP
{
    double V, rV, aV;       // line 258, fields 2-4
    double W, rW, aW;       // line 258, fields 6-8
    double Wd, rWd, aWd;    // line 259, fields 6-8
    double Vso, rVso, aVso; // line 260, fields 2-4
    double Wso, rWso, aWso; // line 260, fields 6-8
};

// Exit-channel nucleus (the residual that the proton scatters from)
const int kResZ = 3, kResA = 8; // 8Li

// ---- exit-channel OMP: REUSES PhysOMP::KoningDelaroche (no duplicated physics) --
// Thin subclass, only to read the parameters the existing class already computed.
// PhysOMP convention: FRESCO-like signs, spin-orbit as L.S product (FRESCO = value/2).
// Conversion to DWUCK4:
//   real/imag VOLUME : sign change
//   SURFACE          : x4
//   SPIN-ORBIT       : sign change and x4, applied to the FRESCO value (L.S / 2)
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

// Ep = proton lab-equivalent energy in the exit channel (see ExitEnergyEquiv)
DwuckOMP ExitOMP(double /*Ex*/, double Ep)
{
    KDView kd(kResZ, kResA, Ep);
    return kd.ToDwuck();
}

// Normalisation passed to normdw4 (line 278):
//   N = (125.2)^2 / 10^3 * (2Jf+1) / (2Ji+1) / (2j+1)
// (2j+1 uses the transferred-neutron j, i.e. kTwoJ.)
double Normalization(int twoJf)
{
    return 125.2 * 125.2 / 1000. * (twoJf + 1.) / (kTwoJi + 1.) / (kTwoJ + 1.);
}

// --------------------------------- helpers -----------------------------------
// Non-relativistic kinematics: proton lab-equivalent energy in the exit channel
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
};

std::string Trim(const std::string& s)
{
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}

// Fortran-like Fw.d, right-justified.  dec==0 prints "3." like your file.
std::string Fmt(double v, int dec, size_t w)
{
    if(std::fabs(v) < 0.5 * std::pow(10., -dec))
        v = 0.; // avoid "-0.000"
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

// Overwrite one field of a template line in place (col is 1-based)
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

// Ex field: last field of the card (cols 73-80), written with explicit sign
// ("-0.981" / "+0.981", like the "-0.00" / "+0.00" of the template).
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

// "3" -> "3-2" (half-integer, '/' not allowed in folder names), 4 -> "2"
std::string SpinLabel(int twoJ)
{
    return (twoJ % 2 == 0) ? std::to_string(twoJ / 2) : std::to_string(twoJ) + "-2";
}
// Transferred-neutron shell as written in the file: n, l letter, j  -> "0p3-2"
std::string ShellLabel()
{
    const char* lett = "spdfghik";
    if(kL < 0 || kL > 7)
        throw std::runtime_error("l out of range for the shell label");
    return std::to_string(kN) + lett[kL] + SpinLabel(kTwoJ);
}
// Final-state label: Jf + parity (Pf = Pi*(-1)^l) -> "2+"
std::string StateLabel(int twoJf)
{
    const int pf = kParityI * ((kL % 2 == 0) ? 1 : -1);
    return SpinLabel(twoJf) + (pf > 0 ? "+" : "-");
}

// ------------------------------ main worker ----------------------------------
// Writes <jfDir>/<Ex tag>/run.sh and returns that folder ('' if the Ex was skipped)
std::string Build(const Lines& tpl, double ex, int twoJf, const std::string& jfDir)
{
    Lines L = tpl;

    // ---- sanity checks: are we really editing the lines we think? -------------
    Require(L, 1, {"NOM="});
    Require(L, kLnQ, {"1.0073", "8.0208"}); // p + 8Li
    Require(L, kLnOMP1, {"+1."});
    Require(L, kLnOMP2, {"+2."});
    Require(L, kLnOMP3, {"-4."});
    Require(L, kLnBound, {"1.0087", "7.0144"}); // n + 7Li
    Require(L, kLnQN, {"50."});
    if(Trim(L.at(kLnNorm - 1)) != "xs_dw4.n")
        throw std::runtime_error("line before the normalisation is not 'xs_dw4.n'");
    if(L.at(kLnTransfer).size() < 12)
        throw std::runtime_error("line 5 is too short");

    // Q(g.s.) and Sn(g.s.) as written in the template (DWUCK adds Ex to both)
    const double qgs = std::atof(L.at(kLnQ).substr(0, 8).c_str());       // e.g. -0.192
    const double sngs = -std::atof(L.at(kLnBound).substr(0, 8).c_str()); // e.g.  2.033

    const std::string tag = Tag(ex);
    const std::string nom = kBase + "_" + tag + "_" + kSuffix;

    // ---- line 1: unique output name ------------------------------------------
    L.at(1) = "NOM=" + nom;

    // ---- line 3: title card, say which state this is --------------------------
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

    // ---- line 5 (30, 1, L, 2J): 3rd field (cols 7-9) = L, 4th (cols 10-12) = 2J ----
    {
        std::string& s5 = L.at(kLnTransfer);
        char b[8];
        std::snprintf(b, sizeof b, "%3d", kL);
        s5.replace(6, 3, b);
        std::snprintf(b, sizeof b, "%3d", kTwoJ);
        s5.replace(9, 3, b);
    }

    // ---- lines 257 / 261: Ex goes in the LAST field (cols 73-80): -Ex and +Ex ---------
    PutEx(L, kLnQ, -ex);     // line 257:  -Ex
    PutEx(L, kLnBound, +ex); // line 261:  +Ex

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

    // ---- quantum-number line: n, l, 2J, 2s (fields 1-4, F8.0) ------------------
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

    // ---- write: one folder per Ex, one script inside ----------------------------
    const std::string exDir = jfDir + "/" + tag;
    if(std::system(("mkdir -p '" + exDir + "'").c_str()) != 0)
        throw std::runtime_error("cannot create " + exDir);
    {
        std::ofstream out((exDir + "/run.sh").c_str());
        if(!out)
            throw std::runtime_error("cannot write " + exDir + "/run.sh");
        for(size_t i = 0; i < L.v.size(); ++i)
            out << L.v[i] << "\n";
    }

    // ---- show what changed (for a quick visual check) -------------------------
    std::printf("=== Ex = %.3f MeV, 2Jf = %d, norm = %.5f  ->  %s/run.sh\n", ex, twoJf, Normalization(twoJf),
                exDir.c_str());
    for(size_t i = 0; i < L.v.size(); ++i)
        if(L.v[i] != tpl.v[i])
            std::printf("%4zu|%s\n", i + 1, L.v[i].c_str());
    return exDir;
}

} // namespace dw

// ROOT entry point (function name = file name).  One Jf, a list of Ex.
void MakeDwuck(int twoJf = dw::kTwoJf, const std::vector<double>& exIn = std::vector<double>(), bool run = true)
{
    std::vector<double> exList = exIn.empty() ? dw::kExList : exIn;
    try
    {
        dw::Lines tpl = dw::ReadTemplate(dw::kTemplate);
        const std::string jfDir = dw::kOutRoot + "/" + dw::ShellLabel() + "_" + dw::StateLabel(twoJf);
        std::cout << "Folder: " << jfDir << std::endl;

        for(size_t i = 0; i < exList.size(); ++i)
        {
            std::string exDir = dw::Build(tpl, exList[i], twoJf, jfDir);
            if(exDir.empty() || !run)
                continue;
            std::cout << "running " << exDir << "/run.sh ..." << std::endl;
            int rc = std::system(("cd '" + exDir + "' && bash run.sh > run.log 2>&1").c_str());
            if(rc != 0)
                std::cerr << "WARNING: run in " << exDir << " returned " << rc << " (see run.log)\n";
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
}