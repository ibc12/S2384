// Twofnr2fort4.C
// ---------------------------------------------------------------------------
// Convierte el potencial tabulado de salida de twofnr (bloque real + bloque
// imaginario) en el fichero fort.4 que lee fresco, reproduciendo lo que hace
// el programa Fortran "transpot".
//
// Ademas guarda, en la misma carpeta que el fort.4 de salida, el fichero
// "r V W" (formato de entrada de transpot) para poder comprobar la conversion.
//
// Uso en ROOT:
//   root -l -b -q 'Twofnr2fort4.C()'                                   // paths por defecto
//   root -l -b -q 'Twofnr2fort4.C("in.dat","out/fort.4")'
//   root -l -b -q 'Twofnr2fort4.C("in.dat","out/fort.4",true)'         // recorta ceros finales
//   root -l -b -q 'Twofnr2fort4.C("in.dat","out/fort.4",false,true)'   // usa ./transpot externo
//
// Modos:
//   useExternalTranspot = false : escribe fort.4 directamente (misma logica y
//                                 mismo formato 1P,6E12.4 que transpot).
//   useExternalTranspot = true  : ejecuta el binario transpot (compilado con
//                                 gfortran) sobre el fichero "r V W" guardado.
//
// trimZeros = true elimina los puntos finales con V = W = 0, tanto del
// fort.4 como del fichero "r V W".
// ---------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace twofnr2fort4_impl {

struct Block {
   std::string title;
   int nr = 0;
   double dr = 0, rfirst = 0;
   std::vector<double> v;
};

// Lee: linea de titulo, linea "nr dr rfirst" y nr valores (en cualquier
// numero de valores por linea).
inline bool readBlock(std::ifstream &f, Block &b)
{
   std::string l;
   if (!std::getline(f, b.title)) return false;
   if (!std::getline(f, l)) return false;
   std::istringstream hs(l);
   if (!(hs >> b.nr >> b.dr >> b.rfirst) || b.nr <= 0) return false;
   b.v.clear();
   b.v.reserve(b.nr);
   while ((int)b.v.size() < b.nr && std::getline(f, l)) {
      std::istringstream ls(l);
      double x;
      while ((int)b.v.size() < b.nr && ls >> x) b.v.push_back(x);
   }
   return (int)b.v.size() == b.nr;
}

inline std::string oneWord(std::string s)
{
   for (auto &c : s)
      if (c == ' ' || c == '\t') c = '_';
   return s;
}

// Carpeta que contiene un path ("a/b/fort.4" -> "a/b"; "fort.4" -> ".")
inline std::string dirOf(const std::string &p)
{
   size_t k = p.find_last_of('/');
   return (k == std::string::npos) ? std::string(".") : p.substr(0, k);
}

} // namespace twofnr2fort4_impl

int Twofnr2fort4(const char *infile = "./ADWA/prior/myAdiabatic/fort.4",
                 const char *outfile = "./ADWA/prior/twocols_changeWithMyAdiabatic/fort.4",
                 bool trimZeros = true,
                 bool useExternalTranspot = false,
                 const char *transpotExe = "./transpot",
                 const char *title = "Adiabatic_KD_8Li",
                 const char *colsName = "twocols.dat")
{
   using namespace twofnr2fort4_impl;

   // ---- 1. Lectura del fichero de twofnr -------------------------------
   std::ifstream fin(infile);
   if (!fin) {
      std::cerr << "ERROR: no puedo abrir " << infile << std::endl;
      return 1;
   }
   Block re, im;
   if (!readBlock(fin, re)) {
      std::cerr << "ERROR: bloque real mal formado en " << infile << std::endl;
      return 2;
   }
   bool cmplx = readBlock(fin, im); // si no hay 2o bloque, potencial real
   if (cmplx && (im.nr != re.nr || im.dr != re.dr || im.rfirst != re.rfirst)) {
      std::cerr << "ERROR: bloques real e imaginario incoherentes" << std::endl;
      return 3;
   }

   int n = re.nr;
   if (trimZeros) // quita puntos finales con V = W = 0
      while (n > 0 && re.v[n - 1] == 0.0 && (!cmplx || im.v[n - 1] == 0.0)) --n;
   if (n < 3) {
      std::cerr << "ERROR: muy pocos puntos (" << n << ")" << std::endl;
      return 4;
   }

   std::cout << "Leido: " << infile << "  nr=" << re.nr << "  dr=" << re.dr
             << "  rfirst=" << re.rfirst << (cmplx ? "  (complejo)" : "  (real)")
             << "  -> puntos escritos: " << n << std::endl;

   // ---- 2. Carpeta de salida y fichero "r V W" de comprobacion ----------
   const std::string outDir = dirOf(outfile);
   std::system(("mkdir -p '" + outDir + "'").c_str());
   const std::string colsPath = outDir + "/" + colsName;
   {
      std::ofstream fc(colsPath);
      if (!fc) {
         std::cerr << "ERROR: no puedo escribir " << colsPath << std::endl;
         return 7;
      }
      // Mismo formato que espera transpot: T/F, titulo (una palabra), r V [W]
      fc << (cmplx ? "T" : "F") << "\n" << oneWord(title) << "\n";
      char buf[160];
      for (int i = 0; i < n; ++i) {
         double r = re.rfirst + i * re.dr;
         if (cmplx) std::snprintf(buf, sizeof buf, "%12.6f %18.7E %18.7E\n", r, re.v[i], im.v[i]);
         else       std::snprintf(buf, sizeof buf, "%12.6f %18.7E\n", r, re.v[i]);
         fc << buf;
      }
   }
   std::cout << "Escrito " << colsPath << " (r V W)" << std::endl;

   // ---- 3a. Modo transpot externo ---------------------------------------
   if (useExternalTranspot) {
      std::string cmd = std::string(transpotExe) + " < '" + colsPath + "' > '" + outfile + "'";
      std::cout << "Ejecutando: " << cmd << std::endl;
      int rc = std::system(cmd.c_str());
      if (rc != 0) {
         std::cerr << "ERROR: fallo al ejecutar " << transpotExe
                   << " (compilalo con: gfortran transpot.f90 -o transpot)" << std::endl;
         return 5;
      }
      std::cout << "Escrito " << outfile << " via transpot externo" << std::endl;
      return 0;
   }

   // ---- 3b. Equivalente nativo de transpot -------------------------------
   // transpot usa dr = r(3)-r(2) y rfirst = r(1); aqui salen del propio header.
   FILE *fo = std::fopen(outfile, "w");
   if (!fo) {
      std::cerr << "ERROR: no puedo escribir " << outfile << std::endl;
      return 6;
   }
   std::fprintf(fo, " #Expressed as fresco %s\n", title);
   std::fprintf(fo, " %d %.6f %.6f\n", n, re.dr, re.rfirst);

   int col = 0; // formato Fortran (1P,6E12.4): 6 numeros por linea
   auto put = [&](double x) {
      std::fprintf(fo, "%12.4E", x);
      if (++col == 6) { std::fputc('\n', fo); col = 0; }
   };
   for (int i = 0; i < n; ++i) {
      put(re.v[i]);              // V_i
      if (cmplx) put(im.v[i]);   // W_i  (intercalados V1 W1 V2 W2 ...)
   }
   if (col) std::fputc('\n', fo);
   std::fclose(fo);

   std::cout << "Escrito " << outfile << std::endl;
   return 0;
}