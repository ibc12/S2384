#include "ActDataManager.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

// Estructura simple para guardar la info de clasificación
struct ReactionInfo
{
    std::string type;
    bool status;
};

// Combina (run, entry) en una única clave de 64 bits para poder usar unordered_map
using Key = long long;
inline Key MakeKey(int run, int entry)
{
    return (static_cast<long long>(run) << 32) | static_cast<unsigned int>(entry);
}

// Lee el CSV y para en cuanto 'type' esté vacío
std::unordered_map<Key, ReactionInfo> ReadReactionMap(const std::string& csvFile)
{
    std::unordered_map<Key, ReactionInfo> reactionMap;
    reactionMap.reserve(10000); // evita rehashes; ajusta si tienes muchas más entradas

    std::ifstream file(csvFile);
    if(!file.is_open())
    {
        std::cerr << "No se pudo abrir " << csvFile << std::endl;
        return reactionMap;
    }

    std::string line;
    std::getline(file, line); // saltar cabecera: run,entry,type,status

    while(std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;

        while(std::getline(ss, field, ','))
            fields.push_back(field);
        if(!line.empty() && line.back() == ',')
            fields.push_back("");

        if(fields.size() < 4 || fields[2].empty())
        {
            std::cout << "Fin de datos útiles (type vacío). Parando lectura.\n";
            break;
        }

        int run = std::stoi(fields[0]);
        int entry = std::stoi(fields[1]);
        std::string type = fields[2];
        bool status = (fields[3] == "True");

        reactionMap[MakeKey(run, entry)] = {type, status};
    }

    return reactionMap;
}

// Devuelve true si el fichero existe y se puede abrir
bool FileExists(const std::string& path)
{
    std::ifstream f(path);
    return f.good();
}

// Construye (o carga desde caché) el dataframe ya clasificado con las
// columnas "type" y "status". Si outputFile ya existe, simplemente lo lee;
// si no, filtra el chain completo con el mapa del CSV y lo guarda para la
// próxima vez.
ROOT::RDF::RNode
GetClassifiedDataFrame(const std::string& csvFile, const std::string& outputFile, const std::string& outputTree)
{
    if(FileExists(outputFile))
    {
        std::cout << "Cache encontrada, cargando dataframe clasificado desde: " << outputFile << std::endl;
        return ROOT::RDF::RNode(ROOT::RDataFrame(outputTree, outputFile));
    }

    std::cout << "No existe cache. Clasificando eventos desde cero...\n";

    // 1) Cargamos la tabla de clasificación
    auto reactionMap = ReadReactionMap(csvFile);
    std::cout << "Entradas cargadas en el mapa: " << reactionMap.size() << std::endl;

    // 2) Chain completo con todos los eventos
    ActRoot::DataManager dataManager {};
    dataManager.ReadDataFile("../../configs/data.conf");

    // For L1 enough with 4 runs (64, 67). For lat sils put at least 10
    dataManager.SetRuns(50, 69);
    auto chain {dataManager.GetChain(ActRoot::ModeType::EMerge)};
    auto df {ROOT::RDataFrame(*chain)};

    // 3) Una sola búsqueda en el mapa por evento
    auto dfInfo = df.Define("info",
                            [reactionMap](ActRoot::MergerData& m) -> const ReactionInfo*
                            {
                                auto key = MakeKey(m.fRun, m.fEntry);
                                auto it = reactionMap.find(key);
                                return (it != reactionMap.end()) ? &it->second : nullptr;
                            },
                            {"MergerData"});

    auto dfFiltered = dfInfo.Filter([](const ReactionInfo* info) { return info != nullptr; }, {"info"});

    auto dfClasificado = dfFiltered.Define("type", [](const ReactionInfo* info) { return info->type; }, {"info"})
                             .Define("status", [](const ReactionInfo* info) { return info->status; }, {"info"});

    // 4) Guardamos SOLO las columnas serializables (no "info", que es un puntero crudo)
    std::cout << "Guardando dataframe clasificado en: " << outputFile << std::endl;
    dfClasificado.Snapshot(outputTree, outputFile, {"MergerData", "type", "status"});

    // 5) Releemos desde disco: nodo limpio, sin depender de reactionMap en memoria
    return ROOT::RDF::RNode(ROOT::RDataFrame(outputTree, outputFile));
}

void FilterEventsByType()
{
    ROOT::EnableImplicitMT();

    const std::string csvFile = "./events_L1_all.csv";
    // Ajusta esta ruta a tu carpeta de outputs real.
    // Si borras este fichero, la próxima ejecución reclasificará desde cero.
    const std::string outputFile = "./Outputs/events_classified.root";
    const std::string outputTree = "classified";

    auto dfClasificado = GetClassifiedDataFrame(csvFile, outputFile, outputTree);

    // A partir de aquí, todo igual que antes: filtras por tipo/status
    auto dfBinaryTrue = dfClasificado.Filter("type == \"Binary\" && status == true");
    auto dfBinaryFalse = dfClasificado.Filter("type == \"Binary\" && status == false");
    auto dfHolesTrue = dfClasificado.Filter("type == \"Holes\"  && status == true");
    auto dfHolesFalse = dfClasificado.Filter("type == \"Holes\"  && status == false");
    auto dfBrokenTrue = dfClasificado.Filter("type == \"Broken\"  && status == true");
    auto dfBrokenFalse = dfClasificado.Filter("type == \"Broken\"  && status == false");

    std::cout << "Binary/True: " << *dfBinaryTrue.Count() << std::endl;
    std::cout << "Binary/False: " << *dfBinaryFalse.Count() << std::endl;
    std::cout << "Holes/True: " << *dfHolesTrue.Count() << std::endl;
    std::cout << "Holes/False: " << *dfHolesFalse.Count() << std::endl;
    std::cout << "Broken/True: " << *dfBrokenTrue.Count() << std::endl;
    std::cout << "Broken/False: " << *dfBrokenFalse.Count() << std::endl;

    // Plot theta and phi for light particle for each type
    auto hBinaryTrue =
        dfBinaryTrue.Histo2D({"theta_phi_binary_true", "Theta vs Phi (Binary True)", 60, 0, 180, 60, -180, 180},
                             "MergerData.fThetaLight", "MergerData.fPhiLight");
    auto hBinaryFalse =
        dfBinaryFalse.Histo2D({"theta_phi_binary_false", "Theta vs Phi (Binary False)", 60, 0, 180, 60, -180, 180},
                              "MergerData.fThetaLight", "MergerData.fPhiLight");
    auto hHolesTrue =
        dfHolesTrue.Histo2D({"theta_phi_holes_true", "Theta vs Phi (Holes True)", 60, 0, 180, 60, -180, 180},
                            "MergerData.fThetaLight", "MergerData.fPhiLight");
    auto hHolesFalse =
        dfHolesFalse.Histo2D({"theta_phi_holes_false", "Theta vs Phi (Holes False)", 60, 0, 180, 60, -180, 180},
                             "MergerData.fThetaLight", "MergerData.fPhiLight");
    auto hBrokenTrue =
        dfBrokenTrue.Histo2D({"theta_phi_broken_true", "Theta vs Phi (Broken True)", 60, 0, 180, 60, -180, 180},
                             "MergerData.fThetaLight", "MergerData.fPhiLight");
    auto hBrokenFalse =
        dfBrokenFalse.Histo2D({"theta_phi_broken_false", "Theta vs Phi (Broken False)", 60, 0, 180, 60, -180, 180},
                              "MergerData.fThetaLight", "MergerData.fPhiLight");

    auto* c = new TCanvas("c", "Theta vs Phi", 1200, 800);
    c->Divide(3, 2);
    c->cd(1);
    hBinaryTrue->DrawClone("COLZ");
    c->cd(4);
    hBinaryFalse->DrawClone("COLZ");
    c->cd(2);
    hHolesTrue->DrawClone("COLZ");
    c->cd(5);
    hHolesFalse->DrawClone("COLZ");
    c->cd(3);
    hBrokenTrue->DrawClone("COLZ");
    c->cd(6);
    hBrokenFalse->DrawClone("COLZ");
}