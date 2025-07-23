#include "ActInputParser.h"
#include "ActModularDetector.h"

void checkActions(){
    ActRoot::InputParser in {"./configs/detector.conf"};
    ActRoot::ModularDetector mod;
    mod.ReadConfiguration(in.GetBlock("Modular"));
    mod.GetParameters()->Print();
}