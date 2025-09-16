#include "ActDataManager.h"

void checkFileExists(){
    ActRoot::DataManager dataman {"./configs/data.conf", ActRoot::ModeType::EReadTPC};
    auto in {dataman.GetInput()};
}