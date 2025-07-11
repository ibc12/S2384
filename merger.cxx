#include "ActDetectorManager.h"
#include "ActMergerDetector.h"
// #include "ActOptions.h"
#include "ActTypes.h"


void merger()
{
    // Detman to init merger
    ActRoot::DetectorManager detman {ActRoot::ModeType::EMerge};
    detman.ReadDetectorFile("./configs/detector.conf");
    // Retrieve merger det from det man
    auto m {detman.GetDetectorAs<ActRoot::MergerDetector>()};
    // m->InitInputData(std::shared_ptr<TTree> tree)
}
