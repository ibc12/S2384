#ifndef Task0_h
#define Task0_h

#include "ActVTask.h"

#include <map>
#include <string>
#include <vector>

namespace ActAlgorithm
{
class Task0 : public VTask
{
public:
    std::map<int, std::vector<std::pair<std::string, int>>> fMap {
        // r0 -> 3,5
        {30, {{"r0", 3}, {"r0", 5}}},
        {31, {{"r0", 3}, {"r0", 5}}},
        {33, {{"r0", 3}, {"r0", 5}}},
        // r0 -> all (0..11)
        {32,
         {{"r0", 0},
          {"r0", 1},
          {"r0", 2},
          {"r0", 3},
          {"r0", 4},
          {"r0", 5},
          {"r0", 6},
          {"r0", 7},
          {"r0", 8},
          {"r0", 9},
          {"r0", 10},
          {"r0", 11}}},
        // r0 y l0 -> all (0..11)
        {34, {{"r0", 0}, {"r0", 1}, {"r0", 2},  {"r0", 3},  {"r0", 4}, {"r0", 5}, {"r0", 6},  {"r0", 7},
              {"r0", 8}, {"r0", 9}, {"r0", 10}, {"r0", 11}, {"l0", 0}, {"l0", 1}, {"l0", 2},  {"l0", 3},
              {"l0", 4}, {"l0", 5}, {"l0", 6},  {"l0", 7},  {"l0", 8}, {"l0", 9}, {"l0", 10}, {"l0", 11}}},
        // f0 -> 5 Already masked its calibration to 0
        // {36, {{"f0", 5}}},
        // {37, {{"f0", 5}}},
        // {38, {{"f0", 5}}},
        // {39, {{"f0", 5}}},
        // {40, {{"f0", 5}}},
        // {41, {{"f0", 5}}},
        // {42, {{"f0", 5}}},
        // {43, {{"f0", 5}}},
        // {44, {{"f0", 5}}},
        // r0 -> 2
        {116, {{"r0", 2}}},
        {117, {{"r0", 2}}}};

    Task0() : VTask("Task0") {}

    bool Run() override;
    void Print() override;
};
} // namespace ActAlgorithm

#endif // !Task0_h