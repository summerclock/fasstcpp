#include "Global.h"

namespace fasst {

// Global singleton instance
std::unique_ptr<GlobalData> g_data = std::make_unique<fasst::GlobalData>();

GlobalData::GlobalData() = default;
GlobalData::~GlobalData() = default;

} // namespace fasst
