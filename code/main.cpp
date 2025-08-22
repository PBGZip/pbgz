#include <iostream>
#include "pbgz_file.h"
#include "log/logger.h"

int main() {
    // Print version string
    PbgzFileHeader header;
    std::string versionStr = header.getVersionStr();
    LOG_STDOUT(LOG_INFO, "PBGZ File Version: %s", versionStr.c_str());
    return 0;
}