#pragma once

#include <string>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

inline void ensureMiscDir() {
#ifdef _WIN32
    _mkdir("misc");
#else
    mkdir("misc", 0755);
#endif
}

inline std::string miscFile(std::string const &name) {
    ensureMiscDir();
    return "misc/" + name;
}
