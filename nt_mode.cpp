#include "nt_mode.h"


NTMode::NTMode() : prefer_traces(true), nauty_directed(true) {}

NTMode::NTMode(const bool prefer_traces, const bool nauty_directed) : 
    prefer_traces(prefer_traces), nauty_directed(nauty_directed) {}
