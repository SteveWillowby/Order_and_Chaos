#include<iostream>

#include "nt_mode.cpp"


int main(void) {
    NTMode a = NTMode();
    NTMode b = NTMode(false, false);

    std::cout<<a.prefer_traces<<a.nauty_directed<<b.prefer_traces<<b.nauty_directed<<std::endl;

    return 0;
};
