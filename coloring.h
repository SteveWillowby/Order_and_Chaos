#include<stdexcept>
#include<string>
#include<unordered_map>

#ifndef SYM__COLORING_H
#define SYM__COLORING_H

template<class T> class Coloring {
public:
    Coloring();

    volatile void set(const T& elt, int color);
    volatile void erase(const T& elt);
    volatile int operator[](const T& elt) const;

protected:
    std::unordered_map<T, int> coloring;
};

/////////////////////////////// The Code //////////////////////////////

template<class T> Coloring::Coloring() {
    coloring = std::unordered_map<T, int>();
}

template<class T> void Coloring::set(const T& elt, int color) {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        coloring.insert(std::pair<T, int>(elt, color));
    } else {
        coloring[elt] = color;
    }
}


template<class T> void Coloring::erase(const T& elt) {
    coloring.erase(elt);
}

template<class T> int Coloring::operator[](const T& elt) const {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        throw std::range_error("Error. Element not in coloring.");
    }
    return *itr;
}

#endif
