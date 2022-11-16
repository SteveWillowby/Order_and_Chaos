#include<stdexcept>
#include<string>
#include<unordered_map>

#ifndef SYM__COLORING_H
#define SYM__COLORING_H

template<class T> class Coloring {
public:
    Coloring();

    virtual void set(const T& elt, int color);
    virtual void erase(const T& elt);
    virtual int operator[](const T& elt) const;

protected:
    std::unordered_map<T, int> coloring;
};

/////////////////////////////// The Code //////////////////////////////

template<class T> Coloring<T>::Coloring() {
    coloring = std::unordered_map<T, int>();
}

template<class T> void Coloring<T>::set(const T& elt, int color) {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        coloring.insert(std::pair<T, int>(elt, color));
    } else {
        coloring[elt] = color;
    }
}


template<class T> void Coloring<T>::erase(const T& elt) {
    coloring.erase(elt);
}

template<class T> int Coloring<T>::operator[](const T& elt) const {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        throw std::range_error("Error. Element not in coloring.");
    }
    return *itr;
}

#endif
