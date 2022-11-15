#include<unordered_map>

#ifndef SYM__COLORING_H
#define SYM__COLORING_H

template<class T> class Coloring {
public:
    Coloring(const std::unordered_map& init_coloring);

    void erase(const T& elt);
    
    int& operator[](const T& elt);

protected:
    std::unordered_map<T, int> coloring;
};





template<class T> void Coloring::erase(const T& elt) {
    coloring.erase(elt);
}

template<class T> int& Coloring::operator[](const T& elt) {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        return *(coloring.insert(std::pair<T, int>(elt, 0)).first);
    }
    return *itr;
}

#endif