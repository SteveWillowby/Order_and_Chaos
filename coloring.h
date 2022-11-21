#include<set>
#include<stdexcept>
#include<string>
#include<unordered_map>
#include<unordered_set>

#ifndef SYM__COLORING_H
#define SYM__COLORING_H

// TODO: Add support for custom T hash function
template<class T> class Coloring {
public:
    Coloring();

    virtual void set(const T& elt, int color);
    virtual void erase(const T& elt);
    virtual int operator[](const T& elt) const;

    virtual const std::set<int>& colors() const;
    virtual const std::unordered_set<T>& cell(int color) const;

protected:
    std::unordered_map<T, int> coloring;
    std::unordered_map<int, std::unordered_set<T>> cells;
    std::set<int> _colors;
};

template<class T> Coloring<T>::Coloring() : Coloring() {
    cells = std::unordered_map<int, std::unordered_set<T>>();
    _colors = std::set<int>();
}

template<class T> void Coloring<T>::set(const T& elt, int color) {
    auto cell_itr = cells.find(color);
    if (cell_itr == cells.end()) {
        cells[color] = std::unordered_set<T>(elt);
        _colors.insert(color);
    } else {
        cell_itr->second.insert(elt);
    }

    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        coloring.insert(std::pair<T, int>(elt, color));
    } else {
        int old_color = itr->second;
        coloring[elt] = color;
        if (cells[old_color].size() == 1) {
            cells.erase(old_color);
            _colors.erase(old_color);
        } else {
            cells[old_color].erase(elt);
        }
    }
}


template<class T> void Coloring<T>::erase(const T& elt) {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        return;
    }
    int color = itr->second;
    coloring.erase(itr);

    cells[color].erase(elt);
    if (cells[color].size() == 0) {
        cells.erase(color);
        _colors.erase(color);
    }
}

template<class T> int Coloring<T>::operator[](const T& elt) const {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        throw std::range_error("Error. Element not in coloring.");
    }
    return *itr;
}

template<class T> const std::set<int>& Coloring<T>::colors() const {
    return _colors;
}

template<class T> const std::unordered_set<T>&
                            Coloring<T>::cell(int color) const {
    auto cell_itr = cells.find(color);
    if (cell_itr == cells.end()) {
        throw std::range_error("Error. Color " + std::to_string(color) +
                               " not in coloring.");
    }
    return cell_itr->second;
}

#endif
