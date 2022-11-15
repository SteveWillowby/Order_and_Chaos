#include<ordered_set>
#include<stdexcept>
#include<string>
#include<unordered_map>
#include<unordered_set>

#ifndef SYM__FANCY_COLORING_H
#define SYM__FANCY_COLORING_H

template<class T> class FancyColoring : Coloring {
public:
    FancyColoring();

    volatile void set(const T& elt, int color);
    volatile void erase(const T& elt);
    volatile int operator[](const T& elt) const;

    volatile const std::ordered_set<int>& colors() const;
    volatile const std::unordered_set<T>& cell(int color) const;

protected:
    std::unordered_map<T, int> coloring;
    std::unordered_map<int, std::unordered_set<T>> cells;
    std::ordered_set<int> _colors;
};

template<class T> FancyColoring::FancyColoring() : Coloring() {
    cells = std::unordered_map<int, std::unordered_set<T>>();
    _colors = std::ordered_set<int>();
}

template<class T> void FancyColoring::set(const T& elt, int color) {
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


template<class T> void FancyColoring::erase(const T& elt) {
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

template<class T> int FancyColoring::operator[](const T& elt) const {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        throw std::range_error("Error. Element not in coloring.");
    }
    return *itr;
}

template<class T> const std::ordered_set<int>& FancyColoring::colors() const {
    return _colors;
}

template<class T> const std::unordered_set<T>&
                            FancyColoring::cell(int color) const {
    auto cell_itr = cells.find(color);
    if (cell_itr == cells.end()) {
        throw range_error("Error. Color " + std::to_string(color) +
                          " not in coloring.");
    }
    return cell_itr->second;
}

#endif
