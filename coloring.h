#include<set>
#include<stdexcept>
#include<string>
#include<unordered_map>
#include<unordered_set>

#ifndef SCHENO__COLORING_H
#define SCHENO__COLORING_H

template<class T, class THash=std::hash<T>> class Coloring {
public:
    Coloring();

    virtual size_t size() const;

    virtual void set(const T& elt, int color);
    virtual void erase(const T& elt);
    virtual int operator[](const T& elt) const;

    virtual const std::set<int>& colors() const;
    virtual const std::unordered_set<T, THash>& cell(int color) const;

protected:
    std::unordered_map<T, int, THash> coloring;
    std::unordered_map<int, std::unordered_set<T, THash>> cells;
    std::set<int> _colors;
};

template<class T, class THash> Coloring<T, THash>::Coloring() {
    cells = std::unordered_map<int, std::unordered_set<T, THash>>();
    _colors = std::set<int>();
}

template<class T, class THash> size_t Coloring<T, THash>::size() const {
    return coloring.size();
}

template<class T, class THash>
                       void Coloring<T, THash>::set(const T& elt, int color) {
    auto cell_itr = cells.find(color);
    if (cell_itr == cells.end()) {
        cells[color] = std::unordered_set<T, THash>({elt});
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


template<class T, class THash>
                       void Coloring<T, THash>::erase(const T& elt) {
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

template<class T, class THash>
                        int Coloring<T, THash>::operator[](const T& elt) const {
    auto itr = coloring.find(elt);
    if (itr == coloring.end()) {
        throw std::range_error("Error. Element not in coloring.");
    }
    return itr->second;
}

template<class T, class THash> const std::set<int>&
                            Coloring<T, THash>::colors() const {
    return _colors;
}

template<class T, class THash> const std::unordered_set<T, THash>&
                            Coloring<T, THash>::cell(int color) const {
    auto cell_itr = cells.find(color);
    if (cell_itr == cells.end()) {
        throw std::range_error("Error. Color " + std::to_string(color) +
                               " not in coloring.");
    }
    return cell_itr->second;
}

#endif
