// This class implements a minimal multimap where the key-value pairs are
//  sorted by key but not by value (thus there can be ties).
//
// The main point of the class is to make accessing and removing a single
//  key-value pair more efficient. With a standard multimap, removing a
//  particular key-value pair is big-O in the number of key-value pairs with
//  the same key. Here is it O(1).


#include<cstddef>
#include<map>
#include<stdexcept>
#include<vector>

// Define this if you want deterministic behavior for testing purposes.
//  When this is defined, the set of values is a set rather than an
//  unordered_set, meaning that whenever a 
#define SCHENO__AUGMENTED_MULTIMAP_FULLY_SORTED

#ifdef SCHENO__AUGMENTED_MULTIMAP_FULLY_SORTED
#include<set>
#define augmented_multimap_values std::set
#else
#include<unordered_set>
#define augmented_multimap_values std::unordered_set
#endif

#ifndef SCHENO__AUGMENTED_MULTIMAP_H
#define SCHENO__AUGMENTED_MULTIMAP_H

template<class K, class V> class KeyValuePair {
public:
    // Used to construct an empty pair. is_none() returns true.
    KeyValuePair();
    // Used to construct a normal pair. is_none() returns false.
    KeyValuePair(const K &key, const V &value);

    bool is_none() const;
    // Throws error when is_none() is true.
    const K &key() const;
    // Throws erros when is_none() is true.
    const V &value() const;

protected:
    const bool _is_none;
    K _key;
    V _value;
};
 

template<class K, class V> class AugmentedMultimap {

public:
    AugmentedMultimap();

    // Returns true iff the key-value pair is new.
    bool insert(const K &key, const V &value);

    // Returns number of key-value pairs that were removed.
    size_t erase(const K &key);
    // Returns true iff the key-value pair was present.
    bool erase(const K &key, const V &value);

    // Returns a key-value pair where the key is one of the smallest keys that
    //  is >= `key`. If no such pair can be found, returns a `is_none` pair.
    KeyValuePair<K, V> lower_bound(const K &key) const;

    bool contains(const K &key) const;
    bool contains(const K &key, const V& value) const;

    // O(n)
    std::vector<std::pair<const K&, const V&>> all_pairs() const;

protected:
    std::map<K, augmented_multimap_values<V>> data;

};

///////////////////////////////// KeyValuePair /////////////////////////////////

template<class K, class V> KeyValuePair<K,V>::KeyValuePair() :
    _is_none(true), _key(K()), _value(V()) {}
template<class K, class V>
    KeyValuePair<K, V>::KeyValuePair(const K &key, const V &value) : 
        _is_none(false), _key(key), _value(value) {}
template<class K, class V> bool KeyValuePair<K, V>::is_none() const {
    return _is_none;
}
template<class K, class V> const K &KeyValuePair<K, V>::key() const {
    if (_is_none) {
        throw std::logic_error(
                "Error! Cannot access key from 'none' KeyValuePair");
    }
    return _key;
}
template<class K, class V> const V &KeyValuePair<K, V>::value() const {
    if (_is_none) {
        throw std::logic_error(
                "Error! Cannot access value from 'none' KeyValuePair");
    }
    return _value;
}

/////////////////////////////// AugmentedMultimap //////////////////////////////

template<class K, class V> AugmentedMultimap<K, V>::AugmentedMultimap() {
    data = std::map<K, augmented_multimap_values<V>>();
}

// Returns true iff the key-value pair is new.
template<class K, class V> bool AugmentedMultimap<K,V>::insert(const K &key,
                                                               const V &value) {
    auto ref = data.find(key);
    if (ref == data.end()) {
        augmented_multimap_values<V> values = augmented_multimap_values<V>();
        values.insert(value);
        data.insert(std::pair<K, augmented_multimap_values<V>>(key, values));
        return true;
    }
    return ref->second.insert(value).second;
}

// Returns number of key-value pairs that were removed.
template<class K, class V> size_t AugmentedMultimap<K, V>::erase(const K &key) {

    auto ref = data.find(key);
    if (ref == data.end()) {
        return 0;
    }

    size_t size = ref->second.size();
    data.erase(ref);
    return size;
}

// Returns true iff the key-value pair was present.
template<class K, class V> bool AugmentedMultimap<K, V>::erase(const K &key,
                                                               const V &value) {
    auto key_ref = data.find(key);
    if (key_ref == data.end()) {
        return false;
    }
    auto value_ref = key_ref->second.find(value);
    if (value_ref == key_ref->second.end()) {
        return false;
    }

    if (key_ref->second.size() == 1) {
        data.erase(key_ref);
    } else {
        key_ref->second.erase(value_ref);
    }
    return true;
}

template<class K, class V> KeyValuePair<K, V>
        AugmentedMultimap<K, V>::lower_bound(const K &key) const {
    const auto &result = data.lower_bound(key);
    if (result == data.end()) {
        return KeyValuePair<K, V>();
    }
    return KeyValuePair<K, V>(result->first, *(result->second.begin()));
}

template<class K, class V> bool
    AugmentedMultimap<K, V>::contains(const K &key) const {
    return data.find(key) != data.end();
}

template<class K, class V> bool
    AugmentedMultimap<K, V>::contains(const K &key, const V& value) const {
    const auto &key_ref = data.find(key);
    return key_ref != data.end() &&
           key_ref->second.find(value) != key_ref->second.end();
}

template<class K, class V>
        std::vector<std::pair<const K&, const V&>>
            AugmentedMultimap<K, V>::all_pairs() const {

    std::vector<std::pair<const K&, const V&>> result = 
                    std::vector<std::pair<const K&, const V&>>();

    for (auto main_itr = data.begin(); main_itr != data.end(); main_itr++) {
        for (auto sub_itr = main_itr->second.begin();
                    sub_itr != main_itr->second.end(); sub_itr++) {
            result.push_back(std::pair<const K&, const V&>(main_itr->first,
                                                           *sub_itr));
        }
    }

    return result;
}

#endif
