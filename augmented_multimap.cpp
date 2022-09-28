#include "augmented_multimap.h"

#include<map>
#include<set>
#include<stdexcept>
#include<unordered_set>

///////////////////////////////// KeyValuePair /////////////////////////////////

template<class K, class V> KeyValuePair<K,V>::KeyValuePair() : _is_none(true) {}
template<class K, class V>
    KeyValuePair<K, V>::KeyValuePair(const K &key, const V &value) : 
        _is_none(false) {

    _key = key;
    _value = value;
}
template<class K, class V> bool KeyValuePair<K, V>::is_none() const {
    return _is_none;
}
template<class K, class V> const K &KeyValuePair<K, V>::key() const {
    if (_is_none) {
        throw std::logic_error(
                "Error! Cannot access key from 'none' KeyValuePair");
    }
    return key;
}
template<class K, class V> const V &KeyValuePair<K, V>::value() const {
    if (_is_none) {
        throw std::logic_error(
                "Error! Cannot access value from 'none' KeyValuePair");
    }
    return value;
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
        data.insert(key, values);
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
