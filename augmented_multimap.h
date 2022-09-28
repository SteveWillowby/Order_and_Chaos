// This class implements a minimal multimap where the key-value pairs are
//  sorted by key but not by value (thus there can be ties).
//
// The main point of the class is to make accessing and removing a single
//  key-value pair more efficient. With a standard multimap, removing a
//  particular key-value pair is big-O in the number of key-value pairs with
//  the same key. Here is it O(1).


#include<map>

// Define this if you want deterministic behavior for testing purposes.
//  When this is defined, the set of values is a set rather than an
//  unordered_set, meaning that whenever a 
#define SYM__AUGMENTED_MULTIMAP_FULLY_SORTED

#ifdef SYM__AUGMENTED_MULTIMAP_FULLY_SORTED
#include<set>
#define augmented_multimap_values std::set
#else
#include<unordered_set>
#define augmented_multimap_values std::unordered_set
#endif

#ifndef SYM__AUGMENTED_MULTIMAP_H
#define SYM__AUGMENTED_MULTIMAP_H

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


protected:
    std::map<K, augmented_multimap_values<V>> data;

};

#endif
