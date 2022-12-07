__CombinatoricUtility::__CombinatoricUtility() {
    // For convenience we say that 0! = 1, and thus log2(0!) = 0.
    log2_factorials = std::vector<double>(2, 0.0);
    // There is no log2 for 0, but oh well. Put a placeholder value.
    log2_s = std::vector<double>(2, 0.0);
    max_access = 1;
}

void __CombinatoricUtility::set_max_access(size_t max_access) {
    this->max_access = max_access;

    log2_factorials = std::vector<double>(max_access + 1, 0.0);
    log2_s = std::vector<double>(max_access + 1, 0.0);

    for (size_t i = 2; i <= max_access; i++) {
        log2_s[i] = std::log2(i);
    }

    // Collect partial log sums so as to minimize risk of floating point
    //  resolution errors.
    std::vector<std::pair<size_t, double>> partial_sums
                                    = std::vector<std::pair<size_t, double>>();

    size_t chunk_size = 8;
    size_t next_sum_needed = 2;

    while (next_sum_needed <= max_access) {
        
    }
}

double __CombinatoricUtility::log2(size_t x) {
    if (x > max_access) {
        throw std::range_error(
                std::string("Error! Each thread must first call ")
                + "__combinatoric_utility.set_max_access(Y) with a value Y >= "
                + "the max possible number of edges for your graph."
                + "  NOTE: This will take O(Y) time.");
    }
    return log2_s[x];
}

double __CombinatoricUtility::log2_factorial(size_t x) {
    if (x > max_access) {
        throw std::range_error(
                std::string("Error! Each thread must first call ")
                + "__combinatoric_utility.set_max_access(Y) with a value Y >= "
                + "the max possible number of edges for your graph."
                + "  NOTE: This will take O(Y) time.");
    }
    return log2_factorials[x];
}

double __CombinatoricUtility::log2_a_choose_b(size_t a, size_t b) {
    if (b > a) {
        throw std::invalid_argument("Error! Cannot do a-choose-b with b > a.");
    } else if (a > max_access) {
        throw std::range_error(
                std::string("Error! Each thread must first call ")
                + "__combinatoric_utility.set_max_access(Y) with a value Y >= "
                + "the max possible number of edges for your graph."
                + "  NOTE: This will take O(Y) time.");
    }
    return log2_factorials[a] - (log2_factorials[b] + log2_factorials[a - b]);
}
