// Funky attempt -- real stuff below

/*
    size_t chunk_size = 2;
    size_t half_chunk_size = 1;
    size_t quarter_chunk_size = 0;

    size_t prev_vec_B = 2;
    size_t prev_vec_A = 0;
    size_t curr_vec = 1;

    while (chunk_size < num_elts) {  // TODO: Update stopping criterion
        // partial_sums[prev_vec][i] is equal to:
        //     SUM_{from j = i - quarter_chunk_size to i} of log2_s[j]
        for (size_t i = 2; i < num_elts; i++) {
            if (i < half_chunk_size) {
                partial_sums[curr_vec][i] =
                              partial_sums[curr_vec][i - half_chunk_size]
                            + partial_sums[prev_vec][i];
            }
            partial_sums[curr_vec][i] = partial_sums[prev_vec][i] +
                                        partial_sums[prev_vec][i - chunk_size];
        }
        quarter_chunk_size = half_chunk_size;
        half_chunk_size = chunk_size;
        chunk_size *= 2;
        prev_vec_ = curr_vec;
        curr_vec = 1 - curr_vec;
    }

    log2_factorials = partial_sums[prev_vec];
*/



__CombinatoricUtility::__CombinatoricUtility() {
    // For convenience we say that 0! = 1, and thus log2(0!) = 0.
    log2_factorials = std::vector<double>(2, 0.0);
    // There is no log2 for 0, but oh well. Put a placeholder value.
    log2_s = std::vector<double>(2, 0.0);
}

double __CombinatoricUtility::log2(size_t x) {
    if (x >= log2_s.size()) {
        for (size_t i = log2_s.size(); i <= x; i++) {
            log2_s.push_back(std::log2(i));
        }
    }
    return log2_s[x];
}

double __CombinatoricUtility::log2_factorial(size_t x) {
    if (x < log2_factorials.size()) {
        return log2_factorials[x];
    }
    double v = log2_factorials[log2_factorials.size() - 1];
    while (x >= log2_factorials.size()) {
        v += log2(log2_factorials.size());
        log2_factorials.push_back(v);
    }
    return log2_factorials[x];
}

double __CombinatoricUtility::log2_a_choose_b(size_t a, size_t b) {
    double v1 = log2_factorial(a);
    return v1 - (log2_factorials[b] + log2_factorials[a - b]);
}
