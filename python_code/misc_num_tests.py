import bigfloat
import math
from prob_calc_utils import log2_a_choose_b

def f(p_plus, p_minus, q, k_plus, k_minus, n, m):
    return k_plus *  (math.log2(q - p_plus) + math.log2(p_plus) + \
                      -math.log2(1 - p_minus) + -math.log2(p_minus - q)) + \
           k_minus * (math.log2(p_minus - q) + math.log2(p_minus) + \
                      -math.log2(1 - p_plus) + -math.log2(q - p_plus)) + \
           -log2_a_choose_b(int((n * (n - 1)) / 2), m + k_minus - k_plus)

if __name__ == "__main__":
    n = 200
    m = 500
    k_plus_orig = 20
    k_minus_orig = 30

    q = m / ((n * (n - 1)) / 2.0)
    p_plus = k_plus_orig / ((n * (n - 1)) / 2.0 - m)
    p_minus = k_minus_orig / float(m)

    assert p_minus > q
    assert q > p_plus

    for k_minus in range(0, 50):
        print(bigfloat.BigFloat(2.0)**f(p_plus, p_minus, q, k_plus_orig, k_minus, n, m))
