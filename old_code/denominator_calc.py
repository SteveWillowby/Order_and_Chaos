import bigfloat as bf
import math
import numpy as np
import scipy

def bigfloat_str(b, k=79):
    s = str(b).split(".")
    if len(s[0]) > k or len(s) == 1:
        return s[0]

    end = ""
    if "e" in s[1]:
        end = s[1][s[1].find("e") - 1:]
        s[1] = s[1][:-len(end)]
    k -= len(end)
    return s[0] + "." + s[1][:(k - min(k,len(s[0])))] + end

def print_bf(b):
    print(bigfloat_str(b))

def undirected_to_directed():
    oeis_undir = [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168, \
                  1018997864, 165091172592, 50502031367952, 29054155657235488, \
                  31426485969804308768, 64001015704527557894928, \
                  245935864153532932683719776, 1787577725145611700547878190848,\
                  24637809253125004524383007491432768]

    print(str([bigfloat_str(bf.log2(x), k=30) + " " for x in oeis_undir]).replace("'", ""))
    print(str([bigfloat_str(bf.log2(x), k=31) for x in oeis_undir]).replace("'", ""))

    oeis_dir = [1, 1, 3, 16, 218, 9608, 1540944, 882033440, 1793359192848, \
                13027956824399552, 341260431952972580352, \
                32522909385055886111197440, 11366745430825400574433894004224, \
                14669085692712929869037096075316220928, \
                70315656615234999521385506555979904091217920]

    print("")
    print(str([bigfloat_str(bf.log2(x), k=30) + " " for x in oeis_dir]).replace("'", ""))
    print(str([bigfloat_str(bf.log2(x), k=31) for x in oeis_dir]).replace("'", ""))

    max_n = min(len(oeis_undir), len(oeis_dir)) - 1

    factors = [bf.BigFloat(oeis_undir[n]) / oeis_dir[n] \
                    for n in range(0, max_n + 1)]

    factors_as_powers_of_n = [bf.log2(factors[n]) / bf.log2(n) \
                                for n in range(0, max_n + 1)]

    for n in range(0, max_n + 1):
        print("")
        print("Factor for %d:" % n)
        print_bf(factors_as_powers_of_n[n])
        print("  Error:")
        print_bf(factors[n] - bf.pow(n, factors_as_powers_of_n[n]))

def fancy_formula(n, directed=False):
    if directed:
        assert "DIRECTED NOT" == "IMPLEMENTED"

    constant = 1.0

    mess = (bf.BigFloat(1.0) + (n*n - n) / bf.pow(2, n - 1) + \
                8 * n * (n - 1) * (n - 2) * (n - 3) * \
                (3*n - 7) * (3*n - 9) / \
                    bf.pow(2, 2 * n))

    needs_constant = bf.BigFloat(n**5) / bf.pow(2, 5 * n / 2.0)

    mess += needs_constant * constant

    assert mess > 1.0

    return mess
    

# TODO: Consider adding support for self-loops
def estimated_num_graphs(n, directed=False):

    if directed:
        E = n * (n - 1)
        assert "DIRECTED NOT" == "IMPLEMENTED"
    else:
        E = (n * (n - 1)) / 2

    log2_n_fact = bf.BigFloat(0)
    for i in range(2, n + 1):
        log2_n_fact += bf.log2(i)

    log2_2_to_E = bf.BigFloat(E)

    log2_big_chunk = log2_2_to_E - log2_n_fact

    mess = fancy_formula(n, directed)

    return bf.exp2(log2_big_chunk) * mess

def beta_part_num_graphs(beta, n_fact_to_beta):
    f = (scipy.special.gammaincc(np.double(1.0 / beta), np.double(1.0)) * \
                bf.gamma(1.0 / beta) - \
         scipy.special.gammaincc(np.double(1.0 / beta),
                                 np.double(n_fact_to_beta)) * \
                bf.gamma(1.0 / beta))
    # print(scipy.special.gammaincc(np.double(1.0 / beta), np.double(1.0)))
    # print(np.double(n_fact_to_beta))  # inf
    # print(scipy.special.gammaincc(np.double(1.0 / beta), np.double(n_fact_to_beta)))  # 0
    return f / beta

def beta_part_num_mat(beta, n_fact_to_beta):
    f = (bf.eint(-n_fact_to_beta) - bf.eint(-1.0))
    return f / beta

def to_minimize(n, beta, directed=False):
    if directed:
        assert "DIRECTED NOT" == "IMPLEMENTED"

    A = factor = fancy_formula(n, directed)

    n_fact_to_beta = 1.0
    for i in range(2, n + 1):
        n_fact_to_beta *= bf.pow(i, beta)

    B = beta_part_num_graphs(beta, n_fact_to_beta)

    C = beta_part_num_mat(beta, n_fact_to_beta)
    # print("A: ")
    # print_bf(A)
    # print("B/C: ")
    # print_bf(B / C)

    return bf.abs(A - (B / C))


def find_beta(n, directed=False):
    if directed:
        assert "DIRECTED NOT" == "IMPLEMENTED"

    low = bf.pow(bf.BigFloat(1.0) / n, 2.0)
    high = bf.exp2(n)
    # mid = bf.sqrt(high * low)
    mid = (high + low) / 2.0

    epsilon_factor = bf.BigFloat(1.0) + bf.min(bf.BigFloat(1.0) / high, 0.00001)

    while high / low > bf.sqrt(epsilon_factor):
        value =        to_minimize(n, mid, directed)
        higher_value = to_minimize(n, mid * epsilon_factor, directed)
        lower_value =  to_minimize(n, mid / epsilon_factor, directed)

        if higher_value < lower_value:
            if value < higher_value:
                return mid
            # Go up
            low = mid
        else:
            if value < lower_value:
                return mid
            # Go down
            high = mid

        # mid = bf.sqrt(high * low)
        mid = (high + low) / 2.0

    print("High:")
    print_bf(high)
    print("Low:")
    print_bf(low)
    print("Error:")
    print_bf(value)
    return mid

def upper_incomplete_gamma(s, z):
    return bf.gamma(s) - lower_incomplete_gamma(s, z)

def lower_incomplete_gamma(s, z, epsilon_factor=(bf.BigFloat(1.0) + (bf.BigFloat(1.0) / 1000))):
    curr_value = bf.pow(z, s) * bf.exp(-z)
    # TODO: Consider finishing

def num_automorphisms_estimate(n, beta, alpha):
    n_fact_to_beta = 1.0
    for i in range(2, n + 1):
        n_fact_to_beta *= bf.pow(i, beta)

    f = (scipy.special.gammaincc(np.double(2.0 / beta), np.double(1.0)) * \
                bf.gamma(2.0 / beta) - \
         scipy.special.gammaincc(np.double(2.0 / beta),
                                 np.double(n_fact_to_beta)) * \
                bf.gamma(2.0 / beta))
    return alpha * f / beta

if __name__ == "__main__":

    MAX_N = 5000
    precision = MAX_N*2*(int(math.log2(MAX_N) + 1))
    emax =  int(MAX_N**1.5)
    emin = -int(MAX_N**1.5)

    bf.setcontext(bf.Context(precision=precision, emax=emax, emin=emin))

    # Don't try to actualize 2^(N^2)

    undirected_to_directed()

    log2_n_fact = bf.BigFloat(0.0)

    f = open("P_values.h", "w")
    f.write("#include<vector>\n\n")
    f.write("// Contains a map n --> X such that in an n-node graph with m\n")
    f.write("//   edges, log2_(1 - p) = P_values[n]\n")
    f.write("std::vector<long double> P_values = {-1.0, -1.0")

    START_N = 50

    for n in range(2, START_N):
        log2_n_fact += bf.log2(n)

    for n in range(START_N, MAX_N + 1):
        log2_n_fact += bf.log2(n)
        big_hunk = (n * (n - 1)) / 2 - log2_n_fact
        big_hunk = bf.exp2(big_hunk)

        B = find_beta(n)
        A = (big_hunk * fancy_formula(n)) / \
                beta_part_num_graphs(B, bf.exp2(B * log2_n_fact))

        num_mat = bf.exp2(n * (n - 1) / 2)
        est_num_mat = bf.exp2(log2_n_fact) * A * \
                      beta_part_num_mat(B, bf.exp2(B * log2_n_fact))

        num_graphs = big_hunk * fancy_formula(n)
        est_num_graphs = A * beta_part_num_graphs(B, bf.exp2(B * log2_n_fact))

        nae = num_automorphisms_estimate(n, B, A)

        if n % 50 == 0:
            print("N = %d" % n)
            print("Meta Error Factors:")
            print_bf(num_mat / est_num_mat)
            print_bf(num_graphs / est_num_graphs)
            print("Num Graphs / Num Automorphisms Estimate")
            print_bf(num_graphs / nae)
            print("fancy_factor:")
            print_bf(fancy_formula(n))
            print("fancy factor reconstruction:")
            print_bf(est_num_graphs / big_hunk)
            print("fancy_factor - (Est Num Graphs / (Est Num Mat / n!))")
            print_bf(fancy_formula(n) -
                     (est_num_graphs / (est_num_mat / bf.exp2(log2_n_fact))))
            print("fancy_factor / (nae / num_graphs)")
            print_bf(fancy_formula(n) / (nae / num_graphs))


        S = -(bf.log2(nae) - log2_n_fact) / ((n * (n - 1)) / 2)
        f.write(",\n%s" % (bigfloat_str(S, k=120)))

    f.write("};\n")
    f.close()

