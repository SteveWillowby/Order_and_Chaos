import bigfloat
import math

# Interface:
#
# log2_automorphisms(directed, has_edge_types, neighbors_collections,
#                    auto_solver_class)
#
#   Returns a float if possible. Might need to return a bigfloat, especially
#       if `auto_solver_class` is not PyNTSession
#
#   Current supported classes for `auto_solver_class` are:
#       * PyNTSession
#       * RAMFriendlyNTSession
#       * Any class that satisfies one of those two interfaces.
#
# log2_factorial(a)
#
#   Returns a float.
#
#
# log2_a_choose_b(a, b)
#
#   Returns a float.


LOG_2_OF_10 = math.log2(10)
FLOAT_INF = float("inf")
BIGFLOAT_INF = bigfloat.BigFloat("inf")

def log2_automorphisms(directed, has_edge_types, neighbors_collections, \
                       auto_solver_class):

    
    auto_solver = auto_solver_class(directed, has_edge_types, \
                                    neighbors_collections)
    num_automorphs = auto_solver.get_num_automorphisms()
    auto_solver.run()
    num_automorphs = num_automorphs.get()
    auto_solver.end_session()
    if type(num_automorphs) is tuple:
        assert len(num_automorphs) == 2
        (a, b) = num_automorphs
        return math.log2(a) + LOG_2_OF_10 * b

    assert type(num_automorphs) is str
    i = 1
    while True:
        # Compute the log-number of automorphisms.
        if i == 1:
            num_automorphisms = float(num_automorphs)
            if num_automorphisms != FLOAT_INF:
                return math.log2(num_automorphisms)

            else:
                print("Not enough bits to store num of automorphisms.")
                print("Retrying with bigfloats.")

        num_automorphisms = bigfloat.BigFloat(num_automorphs)
        if num_automorphisms != BIGFLOAT_INF:
            break
        else:
            print("Still not enough bits to store num of automorphisms.")
            print("Retrying with bigger bigfloats.")

        if i == 1:
            old_context = bigfloat.getcontext()

        bf_context = bigfloat.Context(precision=2000*i, \
                                      emax=1000000*i, emin=-1000000*i)
        bigfloat.setcontext(bf_context)
        i += 1

    if i > 1:
        bigfloat.setcontext(old_context)

    return bigfloat.log2(num_automorphisms)


# A library of log2(a choose b) values keyed by (a, b)
THE_FROZEN_CHOSEN = {}

# A library of log2(a!) values keyed by a
THE_FROZEN_FACTS = [0.0]  # log2(1)

def log2_a_choose_b(a, b):
    assert a >= b
    if a == b or b == 0:
        return 0.0  # log2(1)

    if (a, b) in THE_FROZEN_CHOSEN:
        return THE_FROZEN_CHOSEN[(a, b)]

    if a < len(THE_FROZEN_FACTS):
        value = THE_FROZEN_FACTS[a] - \
                (THE_FROZEN_FACTS[b] + THE_FROZEN_FACTS[a - b])
        THE_FROZEN_CHOSEN[(a, b)] = value
        return value

    # After the first part evaluates, the following values exist.
    value = log2_factorial(a) - \
            (THE_FROZEN_FACTS[b] + THE_FROZEN_FACTS[a - b])
    THE_FROZEN_CHOSEN[(a, b)] = value
    return value

def log2_factorial(a):
    if a < len(THE_FROZEN_FACTS):
        return THE_FROZEN_FACTS[a]

    l = len(THE_FROZEN_FACTS)
    cumulative = THE_FROZEN_FACTS[l - 1]
    while l <= a:
        cumulative += math.log2(l)
        THE_FROZEN_FACTS.append(cumulative)
        l += 1

    return THE_FROZEN_FACTS[a]

if __name__ == "__main__":
    from py_NT_session import PyNTSession
    from ram_friendly_NT_session import RAMFriendlyNTSession

    directed = False
    has_edge_types = False
    neighbors_collections = [set([1, 5]), set([0, 2]), set([1, 3]), set([2, 4]), set([3, 5]), set([4, 0])]
    print("Expecting %s" % bigfloat.log2(12))
    print(log2_automorphisms(directed, has_edge_types, neighbors_collections, PyNTSession))
    print(log2_automorphisms(directed, has_edge_types, neighbors_collections, RAMFriendlyNTSession))

    N = 1000
    neighbors_collections = [set([j for j in range(0, i)] + [j for j in range(i + 1, N)]) \
                                for i in range(0, N)]
    print("Expecting %s" % sum([bigfloat.log2(i) for i in range(1, N + 1)]))
    print(log2_automorphisms(directed, has_edge_types, neighbors_collections, PyNTSession))
    print(log2_automorphisms(directed, has_edge_types, neighbors_collections, RAMFriendlyNTSession))

    print("Expecting %f == %f" % (log2_a_choose_b(12, 5), math.log2(12*11*10*9*8 / (5*4*3*2*1))))
    print("Expecting %f == %f" % (log2_a_choose_b(12, 5), math.log2(12*11*10*9*8 / (5*4*3*2*1))))
    print("Expecting %f == %f" % (log2_a_choose_b(7, 5), math.log2(7*6 / (2*1))))
