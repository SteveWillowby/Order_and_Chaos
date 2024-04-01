#include<cmath>

#include "noise_prob_choice.h"

#include "nt_wrappers/nauty_traces.h"

#include "scoring_function.h"

long double log2_fancy_factor(long double n, long double max_E_no_SL,
                              long double log2_n_fact, bool directed) {

    if (directed && n < 15) {
        std::vector<long double> log2_num_graphs = {0, 0,
            1.58496250072115618145373894395, 4.00000000000000000000000000000,
            7.76818432477692635847878680267, 13.2300204357056340945901830909,
            20.5553830026302805397007004079, 29.7162581119425090062051620927,
            40.7058016139384163498316900063, 53.5324603618554374979747570314,
            68.2094350460413784972197484213, 84.7496586937927632790050761116,
            103.164590177010205079441905893, 123.464118463118692395266911061,
            145.656754038391204022206232903};
        long double log2_ratio =
                        (log2_num_graphs[n] - max_E_no_SL) + log2_n_fact;
        return log2_ratio;
    } else if (!directed && n < 20) {
        std::vector<long double> log2_num_graphs = {0, 0, 1.0, 2.0,
            3.45943161863729725619936304675, 5.08746284125033940825406601081,
            7.28540221886224834185055159820, 10.0279059965698844836272512125,
            13.5917560760229900027396129249, 18.0673293144811357155639053873,
            23.5171522569569572352913536136, 29.9245038813488327373283308491,
            37.2644720253191471328061410100, 45.5214066527967648257388625947,
            54.6895940473113672250797267462, 64.7686147655035370309689021222,
            75.7605128885761279188207305346, 87.6684126016135392619332623025,
            100.495848819277919926324471434, 114.246429206222428663749162609};
        long double log2_ratio =
                        (log2_num_graphs[n] - max_E_no_SL) + log2_n_fact;
        return log2_ratio;
    }

    long double r;
    // These formulae taken from Harary and Palmer's book
    //  "Graphical Enumeration"
    if (directed) {
        r = (n * n - n) / std::exp2l(n - 1)   +
            (n * (n - 1) * (n - 2) * (n - 3) * (3*n - 7)) /
                ((3*n - 9) * std::exp2l(2*n)) +
            (n * n * n * n * n) / std::exp2l(5.0 * n / 2.0);
    } else {
        r = (4 * n * (n - 1)) / std::exp2l(2 * n) +
            (n * (n - 1) * (n - 2) * (n - 3) * (3*n - 7)) /
                ((3*n - 9) * std::exp2l(4*n - 7)) +
            (n * n * n * n * n) / std::exp2l(5 * n);
    }
    return std::log1pl(r) / std::log(2.0);
}

std::vector<long double> log2_noise_probs_fancy_equality(NTSparseGraph& g,
                                      const CombinatoricUtility& comb_util) {

    bool directed = g.directed;
    size_t num_nodes = g.num_nodes();

    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed)) +
            (num_nodes * size_t(g.num_loops() > 0));

    long double log2_n_fact = comb_util.log2_factorial(num_nodes);
    long double max_E = max_possible_edges;
    long double max_E_no_SL = max_E - (num_nodes * double(g.num_loops() > 0));

    long double log2_f_factor = log2_fancy_factor((long double) num_nodes,
                                                  max_E_no_SL, log2_n_fact,
                                                  directed);

    long double power = 2.0;  // raise the un-logged fancy factor to this power
    if (directed && num_nodes <= 12) {
        std::vector<long double> powers = {1.0, 1.0,
          2.70951129135145, 2.85903000604259, 3.4796254752191, 3.58220606015045,
          3.25554047733668, 2.8             , 2.6            , 2.4             ,
          2.3             , 2.2             , 2.1};

        power = powers[num_nodes];

    } else if (!directed && num_nodes <= 16) {
        std::vector<long double> powers = {1.0, 1.0,
          2.0             , 2.26185950714291, 2.48328614770423, 2.88433360919868,
          3.41709692824793, 4.23390912612406, 5.25128327223332, 5.96813211471908,
          2.1 / 0.31      , 1.0 / 0.17      , 0.47 / 0.1      , 0.22 / 0.057    ,
          0.11 / 0.033    , 0.05 / 0.02     , 0.02 / 0.012};
        power = powers[num_nodes];
    }
    log2_f_factor *= power;

    // Equivalent to:
    //  log2_estimated_autos = max_E - (log2_n_fact - log2_f_factor);
    //  log2_1_minus_p = (log2_n_fact - log2_estimated_autos) / max_E
    long double log2_1_minus_p =
                    ((2.0 * log2_n_fact - log2_f_factor) / max_E) - 1.0;

    // NOTE: at as few as 120 nodes, log2_f_factor ceases to be relevant
    // volatile long double check = (2.0 * log2_n_fact - log2_f_factor);
    // check = check - (2.0 * log2_n_fact);
    // if (check == 0.0) {
    //     std::cout<<"Lost in translation."<<std::endl;
    // }

    long double log2_p = std::log2l(1.0 - std::exp2l(log2_1_minus_p));

    return {log2_p, log2_1_minus_p, log2_p, log2_1_minus_p};
}

std::vector<long double> log2_noise_probs_empty_g(NTSparseGraph& g,
                                      const CombinatoricUtility& comb_util) {
    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g, o);

    bool directed = g.directed;
    size_t num_nodes = g.num_nodes();
    size_t num_edges = g.num_edges();

    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed)) +
            (num_nodes * size_t(g.num_loops() > 0));

    long double alpha;
    if (num_edges < (max_possible_edges / 2)) {
        alpha = std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(num_edges));
    } else {
        alpha = std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(max_possible_edges - num_edges));
    }
    long double log2_p_plus = std::log2l(alpha) - std::log2l(1 + alpha);
    long double log2_1_minus_p_plus = -std::log2l(1 + alpha);
    long double log2_p_minus = log2_p_plus;
    long double log2_1_minus_p_minus = log2_1_minus_p_plus;

    return {log2_p_plus, log2_1_minus_p_plus,
            log2_p_minus, log2_1_minus_p_minus};
}

std::vector<long double> log2_noise_probs_empty_g_full(NTSparseGraph& g,
                                      const CombinatoricUtility& comb_util) {
    NautyTracesOptions o;

    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g, o);

    bool directed = g.directed;
    size_t num_nodes = g.num_nodes();
    size_t num_edges = g.num_edges();

    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed));

    long double k1 =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(num_edges));

    long double k2 =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                        (long double)(max_possible_edges - num_edges));

    // Here's hoping the decimal places work OK!
    //  TODO: Double-check the floating point math for safety
    long double log2_p_plus = std::log2l(k1 - (k1 * k2)) -
                              std::log2l(1.0 - (k1 * k2));
    long double log2_1_minus_p_plus = std::log2l(1.0 - k1) -
                                      std::log2l(1.0 - (k1 * k2));

    long double log2_p_minus = std::log2l(k2 - (k1 * k2)) -
                               std::log2l(1.0 - (k1 * k2));
    long double log2_1_minus_p_minus = std::log2l(1.0 - k2) -
                                       std::log2l(1.0 - (k1 * k2));

    return {log2_p_plus, log2_1_minus_p_plus,
            log2_p_minus, log2_1_minus_p_minus};
}
