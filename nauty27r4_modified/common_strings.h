// #define COMMON_STRINGS_FOR_UTILS

#ifndef COMMON_STRINGS_H
#define COMMON_STRINGS_H

static char char_slash[2] = "/";
static char *str_slash = &(char_slash[0]);

static char char_semi_[3] = "; ";
static char *str_semi_ = &(char_semi_[0]);

static char char__orbit[7] = " orbit";
static char *str__orbit = &(char__orbit[0]);

static char char_level_[7] = "level ";
static char *str_level_ = &(char_level_[0]);

static char char_colon__[4] = ":  ";
static char *str_colon__ = &(char_colon__[0]);

static char char__cell[6] = " cell";
static char *str__cell = &(char__cell[0]);

static char char_ssemi_[4] = "s; ";
static char *str_ssemi_ = &(char_ssemi_[0]);

static char char__fixedsemi_index_[15] = " fixed; index ";
static char *str__fixedsemi_index_ = &(char__fixedsemi_index_[0]);

static char char_endl[2] = "\n";
static char *str_endl = &(char_endl[0]);

static char char_endl_[3] = "\n ";
static char *str_endl_ = &(char_endl_[0]);

static char char_endl__[4] = "\n  ";
static char *str_endl__ = &(char_endl__[0]);

static char char_endl___[5] = "\n   ";
static char *str_endl___ = &(char_endl___[0]);

static char char_semiendl[3] = ";\n";
static char *str_semiendl = &(char_semiendl[0]);

static char char_lpar1rparendl[5] = "(1)\n";
static char *str_lpar1rparendl = &(char_lpar1rparendl[0]);

static char char__colon[3] = " :";
static char *str__colon = &(char__colon[0]);


#ifdef COMMON_STRINGS_FOR_UTILS

static char char_none[5] = "none";
static char *str_none = &(char_none[0]);

static char char_twopaths[9] = "twopaths";
static char *str_twopaths = &(char_twopaths[0]);

static char char_unavailable[12] = "unavailable";
static char *str_unavailable = &(char_unavailable[0]);

static char char_adjtriang[10] = "adjtriang";
static char *str_adjtriang = &(char_adjtriang[0]);

static char char_triples[8] = "triples";
static char *str_triples = &(char_triples[0]);

static char char_quadruples[11] = "quadruples";
static char *str_quadruples = &(char_quadruples[0]);

static char char_celltrips[10] = "celltrips";
static char *str_celltrips = &(char_celltrips[0]);

static char char_cellquads[10] = "cellquads";
static char *str_cellquads = &(char_cellquads[0]);

static char char_cellquins[10] = "cellquins";
static char *str_cellquins = &(char_cellquins[0]);

static char char_distances[10] = "distances";
static char *str_distances = &(char_distances[0]);

static char char_distances_sg[13] = "distances_sg";
static char *str_distances_sg = &(char_distances_sg[0]);

static char char_indsets[8] = "indsets";
static char *str_indsets = &(char_indsets[0]);

static char char_cliques[8] = "cliques";
static char *str_cliques = &(char_cliques[0]);

static char char_cellcliq[9] = "cellcliq";
static char *str_cellcliq = &(char_cellcliq[0]);

static char char_cellind[8] = "cellind";
static char *str_cellind = &(char_cellind[0]);

static char char_adjacencies[12] = "adjacencies";
static char *str_adjacencies = &(char_adjacencies[0]);

static char char_adjacencies_sg[15] = "adjacencies_sg";
static char *str_adjacencies_sg = &(char_adjacencies_sg[0]);

static char char_cellfano[9] = "cellfano";
static char *str_cellfano = &(char_cellfano[0]);

static char char_cellfano2[10] = "cellfano2";
static char *str_cellfano2 = &(char_cellfano2[0]);

static char char_refinvar[9] = "refinvar";
static char *str_refinvar = &(char_refinvar[0]);

static char char_empty[1] = "";
static char *str_empty = &(char_empty[0]);

static char char_w[2] = "w";
static char *str_w = &(char_w[0]);

static char char_V[2] = "V";
static char *str_V = &(char_V[0]);

static char char_S[2] = "S";
static char *str_S = &(char_S[0]);

static char char_G[2] = "G";
static char *str_G = &(char_G[0]);

static char char_y[2] = "y";
static char *str_y = &(char_y[0]);

static char char_dollar[2] = "$";
static char *str_dollar = &(char_dollar[0]);

static char char_M[2] = "M";
static char *str_M = &(char_M[0]);

static char char_Mslash[3] = "M/";
static char *str_Mslash = &(char_Mslash[0]);

static char char_l[2] = "l";
static char *str_l = &(char_l[0]);

static char char_copyg_dashI[9] = "copyg -I";
static char *str_copyg_dashI = &(char_copyg_dashI[0]);

static char char_copyg_dashp[9] = "copyg -p";
static char *str_copyg_dashp = &(char_copyg_dashp[0]);

static char char_listg_dashI[9] = "listg -I";
static char *str_listg_dashI = &(char_listg_dashI[0]);

static char char_listg_dashp[9] = "listg -p";
static char *str_listg_dashp = &(char_listg_dashp[0]);

static char char_listg_dasho[9] = "listg -o";
static char *str_listg_dasho = &(char_listg_dasho[0]);

static char char_listg_dashl[9] = "listg -l";
static char *str_listg_dashl = &(char_listg_dashl[0]);

static char char_colondash[3] = ":-";
static char *str_colondash = &(char_colondash[0]);

static char char_stdin[6] = "stdin";
static char *str_stdin = &(char_stdin[0]);

static char char_stdout[7] = "stdout";
static char *str_stdout = &(char_stdout[0]);

static char char_sparse6[12] = ">>sparse6<<";
static char *str_sparse6 = &(char_sparse6[0]);

static char char_digraph6[13] = ">>digraph6<<";
static char *str_digraph6 = &(char_digraph6[0]);

static char char_graph6[11] = ">>graph6<<";
static char *str_graph6 = &(char_graph6[0]);

static char char_planar_code[16] = ">>planar_code<<";
static char *str_planar_code = &(char_planar_code[0]);

static char char_planar_code_le[19] = ">>planar_code le<<";
static char *str_planar_code_le = &(char_planar_code_le[0]);

static char char_planar_code_be[19] = ">>planar_code be<<";
static char *str_planar_code_be = &(char_planar_code_be[0]);

static char char_edge_code[14] = ">>edge_code<<";
static char *str_edge_code = &(char_edge_code[0]);

static char char_labelg_dashC[10] = "labelg -C";
static char *str_labelg_dashC = &(char_labelg_dashC[0]);

static char char_labelg_dashW[10] = "labelg -W";
static char *str_labelg_dashW = &(char_labelg_dashW[0]);

static char char_labelg_dashI[10] = "labelg -I";
static char *str_labelg_dashI = &(char_labelg_dashI[0]);

static char char_labelg_dashi[10] = "labelg -i";
static char *str_labelg_dashi = &(char_labelg_dashi[0]);

static char char_labelg_dashK[10] = "labelg -K";
static char *str_labelg_dashK = &(char_labelg_dashK[0]);

static char char_labelg_dashk[10] = "labelg -k";
static char *str_labelg_dashk = &(char_labelg_dashk[0]);

static char char_labelg_dashM[10] = "labelg -M";
static char *str_labelg_dashM = &(char_labelg_dashM[0]);

static char char_E_dretog_dasho[13] = ">E dretog -o";
static char *str_E_dretog_dasho = &(char_E_dretog_dasho[0]);

static char char_E_dretog_dashn[13] = ">E dretog -n";
static char *str_E_dretog_dashn = &(char_E_dretog_dashn[0]);

static char char_E_amtog_dasho[12] = ">E amtog -o";
static char *str_E_amtog_dasho = &(char_E_amtog_dasho[0]);

static char char_E_amtog_dashn[12] = ">E amtog -n";
static char *str_E_amtog_dashn = &(char_E_amtog_dashn[0]);

static char char_geng_dashd[8] = "geng -d";
static char *str_geng_dashd = &(char_geng_dashd[0]);

static char char_geng_dashD[8] = "geng -D";
static char *str_geng_dashD = &(char_geng_dashD[0]);

static char char_geng_dashx[8] = "geng -x";
static char *str_geng_dashx = &(char_geng_dashx[0]);

static char char_geng_dashX[8] = "geng -X";
static char *str_geng_dashX = &(char_geng_dashX[0]);

#endif

#endif
