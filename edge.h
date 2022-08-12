/* A hashable edge. */

#include<pair>

#ifndef SYM__EDGE_H
#define SYM__EDGE_H

#define EDGE(s, t, d) std::pair<int>((d || s < t ? s : t), (d || s < t ? t : s))

typedef std::pair<int, int> Edge;

#endif
