#include <string.h>
#include "spoa/spoa.hpp"
using namespace std;

static inline char *fasta_to_ccs(vector<string> sequences,
	int l, int m, int n, int g, int e, int q, int c, int v) {

//    vector<string> sequences = {
//        "TGTGTCTGTGTGTATGGCTGTGTCTATGTGTAGCTGTGTATGTGTGGCTGTGTGTATGTGGCTGTGTGTGTCTGTGTCTGTGTCTGTGTGTGTCTGTGTGTGTGTCTCTGTGTG",
//        "TGTGTCTGTGTGTGTATGGCTGTGTCTGTGTGTGGCTGTGTATGTGTGGCTGTGTATGTGGCTGTGTGTCTGTGTCTGTGTGTCTGTGTCTGTTTGTGTGTGTGTGTGTCTC",
//        "TGTGTGTGTGTGTGTGTCTGTGTGTCATATAGCTAATATCTATGCTATACCTCTTGGCTGTGTGTATGTGTCTGTGTGTGTCTGTGTCTGTGTGTCTGTGTGTCTGTTTGTGTC",
//        "TGTGTGTGTGTCTCTGTGTGTGTGTGTCTGTGTGTGTGTGTATGGCAATGCTCTGTGTGTGGCTGTGTATGTGTGGCTGTGTGTATGTGGCTGTGTGTCTGTGTCTG",
//        "TGTGTCTGTGTGTGTCTGTGTGTGTGTGTCTCTGTGTGTGTGTGTCTGTGTGTATGGCT",
//    };

//    int l = 0;
//    int m = 5;
//    int n = -4;
//    int g = -8;
//    int e = -6;

    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(l),
        m, n, g, e, q, c);

    auto graph = spoa::createGraph();

    for (const auto& it: sequences) {
        auto alignment = alignment_engine->align(it, graph);
        graph->add_alignment(alignment, it);
    }

    string consensus = graph->generate_consensus();

	if (v==1) {
		fprintf(stderr, "Consensus (%zu)\n", consensus.size());
		fprintf(stderr, "%s\n", consensus.c_str());

		vector<string> msa;
		graph->generate_multiple_sequence_alignment(msa);

		fprintf(stderr, "Multiple sequence alignment\n");
		for (const auto& it: msa) {
			fprintf(stderr, "%s\n", it.c_str());
		}
	}

	char *cons_str;
	cons_str = new char [consensus.size() + 1];
	strcpy (cons_str, consensus.c_str());
    return cons_str;
}