import argparse
import tempfile
import os
import pandas as pd
import numpy as np
from conga.tcrdist.make_10x_clones_file import make_10x_clones_file
from conga.preprocess import calc_tcrdist_matrix_cpp
# from hello import say_hello_to, parse_charptr_to_py_int


def covepitope_convert_from_10x():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        'fca_file',
        type=str,
        help='The CSV file containing filtered contig annotation',
    )
    parser.add_argument(
        '--organism',
        dest='organism',
        default='human',
        help='human, ',
    )
    parser.add_argument(
        'clones_csv',
        help='The name/path of the output CSV file containing clonotypes',
    )
    parser.add_argument(
        'edges_csv',
        help='The name/path of the output CSV file containing the pairwise TCRDist3 distances',
    )

    args = parser.parse_args()

    temp_clones_file = tempfile.NamedTemporaryFile()

    # parse filtered_contig_annotations.csv into paired clonotypes
    make_10x_clones_file(args.fca_file, args.organism, temp_clones_file.name)
    df = pd.read_csv(temp_clones_file.name, sep='\t')
    df = df.rename(
        columns={
            'clone_id': 'index',
            'subject': 'donor',
            'va_gene': 'v_a_gene',
            'ja_gene': 'j_a_gene',
            'va2_gene': 'va2',
            'ja2_gene': 'ja2',
            'vb_gene': 'v_b_gene',
            'jb_gene': 'j_b_gene',
            'cdr3a': 'cdr3_a_aa',
            'cdr3a_nucseq': 'cdr3_a_nucseq',
            'cdr3a2_nucseq': 'cdr3a2_nt',
            'cdr3b': 'cdr3_b_aa',
            'cdr3b_nucseq': 'cdr3_b_nucseq',
        }
    )
    df.to_csv(args.clones_csv, index=False)

    # tuples of tuples with tcr info
    tcrs = [((l.v_a_gene, l.j_a_gene, l.cdr3_a_aa), (l.v_b_gene, l.j_b_gene, l.cdr3_b_aa)) for l in df.itertuples()]

    D_cpp = calc_tcrdist_matrix_cpp(tcrs, args.organism).astype(np.uint)
    df_edges = pd.DataFrame(D_cpp).reset_index().melt('index').rename(columns={'variable': 'index_end', 'value': 'dist'}).assign(pw_type='pw_both')
    df_edges.to_csv(args.edges_csv, index=False)

    # alternatively we can use the slower python TCRdist calculator.
    # from conga.tcrdist.tcr_distances import TcrDistCalculator
    # tcrdister = TcrDistCalculator(organism)
    # D_py = np.array([tcrdister(x,y) for x in tcrs for y in tcrs]).reshape((len(tcrs), len(tcrs)))
    # np.savetxt("/output/sc5p_v2_hs_PBMC_10k_t_py.dist", D_py, delimiter="\t")
