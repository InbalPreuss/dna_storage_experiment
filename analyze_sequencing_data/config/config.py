import itertools


def build_config():
    shrink_dict_k_mer = {'GCTCCTTACGCGATTCAAGA': 'X1',
                         'CTAATCTGTTGAACGGATCA': 'X2',
                         'TTAGAATGGACCGCACTCAT': 'X3',
                         'CTCTAGTGAGACTCATGGCT': 'X4',
                         'CATACTTATGCAACCTGCGG': 'X5',
                         'ATATACCATCCGAGCGCATG': 'X6',
                         'GGCGGATCAAGCTTGTATCA': 'X7',
                         'CGGTTCGAAGTCAACTGTAC': 'X8',
                         'ACAGCCGTAGTTGTTACATG': 'X9',
                         'TGGTCACCAGGTTCACCGAG': 'X10',
                         'ACGCTAAGATCGATCGTGTA': 'X11',
                         'CAGCGTATCATATGGATCGC': 'X12',
                         'GGACTCTAATCGGAGCTCTG': 'X13',
                         'GGACCATTGTATCCACACAG': 'X14',
                         'CGTGAGCTGTATAACGCATT': 'X15',
                         'GCCTGCGTATCGCATAAGAA': 'X16'}

    shrink_dict_size = len(shrink_dict_k_mer)

    subset_size = 5
    bits_per_z = 12
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    x_combinations = [set(k) for k in k_mer_representative]
    z = itertools.combinations(['Z' + str(i) for i in range(1, len(x_combinations) + 1)], 1)
    z = [i[0] for i in z]

    # z = z[:2 ** bits_per_z]
    z = z
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    # k_mer_representative = list(k_mer_representative)[:2 ** bits_per_z]
    k_mer_representative = list(k_mer_representative)
    k_mer_representative_to_z = dict(zip(k_mer_representative, z))
    z_to_k_mer_representative = dict(zip(z, k_mer_representative))

    config = {
        'payload_pos': [60, 100, 140, 180],
        # 'amount_of_bc': 167,
        'amount_of_bc': 2,
        'design_len': 220,
        'payload_len': 20,
        'universal_len': 20,
        'barcode_len': 20,
        'subset_size': 5,
        'amount_of_payloads': 16,
        'amount_of_cycles': ['c1', 'c2', 'c3', 'c4'],
        'k_mer_representative_to_z': k_mer_representative_to_z,
        'z_to_k_mer_representative': z_to_k_mer_representative,
        'hamming_dist': [0, 1, 2, 3],
        'no_hamming_dist': 0,

        # 167 BC analysis
        'input_file': "analyze_sequencing_data/data/output_prefix.assembled.fastq",
        'results_most_common_file': "analyze_sequencing_data/output/csv/results_most_common.csv",
        'const_design_file': "analyze_sequencing_data/config/design.csv",
        'barcodes_design_file': "analyze_sequencing_data/config/barcodes_design.csv",
        'payload_design_file': "analyze_sequencing_data/config/payload_design.csv",
        'results_good_reads_file': "analyze_sequencing_data/output/csv/results_good_reads.csv",
        'results_good_reads_with_dist_per_cycle_file': "analyze_sequencing_data/output/csv/results_good_reads_with_dist_per_cycle_file.csv",
        'count_reads_for_each_bc_file': "analyze_sequencing_data/output/csv/count_reads_for_each_bc.csv",
        'missing_bcs_file': "analyze_sequencing_data/output/csv/missing_bc.csv",
        'output_csv_folder': "analyze_sequencing_data/output/csv/",
        'output_csv_coupon_collector_folder': "analyze_sequencing_data/output/csv/coupon_collector/",
        'foreach_bc_payload_count_file': "analyze_sequencing_data/output/csv/foreach_bc_payload_count.csv",
        'compare_design_to_experiment_results_output_file': "analyze_sequencing_data/output/csv/compare_design_to_experiment_results.csv",
        'results_table_for_two_bc_file': "analyze_sequencing_data/output/csv/results_table_for_two_bc.csv",
        'output_hist_folder': "analyze_sequencing_data/output/graphs/hist/",
        'output_folder': "analyze_sequencing_data/output/",
        'len_reads_hist_output_file': "analyze_sequencing_data/output/graphs/hist/len_reads_hist.png",
        'hist_coupon_colector_output_file': "analyze_sequencing_data/output/graphs/hist/hist_coupon_colector/",
        'output_graphs_folder': 'analyze_sequencing_data/output/graphs/',
        'output_line_graphs_folder': 'analyze_sequencing_data/output/graphs/line_graphs/',
        'sampling_rate_from_good_reads_graph': 'analyze_sequencing_data/output/graphs/line_graphs/sampling_rate_from_good_reads_graph',
        'output_heatmap_folder': 'analyze_sequencing_data/output/graphs/heatmap/',
        'heatmap_foreach_bc_and_x_count_with_most_common_file':
            "analyze_sequencing_data/output/graphs/heatmap/heatmap_foreach_bc_and_x_count_with_most_common.png",
        'hist_per_bc_file': "analyze_sequencing_data/output/graphs/hist/hist_per_bc",
        'output_hist_coupon_collector_folder': "analyze_sequencing_data/output/graphs/hist/output_hist_coupon_collector_folder/",
        'hist_foreach_bc_read_count_file': "analyze_sequencing_data/output/graphs/hist/hist_foreach_bc_read_count",
        'hist_foreach_read_count_count_bc_file': "analyze_sequencing_data/output/graphs/hist/hist_foreach_read_count_count_bc",
        'hist_foreach_error_count_of_bc_file': "analyze_sequencing_data/output/graphs/hist/hist_foreach_error_count_of_bc",
        'design_simulation_file': 'analyze_sequencing_data/data/sequence_design_file.dna',

        # output all sampling rate
        'output_all_sampling_rate_folder': "analyze_sequencing_data/output_all_sampling_rate/",
        'design_results_only_z_file': 'analyze_sequencing_data/output/csv/design/design_results_only_z.csv',
        'sampling_rate_array': [100, 50, 30, 20, 10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001],
        # 'sampling_rate_array': [0.1, 0.05, 0.01, 0.005, 0.001],

        # Params for the coupon collector problem
        't_list': [1, 2, 3, 4],
        # 'n_list': [500, 1000, 1500],
        'n_list': [500],
        'cycles_list': [[4], [1, 2, 3, 4]],
        'bc_list': [1, 2],
        'hamming_dist_for_count': [1, 2, 3],
        'hamming_dist_to_include_list': [[0], [0, 1, 2, 3]]
    }

    return config
