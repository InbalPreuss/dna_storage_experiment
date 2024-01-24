from analyze_sequencing_data.analyze_fastq import AnalyzeFastqData
from analyze_sequencing_data.config import config

if __name__ == '__main__':

    # Analyze fastq results from sequencing
    config = config.build_config()
    analyze_fastq_data = AnalyzeFastqData(input_file=config['input_file'],
                                          const_design_file=config['const_design_file'],
                                          payload_design_file=config['payload_design_file'],
                                          barcodes_design_file=config['barcodes_design_file'],
                                          len_reads_hist_output_file=config['len_reads_hist_output_file'],
                                          results_good_reads_file=config['results_good_reads_file'],
                                          results_good_reads_with_dist_per_cycle_file=config['results_good_reads_with_dist_per_cycle_file'],
                                          results_most_common_file=config['results_most_common_file'],
                                          design_simulation_file=config['design_simulation_file'],
                                          compare_design_to_experiment_results_output_file=config[
                                              'compare_design_to_experiment_results_output_file'],
                                          foreach_bc_payload_count_file=config['foreach_bc_payload_count_file'],
                                          heatmap_foreach_bc_and_x_count_with_most_common_file=config[
                                              'heatmap_foreach_bc_and_x_count_with_most_common_file'],
                                          count_reads_for_each_bc_file=config['count_reads_for_each_bc_file'],
                                          missing_bcs_file=config['missing_bcs_file'],
                                          output_hist_folder=config['output_hist_folder'],
                                          output_folder=config['output_folder'],
                                          output_graphs_folder=config['output_graphs_folder'],
                                          output_csv_folder=config['output_csv_folder'],
                                          output_heatmap_folder=config['output_heatmap_folder'],
                                          hist_per_bc_file=config['hist_per_bc_file'],
                                          hist_foreach_bc_read_count_file=config['hist_foreach_bc_read_count_file'],
                                          hist_foreach_error_count_of_bc_file=config[
                                              'hist_foreach_error_count_of_bc_file'],
                                          hist_foreach_read_count_count_bc_file=config[
                                              'hist_foreach_read_count_count_bc_file'],
                                          len_reads_to_retrieve=config['design_len'],
                                          amount_of_bc=config['amount_of_bc'],
                                          barcode_len=config['barcode_len'],
                                          subset_size=config['subset_size'],
                                          amount_of_payloads=config['amount_of_payloads'],
                                          z_to_k_mer_representative=config['z_to_k_mer_representative'],
                                          k_mer_representative_to_z=config['k_mer_representative_to_z'],
                                          payload_pos=config['payload_pos'],
                                          sampling_rate_from_good_reads_graph=config['sampling_rate_from_good_reads_graph'],
                                          output_line_graphs_folder=config['output_line_graphs_folder'],
                                          t_list=config['t_list'],
                                          n_list=config['n_list'],
                                          cycles_list=config['cycles_list'],
                                          bc_list=config['bc_list'],
                                          output_csv_coupon_collector_folder=config['output_csv_coupon_collector_folder'],
                                          output_hist_coupon_collector_folder=config['output_hist_coupon_collector_folder'],
                                          amount_of_cycles=config['amount_of_cycles'],
                                          hamming_dist=config['hamming_dist'],
                                          hamming_dist_for_count=config['hamming_dist_for_count'],
                                          hamming_dist_to_include_list=config['hamming_dist_to_include_list'],
                                          results_table_for_two_bc_file=config['results_table_for_two_bc_file'],
                                          output_all_sampling_rate_folder=config['output_all_sampling_rate_folder'],
                                          design_results_only_z_file=config['design_results_only_z_file'],
                                          sampling_rate_array=config['sampling_rate_array'],
                                          )
    analyze_fastq_data.run()
