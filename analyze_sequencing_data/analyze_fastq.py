import ast
import heapq
from pathlib import Path

import copy
from typing import Union, Dict, List, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from matplotlib import pyplot as plt, patches
from pandas import DataFrame
from pandas.io.parsers import TextFileReader
from tqdm import tqdm
import csv

import utilities.utilities as utilities


def compare_with_errors(read, seq, max_dist=3):
    return (sum(r != s for r, s in zip(read, seq)) <= max_dist), sum(r != s for r, s in zip(read, seq))


def verify_const_universal(read, const_design):
    for s_idx, s in const_design.iterrows():
        pos = int(s['Pos'])
        is_no_errors, dist = compare_with_errors(read[pos:pos + 20], s['Seq'])
        if not is_no_errors:
            return False
    return True


def verify_const_universal_and_reverse_complement(read, const_design):
    if verify_const_universal(read, const_design):
        return read
    return SeqIO.SeqRecord("None")


class AnalyzeFastqData:
    def __init__(self, input_file: Union[Path, str],
                 const_design_file: Union[Path, str],
                 payload_design_file: Union[Path, str],
                 barcodes_design_file: Union[Path, str],
                 len_reads_hist_output_file: Union[Path, str],
                 results_good_reads_file: Union[Path, str],
                 results_good_reads_with_dist_per_cycle_file: Union[Path, str],
                 results_most_common_file: Union[Path, str],
                 design_simulation_file: Union[Path, str],
                 compare_design_to_experiment_results_output_file: Union[Path, str],
                 foreach_bc_payload_count_file: Union[Path, str],
                 heatmap_foreach_bc_and_x_count_with_most_common_file: Union[Path, str],
                 count_reads_for_each_bc_file: Union[Path, str],
                 missing_bcs_file: Union[Path, str],
                 output_hist_folder: Union[Path, str],
                 output_folder: Union[Path, str],
                 output_graphs_folder: Union[Path, str],
                 output_csv_folder: Union[Path, str],
                 output_heatmap_folder: Union[Path, str],
                 hist_per_bc_file: Union[Path, str],
                 hist_foreach_bc_read_count_file: Union[Path, str],
                 hist_foreach_error_count_of_bc_file: Union[Path, str],
                 hist_foreach_read_count_count_bc_file: Union[Path, str],
                 len_reads_to_retrieve: int,
                 amount_of_bc: int,
                 barcode_len: int,
                 subset_size: int,
                 amount_of_payloads: int,
                 z_to_k_mer_representative: Dict,
                 k_mer_representative_to_z: Dict,
                 payload_pos: List,
                 sampling_rate_from_good_reads_graph: Union[Path, str],
                 output_line_graphs_folder: Union[Path, str],
                 t_list: List,
                 n_list: List,
                 cycles_list: List[List],
                 bc_list: List,
                 output_csv_coupon_collector_folder: Union[Path, str],
                 output_hist_coupon_collector_folder: Union[Path, str],
                 amount_of_cycles: List
                 ):
        self.input_file = input_file
        self.const_design_file = const_design_file
        self.payload_design_file = payload_design_file
        self.barcodes_design_file = barcodes_design_file
        self.results_good_reads_file = results_good_reads_file
        self.results_good_reads_with_dist_per_cycle_file = results_good_reads_with_dist_per_cycle_file
        self.len_reads_hist_output_file = len_reads_hist_output_file
        self.results_most_common_file = results_most_common_file
        self.design_simulation_file = design_simulation_file
        self.compare_design_to_experiment_results_output_file = compare_design_to_experiment_results_output_file
        self.foreach_bc_payload_count_file = foreach_bc_payload_count_file
        self.heatmap_foreach_bc_and_x_count_with_most_common_file = heatmap_foreach_bc_and_x_count_with_most_common_file
        self.count_reads_for_each_bc_file = count_reads_for_each_bc_file
        self.missing_bcs_file = missing_bcs_file
        self.output_hist_folder = output_hist_folder
        self.output_folder = output_folder
        self.output_csv_folder = output_csv_folder
        self.output_graphs_folder = output_graphs_folder
        self.output_heatmap_folder = output_heatmap_folder
        self.hist_per_bc_file = hist_per_bc_file
        self.hist_foreach_bc_read_count_file = hist_foreach_bc_read_count_file
        self.hist_foreach_error_count_of_bc_file = hist_foreach_error_count_of_bc_file
        self.hist_foreach_read_count_count_bc_file = hist_foreach_read_count_count_bc_file
        self.len_reads_to_retrieve = len_reads_to_retrieve
        self.amount_of_bc = amount_of_bc
        self.barcode_len = barcode_len
        self.subset_size = subset_size
        self.amount_of_payloads = amount_of_payloads
        self.z_to_k_mer_representative = z_to_k_mer_representative
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.payload_pos = payload_pos
        self.sampling_rate_from_good_reads_graph = sampling_rate_from_good_reads_graph
        self.output_line_graphs_folder = output_line_graphs_folder
        self.t_list = t_list
        self.n_list = n_list
        self.cycles_list = cycles_list
        self.bc_list = bc_list
        self.output_csv_coupon_collector_folder = output_csv_coupon_collector_folder
        self.output_hist_coupon_collector_folder = output_hist_coupon_collector_folder
        self.amount_of_cycles = amount_of_cycles

    # Verify universal

    def identify_oligo(self, read, payload_design, barcodes_design):
        bc = self.identify_barcode(read, barcodes_design)
        payload, res_dist = self.identify_payload(read, payload_design)
        payload.insert(0, bc)
        return payload, res_dist

    # Identify which payload in each position
    def identify_payload(self, read, payload_design):
        res = [0, ] * len(self.payload_pos)
        res_dist = [0, ] * len(self.payload_pos)
        for p_idx, pos in enumerate(self.payload_pos):
            for s_idx, s in payload_design.iterrows():
                is_no_errors, dist = compare_with_errors(read[pos:pos + 20], s['Seq'])
                if is_no_errors:
                    res[p_idx] = s_idx
                    res_dist[p_idx] = dist
                    break
        return res, res_dist

    # Identify which barcode in each position
    def identify_barcode(self, read, barcodes_design):
        pos = self.barcode_len
        for s_idx, s in barcodes_design.iterrows():
            is_no_errors, dist = compare_with_errors(read[pos:pos + 20], s['Seq'])
            if is_no_errors:
                return s_idx
        return 0

    def open_fastq(self) -> List[str]:
        with open(self.input_file, 'r') as inf:
            reads = list(SeqIO.parse(inf, 'fastq'))
        return reads

    def reads_len_hist(self, reads: List[str]) -> None:
        len_reads = len(reads)
        print(f'reads = {len_reads}')
        lens = [len(r) for r in reads]
        plt.hist(lens, bins=50)
        plt.xlabel('length')
        plt.ylabel('# reads')
        plt.savefig(self.len_reads_hist_output_file)
        plt.close()

    def upload_design(self) -> Tuple[Union[TextFileReader, DataFrame],
                                     Union[TextFileReader, DataFrame],
                                     Union[TextFileReader, DataFrame]]:
        const_design_pd = pd.read_csv(self.const_design_file, index_col=0, header=0)
        payload_design_pd = pd.read_csv(self.payload_design_file, index_col=0, header=0)
        barcodes_design_pd = pd.read_csv(self.barcodes_design_file, index_col=0, header=0)

        return const_design_pd, payload_design_pd, barcodes_design_pd

    def retrieve_reads_in_specific_len(self, reads: List[str]) -> List[str]:
        reads = [r for r in reads if len(r) == self.len_reads_to_retrieve]
        print(len(reads))
        return reads

    def reads_results_to_csv(self, reads: List[str],
                             const_design: pd.DataFrame,
                             payload_design: pd.DataFrame,
                             barcodes_design: pd.DataFrame) -> None:
        res = list()
        res_with_dist = list()
        failed = 0
        read_idx = 0
        none = SeqIO.SeqRecord("None")

        with open(self.results_good_reads_file, "ab") as f:
            cols_names = [['bc', 'c1', 'c2', 'c3', 'c4']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        with open(self.results_good_reads_with_dist_per_cycle_file, "ab") as f:
            cols_names = [['bc', 'c1', 'c2', 'c3', 'c4', 'd1', 'd2', 'd3', 'd4']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        for read_idx, read in enumerate(reads):
            if read_idx % 1000 == 999:
                print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
                with open(self.results_good_reads_file, "ab") as f:
                    np.savetxt(f, res, fmt='%i', delimiter=",")
                with open(self.results_good_reads_with_dist_per_cycle_file, "ab") as f:
                    np.savetxt(f, res_with_dist, fmt='%i', delimiter=",")
                res = list()
                res_with_dist = list()
            read = verify_const_universal_and_reverse_complement(read, const_design)
            if read.seq == none.seq:
                failed += 1
                continue

            # Results for good reads
            payload, res_dist = self.identify_oligo(read, payload_design, barcodes_design)
            res.append(payload)

            payload_copy = copy.deepcopy(payload)
            payload_copy.extend(res_dist)
            res_with_dist.append(payload_copy)

            if res[-1].__contains__(0):
                failed += 1

        print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
        with open(self.results_good_reads_file, "ab") as f:
            np.savetxt(f, res, fmt='%i', delimiter=",")
        with open(self.results_good_reads_with_dist_per_cycle_file, "ab") as f:
            np.savetxt(f, res_with_dist, fmt='%i', delimiter=",")

    def compare_payloads_of_two_bc(self, bc1, bc2):
        df_foreach_bc_payload_count = pd.read_csv(self.foreach_bc_payload_count_file)
        for cycle in self.amount_of_cycles:
            bc1_cycle_i = eval(df_foreach_bc_payload_count[str(bc1)][0])[cycle]
            bc2_cycle_i = eval(df_foreach_bc_payload_count[str(bc2)][0])[cycle]

            # Normalize values in the dictionaries
            total1 = sum(bc1_cycle_i.values())
            bc1_cycle_i = {k: v / total1 for k, v in bc1_cycle_i.items()}

            total2 = sum(bc2_cycle_i.values())
            bc2_cycle_i = {k: v / total2 for k, v in bc2_cycle_i.items()}

            plt.bar(bc1_cycle_i.keys(), bc1_cycle_i.values(), color='#4169E1')
            plt.bar(bc2_cycle_i.keys(), bc2_cycle_i.values(), color='gray', alpha=0.7)
            amount_of_payloads = self.amount_of_payloads + 1
            plt.xticks(range(amount_of_payloads))
            plt.xlim(0.5, amount_of_payloads)
            plt.ylim(0, 0.5)
            plt.title('bc=[' + str(bc1) + ',' + str(bc2) + '], cycle=' + cycle)
            utilities.is_dir_exists(self.hist_per_bc_file)
            plt.savefig(self.hist_per_bc_file + '/bc=[' + str(bc1) + ',' + str(bc2) + ']_cycle=' + cycle + '_hist.png')
            plt.close()

    def hist_per_bc(self, dict_bc_payload, bc, payload):
        plt.bar(dict_bc_payload.keys(), dict_bc_payload.values())
        amount_of_payloads = self.amount_of_payloads + 1
        plt.xticks(range(amount_of_payloads))
        plt.xlim(0.5, amount_of_payloads)
        plt.title('bc=' + str(bc) + 'payload=' + payload)
        utilities.is_dir_exists(self.hist_per_bc_file)
        plt.savefig(self.hist_per_bc_file + '/bc=' + str(bc) + '_payload=' + payload + '_hist.png')
        plt.close()

    def analyze_results_good_reads(self) -> Dict:
        # read your csv to a dataframe

        df = pd.read_csv(self.results_good_reads_file)
        dict_bc = {}

        amount_of_bc = self.amount_of_bc + 1
        for i in range(amount_of_bc):
            dict_bc_i = {}
            for payload in self.amount_of_cycles:
                dict_p = {
                    payload: {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0,
                              14: 0,
                              15: 0, 16: 0}}
                dict_bc_i.update(dict_p)
            dict_bc[i] = dict_bc_i
        df_bc = df.sort_values(["bc", "c1", "c2", "c3", "c4"], ascending=True, key=np.sin)

        df_bc.fillna(0)
        for row_idx, row in tqdm(df_bc.iterrows()):
            for location_payload, payload in row[1:].items():
                try:
                    dict_bc[int(row['bc'])][location_payload][int(payload)] += 1

                except:
                    print(f"An exception occurred {row_idx}, row={row}")
                    continue

                dict_bc[int(row['bc'])][location_payload][int(payload)] += 1

        utilities.write_dict_to_csv(dict_bc, self.foreach_bc_payload_count_file)

        return dict_bc

    def most_common_for_each_bc(self, dict_bc: Dict) -> None:
        dict_bc_most_common = {}
        amount_of_bc = self.amount_of_bc + 1

        # init dict_bc_most_common_i
        for i in range(amount_of_bc):
            dict_bc_most_common_i = {}
            for payload in self.amount_of_cycles:
                dict_p = {payload: []}
                dict_bc_most_common_i.update(dict_p)
            dict_bc_most_common[i] = dict_bc_most_common_i

        # find the config['subset_size'] most common
        for bc_i in range(1, amount_of_bc):
            for payload in self.amount_of_cycles:
                del dict_bc[bc_i][payload][0]

                most_common = heapq.nlargest(self.subset_size, dict_bc[bc_i][payload],
                                             key=dict_bc[bc_i][payload].get)
                dict_bc_most_common[bc_i][payload] = most_common

                print(f'bc = {bc_i}, cycle = {payload}, 5 most common = {most_common}')
                self.hist_per_bc(dict_bc_payload=dict_bc[bc_i][payload], bc=bc_i, payload=payload)

        utilities.write_dict_to_csv(dict_bc_most_common, self.results_most_common_file)

    def convert_most_common_to_letters_in_new_alphabet(self) -> Dict:
        df = pd.read_csv(self.results_most_common_file, index_col=0)
        dict_most_common = df.to_dict("list")
        result_payload = {}
        print(dict_most_common)

        for bc_i in range(1, (self.amount_of_bc + 1)):
            result_payload[bc_i] = []
            for payload in self.amount_of_cycles:
                bc_i_str = str(bc_i)
                d = ast.literal_eval(dict_most_common[bc_i_str][0])
                d = d[payload]
                dict_convert_to_x = []
                for i in d:
                    dict_convert_to_x.append('X' + str(i))
                dict_convert_to_x = tuple(utilities.sorted_human(dict_convert_to_x))

                try:
                    z = self.k_mer_representative_to_z[dict_convert_to_x]
                except KeyError:
                    z = 'Z1'  # The X tuple is out of range
                result_payload[bc_i].append({z: dict_convert_to_x})

        return result_payload

    def compare_most_common_to_design(self, result_payload: Dict) -> None:
        with open(self.compare_design_to_experiment_results_output_file, "ab") as f:
            cols_names = [
                ['bc', 'design_z', 'design_x', 'experiment_results_z', 'experiment_results_x', 'is_retrieved_all_z',
                 'count_all_z_mm']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        with open(self.design_simulation_file, 'r') as inf:
            for bc_i_idx, bc_i_design in enumerate(inf):
                count_all_z_mm = 0

                array_bc_i_design = bc_i_design.split("\n")
                array_bc_i_design = array_bc_i_design[0].split(",")
                for z_i in zip(array_bc_i_design[1:self.subset_size], result_payload[bc_i_idx + 1]):
                    if z_i[0] != list(z_i[1])[0]:
                        count_all_z_mm += 1
                with open(self.compare_design_to_experiment_results_output_file, "ab") as f:
                    bc = bc_i_idx + 1
                    design_z = " ".join(array_bc_i_design[1:self.subset_size])
                    design_x = str(
                        {z: self.z_to_k_mer_representative[z] for z in array_bc_i_design[1:self.subset_size]}).replace(
                        ",",
                        " ")
                    experiment_results_z = " ".join([list(z.keys())[0] for z in result_payload[bc_i_idx + 1]])
                    experiment_results_x = str(result_payload[bc_i_idx + 1]).replace(",", " ")
                    is_retrieved_all_z = (count_all_z_mm == 0)
                    count_all_z_mm = count_all_z_mm
                    cols = [[bc, design_z, design_x, experiment_results_z, experiment_results_x, is_retrieved_all_z,
                             count_all_z_mm]]
                    np.savetxt(f, cols, fmt='%s', delimiter=",")

    def parse_design(self, design_str):
        # Remove curly braces
        design_str = design_str[1:-1]
        # Split into individual entries
        entries = design_str.split(")  ")
        design_dict = {}
        for entry in entries:
            # Split the key from the values
            key, values_str = entry.split(": ")
            # Remove quotation marks from the key
            key = key.strip("'")
            # Split the values string into individual values and remove the 'X' prefix
            values_str = values_str.replace("(", "").replace(")", "").replace("X", "").replace("  ", ",").replace("'",
                                                                                                                  "").replace(
                " ", "").split(",")
            values = {int(number) for number in values_str}
            # Add the key-value pair to the dictionary
            design_dict[key] = values
        return design_dict

    def sample_till_K_unique(self, df: pd.DataFrame, column_idx: List, K: int, t: int):
        unique_values = {}
        sample_rows = []
        sampled_indices = []
        column = 'c' + str(column_idx)
        skipped_zero_count = 0

        while len([k for k, v in unique_values.items() if v >= t]) < K:
            sample = df.drop(sampled_indices).sample()
            sampled_indices.append(sample.index[0])
            value = sample.iloc[0][column]
            if value != 0:
                unique_values[value] = unique_values.get(value, 0) + 1
            else:
                skipped_zero_count += 1
            sample_rows.append(sample)

        # Select K unique values that each appeared at least t times
        selected_values = set([k for k, v in unique_values.items() if v >= t][:K])

        # Create a DataFrame with only the selected values
        result_df = pd.concat([row for row in sample_rows if row.iloc[0][column] in unique_values])

        return result_df, selected_values, unique_values, skipped_zero_count

    def random_sample_and_compare(self, data_df, design_df, n, K, bc, cycles_list, t):
        sampling_results = []
        design_x = design_df.loc[design_df['bc'] == bc, 'design_x'].iloc[0]

        for _ in range(n):
            sampling_results_per_cycles_option = []
            for ci in cycles_list:
                sampled_df, selected_values, unique_values, skipped_zero_count = self.sample_till_K_unique(
                    data_df.loc[data_df['bc'] == bc],
                    ci, K, t)
                is_design_equal_sampled_data = list(design_x.values())[ci - 1] == selected_values
                sampling_results_per_cycles_option.append((len(sampled_df), is_design_equal_sampled_data,
                                                           list(design_x.values())[ci - 1], selected_values, ci,
                                                           unique_values))

            total_samples = sum(item[0] for item in sampling_results_per_cycles_option)
            status = all(item[1] for item in sampling_results_per_cycles_option)
            design = [item[2] for item in sampling_results_per_cycles_option]
            selected_values = [item[3] for item in sampling_results_per_cycles_option]
            ci = [item[4] for item in sampling_results_per_cycles_option]
            unique_values = [item[5] for item in sampling_results_per_cycles_option]
            sampling_results.append([total_samples, status, design, selected_values, ci, unique_values])

        return sampling_results, skipped_zero_count

    def plot_results(self, results):
        xs = [r[0] for r in results]
        ys = [r[1] for r in results]
        plt.scatter(xs, ys)
        plt.xlabel('Number of samples')
        plt.ylabel('Match with design')
        plt.show()

    def the_coupon_collector_problem_with_t(self):
        K = self.subset_size

        # read csv data
        data_df = pd.read_csv(self.results_good_reads_file)

        # convert design file to dataframe
        design_df = pd.read_csv(self.compare_design_to_experiment_results_output_file)
        design_df['design_x'] = design_df['design_x'].apply(self.parse_design)

        for bc in self.bc_list:
            for cycles in list(self.cycles_list):
                for t in self.t_list:
                    for n in self.n_list:
                        # perform the operation
                        results, skipped_zero_count = self.random_sample_and_compare(data_df, design_df, n=n, K=K,
                                                                                     bc=bc,
                                                                                     cycles_list=cycles, t=t)

                        # plot the results
                        self.hist_of_the_coupon_collector_problem_with_t(results, K, t, n, cycles, bc,
                                                                         skipped_zero_count)

    def hist_of_the_coupon_collector_problem_with_t(self, results, K, t, n, cycles_list, bc, skipped_zero_count,
                                                    is_histogram_of_true_n_false=True):

        xs_results = []
        if is_histogram_of_true_n_false:
            xs_results.append([[r[0] for r in results], 'True And False'])
        else:
            xs_results.append([[r[0] for r in results if r[1] is True], True])
            xs_results.append([[r[0] for r in results if r[1] is False], False])

        for xs, status in xs_results:
            fig, ax = plt.subplots(figsize=(10, 8))

            hist_name_file = ['bc=' + str(bc)
                , 'K=' + str(K)
                , 't=' + str(t)
                , 'n=' + str(n)
                , 'cycles_list=' + str(cycles_list)]

            # Expectation calculation
            expectation_value = sum(xs) / len(xs) if xs else 0
            hist_name_file.extend(['expectation=' + str(round(expectation_value)),
                                   'skipped_zero_count=' + str(skipped_zero_count),
                                   'status=' + str(status)])

            # Plot histogram
            ax.hist(xs, bins=10, alpha=0.5, label='Number of samples')

            ax.set_xlabel('# Reads')
            ax.set_ylabel('Frequency')
            ax.legend(loc='upper right')
            ax.set_title("Sampling Rate")
            # plt.figtext(0.85, 0.85, "\n".join(hist_name_file), fontsize=8)
            # ax.annotate("\n".join(hist_name_file), xy=(0.5, 0.5), xytext=(0, 0),
            #             textcoords='figure fraction', fontsize=10,
            #             bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=2))
            plt.subplots_adjust(right=0.8)  # Adjust the bottom leaving space for your description
            fig.text(0.82, 0.5, "\n".join(hist_name_file), ha='left', va='center',
                     fontsize=10, bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=2))

            fig.savefig(self.output_hist_coupon_collector_folder + "_".join(hist_name_file))
            plt.close(fig)
            print(hist_name_file)

        file_name = self.output_csv_coupon_collector_folder + "_".join(hist_name_file) + '.csv'

        headers = ["Sample", "Status", "Design", "Sample", "Cycle", "Unique values"]

        # Write to CSV
        with open(file_name, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            for csv_row in results:
                writer.writerow(csv_row)

        print(f'Data written to {file_name}')

    def create_heatmap_with_rectangles_on_most_common(self, dict_foreach_bc_and_x_count_all_cycles_matrix: Dict,
                                                      ax) -> None:

        # remove col 0 and row 0
        del dict_foreach_bc_and_x_count_all_cycles_matrix[0]
        dict_foreach_bc_and_x_count_all_cycles_matrix.drop(index=['0'])

        # for each bc, normalize the data to 1
        dict_foreach_bc_and_x_count_all_cycles_matrix = dict_foreach_bc_and_x_count_all_cycles_matrix.div(
            dict_foreach_bc_and_x_count_all_cycles_matrix.sum(axis=1), axis=0)

        # fill nan cells with 0
        dict_foreach_bc_and_x_count_all_cycles_matrix = dict_foreach_bc_and_x_count_all_cycles_matrix.fillna(0)

        # create heatmap
        sns.heatmap(dict_foreach_bc_and_x_count_all_cycles_matrix.T)
        plt.xlabel("payloads")
        plt.ylabel("bc")
        # Add lines to distinguish payloads at each cycle
        ax.hlines(
            [self.amount_of_payloads, (self.amount_of_payloads * 2), (self.amount_of_payloads * 3)],
            *ax.get_xlim())

        df_results_most_common = pd.read_csv(self.results_most_common_file)
        del df_results_most_common['0']
        for bc_idx, cycles in df_results_most_common.items():
            if bc_idx == '0':
                continue
            dict_cycles = ast.literal_eval(cycles[0])
            cycle_idx = 1
            for cycle_name, cycle_data in dict_cycles.items():
                for one_of_most_common in cycle_data:
                    ax.add_patch(
                        patches.Rectangle(
                            (int(bc_idx),
                             ((self.amount_of_payloads * (cycle_idx - 1)) + int(one_of_most_common)) - 1),
                            1.0,
                            1.0,
                            edgecolor='red',
                            fill=False,
                            lw=0.5
                        ))
                cycle_idx += 1

        plt.savefig(self.heatmap_foreach_bc_and_x_count_with_most_common_file, dpi=400)
        plt.close()

    def heatmap_foreach_bc_and_x_count_with_most_common(self) -> None:
        df = pd.read_csv(self.foreach_bc_payload_count_file)
        dict_foreach_bc_and_x_count_str = df.to_dict("list")

        dict_foreach_bc_and_x_count_split_to_cycles = {}
        dict_foreach_bc_and_x_count_all_cycles = {}
        for bc_i_str, payloads in dict_foreach_bc_and_x_count_str.items():
            dict_foreach_bc_and_x_count_all_cycles[bc_i_str] = {}
            dict_foreach_bc_and_x_count_split_to_cycles[bc_i_str] = ast.literal_eval(
                dict_foreach_bc_and_x_count_str[bc_i_str][0])

            cycle_idx = 1
            for cycle_name, cycle_data in dict_foreach_bc_and_x_count_split_to_cycles[bc_i_str].items():
                cycle_data_new = cycle_data
                if cycle_name != 'c1':
                    cycle_idx += 1
                    cycle_data_new = {}
                    for payload_idx, payload_count in cycle_data.items():
                        if payload_idx == 0:
                            continue
                        payload_new_x = (self.amount_of_payloads * (cycle_idx - 1)) + payload_idx
                        cycle_data_new[payload_new_x] = payload_count

                dict_foreach_bc_and_x_count_all_cycles[bc_i_str].update(cycle_data_new)

        fig, ax = plt.subplots(figsize=(10, 7))

        # make the matrix bcs x payloads
        dict_foreach_bc_and_x_count_all_cycles_matrix = pd.DataFrame(dict_foreach_bc_and_x_count_all_cycles).T.fillna(0)

        self.create_heatmap_with_rectangles_on_most_common(dict_foreach_bc_and_x_count_all_cycles_matrix=
                                                           dict_foreach_bc_and_x_count_all_cycles_matrix, ax=ax)

    def find_most_common(self) -> None:
        dict_bc = self.analyze_results_good_reads()
        self.most_common_for_each_bc(dict_bc=dict_bc)
        result_payload = self.convert_most_common_to_letters_in_new_alphabet()
        self.compare_most_common_to_design(result_payload=result_payload)

        self.heatmap_foreach_bc_and_x_count_with_most_common()

    def missing_bc_to_csv(self, dict_append_missing_bc):
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        ser_append_missing_bc.to_csv(self.missing_bcs_file, mode='a', header=False)

    def create_sampling_rate_from_good_reads_graph(self) -> None:
        """
        This function takes in a list of integers and creates a graph with the x-axis being the sampling rate
        and the y-axis being the count of different values in arr[0].
        The sampling rate will be in increments of 10% starting from 0% until 100%.
        """

        df_good_reads = pd.read_csv(self.results_good_reads_file)

        sampling_rates = [i / 10 for i in range(11)]  # create a list of sampling rates from 0% to 100%
        counts = []  # list to store counts of different values in arr[0]
        # Iterate over the sampling rates
        for sampling_rate in sampling_rates:
            # Sample the dataframe with the current sampling rate
            sampled_df = df_good_reads.sample(frac=sampling_rate)

            # Count the number of different values in the "bc" column
            count = sampled_df["bc"].nunique()

            # Add the count to the list
            counts.append(count)

        # Plot the graph
        plt.plot(sampling_rates, counts)
        plt.title("sampling_rate_graph")
        plt.xlabel('Sampling Rate %')
        plt.ylabel('Count of Different Values')
        plt.savefig(self.sampling_rate_from_good_reads_graph)
        plt.close()

    def for_each_bc_count_reads_read_to_csv(self, output_file: Union[Path, str]) -> pd.DataFrame:
        with open(output_file, "ab") as f:
            cols_names = [['bc', 'count']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        df = pd.read_csv(self.results_good_reads_file)
        forth_column = df.iloc[:, 0]
        counts = forth_column.value_counts()
        count_sorted = counts.sort_index()

        dict_append_missing_bc = {}
        idx = 0
        for bc, count in count_sorted.iteritems():
            if idx != bc:
                for idx in range(idx, bc):
                    dict_append_missing_bc[idx] = 0
                    idx += 1
            idx += 1

        self.missing_bc_to_csv(dict_append_missing_bc)
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        count_sorted = count_sorted.append(ser_append_missing_bc)
        count_sorted_with_missing_bc = count_sorted.sort_index()
        count_sorted_with_missing_bc.to_csv(output_file, mode='a', header=False)

        return count_sorted_with_missing_bc

    def hist_foreach_bc_read_count(self, csv_output_file: Union[Path, str]) -> None:
        count_sorted_with_missing_bc = pd.read_csv(csv_output_file)
        x = count_sorted_with_missing_bc["count"].index
        y = count_sorted_with_missing_bc["count"]
        plt.bar(x, y, align='center')
        plt.xlabel('bc')
        plt.ylabel('count')
        plt.savefig(self.hist_foreach_bc_read_count_file)
        plt.close()

    def hist_foreach_error_count_of_bc(self) -> None:
        count_sorted_with_missing_bc = pd.read_csv(self.compare_design_to_experiment_results_output_file)
        x = count_sorted_with_missing_bc["count_all_z_mm"]
        plt.hist(x, bins=5)
        plt.xlabel('errors in bc')
        plt.ylabel('count')
        plt.savefig(self.hist_foreach_error_count_of_bc_file)
        plt.close()

    def hist_foreach_read_count_count_bc(self, csv_output_file: Union[Path, str]) -> None:
        count_sorted_with_missing_bc = pd.read_csv(csv_output_file)
        x = count_sorted_with_missing_bc["count"]
        plt.hist(x, bins=20)
        plt.xlabel('# reads')
        plt.ylabel('bc')
        plt.savefig(self.hist_foreach_read_count_count_bc_file)
        plt.close()

    def for_each_bc_count_reads_read(self, csv_output_file: Union[Path, str]) -> None:
        self.for_each_bc_count_reads_read_to_csv(csv_output_file)
        self.hist_foreach_bc_read_count(csv_output_file)
        self.hist_foreach_read_count_count_bc(csv_output_file=csv_output_file)
        self.hist_foreach_error_count_of_bc()

    def create_folders(self) -> None:
        utilities.is_dir_exists(self.output_hist_folder)
        utilities.is_dir_exists(self.output_folder)
        utilities.is_dir_exists(self.output_heatmap_folder)
        utilities.is_dir_exists(self.output_graphs_folder)
        utilities.is_dir_exists(self.output_csv_folder)
        utilities.is_dir_exists(self.output_line_graphs_folder)
        utilities.is_dir_exists(self.output_csv_coupon_collector_folder)
        utilities.is_dir_exists(self.output_hist_coupon_collector_folder)

    def run(self):
        self.create_folders()

        # upload design
        const_design_pd, payload_design_pd, barcodes_design_pd = self.upload_design()

        # reads
        reads = self.open_fastq()

        # reads len showed in histogram
        self.reads_len_hist(reads=reads)

        # good reads with len  220
        good_reads = self.retrieve_reads_in_specific_len(reads=reads)

        # Write the good reads with len 220 to results_good_reads.csv
        self.reads_results_to_csv(reads=good_reads,
                                  const_design=const_design_pd,
                                  payload_design=payload_design_pd,
                                  barcodes_design=barcodes_design_pd)

        # Find most common for each bc and for every cycle in that bc in results of good reads
        self.find_most_common()

        # For each bc count amount of reads sequenced
        self.for_each_bc_count_reads_read(csv_output_file=self.count_reads_for_each_bc_file)

        # Create graph with sampling rate
        self.create_sampling_rate_from_good_reads_graph()

        # Compare Payloads_of_two_bc
        self.compare_payloads_of_two_bc(bc1=1, bc2=2)

        # We want to find
        self.the_coupon_collector_problem_with_t()
