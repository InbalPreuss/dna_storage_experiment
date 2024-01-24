import ast
import collections
import heapq
import os
import shutil
from pathlib import Path

import copy
from typing import Union, Dict, List, Tuple, Any

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
                 amount_of_cycles: List,
                 hamming_dist: List,
                 hamming_dist_for_count: List,
                 hamming_dist_to_include_list: List,
                 results_table_for_two_bc_file: Union[Path, str],
                 output_all_sampling_rate_folder: Union[Path, str],
                 design_results_only_z_file: Union[Path, str],
                 sampling_rate_array: List,
                 ):
        self.output_all_sampling_rate_folder = output_all_sampling_rate_folder
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
        self.hamming_dist = hamming_dist
        self.hamming_dist_for_count = hamming_dist_for_count
        self.hamming_dist_to_include_list = hamming_dist_to_include_list
        self.results_table_for_two_bc_file = results_table_for_two_bc_file
        self.design_results_only_z_file = design_results_only_z_file
        self.sampling_rate_array = sampling_rate_array

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

            plt.bar(bc1_cycle_i.keys(), bc1_cycle_i.values(), color='#4169E1', label='BC1')
            plt.bar(bc2_cycle_i.keys(), bc2_cycle_i.values(), color='gray', alpha=0.7, label='BC2')
            amount_of_payloads = self.amount_of_payloads + 1
            plt.xticks(range(amount_of_payloads))
            plt.xlabel("Payloads")
            plt.ylabel("Fraction of Reads")
            plt.xlim(0.5, amount_of_payloads)
            plt.ylim(0, 0.5)
            plt.legend()
            plt.title('BC=[' + str(bc1) + ',' + str(bc2) + '], Cycle=' + cycle)
            utilities.is_dir_exists(self.hist_per_bc_file)
            plt.savefig(self.hist_per_bc_file + '/BC=[' + str(bc1) + ',' + str(bc2) + '] Cycle=' + cycle + '_hist.svg')
            plt.close()

    def create_results_table_for_two_bc(self, bc1, bc2):
        # 1. Initialize the multi-level columns DataFrame
        columns = pd.MultiIndex.from_product([['bc' + str(bc1), 'bc' + str(bc2)], ['c1', 'c2', 'c3', 'c4']],
                                             names=['barcode', 'cycle'])
        df = pd.DataFrame(columns=columns, index=['X' + str(i) for i in range(1, 17)])
        # Assuming you already have this:
        df_foreach_bc_payload_count = pd.read_csv(self.foreach_bc_payload_count_file)
        # 2. Process the loaded CSV file for bc values for each cycle
        for cycle in self.amount_of_cycles:
            bc1_cycle_i = eval(df_foreach_bc_payload_count[str(bc1)][0])[cycle]
            bc2_cycle_i = eval(df_foreach_bc_payload_count[str(bc2)][0])[cycle]
            # Normalize values in the dictionaries
            total1 = sum(bc1_cycle_i.values())
            bc1_cycle_i = {k: v / total1 for k, v in bc1_cycle_i.items()}
            total2 = sum(bc2_cycle_i.values())
            bc2_cycle_i = {k: v / total2 for k, v in bc2_cycle_i.items()}
            # 3. Insert values into the DataFrame. Here's an example:
            for key in bc1_cycle_i:
                df.loc['X' + str(key), ('bc1', cycle)] = bc1_cycle_i[key]
            for key in bc2_cycle_i:
                df.loc['X' + str(key), ('bc2', cycle)] = bc2_cycle_i[key]

        df = df.drop(['X0'])
        df.to_csv(self.results_table_for_two_bc_file)
        print(df)

    def hist_per_bc(self, dict_bc_payload, bc, payload):
        # Normalize dict_bc_payload if possible
        total = sum(dict_bc_payload.values())
        if total != 0:
            dict_bc_payload = {k: v / total for k, v in dict_bc_payload.items()}

        plt.bar(dict_bc_payload.keys(), dict_bc_payload.values())
        amount_of_payloads = self.amount_of_payloads + 1
        plt.xticks(range(amount_of_payloads))
        plt.xlim(0.5, amount_of_payloads)
        plt.ylim(0, 0.5)
        plt.xlabel("Payloads")
        plt.ylabel("Fraction of Reads")
        plt.title('BC=' + str(bc) + 'Cycle=' + payload)
        utilities.is_dir_exists(self.hist_per_bc_file)
        plt.savefig(self.hist_per_bc_file + '/bc=' + str(bc) + '_cycle=' + payload + '_hist.svg')
        plt.close()

    def analyze_results_good_reads(self, input_csv_path: Union[Path, str],
                                   dict_to_csv: Union[Path, str]) -> Dict:
        # read your csv to a dataframe

        df = pd.read_csv(input_csv_path)
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

        utilities.write_dict_to_csv(dict_bc, dict_to_csv)

        return dict_bc

    def most_common_for_each_bc(self, dict_bc: Dict, dict_to_csv_path: Union[Path, str]) -> None:
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

                print(f'Bc = {bc_i}, Cycle = {payload}, 5 most common = {most_common}')
                self.hist_per_bc(dict_bc_payload=dict_bc[bc_i][payload], bc=bc_i, payload=payload)

        utilities.write_dict_to_csv(dict=dict_bc_most_common, csv_path=dict_to_csv_path)

    def convert_most_common_to_letters_in_new_alphabet(self, results_most_common_file: Union[Path, str]) -> tuple[
        dict[int, list[Any]], dict[int, list[Any]]]:
        df = pd.read_csv(results_most_common_file, index_col=0)
        dict_most_common = df.to_dict("list")
        result_payload = {}
        result_only_z = {}
        print(dict_most_common)

        for bc_i in range(1, (self.amount_of_bc + 1)):
            result_payload[bc_i] = []
            result_only_z[bc_i] = []
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
                result_only_z[bc_i].append(z)

        return result_payload, result_only_z

    def compare_most_common_to_design(self, result_payload: Dict,
                                      compare_design_to_experiment_results_output_file: Union[Path, str],
                                      design_simulation_file: Union[Path, str]) -> None:
        with open(compare_design_to_experiment_results_output_file, "ab") as f:
            cols_names = [
                ['bc', 'design_z', 'design_x', 'experiment_results_z', 'experiment_results_x', 'is_retrieved_all_z',
                 'count_all_z_mm']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        with open(design_simulation_file, 'r') as inf:
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
            try:
                key, values_str = entry.split(": ")
            except:
                print(f'entry={entry}')
            # Remove quotation marks from the key
            try:
                key = key.strip("'")
            except:
                print(key)
            # Split the values string into individual values and remove the 'X' prefix
            values_str = values_str.replace("(", "").replace(")", "").replace("X", "").replace("  ", ",").replace("'",
                                                                                                                  "").replace(
                " ", "").split(",")
            values = {int(number) for number in values_str}
            # Add the key-value pair to the dictionary
            design_dict[key] = values
        return design_dict

    def sample_till_K_unique(self, df: pd.DataFrame, column_idx: List, K: int, t: int, hamming_dist_to_include: List):
        unique_values = {}
        sample_rows = []
        sampled_indices = []
        column_cycle_name = 'c' + str(column_idx)
        column_dist_name = 'd' + str(column_idx)
        skipped_zero_count = 0
        hamming_dist_bigger_then_0 = 0

        while len([k for k, v in unique_values.items() if v >= t]) < K:
            # Choose sample from file
            sample_read = df.drop(sampled_indices).sample()
            # Find the dist value
            dist_value = sample_read.iloc[0][column_dist_name]

            value = sample_read.iloc[0][column_cycle_name]

            if dist_value not in hamming_dist_to_include or value == 0:
                skipped_zero_count += 1

            if dist_value in self.hamming_dist_for_count:
                hamming_dist_bigger_then_0 += 1

            sampled_indices.append(sample_read.index[0])
            if value != 0:
                unique_values[value] = unique_values.get(value, 0) + 1

            sample_rows.append(sample_read)

        # Select K unique values that each appeared at least t times
        selected_values = set([k for k, v in unique_values.items() if v >= t][:K])

        # Create a DataFrame with only the selected values
        result_df = pd.concat([row for row in sample_rows if row.iloc[0][column_cycle_name] in unique_values])

        return result_df, selected_values, unique_values, skipped_zero_count, hamming_dist_bigger_then_0

    def random_sample_and_compare(self, data_df: pd.DataFrame, design_df: pd.DataFrame, n: int, K: int, bc: int,
                                  cycles_list: List, t: int, hamming_dist_to_include: List):

        total_samples_list = []
        status_list = []
        design_list = []
        selected_values_list = []
        ci_list = []
        unique_values_list = []
        hamming_dist_bigger_then_0_list = []
        skipped_zero_count_list = []
        design_x = design_df.loc[design_df['bc'] == bc, 'design_x'].iloc[0]

        for _ in range(n):
            sampled_count_all_cycles = 0
            status_all_cycles_bool = True
            design_list_all_cycles = []
            selected_values_all_cycles = []
            unique_values_all_cycles = []
            hamming_dist_bigger_then_0_count_all_cycles = 0
            skipped_zero_count_all_cycles = 0
            for ci in cycles_list:
                sampled_df, selected_values, unique_values, skipped_zero_count, hamming_dist_bigger_then_0 = self.sample_till_K_unique(
                    data_df.loc[data_df['bc'] == bc],
                    ci, K, t, hamming_dist_to_include)
                is_design_equal_sampled_data = list(design_x.values())[ci - 1] == selected_values

                sampled_count_all_cycles += len(sampled_df)
                status_all_cycles_bool = status_all_cycles_bool and is_design_equal_sampled_data
                hamming_dist_bigger_then_0_count_all_cycles += hamming_dist_bigger_then_0
                skipped_zero_count_all_cycles += skipped_zero_count
                design_list_all_cycles.extend([list(design_x.values())[ci - 1]])
                selected_values_all_cycles.extend([selected_values])
                unique_values_all_cycles.extend([list(unique_values.keys())])

            total_samples_list.append([sampled_count_all_cycles])
            status_list.append(status_all_cycles_bool)
            design_list.extend([design_list_all_cycles])
            selected_values_list.extend(selected_values)
            ci_list.append(cycles_list)
            unique_values_list.extend([unique_values_all_cycles])
            hamming_dist_bigger_then_0_list.append(hamming_dist_bigger_then_0_count_all_cycles)
            skipped_zero_count_list.append(skipped_zero_count_all_cycles)

        return total_samples_list, status_list, design_list, \
               selected_values_list, ci_list, unique_values_list, \
               hamming_dist_bigger_then_0_list, skipped_zero_count_list

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
        data_df = pd.read_csv(self.results_good_reads_with_dist_per_cycle_file)

        # convert design file to dataframe
        design_df = pd.read_csv(self.compare_design_to_experiment_results_output_file)
        design_df['design_x'] = design_df['design_x'].apply(self.parse_design)

        for bc in self.bc_list:
            for cycles in list(self.cycles_list):
                for t in self.t_list:
                    for n in self.n_list:
                        for hamming_dist_to_include in self.hamming_dist_to_include_list:
                            # perform the operation
                            total_samples_list, status_list, design_list, \
                            selected_values_list, ci_list, unique_values_list, \
                            hamming_dist_bigger_then_0_list, skipped_zero_count_list \
                                = self.random_sample_and_compare(data_df, design_df, n=n, K=K,
                                                                 bc=bc,
                                                                 cycles_list=cycles, t=t,
                                                                 hamming_dist_to_include=hamming_dist_to_include)

                            # plot the results
                            self.hist_of_the_coupon_collector_problem_with_t(K, t, n, cycles, bc, total_samples_list,
                                                                             status_list, design_list,
                                                                             selected_values_list,
                                                                             ci_list, unique_values_list,
                                                                             hamming_dist_bigger_then_0_list,
                                                                             skipped_zero_count_list,
                                                                             hamming_dist_to_include)

    def hist_of_the_coupon_collector_problem_with_t(self, K, t, n, cycles_list, bc,
                                                    total_samples_list, status_list, design_list,
                                                    selected_values_list, ci_list, unique_values_list,
                                                    hamming_dist_bigger_then_0_list, skipped_zero_count_list,
                                                    hamming_dist_to_include,
                                                    is_histogram_of_true_n_false=True) -> str:

        xs_results = []
        results = list(zip(total_samples_list, status_list, design_list,
                           selected_values_list, ci_list, unique_values_list,
                           hamming_dist_bigger_then_0_list, skipped_zero_count_list))
        if is_histogram_of_true_n_false:
            xs_results.append([[r[0] for r in results], 'T&F'])

            true_false_values = [item[1] for item in results]

            # Count the number of True and False values
            num_true = true_false_values.count(True)
            num_false = true_false_values.count(False)

            # Calculate the ratio
            true_false_ratio = num_true / (num_true + num_false)
        else:
            xs_results.append([[r[0] for r in results if r[1] is True], True])
            xs_results.append([[r[0] for r in results if r[1] is False], False])

        for xs, status in xs_results:
            plt.rc('xtick', labelsize=16)
            plt.rc('ytick', labelsize=16)
            fig, ax = plt.subplots(figsize=(20, 16))

            hist_name_file = ['bc=' + str(bc)
                , 'K=' + str(K)
                , 't=' + str(t)
                , 'n=' + str(n)
                , 'cycles=' + str(cycles_list)
                , 'hamming2include=' + str(hamming_dist_to_include)]

            hist_description_file = copy.deepcopy(hist_name_file)

            # Expectation calculation
            xs_flattened_list = [num for sublist in xs for num in sublist]
            expectation_value = sum(xs_flattened_list) / len(xs) if xs else 0
            hamming_dist_bigger_then_0_sum = sum(hamming_dist_bigger_then_0_list)
            hist_description_file.extend(['E=' + str(round(expectation_value, 2)),
                                          'skipped0=' + str(sum(skipped_zero_count_list)),
                                          'status=' + str(status),
                                          'hammingBiggerThen0=' + str(hamming_dist_bigger_then_0_sum),
                                          'T/(T+F)=' + str(true_false_ratio)])

            # Plot histogram
            ax.hist(xs_flattened_list, bins=20, alpha=0.5, label='Number of samples')
            ax.set_xlabel('# Reads', fontsize=20)
            ax.set_ylabel('Frequency', fontsize=20)
            ax.legend(loc='upper right', fontsize=18)
            ax.set_title("Sampling Rate", fontsize=22)

            plt.subplots_adjust(right=0.8)  # Adjust the bottom leaving space for your description
            fig.text(0.82, 0.5, "\n".join(hist_description_file), ha='left', va='center',
                     fontsize=16, bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=2))

            fig.savefig(self.output_hist_coupon_collector_folder + "".join(hist_name_file) + str('.svg'))
            plt.close(fig)
            print(hist_name_file)

        file_name = self.output_csv_coupon_collector_folder + "_".join(hist_name_file) + '.csv'

        headers = ["Total Sample", "Status", "Design", "Sample", "Cycle", "Unique values", "Hamming Dist bigger then 0",
                   "Skipped zero count"]

        # Write to CSV
        with open(file_name, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            for csv_row in results:
                writer.writerow(csv_row)

        print(f'Data written to {file_name}')

    def create_heatmap_with_rectangles_on_most_common(self, dict_foreach_bc_and_x_count_all_cycles_matrix: Dict,
                                                      heatmap_foreach_bc_and_x_count_with_most_common_file: Union[
                                                          Path, str], ax) -> None:

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

        plt.savefig(heatmap_foreach_bc_and_x_count_with_most_common_file, dpi=400)
        plt.close()

    def heatmap_foreach_bc_and_x_count_with_most_common(self,
                                                        heatmap_foreach_bc_and_x_count_with_most_common_file: Union[
                                                            Path, str]) -> None:
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
                                                           dict_foreach_bc_and_x_count_all_cycles_matrix,
                                                           heatmap_foreach_bc_and_x_count_with_most_common_file=heatmap_foreach_bc_and_x_count_with_most_common_file,
                                                           ax=ax)

    def find_most_common(self, input_csv_path: Union[Path, str],
                         foreach_bc_payload_count_file_dict_to_csv: Union[Path, str],
                         most_common_dict_to_csv_path: Union[Path, str],
                         compare_design_to_experiment_results_output_file: Union[Path, str],
                         design_simulation_file: Union[Path, str],
                         heatmap_foreach_bc_and_x_count_with_most_common_file: Union[Path, str]) -> None:
        dict_bc = self.analyze_results_good_reads(input_csv_path=input_csv_path,
                                                  dict_to_csv=foreach_bc_payload_count_file_dict_to_csv)
        self.most_common_for_each_bc(dict_bc=dict_bc, dict_to_csv_path=most_common_dict_to_csv_path)
        result_payload, result_only_z = self.convert_most_common_to_letters_in_new_alphabet(
            results_most_common_file=most_common_dict_to_csv_path)
        self.compare_most_common_to_design(result_payload=result_payload,
                                           compare_design_to_experiment_results_output_file=compare_design_to_experiment_results_output_file,
                                           design_simulation_file=design_simulation_file)

        self.heatmap_foreach_bc_and_x_count_with_most_common(
            heatmap_foreach_bc_and_x_count_with_most_common_file=heatmap_foreach_bc_and_x_count_with_most_common_file)

    def missing_bc_to_csv(self, dict_append_missing_bc):
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        ser_append_missing_bc.to_csv(self.missing_bcs_file, mode='a', header=False)

    def create_sampling_rate_from_good_reads_graph(self, input_csv_path: Path) -> None:
        """
        This function takes in a list of integers and creates a graph with the x-axis being the sampling rate
        and the y-axis being the count of different values in arr[0].
        The sampling rate will be in increments of 10% starting from 0% until 100%.
        """

        df_good_reads = pd.read_csv(input_csv_path)

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

    def for_each_bc_count_reads_read_to_csv(self, output_file: Union[Path, str], input_csv_path: Path) -> pd.DataFrame:
        with open(output_file, "ab") as f:
            cols_names = [['bc', 'count']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        df = pd.read_csv(input_csv_path)
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

    def move_output_to_sampling_percentage_directory(self, new_output_dir: str):
        src_dir = self.output_folder
        dst_dir = new_output_dir
        # Create the destination directory if it does not exist
        if not os.path.exists(dst_dir):
            try:
                os.makedirs(dst_dir)
            except PermissionError:
                print(f"Unable to create directory {dst_dir}: Permission denied")
                return

        # Change the permissions of the destination directory
        try:
            os.chmod(dst_dir, 0o777)
        except PermissionError:
            print(f"Unable to change permissions for directory {dst_dir}: Permission denied")
            return

        # Move the output files to the destination directory
        for item in os.listdir(src_dir):
            src_path = os.path.join(src_dir, item)
            dst_path = os.path.join(dst_dir, item)
            dst_path = dst_dir
            if os.path.isfile(src_path):
                try:
                    # Check if a file with the same name already exists in the destination folder
                    if os.path.exists(dst_path):
                        # Rename the source file
                        dst_path = os.path.join(dst_dir,
                                                f"{os.path.splitext(item)[0]}_{len(os.listdir(dst_dir))}{os.path.splitext(item)[1]}")
                    shutil.copy2(src_path, dst_path)
                except PermissionError:
                    print(f"Unable to copy file {src_path} to {dst_path}: Permission denied")
            elif os.path.isdir(src_path):
                try:
                    # Check if the destination path already exists and is a directory
                    if os.path.exists(dst_path) and os.path.isdir(dst_path):
                        dst_path = os.path.join(dst_path, os.path.basename(src_path))
                    os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                    shutil.copytree(src_path, dst_path)
                except PermissionError:
                    print(f"Unable to copy directory {src_path} to {dst_path}: Permission denied")
            else:
                print(f"Skipping {src_path} as it is not a file or directory")

        # Change the permissions of the destination directory back to its original value
        try:
            os.chmod(dst_dir, 0o755)
        except PermissionError:
            print(f"Unable to change permissions for directory {dst_dir}: Permission denied")
    def graph_of_all_sampling_most_common_x_to_each_bc(self):
        all_results_count = {}
        x_axis = []
        df_design_results_only_z_file = pd.read_csv(self.design_results_only_z_file, index_col=0)

        for folder_name in os.listdir(self.output_all_sampling_rate_folder):
            # folder_path = os.path.join(self.output_all_sampling_rate_folder, folder_name, "csv")
            if float(folder_name.split('_')[1]) not in [0.01, 0.005, 0.001]:
                continue
            sampling_percentage_folder = os.path.join(self.output_all_sampling_rate_folder, folder_name)
            if os.path.isdir(sampling_percentage_folder):
                sampling_percentage = float(folder_name.split("_")[-1])
                all_results_count[sampling_percentage] = {}
                if sampling_percentage not in self.sampling_rate_array:
                    continue

                for iteration_folder in os.listdir(sampling_percentage_folder):
                    # sampling_percentage_folder = "output_all_sampling_rate/output_100"
                    if iteration_folder.startswith("iter_"):
                        all_results_count[sampling_percentage][iteration_folder] = []
                        iter_path = os.path.join(sampling_percentage_folder, iteration_folder)
                        csv_file = os.path.join(iter_path, "csv/foreach_bc_payload_count.csv")

                        if os.path.exists(csv_file):
                            data = pd.read_csv(csv_file, header=None)

                # csv_file = os.path.join(folder_path, "foreach_bc_payload_count.csv")

                            for bc_i, item in data.items():
                                x_count_to_each_bc = []
                                if bc_i == 0 or bc_i not in self.bc_list:
                                    continue
                                # Parse the dictionary from the cell
                                dict_data = ast.literal_eval(item[1])

                                for cycle_i, (cycle_name, payload_count_dict) in enumerate(dict_data.items()):
                                    # Create a Counter from the dictionary
                                    counter = collections.Counter(payload_count_dict)
                                    # Remove if the payload is 0
                                    del counter[0]
                                    # Get the 4 most common items
                                    most_common_payloads = counter.most_common(self.subset_size)
                                    # Filter the payloads that have 0 appearances.
                                    filtered_common = [item for item, count in most_common_payloads if count > 0]
                                    x_count_to_each_bc.append(self.compare_amount_of_x_to_design(set_most_common_x_to_each_bc=set(filtered_common),
                                                                       z=df_design_results_only_z_file[str(bc_i)][cycle_i]))

                                all_results_count[sampling_percentage][iteration_folder].append(x_count_to_each_bc)

        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='errorbar', x_axis_type='sampling_rate')
        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='boxplot', x_axis_type='sampling_rate')
        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='errorbar', x_axis_type='bc_avg_coverage')
        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='boxplot', x_axis_type='bc_avg_coverage')


    def compare_amount_of_x_to_design(self, set_most_common_x_to_each_bc, z):
        z_to_x_tuple = self.z_to_k_mer_representative[z]
        set_of_x_vals = set([int(x[1:]) for x in z_to_x_tuple])

        # Get common elements
        common_elements = set_most_common_x_to_each_bc.intersection(set_of_x_vals)

        # Get count of common elements
        count = len(common_elements)

        return count

    def calculate_average_coverage(self, csv_file: Union[Path, str]):
            bc_count_df = pd.read_csv(csv_file)

            # filter the DataFrame
            filtered_df = bc_count_df[bc_count_df['bc'].isin(self.bc_list)]

            # calculate the average of the 'count' column
            avg_count = np.mean(filtered_df['count'])

            return avg_count

    def create_graph_avg_x_to_each_sampling(self, all_counts, graph_type, x_axis_type):
        # x_axis = []
        data = {}
        data_avg = {}

        for sampling_percentage, iterations_data in all_counts.items():
            data[sampling_percentage] = []

            # Determine the number of elements in each iteration (assuming equal for all iterations)
            num_elements = len(
                iterations_data[next(iter(iterations_data))])  # Get the length of one of the iteration lists

            # Iterate over the number of elements
            for i in range(num_elements):
                temp_data = []

                # Iterate over each iteration
                for iteration in iterations_data:
                    # Append the i-th element of each iteration
                    # Check if the i-th element exists in this iteration
                    if i < len(iterations_data[iteration]):
                        temp_data.extend(iterations_data[iteration][i])
                data[sampling_percentage].append(temp_data)  # Extract the value from the list

        for folder_name in os.listdir(self.output_all_sampling_rate_folder):
            if float(folder_name.split('_')[1]) not in [0.01, 0.005, 0.001]:
                continue
            x_axis = []
            avg_count_per_iter = []
            folder_path = os.path.join(self.output_all_sampling_rate_folder, folder_name)
            sampling_percentage = float(folder_name.split("_")[-1])
            # Add the sampling percentage and counts to the lists
            if x_axis_type == 'sampling_rate':
                x_axis.extend([sampling_percentage] * len(all_counts[sampling_percentage]))
            elif x_axis_type == 'bc_avg_coverage':
                    for iter_i, iterations in all_counts[sampling_percentage].items():
                        csv_file_name = self.count_reads_for_each_bc_file.replace("analyze_sequencing_data/output/", "")
                        csv_file = os.path.join(folder_path,iter_i ,csv_file_name)
                        try:
                            avg_count = round(self.calculate_average_coverage(csv_file=csv_file))
                        except:
                            print(csv_file)
                        avg_count_per_iter.append([avg_count])
                    x_axis_avg = round(sum([item[0] for item in avg_count_per_iter]) / len(avg_count_per_iter))
                    x_axis.extend([x_axis_avg] * self.amount_of_bc)

                    # Group data points by sampling percentage
                    for bc in range(self.amount_of_bc):
                        data_avg_temp = []
                        for iter_i, iterations in all_counts[sampling_percentage].items():

                            if x_axis_avg not in data_avg:
                                data_avg[x_axis_avg] = []
                            data_avg_temp.extend(all_counts[sampling_percentage][iter_i][bc])
                        data_avg[x_axis_avg].append(data_avg_temp)

        if x_axis_type == 'bc_avg_coverage':
            data = data_avg
        # Calculate mean and standard deviation for each group
        means = []
        std_devs = []
        sorted_keys = sorted(data.keys())
        for key in sorted_keys:
            mean = np.mean(data[key])
            std_dev = np.std(data[key])
            means.append(mean)
            std_devs.append(std_dev)

        # Create a saturation analysis plot
        if graph_type == 'errorbar':
            plt.errorbar(sorted_keys, means, yerr=std_devs, capsize=5, fmt='o-', markersize=8)
            plt.grid(True)
        elif graph_type == 'boxplot':
            # plt.boxplot(sorted_keys, labels=sorted_keys)
            flattened_data = [[item for sublist in lst for item in sublist] for lst in
                              [data[key] for key in sorted_keys]]
            plt.boxplot(flattened_data)
            sorted_keys_xticks_str_arr = [str(round(element,2)) for element in sorted_keys]
            positions = range(1, len(sorted_keys) + 1)
            plt.xticks(positions, sorted_keys_xticks_str_arr, rotation=45)


        # plt.ylim(-0.5, 4.5)
        plt.ylabel('Number Of Observed $s_i$')
        # plt.title('Saturation Analysis of Design X and Experiment Results X')
        plt.subplots_adjust(bottom=0.2)
        if x_axis_type == 'sampling_rate':
            plt.xlabel('Sampling Percentage')
            plt.savefig(self.output_line_graphs_folder + 'sampling_rate_' + graph_type+'.svg')
        elif x_axis_type == 'bc_avg_coverage':
            plt.xlabel('Average Reads Per Strand')
            plt.savefig(self.output_line_graphs_folder + 'bc_avg_coverage_' + graph_type+'.svg')
        plt.show()
        plt.close()

    def save_output_into_separate_folder(self, sampling_percentage: float, iteration: int):
        output_folder = f"analyze_sequencing_data/{self.output_all_sampling_rate_folder}output_{sampling_percentage}/iter_{iteration}"
        os.makedirs(output_folder, exist_ok=True)
        self.move_output_to_sampling_percentage_directory(new_output_dir=output_folder)
        return output_folder

    def for_each_bc_count_reads_read(self, csv_output_file: Union[Path, str], input_csv_path: Path) -> None:
        self.for_each_bc_count_reads_read_to_csv(output_file=csv_output_file, input_csv_path=input_csv_path)
        self.hist_foreach_bc_read_count(csv_output_file)
        self.hist_foreach_read_count_count_bc(csv_output_file=csv_output_file)
        self.hist_foreach_error_count_of_bc()

    def create_new_input_csv_according_to_sampling_rate(self, sampling_rate: float,
                                                        input_csv_path: str,
                                                        iteration: int) -> str:
        # Set the chunk size and sampling fraction
        chunk_size = 100000  # 100,000 rows per chunk
        sampling_fraction = float(sampling_rate / 100)

        # Create a directory to store the sampled CSV files
        output_sampling_rate_dir = 'analyze_sequencing_data/'+self.output_all_sampling_rate_folder + f'output_{sampling_rate}/iter_{iteration}/'
        os.makedirs(output_sampling_rate_dir, exist_ok=True)

        # Iterate over the CSV file in chunks
        sampled_filepaths = []
        for i, chunk in enumerate(pd.read_csv(input_csv_path, chunksize=chunk_size)):
            # Randomly sample rows from the chunk with the specified fraction
            sampled_chunk = chunk.sample(frac=sampling_fraction)
            # Write the sampled chunk to a new CSV file
            chunk_filename = f'sampled_chunk_{i}_iter_{iteration}.csv'
            chunk_filepath = os.path.join(output_sampling_rate_dir, chunk_filename)
            sampled_chunk.to_csv(chunk_filepath, index=False)
            sampled_filepaths.append(chunk_filepath)

        # Concatenate the sampled files into a single output file
        output_filepath = os.path.join(output_sampling_rate_dir,
                                       f'sampled_{sampling_rate}_iter_{iteration}.csv')
        write_header = True
        with open(output_filepath, 'w+') as output_file:
            for filepath in sampled_filepaths:
                with open(filepath, 'r') as input_file:
                    if write_header:
                        output_file.write(input_file.readline())  # Write the header
                        write_header = False
                    else:
                        input_file.readline()  # Skip the header

                    # Copy the rest of the file content
                    shutil.copyfileobj(input_file, output_file)
                os.remove(filepath)

        # Change the permissions of the output file to 0o777
        os.chmod(output_filepath, 0o777)

        return output_filepath

    def create_folders(self) -> None:
        utilities.is_dir_exists(self.output_hist_folder)
        utilities.is_dir_exists(self.output_folder)
        utilities.is_dir_exists(self.output_heatmap_folder)
        utilities.is_dir_exists(self.output_graphs_folder)
        utilities.is_dir_exists(self.output_csv_folder)
        utilities.is_dir_exists(self.output_line_graphs_folder)
        utilities.is_dir_exists(self.output_csv_coupon_collector_folder)
        utilities.is_dir_exists(self.output_hist_coupon_collector_folder)

    def delete_csv_files_from_folder(self):
        for item in os.listdir(self.output_csv_folder):
            if item.endswith(".csv"):
                file_path = os.path.join(self.output_csv_folder, item)
                try:
                    os.remove(file_path)
                    print(f"Deleted {file_path}")
                except OSError as e:
                    print(f"Error deleting {file_path}: {e}")

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

        input_csv_path = self.results_good_reads_file
        for sampling_percentage in self.sampling_rate_array:
            # if sampling_percentage == 100:
            #     input_csv_path_sampling_percentage_100 = 'analyze_sequencing_data/output_all_sampling_rate/output_100/iter_0/csv/results_good_reads.csv'
            #     continue
            for iteration in range(10):
                if sampling_percentage != 100:
                    input_csv_path = self.create_new_input_csv_according_to_sampling_rate(
                        sampling_rate=sampling_percentage,
                        input_csv_path=input_csv_path_sampling_percentage_100,
                        iteration=iteration)
                # Find most common for each bc and for every cycle in that bc in results of good reads
                self.find_most_common(input_csv_path=input_csv_path,
                                      foreach_bc_payload_count_file_dict_to_csv=self.foreach_bc_payload_count_file,
                                      most_common_dict_to_csv_path=self.results_most_common_file,
                                      compare_design_to_experiment_results_output_file=self.compare_design_to_experiment_results_output_file,
                                      design_simulation_file=self.design_simulation_file,
                                      heatmap_foreach_bc_and_x_count_with_most_common_file=self.heatmap_foreach_bc_and_x_count_with_most_common_file)

                # For each bc count amount of reads sequenced
                self.for_each_bc_count_reads_read(csv_output_file=self.count_reads_for_each_bc_file, input_csv_path=input_csv_path)

                # Create graph with sampling rate
                self.create_sampling_rate_from_good_reads_graph(input_csv_path=input_csv_path)

                # Compare Payloads_of_two_bc
                self.compare_payloads_of_two_bc(bc1=1, bc2=2)

                # Create a summary of the results for BC1 and BC2 that will show the fraction of each k-mer in each cycle to each BC
                self.create_results_table_for_two_bc(bc1=1, bc2=2)

                # We want to find
                self.the_coupon_collector_problem_with_t()

                # Move the data into the designated folder for each sampling rate
                output_folder = self.save_output_into_separate_folder(sampling_percentage=sampling_percentage,
                                                                      iteration=iteration)
                if sampling_percentage == 100:
                    input_csv_path_sampling_percentage_100 = output_folder + '/csv/results_good_reads.csv'

            # Delete all csvs from the output folder
            self.delete_csv_files_from_folder()

        self.graph_of_all_sampling_most_common_x_to_each_bc()
