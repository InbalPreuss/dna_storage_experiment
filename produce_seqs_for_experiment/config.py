config_1 = {
    'number_of_cycles': 1,
    'data_in_kb': 1,
    'bits_per_synthesis_cycle': 12,
    'subset_size': 5,
    'num_of_building_blocks': 16,
    'sequence_len': 20,
    'num_hamming_distance': 10,
    'forbidden_x_letters_in_a_row': 3,
    'each_letter_has_to_be_x_percent_occurrences_from': 0.18,
    'each_letter_has_to_be_x_percent_occurrences_to': 0.32,
    'stick_ends_overhang': 8
}
config_2 = {
    'number_of_cycles': 2,
    'data_in_kb': 1,
    'bits_per_synthesis_cycle': 12,
    'subset_size': 5,
    'num_of_building_blocks': 16,
    'sequence_len': 20,
    'num_hamming_distance': 11,
    'forbidden_x_letters_in_a_row': 3,
    'each_letter_has_to_be_x_percent_occurrences_from': 0.18,
    'each_letter_has_to_be_x_percent_occurrences_to': 0.32,
    'stick_ends_overhang': 8
}


config_3 = {
    'number_of_cycles': 3,
    'data_in_kb': 1,
    'bits_per_synthesis_cycle': 12,
    'subset_size': 5,
    'num_of_building_blocks': 16,
    'sequence_len': 20,
    'num_hamming_distance': 11,
    'forbidden_x_letters_in_a_row': 3,
    'each_letter_has_to_be_x_percent_occurrences_from': 0.18,
    'each_letter_has_to_be_x_percent_occurrences_to': 0.32,
    'stick_ends_overhang': 8
}

config_4 = {
    'number_of_cycles': 4,
    'data_in_kb': 1,
    'bits_per_synthesis_cycle': 12,
    'subset_size': 5,
    'num_of_building_blocks': 16,
    'sequence_len': 20,
    'num_hamming_distance': 12,
    'forbidden_x_letters_in_a_row': 3,
    'each_letter_has_to_be_x_percent_occurrences_from': 0.18,
    'each_letter_has_to_be_x_percent_occurrences_to': 0.32,
    'stick_ends_overhang': 8
}


def get_configs():
    return config_1, config_2, config_3, config_4
