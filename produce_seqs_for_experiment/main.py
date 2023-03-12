from produce_seqs_for_experiment.produce_seqs import ProduceSeqs
from config import get_configs

if __name__ == '__main__':
    configs = get_configs()
    for config in configs:
        produce_seqs = ProduceSeqs(number_of_cycles=config['number_of_cycles'],
                                   data_in_kb=config['data_in_kb'],
                                   bits_per_synthesis_cycle=config['bits_per_synthesis_cycle'],
                                   sequence_len=config['sequence_len'],
                                   num_of_building_blocks=config['num_of_building_blocks'],
                                   num_hamming_distance=config['num_hamming_distance'],
                                   forbidden_x_letters_in_a_row=config['forbidden_x_letters_in_a_row'],
                                   each_letter_has_to_be_x_percent_occurrences_from=config[
                                       'each_letter_has_to_be_x_percent_occurrences_from'],
                                   each_letter_has_to_be_x_percent_occurrences_to=config[
                                       'each_letter_has_to_be_x_percent_occurrences_to'],
                                   stick_ends_overhang=config['stick_ends_overhang']
                                   )

        produce_seqs.run()



