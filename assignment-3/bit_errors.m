function num_of_bit_errors = bit_errors(est_bit_seq,bit_seq)

num_of_bit_errors = sum(est_bit_seq(:) ~= bit_seq(:));

end
