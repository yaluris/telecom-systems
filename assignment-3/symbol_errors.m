function num_of_symbol_errors = symbol_errors(est_X,X)

num_of_symbol_errors = sum(est_X(:) ~= X(:));

end
