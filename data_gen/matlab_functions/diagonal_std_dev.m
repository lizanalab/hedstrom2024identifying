function diagonal_devs = diagonal_std_dev(matrix, min_elements_diag)

    [num_rows, num_columns] = size(matrix);
    diagonal_devs = zeros(1, num_columns);

    for i = 0:num_columns
        diagonal_elements = diag(matrix, i);
        if (num_rows - i + 1) < min_elements_diag || sum(diagonal_elements > 0) == 0
            break;
        end
        diagonal_devs(i+1) = std(diagonal_elements(~isinf(diagonal_elements)));
    end
end
