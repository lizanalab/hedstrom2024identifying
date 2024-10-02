function [matrix_without_centromere, kept_rows, kept_cols]  = delete_centromere(matrix)

    row_sums = sum(matrix, 2);
    col_sums = sum(matrix, 1);

    kept_rows = row_sums ~= 0;
    kept_cols = col_sums ~= 0;
    
    matrix_without_centromere = matrix(kept_rows, kept_cols);
end