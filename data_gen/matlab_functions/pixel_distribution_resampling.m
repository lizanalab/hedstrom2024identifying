function resamp_matrix = pixel_distribution_resampling(matrix, diagonal_devs, min_elements_diag)

    [num_rows, num_columns] = size(matrix);
    resamp_matrix = zeros(num_rows, num_columns);
    
    for i = 1:num_rows
        for j = i:num_columns
            diagonal_index = j - i;
            if diagonal_index > num_rows - min_elements_diag
                break;
            end
           
            if isinf(matrix(j,i))
                rand_pixel = 0;
            else
                rand_pixel = exp(normrnd(matrix(j, i), diagonal_devs(diagonal_index+1), 1, 1)); % normrnd takes std, not var
                rand_pixel = round(rand_pixel); % to keep the values as discrete counts
                if rand_pixel == 0
                    rand_pixel = 1; % to avoid zero cols/rows
                end
            end
            
            resamp_matrix(i, j) = rand_pixel;
            resamp_matrix(j, i) = rand_pixel;
        end
    end
end
