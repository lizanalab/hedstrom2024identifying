function matrix = get_matrix_from_raw_data(data_path, resolution)

    % Read the data from the file
    data = importdata(data_path);

    % Get the indices and values
    indices = data(:, 1:2)./resolution + 1;  % First two columns are indices
    values = data(:, 3);     % Third column is the value

    % Determine the size of the matrix
    num_rows = max(indices(:, 1));
    num_cols = max(indices(:, 2));
    
    
    % Initialize the matrix with zeros
    matrix = zeros(num_rows, num_cols);

    % Populate the matrix with the values
    for i = 1:size(indices, 1)
        row = indices(i, 1);
        col = indices(i, 2);
        value = values(i);
        matrix(row, col) = value;
        matrix(col,row) = value;
    end

end