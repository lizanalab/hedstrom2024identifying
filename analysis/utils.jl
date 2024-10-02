using Random
using Distributions

# Function to calculate the standard deviation for the diagonals of a given matrix.
# Arguments:
#   - matrix: a 2D array (matrix) from which the diagonal elements will be extracted.
# Returns:
#   - diagonal_devs: a vector containing the standard deviations of each diagonal of the matrix.
function diagonal_std_dev(matrix)
    _, num_columns = size(matrix)  # Get the number of columns in the matrix.
    diagonal_devs = zeros(Float64, 1, num_columns)  # Initialize the diagonal deviations array.

    # Loop through each diagonal of the matrix.
    for i in 1:num_columns
        diagonal_elements = diag(matrix, i)  # Extract the diagonal elements.
        
        # If all diagonal elements are non-positive, break out of the loop.
        if all(diagonal_elements .<= 0)
            break
        end
        
        # Filter out infinite values and compute the standard deviation.
        diagonal_devs[i] = std(filter(!isinf, diagonal_elements))
    end

    return diagonal_devs
end


# Function to perform resampling of pixel values based on a normal distribution centered on the original matrix values.
# Arguments:
#   - matrix: a 2D array (matrix) where each value represents a pixel value.
#   - diagonal_devs: a vector of standard deviations for each diagonal, calculated from the `diagonal_std_dev` function.
# Returns:
#   - resamp_matrix: a new resampled matrix where pixel values are randomly generated based on a normal distribution.
function pixel_distribution_resampling(matrix, diagonal_devs)
    num_rows, num_columns = size(matrix)  # Get the number of rows and columns in the matrix.
    resamp_matrix = zeros(Float64, num_rows, num_columns)  # Initialize the resampled matrix.

    # Loop through each element in the upper triangular portion of the matrix.
    for i in 1:num_rows
        for j in i:num_columns
            diagonal_index = j - i + 1  # Determine the corresponding diagonal index.

            # If the matrix element is infinite, set a random pixel value of 0.
            if isinf(matrix[j, i])
                rand_pixel = 0
            else
                # Sample a random pixel value from a normal distribution centered around the matrix element.
                rand_pixel = exp(rand(Normal(matrix[j, i], diagonal_devs[diagonal_index])))
                rand_pixel = round(rand_pixel)
                
                # Ensure that the pixel value is not 0; if it is, set it to 1.
                if rand_pixel == 0
                    rand_pixel = 1
                end
            end

            # Set the resampled values for both the current element and its symmetric counterpart.
            resamp_matrix[i, j] = round(rand_pixel)
            resamp_matrix[j, i] = round(rand_pixel)
        end
    end

    return resamp_matrix
end
