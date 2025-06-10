function generate_communities(chr_num, resolution, gamma, delta, a, num_samp)

    addpath("matlab_functions", "genlouvain");
    base_path = "data/";
    in_path = base_path + "chr" + chr_num + "_" + resolution/1000 + "kb.RAWobserved";

    % when are we ignoring diagonals due to poor statistics
    min_elements_diag = 1;

    % Prepare data
    matrix_raw = get_matrix_from_raw_data(in_path,resolution); % get the cmap data as matrix
    [N, ~] = size(matrix_raw);
    template_matrix = zeros(N);

    [matrix_raw, kept_rows, kept_cols] = delete_centromere(matrix_raw); % get rid of centromere
    matrix_kr = kr_norm(matrix_raw); % kr normalize experimental cmap
    
    % Get community from original map
    template_matrix(kept_rows, kept_cols) = matrix_kr;
    B = heirarchical_domain_model_ana(template_matrix, gamma, a); % background heirarchical domain null model
    commun_exp = iterated_genlouvain(B, 10000, 0, 1, 1); % community search in experimental cmap
    template_matrix(:) = 0;

    n_nodes = numel(commun_exp);

    % Deviations along diag
    log_std_devs = diagonal_std_dev(log(matrix_raw), min_elements_diag); % get the standard deviation per distance 
    log_scaled_std = delta*log_std_devs;

    % alike_communities = cell(num_samp,1);
    communities = cell(num_samp, 1);

    out_filename = "output/chr" + chr_num + "_" + resolution/1000 + "kb_gamma" + gamma + "_delta" + delta + "_a" + a + ".csv";
    fprintf("Running for " + out_filename + "\n");
    for i = 1:num_samp
        % Generate noisy map
        resamp_matrix = pixel_distribution_resampling(log(matrix_raw), log_scaled_std, min_elements_diag); % perform the resampling following the log-normal distrib
        matrix_kr = kr_norm(resamp_matrix); % kr-normalize the newly generated cmap

        % Run GenLouvain
        template_matrix(kept_rows, kept_cols) = matrix_kr;
        B = heirarchical_domain_model_ana(template_matrix, gamma, a); % obtain the background from the Heirarchical domain null model
        commun = iterated_genlouvain(B, 10000, 0, 1, 1, commun_exp); % detect communities using GenLouvain algorithm
        template_matrix(:) = 0;

        communities{i} = commun;
    end

    fid = fopen(out_filename, "w");
    fprintf(fid, "samp\t" + strjoin(arrayfun(@num2str, 1:n_nodes, 'UniformOutput', false), "\t") + "\n");
    fprintf(fid, "0" + "\t" + strjoin(arrayfun(@num2str, commun_exp, 'UniformOutput', false), "\t") + "\n");
    for i = 1:num_samp
        fprintf(fid, i + "\t" + strjoin(arrayfun(@num2str, communities{i}, 'UniformOutput', false), "\t") + "\n");
    end
    fclose(fid);
end
