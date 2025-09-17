%% Evaluation of Background Discounting hypothesis

% This is the code for generating Fig. 5 in the main paper.
% If you use this code, please consider citing the following paper:
% D. Kim, M. Ashraf, A. Chapiro, and R. K. Mantiuk, 
% 'Supra-threshold Contrast Perception in Augmented Reality'
% Conference Proceeding, Siggraph Asia 2025

clear all; close all;

Label_FontSize = 22;
Legend_FontSize = 16;
FontName = 'Linux Biolinum';

current_dir = pwd;
[parent_dir, ~, ~] = fileparts(current_dir);

bkg_labels = {'single-plane', 'dual-flat', 'dual-noise-static', 'dual-noise-dynamic', 'dual-noise-pinhole'};
bkg_index_labels = {'SP','DF','DNS','DND','DNP'};
plot_bkg_idx = [2:5];
% 'lum_total' average regardless of lum  
% 'lum_dark': 30 vs 15 
% 'lum_bright': 30 vs 60
mode = 'lum_total'; 
plot_ind = true;
ind_calib = false;
add_stat_sig = false;
ref_calib = true;
error_mode = 'ci';
alpha = 0.95;

average_csv_file = fullfile('data/exp1', 'exp1_average_data.csv');
D_average = readtable(average_csv_file);

ind_csv_file = fullfile('data/exp1', 'exp1_individual_data.csv');
D_ind = readtable(ind_csv_file);

%% Response scaling method
sub_folder = 'exp1';
savepath = fullfile('result/', sub_folder);
mkdir(savepath);

xlim_val = [min(plot_bkg_idx)-0.5, max(plot_bkg_idx)+0.5];
phys_ylim = [0.15 0.6];
df_ylim = [0.2 2.5];
fig_pos = [100, 100, 500, 600];
ytick_val = [0.2 0.3 0.5];
xlabel_string = 'Condition';
ylabel_string = 'Matched test contrast';

bkg_types_vec = unique(D_average.bkg_type);
c_phys_vec = unique(D_average.c_phys);
sf_vec = unique(D_average.cpd);
obs_vec = unique(D_ind.obs);

num_bkg_type = numel(bkg_types_vec);
num_c_phys = numel(c_phys_vec);
num_sf = numel(sf_vec);
num_obs = numel(obs_vec);
COLORs_bkg = lines(numel(unique(D_average.bkg_type)));

total_mean_log_c_phys = zeros(num_c_phys, num_bkg_type);
total_se_log_c_phys = zeros(num_c_phys, num_bkg_type);

% Desired order
desired_order = {
    'VR-flat-normal-bkg-same-bkg-contrast-1-opposite-static';
    'flat-normal-bkg-same-bkg-contrast-1-opposite-static';
    'optics-normal-bkg-same-bkg-contrast-1-opposite-static';
    'optics-normal-bkg-same-bkg-contrast-1-opposite-slow-dynamic';
    'noise-normal-bkg-same-bkg-contrast-1-opposite-static'};

% Create a lookup table for the desired order
[~, bkg_sort_idx] = ismember(bkg_types_vec, desired_order);

% Define the tiled layout
merged_fig = figure(); % Create a new figure
t = tiledlayout(1, num_c_phys, 'TileSpacing', 'compact', 'Padding', 'compact'); % One column, num_c_phys rows
merged_fig.Position = [100 100 1000 600];
% Iterate through each condition to plot
for c = 1:num_c_phys
    nexttile; % Move to the next tile in the layout

    % Extract contrast value
    c_phys = c_phys_vec(c);
    
    % Prepare the data for the current condition
    for b = plot_bkg_idx
        b_idx = bkg_sort_idx(b);
        idx = find((D_average.c_phys == c_phys) .* (strcmp(D_average.bkg_type, string(bkg_types_vec{b_idx}))));

        if ref_calib
            ref_idx = find((D_average.c_phys == c_phys) .* (strcmp(D_average.bkg_type, string(bkg_types_vec{1}))));
            ref_logmean_c_phys = 10.^D_average.mean_log_c_phys_meas(ref_idx);
        else
            ref_logmean_c_phys = c_phys;
        end

        logmean_c_phys = 10.^D_average.mean_log_c_phys_meas(idx) -  (ref_logmean_c_phys-c_phys);
        se_log_c_phys = 10.^(D_average.mean_log_c_phys_meas(idx)) .* ...
                        log(10) .* D_average.std_log_c_phys_meas(idx) ./ sqrt(D_average.N(idx));

        total_mean_log_c_phys(c,b) = logmean_c_phys;
        total_se_log_c_phys(c,b) = se_log_c_phys;

        if strcmp(error_mode, 'se')
            error_values = 10.^(D_average.mean_log_c_phys_meas(idx)) .* ...
                log(10) .* D_average.std_log_c_phys_meas(idx) ./ sqrt(D_average.N(idx));
        elseif strcmp(error_mode, 'ci')
            c_lin = 10.^D_average.mean_log_c_phys_meas(idx);
            ci_lower = c_lin - (c_lin .* 10.^(-D_average.std_log_c_phys_meas(idx) * 1.96./sqrt(D_average.N(idx))));
            ci_upper = (c_lin .* 10.^( D_average.std_log_c_phys_meas(idx)* 1.96./sqrt(D_average.N(idx)) )) - c_lin;
            error_values = [ci_lower ci_upper];
        else
            error('Invalid error_mode selected. Choose ''se'' or ''ci''.');
        end

        bar(b, logmean_c_phys, 'FaceColor', COLORs_bkg(b, :), 'BarWidth', 0.6, 'DisplayName', [bkg_labels{b} ' (' bkg_index_labels{b} ')']);
        hold on;
        % Add error bars
        if strcmp(error_mode, 'se')
            errorbar(b, total_mean_log_c_phys(c, b), error_values, ...
                'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
        elseif strcmp(error_mode, 'ci')
            errorbar(b, total_mean_log_c_phys(c, b), error_values(1), error_values(2), ...
                'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');
        end

        [string(c_phys) string(bkg_types_vec{b_idx}) logmean_c_phys]
    end

    % Plot individual data points
    if plot_ind
        for i = 1:numel(unique(D_ind.obs))
            obs_name = obs_vec(i);
            for b = plot_bkg_idx
                b_idx = bkg_sort_idx(b);
                idx = find((D_ind.c_phys == c_phys) .* strcmp(D_ind.obs, string(obs_name)) .* ...
                           strcmp(D_ind.bkg_type, string(bkg_types_vec(b_idx))));
                N_dp = nnz(idx);
                
                if ref_calib
                    ref_idx = find((D_ind.c_phys == c_phys) .* (strcmp(D_ind.bkg_type, string(bkg_types_vec{1}))));
                    ref_logmean_c_phys = mean(10.^D_ind.mean_log_c_phys_meas(ref_idx));
                else
                    ref_logmean_c_phys = c_phys;
                end

                if ind_calib
                    yy = 10.^(D_ind.median_log_c_phys_meas(idx))- (ref_logmean_c_phys-c_phys);
                else
                    yy = 10.^(D_ind.median_log_c_phys_meas(idx))- (ref_logmean_c_phys-c_phys);
                end
                scatter(b+(rand(N_dp,1)-0.5)*0.2,yy, ...
                    50, 'filled', 'MarkerFaceColor', 'k', ...
                    'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
            end
        end
    end

    if c==1
        legend_handle = legend('show', 'Location', 'northwest');
        set(legend_handle, 'FontName', FontName, 'FontSize', Legend_FontSize, 'Interpreter', 'tex');
    else
        text(min(plot_bkg_idx), 0.5 * 1.07, ...
            '$c_{\mathrm{disp,max}}=\frac{Y_{\mathrm{FG}}}{Y_{\mathrm{FG}} + Y_{\mathrm{BG}}}=0.5$', ...
            'FontSize', Label_FontSize, 'Color', 'r', 'HorizontalAlignment', 'left', ...
            'Interpreter', 'latex');
    end

    if add_stat_sig
        %        Find the index for "DNP"
        dnp_idx = find(strcmp(bkg_index_labels, 'DNP'));
        dnp_y = total_mean_log_c_phys(c, dnp_idx);

        % Add an asterisk above the "DNP" bar
        text(dnp_idx, dnp_y-0.01, '***', 'FontSize', Label_FontSize + 5, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontName', FontName, 'Color', COLORs_bkg(dnp_idx,:), 'Interpreter', 'tex');
        % if c==1
        % 
        % end
    end

    % Set plot properties
    set(gca, 'YScale', 'log','FontSize',Legend_FontSize,'FontName',FontName); % Log scale for y-axis
    xlabel(xlabel_string, 'FontSize', Label_FontSize, 'FontName', FontName,'Interpreter', 'tex');
    if c==1
        ylabel(ylabel_string, 'FontSize', Label_FontSize, 'FontName', FontName, 'Interpreter', 'tex');
    else
        ylabel([]);
    end
    grid on;
    xlim(xlim_val);
    ylim(phys_ylim);
    yticks(ytick_val);
    yline(0.5, '--', 'HandleVisibility', 'off', 'Color', 'r', 'LineWidth', 1.5);
    yline(c_phys, '--', 'HandleVisibility', 'off', 'Color', 'k', 'LineWidth', 1.5);
    xticks(1:num_bkg_type);
    xticklabels(bkg_index_labels);
    title(sprintf('Reference contrast: %.1f (N=%d)', c_phys, num_obs), 'FontSize', Label_FontSize, 'FontName', FontName, 'Interpreter', 'tex');
end

if ind_calib
    output_pdf = fullfile(savepath, 'exp1_result_ind_calib.pdf');
    output_png = fullfile(savepath, 'exp1_result_ind_calib.png');
else
    output_pdf = fullfile(savepath, 'exp1_result.pdf');
    output_png = fullfile(savepath, 'exp1_result.png');
end
exportgraphics(t, output_pdf, 'ContentType', 'vector');
exportgraphics(t, output_png);
disp(['Merged PDF/PNG saved']);