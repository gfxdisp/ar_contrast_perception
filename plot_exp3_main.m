%% Generalization to complex images in AR

% This is the code for generating Fig. 8 in the main paper.
% If you use this code, please consider citing the following paper:
% D. Kim, M. Ashraf, A. Chapiro, and R. K. Mantiuk, 
% 'Supra-threshold Contrast Perception in Augmented Reality'
% Conference Proceeding, Siggraph Asia 2025

clear all; close all;

savepath = './result/exp3/';
mkdir(savepath)

Font = 'Linux Biolinum';
FontSize= 16;

binocular = false;
if binocular ==true
    scale = 1;
    exp_mode = 'Stereoscopic';
else
    scale = sqrt(2);
    exp_mode = 'Haploscopic';
end

Y_bg = logspace( -2, 2, 4*2+1 );
Y_bg_samples = 4*10+1;

Y_bg = reshape( Y_bg, [1 1 length(Y_bg)] );
Y_dmax = logspace( -1, log10(5000), 1024 ); % AR display peak luminance levels to consider
x = (0:-1:-10)'; % log-2 relative luminance levels of the display

if ~exist( 'CSF_castleCSF', 'class' )
    addpath( fullfile( pwd,  '..', '..', 'castleCSF', 'matlab' ) );
end

csf = CSF_castleCSF();
c_src = logspace( -2, log10(0.9), 32 );  % contrast levels to consider (for Kulikowski's model)
c_src = reshape( c_src, [1 1 1 length(c_src)] );

s_freq_vec = logspace(log10(1), log10(30), 31);
E = readtable( './data/exp3/exp3_raw_data.csv' );

E.log_Y_fg = log2(E.Y_fg);
% Mean by observer first
Eg = grpstats( E, { 'observer','Y_bg', 'Y_fg_ref' }, { 'mean'}, 'DataVars', 'log_Y_fg' );
% and get the ci based on the mean value per observer
Es = grpstats( Eg, {'Y_bg','Y_fg_ref'}, {'mean', 'meanci'},'DataVars','mean_log_Y_fg');

app_cont_model_vec = {'Kulikowski'};
sf_mode_vec = {'natural'};
loss_type_vec = {'log-l1'};
epsilon = 1e-12;

result_save =true;
ab_save = true;

Y_rr = [200 800]; %200:100:1000;
COLORs = lines(length(Y_rr));
img_res = [1080 1920];
fov = 5; ppd = 93;

for aa = 1:length(app_cont_model_vec)
    for sfm = 1:length(sf_mode_vec)
        for ll = 1:length(loss_type_vec)
            loss_type = loss_type_vec{ll};
            app_cont_model = app_cont_model_vec{aa}; % 'Peli', 'Kulikowski', 'Physical'
            sf_mode= sf_mode_vec{sfm}; % 'natural', 'image', 'image_peak', 'image_2cpd', 'image_4cpd', 'image_8cpd', 'cvvdp'

            figure();
            html_change_figure_print_size( gcf, 16, 12 );
            errors = [];
            for rr=1:length(Y_rr)
                Y_smax = Y_rr(rr);  % The peak (white) luminance of the source content
                Y_src_bg = 100;
                Y_src = Y_smax * 2.^x;  % Source content luminance levels

                %% Approximate the image frequency spectrum
                if contains(sf_mode, 'natural')
                    weight = zeros(1,numel(s_freq_vec)-1);
                    % weight(1) = 1e-1;
                    if contains(sf_mode, 'natural')
                        alpha = 1.8;  % Following natural image stats (This already considers energy)
                        for ss = 2:numel(s_freq_vec)
                            f1 = s_freq_vec(ss-1);
                            f2 = s_freq_vec(ss);
                            if alpha ~= 1
                                weight(ss-1) = (f2^(1-alpha) - f1^(1-alpha));
                            else
                                weight(ss-1) = log(f2) - log(f1);
                            end
                        end
                    end
                    weight = weight/sum(weight);

                    if contains(sf_mode, 'peak')
                        weight(weight==max(weight)) = 1;
                        weight(weight<1) =0;
                    elseif contains(sf_mode, 'cpd')
                        tmp_str = split(sf_mode, '_');
                        num_str = tmp_str{end};
                        sf_val = sscanf(num_str, '%d');
                        sf_idx = find(abs(s_freq_vec-sf_val) == min(abs(s_freq_vec-sf_val)));
                        weight(sf_idx) = 1;
                        weight(weight<1)=0;
                    end

                    N = numel(Y_bg);
                    Y_ab = zeros(numel(Y_rr),N);
                    s_freq = reshape(s_freq_vec(1:end-1), [1 1 1 1 length(s_freq_vec)-1]);
                    csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, ...
                        'luminance', Y_src + Y_src_bg, 'area', 1, 'eccentricity', 0 );
                    S_src = csf.sensitivity(csf_pars);

                    % Display luminance levels, no ambient light
                    Y_fg = Y_dmax .* 2.^x;
                    Y_dst = Y_fg + Y_bg;    % With ambient light

                    % Sensitivity for the display
                    csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, ...
                        'luminance', Y_dst, 'area', 1, 'eccentricity', 0 );
                    S_dst = csf.sensitivity(csf_pars);


                    % Physical contrast
                    c_phys = c_src .* Y_fg ./ (Y_fg + Y_bg); %% Test
                    c_src_after = c_src .* Y_src ./ (Y_src + Y_src_bg);

                    % Perceptual contrast based on selected model
                    switch app_cont_model
                        case 'Kulikowski'
                            c_perc_dst = max(repmat(c_phys, [1 1 1 1 length(s_freq_vec)-1]) - scale * 1 ./ repmat(S_dst, [1 1 1 length(c_src) 1]), 0);
                            c_perc_src = max(repmat(c_src_after, [1 1 1 1 length(s_freq_vec)-1]) - scale * 1 ./ repmat(S_src, [1 1 1 length(c_src) 1]), 0);
                        case 'Peli'
                            c_perc_dst = repmat(c_phys, [1 1 1 1 length(s_freq_vec)-1]) .* repmat(S_dst, [1 1 1 length(c_src) 1]) / scale;
                            c_perc_src = repmat(c_src_after, [1 1 1 1 length(s_freq_vec)-1]) .* repmat(S_src, [1 1 1 length(c_src) 1]) / scale;
                        case 'Physical'
                            c_perc_dst = repmat(c_phys, [1 1 1 1 length(s_freq_vec)-1]);
                            c_perc_src = repmat(c_src_after, [1 1 1 1 length(s_freq_vec)-1]);
                    end

                    % Dimension
                    % 1: Luminance level
                    % 2: Y_test_dmax sweep
                    % 3: Y_test_bg sweep
                    % 4: contrast levels
                    % 5: spatial frequency

                    % Contrast ratio - See Eq. (11)
                    W = repmat(reshape(weight, [1 1 1 1 length(s_freq_vec)-1]), [numel(Y_src) numel(Y_dmax) N length(c_src) 1]);

                    if strcmp(loss_type, 'log-l1')
                        log_diff = log10(c_perc_dst + epsilon) - log10(c_perc_src + epsilon);
                        c_loss = squeeze(mean(abs(W .* log_diff), [1 4 5]));
                    elseif strcmp(loss_type, 'log-l2')
                        log_diff = log10(c_perc_dst + epsilon) - log10(c_perc_src + epsilon);
                        c_loss = squeeze(mean((W .* log_diff).^2, [1 4 5]));
                    end
                    % errors = [];

                    Y_ref = Y_rr(rr);
                    D = abs(Y_dmax - Y_ref);
                    i_ref = find(D == min(D), 1);

                    % c_loss_ref = c_loss(i_ref, end); %% Y_dmax: Y_ref & Y_bg: 100 nits
                    for kk = 1:N
                        D = c_loss(:, kk);
                        ind = find(D == min(D), 1);
                        Y_ab(rr, kk) = Y_dmax(ind);
                    end
                    Y_ab_sum = Y_ab;

                    xx_raw = squeeze(Y_bg);
                    xx = logspace(log10(xx_raw(1)), log10(xx_raw(end)), Y_bg_samples);

                    pred = Y_ab_sum(rr,:);
                    ft = fit(xx_raw, pred', 'smoothingspline');
                    pred_interp = transpose(feval(ft, xx));

                    if rr ==1
                        plot( Y_src_bg, Y_rr(rr), 's', 'Color', COLORs(rr,:), ...
                            'MarkerFaceColor', COLORs(rr,:), 'MarkerSize', 8,'HandleVisibility','on');
                        hold on;
                        plot( xx, pred_interp, '-', 'Color', COLORs(rr,:) ,'LineWidth', 1,'HandleVisibility','on');
                        hold on;
                        plot( xx, Y_ref/xx(end).*xx, '--', 'Color', COLORs(rr,:),'LineWidth', 1,'HandleVisibility','on' );
                    else
                        plot( Y_src_bg, Y_rr(rr), 's', 'Color', COLORs(rr,:), ...
                            'MarkerFaceColor', COLORs(rr,:), 'MarkerSize', 8,'HandleVisibility','off');
                        hold on;
                        plot( xx, pred_interp, '-', 'Color', COLORs(rr,:) ,'LineWidth', 1,'HandleVisibility','off');
                        hold on;
                        plot( xx, Y_ref/xx(end).*xx, '--', 'Color', COLORs(rr,:),'LineWidth', 1,'HandleVisibility','off' );
                    end

                    Ess = Es( Es.Y_fg_ref==Y_ref,:);
                    ci = 2.^Ess.meanci_mean_log_Y_fg - 2.^Ess.mean_mean_log_Y_fg;
                    meas = 2.^Ess.mean_mean_log_Y_fg;
                    if rr ==1
                        errorbar( Ess.Y_bg, meas, ci(:,1), ci(:,2), 'o', 'Color', COLORs(rr,:),'LineWidth', 1,'HandleVisibility','on' );
                    else
                        errorbar( Ess.Y_bg, meas, ci(:,1), ci(:,2), 'o', 'Color', COLORs(rr,:),'LineWidth', 1,'HandleVisibility','off' );
                    end

                    legend({'Reference','Prediction', 'Physical', 'Measurements'}, 'Location', 'southeast');
                    [common_elements, indices] = ismember(Ess.Y_bg, xx_raw);
                    errors = [errors; abs(log10(pred(indices)') - log10(meas)).^2 ];

                end
                error = mean(errors);

                set_axis_tick_label( 'x', 'luminance', Y_bg );
                xlabel( 'Background luminance [cd/m^2]' );
                set_axis_tick_label( 'y', 'luminance', Y_dmax );
                ylabel( 'Display peak luminance [cd/m^2]' );

                grid on;
                hold on;
                set(gca, 'FontName', Font, 'FontSize', FontSize); % Axis font
                set(findall(gcf, '-property', 'FontName'), 'FontName', Font);
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', FontSize);

            end
            if result_save
                file_str ='exp3.png';
                saveas(gca,[savepath file_str]);
                file_str ='exp3.pdf';
                exportgraphics(gca, [savepath file_str], 'ContentType', 'vector');
            end
        end
    end
end