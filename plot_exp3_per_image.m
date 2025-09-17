%% Generalization to complex images in AR per image

% This is the code for generating Figures in supp.
% If you use this code, please consider citing the following paper:
% D. Kim, M. Ashraf, A. Chapiro, and R. K. Mantiuk, 
% 'Supra-threshold Contrast Perception in Augmented Reality'
% Conference Proceeding, Siggraph Asia 2025

%% AR contrast paper condition: ('natural', 'Kulikowski', 'log-l1')

%% Power spectrum
% 'natural': power spectrum of natural image statistics
% 'image': power spectrum of individual image

%% Supra-threshold contrast model
% 'Peli','Kulikowski', 'Physical'

%% Operation & loss type
% {operation}-{error_dist}
% operation: log, linear, sqrt
% error_dist: 'l1', 'lp' (p=1.5), 'l2'

clear all; close all;

savepath = './result/exp3/';
imgpath = './image/';
mkdir(savepath)

Font = 'Linux Biolinum';
FontSize= 14;

binocular = false;
if binocular ==true
    scale = 1;
    exp_mode = 'Stereoscopic';
else
    scale = sqrt(2);
    exp_mode = 'Haploscopic';
end

Y_bg = logspace( -2, 2, 4*2+1 );
Y_bg = reshape( Y_bg, [1 1 length(Y_bg)] );
Y_dmax = logspace( -1, log10(5000), 1024 ); % AR display peak luminance levels to consider
x = (0:-1:-10)'; % log-2 relative luminance levels of the display
Y_bg_samples = 4*10+1;

%% Add path of castleCSF
if ~exist( 'CSF_castleCSF', 'class' )
    addpath( fullfile( pwd,  '..', '..', 'castleCSF', 'matlab' ) );
end

csf = CSF_castleCSF();
c_src = logspace( -2, log10(0.9), 32 );  % contrast levels to consider (for Kulikowski's model)
c_src = reshape( c_src, [1 1 1 length(c_src)] );

s_freq_vec = logspace(log10(1), log10(30), 31);
E = readtable('data/exp3/exp3_raw_data.csv');

E.log_Y_fg = log2(E.Y_fg);
% Mean by observer first
Eg = grpstats( E, { 'image','observer','Y_bg', 'Y_fg_ref' }, { 'mean'}, 'DataVars', 'log_Y_fg' );
% and get the ci based on the mean value per observer

app_cont_model_vec = {'Kulikowski'}; % 'Peli', 'Kulikowski', 'Physical'
sf_mode_vec = {'natural'}; % 'natural', 'image'
loss_type_vec = {'log-l1'}; % 'log-l1', 'log-lp', 'log-l2',  'linear-l1', 'linear-lp', 'linear-l2', 'sqrt-l1', 'sqrt-lp', 'sqrt-l2'

epsilon = 1e-12;
p=1.5;

result_save =true;
ab_save = false;

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

            img_vec = unique(Eg.image);
            img_error = zeros(numel(img_vec), 1);

            for i =1:numel(img_vec)
                img_str = string(img_vec(i));
                Es = Eg(Eg.image == img_str,:);
                Es = grpstats( Es, { 'image', 'Y_bg','Y_fg_ref'}, {'mean', 'meanci'},'DataVars','mean_log_Y_fg');
                figure();
                errors = [];
                N = numel(Y_bg);
                Y_ab = zeros(numel(Y_rr),N);
                for rr=1:length(Y_rr)

                    Y_smax = Y_rr(rr);  % The peak (white) luminance of the source content
                    Y_src_bg = 100;
                    Y_src = Y_smax * 2.^x;  % Source content luminance levels

                    %% Approximate the image frequency spectrum
                    if contains(sf_mode, 'natural') || contains(sf_mode, 'image')
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
                        elseif contains(sf_mode, 'image')
                            weight = zeros(1,numel(s_freq_vec)-1);

                            files = dir(imgpath);
                            fileNames = {files.name};

                            % Filter the file names to find those that contain the search string
                            matchedFiles = fileNames(contains(fileNames, img_str));
                            img = cm_rgb2xyz(srgb2lin(im2double(imread(strjoin([imgpath string(matchedFiles)],'')))));
                            img = img(:,:,2);
                            img = imresize(img, img_res);
                            F = fftshift(fft2(img));
                            fx = linspace(-ppd/2,ppd/2,img_res(2));
                            fy = linspace(-ppd/2,ppd/2,img_res(1));
                            [FX FY] = meshgrid(fx,fy);

                            FR = sqrt(FX.^2+FY.^2);

                            for ss = 2:numel(s_freq_vec)
                                f1 = s_freq_vec(ss-1);
                                f2 = s_freq_vec(ss);
                                fm = (abs(FR)>=f1).*(abs(FR)<f2);
                                weight(ss-1) = sum(sum(abs(F.*fm).^2));
                            end
                            total_energy = sum(sum(abs(F(:)).^2));
                            weight = weight / total_energy;
                        elseif contains(sf_mode, 'default')
                            weight = ones(1,numel(s_freq_vec)-1);
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

                        W = repmat(reshape(weight, [1 1 1 1 length(s_freq_vec)-1]), [numel(Y_src) numel(Y_dmax) N length(c_src) 1]);

                        if strcmp(loss_type, 'log-l1')
                            log_diff = log2(c_perc_dst + epsilon) - log2(c_perc_src + epsilon);
                            c_loss = squeeze(mean(abs(W .* log_diff), [1 4 5]));
                        elseif strcmp(loss_type, 'log-l2')
                            log_diff = log10(c_perc_dst + epsilon) - log10(c_perc_src + epsilon);
                            c_loss = squeeze(mean((W .* log_diff).^2, [1 4 5]));
                       elseif startsWith(loss_type, 'log-lp')
                            log_diff = log10(c_perc_dst + epsilon) - log10(c_perc_src + epsilon);
                            c_loss = squeeze(mean((W .* abs(log_diff)).^p, [1 4 5])).^(1/p);     
                        elseif strcmp(loss_type, 'linear-l1')
                            diff = c_perc_dst - c_perc_src ;
                            c_loss = squeeze(mean(abs(W .* diff), [1 4 5]));
                        elseif strcmp(loss_type, 'linear-l2')
                            diff = c_perc_dst - c_perc_src ;
                            c_loss = squeeze(mean((W .* diff).^2, [1 4 5]));
                        elseif startsWith(loss_type, 'linear-lp')
                            diff = c_perc_dst - c_perc_src ;
                            c_loss = squeeze(mean((W .* abs(diff)).^p, [1 4 5])).^(1/p);
                        elseif strcmp(loss_type, 'sqrt-l1')
                            diff = sqrt(c_perc_dst) - sqrt(c_perc_src) ;
                            c_loss = squeeze(mean(abs(W .* diff), [1 4 5]));
                        elseif strcmp(loss_type, 'sqrt-l2')
                            diff = sqrt(c_perc_dst) - sqrt(c_perc_src) ;
                            c_loss = squeeze(mean((W .* diff).^2, [1 4 5]));
                        elseif startsWith(loss_type, 'sqrt-lp')
                            diff = sqrt(c_perc_dst) - sqrt(c_perc_src) ;
                            c_loss = squeeze(mean((W .* abs(diff)).^p, [1 4 5])).^(1/p);

                        end

                        Y_ref = Y_rr(rr);
                        D = abs(Y_dmax - Y_ref);
                        i_ref = find(D == min(D), 1);
                        for kk = 1:N
                            D = c_loss(:, kk);
                            ind = find(D == min(D), 1);
                            Y_ab(rr, kk) = Y_dmax(ind);
                        end
                    end
                    Y_ab_sum = Y_ab;
                    Y_ref = Y_rr(rr);
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
                    % pred(indices)
                end
                if ab_save && contains(sf_mode, 'natural') && strcmp(app_cont_model, 'Kulikowski') && (binocular==false)
                    % Y_bg_vec = Y_bg;
                    % cr_vec = Y_ab_sum;
                    xx_raw = squeeze(Y_bg);  % original 5 points
                    Y_bg_vec = logspace(log10(xx_raw(1)), log10(xx_raw(end)), 41);  % interpolate to 41 points
                    cr_vec = zeros(length(Y_rr), numel(Y_bg_vec));
                    for rr=1:length(Y_rr)
                       Y_ab_tmp = Y_ab_sum(rr,:);
                       ft = fit(xx_raw, Y_ab_tmp', 'smoothingspline');
                       cr_vec(rr,:) = feval(ft, Y_bg_vec);

                    end
                    % ft = fit(xx_raw, pred', 'smoothingspline');
                    % pred = feval(ft, xx);


                    save_filename = sprintf("%s_%s.mat",app_cont_model, sf_mode);
                    save(save_filename, 'Y_bg_vec', 'cr_vec', 'Y_rr');
                end


                error = mean(errors);
                img_error(i) = error;
                set_axis_tick_label( 'x', 'luminance', Y_bg );
                xlabel( 'Background luminance [cd/m^2]' );
                set_axis_tick_label( 'y', 'luminance', Y_dmax );
                ylabel( 'Display peak luminance [cd/m^2]' );
                if contains(sf_mode, 'natural')
                    title_str =sprintf( ' %s / RMSE: %.3f ', img_str, sqrt(error));
                    % title_str =sprintf( ' %s / %s / $1/f^{%.1f}$/ %s  / RMSE(log10):%.3f ', app_cont_model, loss_type, alpha, img_str, sqrt(error));
                elseif contains(sf_mode, 'image')
                    title_str =sprintf( ' %s / %s / %s / %s  / RMSE(log10):%.3f ', app_cont_model, loss_type, sf_mode, img_str, sqrt(error));
                end

                title( title_str , 'Interpreter','latex');
                grid on;

                set(gca, 'FontName', Font, 'FontSize', FontSize); % Axis font
                set(findall(gcf, '-property', 'FontName'), 'FontName', Font);
                set(findall(gcf, '-property', 'FontSize'), 'FontSize', FontSize);

                if result_save
                    file_str =sprintf( '%s_%s_%s_weight_%s_%s_%d.png', exp_mode, loss_type, app_cont_model, sf_mode, img_str, length(unique(E.observer)));
                    saveas(gca,[savepath file_str]);
                    file_str =sprintf( '%s_%s_%s_weight_%s_%s_%d.pdf', exp_mode, loss_type, app_cont_model, sf_mode, img_str, length(unique(E.observer)));
                    exportgraphics(gca, [savepath file_str], 'ContentType', 'vector');
                end
            end
            close all;
            sprintf('%s / %s / %s / %s /  Total error RMSE(log10(pred) - log10(measured)): %.4f', exp_mode, loss_type, app_cont_model, sf_mode, sqrt(mean(img_error)))
        end
    end
end