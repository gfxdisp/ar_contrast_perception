%% Plot equi-perceived-contrast lines for autobrightness

% This is the code for generating Fig. 9 in the main paper.
% If you use this code, please consider citing the following paper:
% D. Kim, M. Ashraf, A. Chapiro, and R. K. Mantiuk, 
% 'Supra-threshold Contrast Perception in Augmented Reality'
% Conference Proceeding, Siggraph Asia 2025

figure(1);
html_change_figure_print_size( gcf, 16, 12 );
clf;
ddm = @(x) 5.464 .* (pi ./ 0.05 .* x) .^ 0.4922 + 41;
Font = 'Linux Biolinum';
FontSize= 16;

% Candidate test luminance
Y_test_cand = logspace( 0, 3.8, 256 );
Y_test_cand_5d = reshape( Y_test_cand, 1, 1, 1, 1, [] ); % 5th dimension

Y_bg_test = logspace( -2, 2, 8 );
Y_bg_test_4d = reshape( Y_bg_test, 1, 1, 1, [] ); % 4th dimension

trans_min = 0.2;
trans_max = 0.8;

T = min( (log10(Y_bg_test_4d))/2 * (trans_min - trans_max) + trans_max, trans_max);
cont_ratio = [64 16 4];

COLORs = lines(length(cont_ratio));

for rr=1:length(cont_ratio)

    Y_black_ref = 0.5;
    Y_ref = Y_black_ref * cont_ratio(rr);
    disp_dr = log2(cont_ratio(rr));

    for tt=1:3
        if tt==1 
            L = display_contrast_loss( Y_test_cand_5d, Y_bg_test_4d, Y_ref, 0, disp_dr );
            lstyle = '-';
        elseif tt==3
            L = display_contrast_loss( Y_test_cand_5d, T.*Y_bg_test_4d, Y_ref, 0, disp_dr );
            lstyle = '--';
        elseif tt==2
            L = abs(ddm(Y_bg_test_4d) - Y_test_cand_5d);
            lstyle = '-.';
        end


        L_p = squeeze(L);

        Y_test = zeros(size(Y_bg_test));
        for kk=1:length(Y_bg_test)
            ind = find( L_p(kk,:) == min(L_p(kk,:)), 1 );
            Y_test(kk) = Y_test_cand(ind);
        end

        if rr ==1
            hc = plot( Y_bg_test, Y_test, lstyle, 'LineWidth', 1.2, 'Color', COLORs(rr,:),'HandleVisibility','on');
        else
            if tt~=2
                hc = plot( Y_bg_test, Y_test, lstyle, 'LineWidth', 1.2, 'Color', COLORs(rr,:),'HandleVisibility','off');
            end
        end

        hold on
        if tt==1
            text( Y_bg_test(2), Y_test(2)*1.05, sprintf( '%d:1', cont_ratio(rr)), 'Color', COLORs(rr,:), 'VerticalAlignment', 'bottom', 'FontSize', 14 );
        end
        drawnow;
    end
end
legend({'Perceptual AB (ours)',  'DDM [Hou et al., 2021]','Ambient dimming'});
legend( 'Location', 'northwest');

set_axis_tick_label( 'x', 'luminance', Y_bg_test );
xlabel( 'Background luminance [cd/m^2]' );
set_axis_tick_label( 'y', 'luminance', Y_test_cand );
ylabel( 'Display peak luminance [cd/m^2]' );
ylim( [Y_test_cand(1) Y_test_cand(end)] );

yyaxis right
P = display_power( Y_test_cand' )/1e3;
ylim( [Y_test_cand(1) Y_test_cand(end)] );
set( gca, 'yscale', 'log' );

power_ticks = [0.4 0.5 1 2 4 8];

Y_tick = interp1( P, Y_test_cand, power_ticks );

set( gca, "YTick", Y_tick );
set( gca, "YTickLabel", num2str(power_ticks') );
set( gca, 'YMinorTick', 'off' );
set( gca, 'YColor', [0.5 0.5 0]);
ylabel( 'Power consumption [W]' );

set(gca, 'FontName', Font, 'FontSize', FontSize); % Axis font
set(findall(gcf, '-property', 'FontName'), 'FontName', Font);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', FontSize);

exportgraphics( gcf, "result/autobrightness/app_autobrightness.pdf" );