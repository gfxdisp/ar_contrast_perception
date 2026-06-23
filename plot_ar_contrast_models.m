FontName = 'Linux Biolinum';
FontSize = 16;
Label_FontSize = 14;
Legend_FontSize = 14;
Y_amb = logspace( -2, 2 );
Y_dmax = 100;

if ~exist( 'CSF_castleCSF', 'class' )
    addpath( fullfile( pwd,  '..', 'castleCSF', 'matlab' ) );
end

csf = CSF_castleCSF();

Y_mean = Y_dmax * 2.^(0:-3:-10);

app_cont_model = 'Peli';
%app_cont_model = 'Kulikowski';


html_change_figure_print_size( gcf, 16, 12 );
clf;

K = length(Y_mean);
COLORs = lines(K);

c_base = 0.2;
s_freq = 4;

hh = [];
pp = 1;
for kk=1:K % For each luminance

    Y = Y_mean(kk);
    % The parameters of the CSF
    csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, 'luminance', Y, 'area', 1, 'eccentricity', 0 );
    S_ref = csf.sensitivity( csf_pars );

    csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, 'luminance', (Y+Y_amb)', 'area', 1, 'eccentricity', 0 );
    S_test = csf.sensitivity( csf_pars )';

    c_phys = c_base * Y./(Y+Y_amb);
    c_phys(c_phys<1./S_test) = nan;

    label = 'Physical contrast';
    ph = plot( Y_amb, c_phys, '-', 'Color', COLORs(kk,:), 'DisplayName', label  );
    if kk==1
        hh(pp) = ph;
        pp = pp+1;
    end
    hold on;

    ph = plot( Y_amb, 1./S_test, ':', 'Color', COLORs(kk,:), 'DisplayName', 'Detection threshold' );
    if kk==1
        hh(pp) = ph;
        pp = pp+1;
    end

    %if cc==1 && ff==1 && gg==1
        ind = ceil(0.8*(length(c_phys) - nnz(isnan(c_phys))));
        text( Y_amb(ind), c_phys(ind), sprintf( 'Y_{FG}=%.3g [cd/m^2]', Y), ...
            'Color', COLORs(kk,:), 'HorizontalAlignment', 'right','FontName', FontName, 'FontSize', Label_FontSize );
    %end

    LSs = { '--', '-.' };
    LABELs = { 'Peli''s model', 'Kulikowski''s model' };
    for mm=1:2

        switch mm
            case 1 %'Peli'
                c_match = S_test .* c_phys ./ S_ref;
            case 2 %'Kulikowski'
                c_match = (c_phys - 1./S_test) + 1./S_ref;
        end

        % We need to remove the points for which we cannot match with the
        % reference
        c_match(c_match<1/S_ref) = nan;

        ph = plot( Y_amb, c_match, LSs{mm}, 'Color', COLORs(kk,:), 'DisplayName', LABELs{mm} );
        if kk==1
            hh(pp) = ph;
            pp = pp+1;
        end
    end
end

hold off;

set_axis_tick_label( 'x', 'luminance', Y_amb );
set_axis_tick_label( 'y', 'contrast', [0.001 1] );

ax = gca;
ax.XAxis.FontName  = FontName;
ax.XAxis.FontSize  = Label_FontSize;

ax.YAxis.FontName  = FontName;
ax.YAxis.FontSize  = Label_FontSize;


ylabel( 'Contrast' , 'FontSize', Label_FontSize, 'FontName', FontName, 'Interpreter', 'tex');
xlabel( 'Background luminance Y_{BG} [cd/m^2]' ,'FontSize', Label_FontSize, 'FontName', FontName, 'Interpreter', 'tex');

legend_handle = legend( hh, 'Location', 'Southwest' );
set(legend_handle, 'FontName', FontName, 'FontSize', Legend_FontSize, 'Interpreter', 'tex');

ylim( [0.01 0.25]);
xlim( [Y_amb(1) Y_amb(end)] );

exportgraphics( gcf, 'result/model/ar_contrast_models.pdf' );