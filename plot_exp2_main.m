%% Evaluation of Supra-threshold Contrast Matching

% This is the code for generating Fig. 7 in the main paper.
% If you use this code, please consider citing the following paper:
% D. Kim, M. Ashraf, A. Chapiro, and R. K. Mantiuk, 
% 'Supra-threshold Contrast Perception in Augmented Reality'
% Conference Proceeding, Siggraph Asia 2025

D = readtable('data/exp2/exp2_average_data.csv');
FontName = 'Linux Biolinum';
FontSize = 18;
Label_FontSize = 16;

Y_amb = logspace( -0.2, 3 );
Y_dmax = 100;

if ~exist( 'CSF_castleCSF', 'class' )
    addpath( fullfile( pwd,  '..', 'castleCSF', 'matlab' ) );
end

csf = CSF_castleCSF();

%c_base = 0.3;

Y_mean = unique(D.Y_fg);

%s_freq = 2;

% app_cont_model = 'Peli';
app_cont_model = 'Kulikowski';

f = figure(1);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7*3, 3.5*3]); 
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [7*3, 3.5*3]); 
clf;

K = length(Y_mean);
COLORs = lines(K);
hh = [];
LineWidth = 1;

c_base_achrom = [0.3 0.8];
txt_pos = [0.007 0.007*log10(10*5) 0.007*log10(10*5)^2];
txt_idx = 1;

coldir_labs = { 'achromatic', 'red-green', 'yellow-violet' };

for cc=1:3 % for each colour direction

    Dcc = D(D.col_direction==cc,:);

    SFREQs = unique(Dcc.cpd);
    CBASEs = unique(Dcc.c_base);


    for ff=1:length(SFREQs) % for each frequency

        for gg=1:length(CBASEs) % for each base contrast

            if cc==1
                col = ff;
            else
                col = 1+cc;
            end
            subplot( 2, 4, col+(gg-1)*4 );

            c_base = CBASEs(gg);
            s_freq = SFREQs(ff);

            pp = 1;
            for kk=1:K % For each luminance

                Y = Y_mean(kk);
                % The parameters of the CSF
                dkl_dirs = eye(3);
                lms_delta = dkl2lms_d65( dkl_dirs(cc,:) );
                csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, 'luminance', Y, 'area', 1, 'eccentricity', 0, 'lms_delta', lms_delta );
                S_ref = csf.sensitivity( csf_pars );

                csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, 'luminance', (Y+Y_amb)', 'area', 1, 'eccentricity', 0, 'lms_delta', lms_delta );
                S_test = csf.sensitivity( csf_pars )';

                c_phys = c_base * Y./(Y+Y_amb);
                c_phys(c_phys<1./S_test) = nan;

                label = 'Physical contrast';
                ph = plot( Y_amb, c_phys, '-', 'LineWidth', LineWidth, 'Color', COLORs(kk,:), 'DisplayName', label  );
                if kk==1
                    hh(pp) = ph;
                    pp = pp+1;
                end
                hold on;

                if cc==1 && ff==1 && gg==1
                    ind = ceil(0.8*(length(c_phys) - nnz(isnan(c_phys))));
                    % text( Y_amb(ind), c_phys(ind), sprintf( 'Y_{ref}=%g [cd/m^2]', Y), ...
                    %     'Color', COLORs(kk,:), 'HorizontalAlignment', 'right' );
                    text( 1, txt_pos(txt_idx), sprintf( 'Y_{FG}=%g [cd/m^2]', Y), ...
                        'Color', COLORs(kk,:), 'HorizontalAlignment', 'left','Interpreter','tex','FontName',FontName,'FontSize', Label_FontSize );
                    txt_idx= txt_idx+1;
                end

                LSs = { '--', ':' };
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

                    ph = plot( Y_amb, c_match, LSs{mm}, 'LineWidth', LineWidth,'Color', COLORs(kk,:), 'DisplayName', LABELs{mm} );
                    if kk==1
                        hh(pp) = ph;
                        pp = pp+1;
                    end
                end

                Ds = Dcc(Dcc.c_base==c_base & Dcc.cpd==s_freq & Dcc.Y_fg==Y,:);
                
                c_lin = 10.^Ds.mean_log_c_phys_meas;
                err_neg = c_lin - (c_lin .* 10.^(-Ds.std_log_c_phys_measured * 1.96./sqrt(Ds.N)));
                err_pos = (c_lin .* 10.^( Ds.std_log_c_phys_measured* 1.96./sqrt(Ds.N) )) - c_lin;
                
                errorbar( Ds.Y_ambient, c_lin, err_neg, err_pos,...
                    'o', 'Color','k','MarkerFaceColor', COLORs(kk,:), 'MarkerEdgeColor', 'k' );

                %plot( Ds.Y_ambient, Ds.mean_c_phys_measured, 'o', 'MarkerFaceColor', COLORs(kk,:), 'MarkerEdgeColor', 'k' );
                set(gca, 'FontName', FontName, 'FontSize', FontSize);
            end

            hold off;

            set_axis_tick_label( 'x', 'luminance', Y_amb );
            set_axis_tick_label( 'y', 'large_contrast', [0.001 1] );
            % yticks([log2(0.01) log2(0.1) log2(1)]);

            if cc==1 && ff==1
                ylabel( 'Matched test contrast', 'FontName',FontName, 'Interpreter','tex' );
            else
                ylabel( [] );
            end
            if gg==2
                xlabel( 'Background luminance Y_{BG} [cd/m^2]', 'FontName',FontName, 'Interpreter','tex');
            else
                xlabel( [] );
            end
            if cc==1 && ff==1 && gg==2
                legend( hh, 'Location', 'Best' , 'Interpreter','tex');
            end
            ylim( [0.005 1]);
            title( { coldir_labs{cc}, sprintf( 'c_{ref}=%g; freq=%g [cpd]', c_base, s_freq ) }, 'FontWeight', 'normal','Interpreter','tex');
            xlim( [Y_amb(1) Y_amb(end)] );
            grid on;


        end
    end
end
% resolution = '-r300'; 
exportgraphics(gcf, 'result/exp2/exp2_result.pdf', 'ContentType', 'vector', 'Resolution', 600);
exportgraphics(gcf,'result/exp2/exp2_result.png','Resolution',300)
