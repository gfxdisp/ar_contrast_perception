function L = display_contrast_loss( Y_dmax_test, Y_bg_test, Y_dmax_ref, Y_bg_ref, disp_dr )
% Reproduces Eq. (12) in the paper

c_src = logspace( -2, log10(0.9), 32 );  % contrast levels to consider (for Kulikowski's model)
f_min = 1;
f_max = 30;
alpha = 0.8;
s_freq = (linspace( 0, 1, 32 ) * (f_max^-alpha-f_min^-alpha) + f_min^-alpha).^(-1/alpha);

y_rel = logspace( -disp_dr, 0, 64 );

% dimensions: 
% 1 - s_freq
% 2 - luminance
% 3 - contrast

s_freq = s_freq';
Y_test = Y_dmax_test .* y_rel;
c_src = reshape( c_src, 1, 1, [] );

cp_test = c_perc( Y_test, Y_bg_test, c_src, s_freq );

Y_ref = Y_dmax_ref.* y_rel;
cp_ref = c_perc( Y_ref, Y_bg_ref, c_src, s_freq );

epsilon = 1e-12;

L = sum( abs( log(cp_test+epsilon) - log(cp_ref+epsilon) ), [1 2 3] );

end

function cp = c_perc(Y, Y_bg, c, s_freq)

if ~exist( 'CSF_castleCSF', 'class' )
    addpath( fullfile( pwd,  '..', '..', 'castleCSF', 'matlab' ) );
end

csf = CSF_castleCSF();

csf_pars = struct( 's_frequency', s_freq, 't_frequency', 0, 'orientation', 0, ...
    'luminance', Y + Y_bg, 'area', 1, 'eccentricity', 0 );
t = 1./csf.sensitivity(csf_pars);

c_disp = c.*(Y./(Y+Y_bg));

cp = max( c_disp - t, 0 );

end