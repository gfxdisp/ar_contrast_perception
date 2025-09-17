function xyz = cm_rgb2xyz(rgb, rgb_space)
%RGB2XY  Convert RGB to CIE XYZ 1931 2-deg color space. Choose from several RGB
%        variants
%
% xyz = rgb2xyz(rgb)
% xyz = rgb2xyz(rgb, rgb_space)
%
% "rgb" could be an image [width x height x 3] or [rows x 3] matrix of colour
%     vectors
% "rgb_space" can be one of: 'Adobe', 'NTSC', 'rec709', 'rec2020'
%     The default is 'rec709'


if ~exist( 'rgb_space', 'var' )
    rgb_space = 'rec709';
end

switch rgb_space
    case 'Adobe'
        M_rgbxyz = [
            0.576700  0.297361  0.0270328
            0.185556  0.627355  0.0706879
            0.188212  0.0752847 0.991248];
    case 'NTSC'
        M_rgbxyz = [
            0.6068909  0.1735011  0.2003480
            0.2989164  0.5865990  0.1144845
            -0.0000000  0.0660957  1.1162243];
    case { 'sRGB', 'rec709' }
        M_rgbxyz = [
            0.4124 0.3576 0.1805;
            0.2126 0.7152 0.0722;
            0.0193 0.1192 0.9505];
    case 'rec2020'                   
        % Source: https://colour.readthedocs.io/en/v0.3.7/colour.models.dataset.rec_2020.html
        M_rgbxyz = ( [0.636953507, 0.144619185, 0.168855854;
               0.262698339, 0.678008766, 0.0592928953;
               4.99407097e-17, 0.0280731358, 1.06082723] );        
    otherwise
        error( 'Unknown RGB colour space' );
end

xyz = cm_colorspace_transform( rgb, M_rgbxyz );

end
