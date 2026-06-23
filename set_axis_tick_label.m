function set_axis_tick_label( ax_name, unit, range )
% Set axis ticks, scale and label for common units, such as cpd, or cd/m^2
%
% set_axis_tick_label( ax_name, unit, range )
%
% ax_name - 'x' for x-axis, 'y' for y-axis
% unit - 'luminance', 'sensitivity', or 'frequemncy'
% range - a vector defining the range. The first and last elements are used
%         as min and max values. 

switch unit

    case { "l", "luminance" }

        f_min = floor( log10(range(1)) );
        f_max = ceil( log10(range(end) ) );

        ticks = 10.^(f_min:f_max);

        xtick_labs = num2str( ticks' );
        axis_label = 'Luminance [cd/m^2]';                
        scale = 'log';
    
    case { "s", "sensitivity" }

        f_min = floor( log10(range(1)) );
        f_max = ceil( log10(range(end) ) );

        ticks = 10.^(f_min:f_max);

        xtick_labs = cell(length(ticks),1);

        for kk=1:length(ticks)
            if ticks(kk) < 0.01 || ticks(kk) > 1000
                xtick_labs{kk} = sprintf( '10^{%d}', round(log10(ticks(kk))) );
            else
                xtick_labs{kk} = num2str( ticks(kk) );
            end
        end
        axis_label = 'Sensitivity';        
        scale = 'log';

    case { "contrast", "Michelson", "large_contrast" } % Michelson contrast

        f_min = floor( log10(range(1)) );
        f_max = ceil( log10(range(end) ) );

%         if unit == "large_contrast"
            step = 1;
%         else
%             step = 0.2;
%         end
        ticks = 10.^(f_min:step:f_max);

        xtick_labs = num2str( ticks' );
        axis_label = 'Contrast';        
        scale = 'log';
        
    case { 'cpd', 'frequency' }
        
        f_min = floor( log2(range(1)) );
        f_max = ceil( log2(range(end) ) );

        ticks = 2.^(f_min:f_max);

        xtick_labs = cell(length(ticks),1);
        for kk=1:length(ticks)
            if ticks(kk) < 0.5
                xtick_labs{kk} = sprintf( '2^{%d}', round(log2(ticks(kk))) );
            else
                xtick_labs{kk} = num2str( ticks(kk) );
            end
        end        
        axis_label = 'Spatial frequency [cpd]';
        scale = 'log';
        
    case { 'deg', 'sigma' }
        
        f_min = floor( log2(range(1)) );
        f_max = ceil( log2(range(end) ) );
        
        ticks = 2.^(f_min:f_max);
        
        xtick_labs = cell(length(ticks),1);
        for kk=1:length(ticks)
            if ticks(kk) < 0.5
                xtick_labs{kk} = sprintf( '2^{%d}', round(log2(ticks(kk))) );
            else
                xtick_labs{kk} = num2str( ticks(kk) );
            end
        end
        axis_label = 'Sigma [deg]';
        scale = 'log';
        
    case { 'lambda' }
        
        f_min = floor( log2(range(1)) );
        f_max = ceil( log2(range(end) ) );
        
        ticks = 2.^(f_min:f_max);
        
        xtick_labs = cell(length(ticks),1);
        for kk=1:length(ticks)
            if ticks(kk) < 0.5
                xtick_labs{kk} = sprintf( '2^{%d}', round(log2(ticks(kk))) );
            else
                xtick_labs{kk} = num2str( ticks(kk) );
            end
        end
        axis_label = 'Lambda [n cycles]';
        scale = 'log';
    
    case { 'tf', 'Hz', 'temporal_frequency'}
        f_min = floor( log2(range(1)) );
        f_max = ceil( log2(range(end) ) );
        
        ticks = 2.^(f_min:f_max);
        
        xtick_labs = cell(length(ticks),1);
        for kk=1:length(ticks)
            if ticks(kk) < 0.5
                xtick_labs{kk} = sprintf( '2^{%d}', round(log2(ticks(kk))) );
            else
                xtick_labs{kk} = num2str( ticks(kk) );
            end
        end
        axis_label = 'Temporal frequency [Hz]';
        scale = 'log';
        
    otherwise 
        error( 'Unknown unit "%s"', unit );
    
end

switch ax_name
    case 'x'
        set( gca, 'XTick', ticks );
        set( gca, 'XTickLabel', xtick_labs );
        set( gca, 'XScale', scale );
        xlabel( axis_label );
    case 'y'
        set( gca, 'YTick', ticks );
        set( gca, 'YTickLabel', xtick_labs );
        ylabel( axis_label );
        set( gca, 'YScale', scale );
    case 'z'
        set( gca, 'ZTick', ticks );
        set( gca, 'ZTickLabel', xtick_labs );
        zlabel( axis_label );
        set( gca, 'ZScale', scale );
    otherwise
        error( 'ax_name must be either "x", "y" or "z"' );
end
       

end