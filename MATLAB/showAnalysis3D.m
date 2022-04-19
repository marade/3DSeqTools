%% showAnalysis3D( data, p_sig )
% This function visulaiizes the 3D seq analysis
% data -- analyzed 3D-seq data set
% p_sig -- significance level for including peaks
function data = showAnalysis3D( data, p_sig )

%% Set flags for data visualization here
show_p_values     = true;
show_single_peaks = true;
show_im_fig       = true;

%% Set significance level for peak inclusion
if ~exist( 'p_sig', 'var' ) || isempty( p_sig )
    if isfield( data, 'p_sig')
        p_sig = data.p_sig;
    else
        p_sig = 1e-2;
    end
end

data.p_sig = p_sig;

%% Make figures directory 'Figs'

fig_dirname = 'Figs';

if ~exist( fig_dirname, 'dir' );
    mkdir( fig_dirname );
end


%% get conversion fraction and positions out
xx = data.xdec;
DD = data.Ddec;
D__ = data.D;
x__ = data.x;

data = makeModelMeans3D( data );
du   = data.u_l;

%% Draw main model figure

% size up the figure for the Decimated data first
figure(1);
doPlot();
doPageFormat( [5,3] );
print( '-dpdf', [fig_dirname,filesep,'Full.pdf'] );

% %% Make tabs: High
% tab = table(  x__', D__', data.u_h',...
%      'VariableNames',{'Positions','Frequencies', 'Model'} );
% filename = [fig_dirname,filesep,'Full_high','.xls']
% writetable(tab,filename,'Sheet',1);

                
%% Make tabs: Low
tab = table(  xx', DD', data.u_l',...
     'VariableNames',{'Positions','Frequencies', 'Model'} );
filename = [fig_dirname,filesep,'Full_low','.xls']
writetable(tab,filename,'Sheet',1);



width        = drill(data.Tvec, '.L0')';
height       = drill(data.Tvec, '.I0')';
position     = drill(data.Tvec, '.x0')';
log10pval    = drill(data.Tvec, '.lpval')';
position_err = drill(data.Tvec, '.fisherI(1,1)')';
a0 = drill(data.Tvec, '.a0')';

tab = table( log10pval, position, position_err, height, width, a0,...
    'VariableNames',{'log10pval', 'position', 'position_err', 'height', 'width', 'a0'});
filename = [fig_dirname,filesep,'Model','.xls']
writetable(tab,filename,'Sheet',1);



%% Show single peak figures
if show_single_peaks
    fig_num = 100;
    if isfield( data, 'Tvec' )
        for jj = 1:numel( data.Tvec )
            
            That = data.Tvec{jj};
            pval = That.pval;
            
            if pval < data.p_sig
                figure(fig_num);
                doPlot();
                
                dX = 2e4;
                xlim( [-dX,dX] + That.x0 );
                drawnow;
                
                doPageFormat( [5,3] );
                print( '-dpdf', [fig_dirname,filesep,'Peak_',...
                    num2str( round(That.x0), '%07d'),'.pdf'] );
                
                
                %% Make tabs: High
                ranger = and( x__> (That.x0-dX), x__< (That.x0+dX) );
                tab = table(  x__(ranger)', D__(ranger)',data.u_h(ranger)',...
                    'VariableNames',{'Positions','Frequencies', 'Model'} );
                filename = [fig_dirname,filesep,'Peak_high_',...
                    num2str( round(That.x0), '%07d'),'.xls'];
                writetable(tab,filename,'Sheet',1);
                
                %% Make tabs: Low
                ranger = and( xx> (That.x0-dX), xx< (That.x0+dX) );
                tab = table(  xx(ranger)', DD(ranger)', data.u_l(ranger)',...
                    'VariableNames',{'Positions','Frequencies', 'Model'} );
                filename = [fig_dirname,filesep,'Peak_low_',...
                    num2str( round(That.x0), '%07d'),'.xls'];
                writetable(tab,filename,'Sheet',1);
                
                
                fig_num = fig_num + 1;
            end
            
        end
    end
end


%% plot function
    function doPlot
        clf;
        plot( xx, DD, '.-b' );
        ylim_vals = ylim;
        ylim_vals(2) = 3*ylim_vals(2);
        
        % Now clear the figure and start again.
        clf;
        % full resolution data in light blue.
        plot( x__, D__, '.-', 'Color', [.75,.75,1] );
        hold on;
        
        % decimated data in dark blue.
        plot( xx, DD, '.-b' );
        ylim( ylim_vals );
        
        plot( xx, du, 'r-', 'LineWidth', 1 );
        
        % construct the model mean
        if isfield( data, 'Tvec' )
            for ii = 1:numel( data.Tvec )
                
                That_ = data.Tvec{ii};
                
                pval = That_.pval;
                lp = makelp3D( That_.lam, data.datap );
                
                umax = That_.umax;
                xmax = That_.xmax;
                
                if show_p_values
                    if pval < data.p_sig
                        % If the peak is significant, include it
                        if pval>1e-3
                            text( xmax, umax, ...
                                ['  ',addComma(round(That_.x0)), ' bp', 10,...
                                '  p = ',num2str( pval, '%1.0e' )],...
                                'HorizontalAlignment', 'left' );
                            plot( xmax, umax,'r.' );
                        else
                            text( xmax, umax, ...
                                ['  ', addComma(round(That_.x0)), ' bp', 10,...
                                '  log_{10}p = ',num2str( lp, '%.1f' )],...
                                'HorizontalAlignment', 'left' );
                            plot( xmax, umax,'r.' );
                            
                        end
                    end
                end
            end
        end
        
        legend( {'High res','Low res','Model'} );
        
        ylabel( 'Allele frequency: r' );
        xlabel( 'Genomic position: x (bp)' );
        
        xlim( [0,xx(end)] );
    end
end


function numOut = addComma(numIn)
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out
end

