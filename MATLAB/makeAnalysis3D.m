function data = makeAnalysis3D( data )

%% set up the defaults
% what is the maximum number of peaks to be fit?
Nloop = data.max_peak_num;

% Tvec is the vector of peak parameters
Tvec = {};

% copy decimated genomic positions into xx
xx = data.xdec;

% decimated conversion frequencies
DD = data.Ddec;

% grab default peak parameters
T0 = data.T0;

% get the background conversion level
u0 = data.T0.u0;

% get the length-scale of the peak
L0 = data.T0.L0;
       
%% show raw data first
figure(4);
clf;
plot( xx, DD, '.-' );
hold on;

xlim( [0,xx(end)] );
ylabel( 'Allele frequency' );
xlabel( 'Genomic position: x (bp)' );

% make 
[T_, mll_vec, mll] = makeMLE( DD, xx, T0);


u_h = data.T0.u0 + data.x*0;
u_l = data.T0.u0 + data.xdec*0;

% Loop through all the putative peaks stopping at either the total possible 
% number numel(T_) or the specified max number requested Nloop, which ever
% is smaller.
for ii = 1:(min([numel(T_),Nloop]))
    
    % Set the current peak params to the ii-th peak parameters found: 
    That = T_(ii);
    
    % let the minus-log-likelihood to the ii-th value
    mll_min = mll_vec(ii);
    
    % set the mean level of the null hypothesis to the current model
    u_l0 = u_l;
    
    % optimize all the parameter values carefully, first on the decimanted
    % data and then on the full resultion data.
    That = makeOpti(data, That, DD, u_h, u_l );
    
    % use these optimized parameter values to make a new model for the mean
    % which contains ONLY the current peak.
    uhat = u0 + makedu(xx, That);
    
    % figure out what position to label as the peak (ind)
    [umax,ind] = max(uhat);
    That.umax = umax; % conversion ratio at top of the peak.
    That.xmax = xx(ind); % approx genomic position of the max
    
    % Now make a model that adds the current peak to already accepted peaks
    % to generate a mean profile for the alternative hypothesis.
    u_h = u_h + makedu(data.x, That); % high-res copy
    u_l = u_l + makedu(data.xdec, That); % low-res copy
    
    % Make minus-log-likelihood for alternative hypothesis (1) and null
    % hypothesis (0)
    mll_1 = makeMLL( DD, u_l-That.u0, That.u0, That.s );
    mll_0 = makeMLL( DD, u_l0-That.u0, That.u0, That.s );
    
    % difference in minus-log-likelihoods
    mll = mll_1-mll_0;
    mll_min = mll;
    
    % test statistic is log-likelihood difference.
    That.lam  = -mll_min;
    
    % convert the test statistic to a p-val
    pval = makeP(-mll_min, data.datap);
    
    % if the p value is too small
    if pval == 0
        % use the log10 p values
        lpval = makelp3D(-mll_min, data.datap);
    else
        % if it is large enough, just log10 it directly.
        lpval = log10( pval );
    end
    
    % Write the p value to the screen
    disp( "P value is " );
    disp( pval );
    
    % store the p value in the model.
    That.pval = pval;
    That.lpval = lpval;
    
    % write the test statistic to the screen.
    disp( 'Test statistic is' );
    disp(-mll_min);
    
    % store the optimized parameter values in the vector 
    Tvec{ii} = That;
        
end

% extract a vector of p values
lpval = drill(Tvec,'.lpval' );

% sort the p values and then reorder the vector of peak parameters from
% smallest to largest p value
[~,ind] = sort( lpval );
data.Tvec = Tvec(ind);

data = makeWeight( data );

end

%% [Thatp,lam] = makeOpti(data, That, DD, u_h, u_l)
% this function does a full optimization of the peak parameters, first on
% the low-resolution (decimated) data and then on the full resolution data.
function [Thatp,lam] = makeOpti(data, That, DD, u_h, u_l)

% That is theta hat i.e. the MLEs. 

%% First we will perform the optimization of the decimated (i.e. low
% resolution data.

% peak position
x0 = That.x0;

% low resolution genomic positions
x  = data.xdec;

% only calculate the mll over a finite range 32x the default scale L0
DL = 32*That.L0;

% minus (m) and plus (p) side of the range.
xm = x0-DL;
xp = x0+DL;

% make a flag for the positions that are in this range
ind = and(x>xm,x<xp);

% make the shorted version of the postions and conversion ratios
x_ = x(ind);
D_ = DD(ind);

% make a three vector for optimization
T0 = [That.x0,... % position of peak
      That.L0,... % genomic scale of peak
      That.a0];   % shape paraemter of peak
% Note that I0 is optimized analytically.

% That are the start values
% Thatp are the post optimization values
Thatp = That;
u0    = That.u0;  % mean of signal
s     = That.s;   % std of signal

% Make null hypothesis mean on just the local range
uu = u_l(ind);

% do the least squares optimization
T_fit = lsqnonlin( @fitter, T0 );

% Load the fit values into Thatp
Thatp.x0 = T_fit(1);
Thatp.L0 = T_fit(2);

Thatp.a0 = T_fit(3);
Thatp.a0 = 1;

%% Next we will perform the optimization on the high resolution data.

% copy the fit values into That
That = Thatp;

% get full resolution genomic positions
x  = data.x;

% get only the converion ratios and positions in the right range
ind = and(x>xm,x<xp);

x_ = x(ind);
D_ = data.D(ind);

% Make the three vector for fitting
T0 = [That.x0,That.L0,That.a0];

Thatp = That;
u0    = That.u0; 
s     = That.s; 

% Make the high-resolution null hypothesis mean over the right range
uu = u_h(ind);

% Do the least-squares optimization and keep the jacobian 
[T_fit,~,~,~,~,~,jacobian]  = lsqnonlin( @fitter, T0 );

% Copy the fit values into Thatp
Thatp.x0 = T_fit(1);
Thatp.L0 = T_fit(2);
Thatp.a0 = T_fit(3);

% keep the jacobian
Thatp.jac= jacobian;

% calculate the Fisher information
tmp = jacobian'* jacobian/s^2;
Thatp.fisher =  full(tmp(1:2,1:2));
Thatp.fisherI = inv(Thatp.fisher);


du  = makedu( x_, Thatp );
[~, Ihat] = makeMLL( D_, du, u0, s);

% Now extract the amplitude value using the algerbraic trick.
Thatp.I0 = Thatp.I0*Ihat;

    % The fit function
    % lsqnonlin minimizes the sum of the squares of dy
    function dy = fitter( TT )
        
           a = TT(3);
           
           if a<.5
               a=.5;
           elseif a>2
               a = 2;
           end
        
            Thatp.x0 = TT(1);
            Thatp.L0 = TT(2);
            Thatp.a0 = a;
            
            

            du  = makedu( x_, Thatp );

            dy = D_-du-uu;
            
            % show the fitting as it happens
            figure(55);
            clf;
            plot( x_, D_, '.-b' );
            hold on;
            plot( x_, du+uu, '.-r' );
            drawnow;
            


    end
end

%% makeMLE( I, x__, T)
% This function generates mll for a peak located in every position along
% the genome (fitting the amplitude only) in the interest of speed.
function [T_, mll_min, mll] = makeMLE( I, x__, T)

xx = x__;

[T_, mll_min, mll] = makeMLE_( I, x__, T);

%% Show the raw data versus the log-likelihood 
% of the peak as a function of position. Peaks in the raw data correspond
% with peaks in the log-likhood.
figure(57);
clf;
subplot(2,1,1);
plot( x__, I, '.-' );
ylabel( 'Conversion Ratio' );
xlabel( 'Genomic Position: x (bp)' );
xlim( [1,x__(end)]);

subplot(2,1,2);
plot( x__, -mll,'.-');
ylabel( 'Log Likelihood' );
xlabel( 'Genomic Position of Peak: x (bp)' );
xlim( [1,x__(end)]);

%% identification of peaks
% since the peaks have non-trivial extent, we need to identify the local
% maxima in the log-likelihoods
mmll = -mll;

% make a flag that is true when a pixel is larger than its neighbors
local_max = [false,...
    and(mmll(2:end-1)>mmll(1:end-2),...
        mmll(2:end-1)>mmll(3:end)),false];

% make a list of these local maxes
ind = find(local_max);

% and then sort by the log-likelihoods
[mmll_ind,indind] = sort(mmll(ind),'descend');

% cut out all of those below a critical value
nnn = sum(mmll_ind>5);

% resort the indices in order of size of the ll
ind = ind(indind);
% and cut them off with the cutoff
ind = ind(1:nnn);

% grab the mll min values 
mll_min = -mmll_ind(1:nnn);


%% Take out neighbors
% now that we have an ordered list, we start at the top and zero out
% anyhing within 400 bp of a larger peak
ind0 = ind;
for ii = 1:numel(ind)
   
    x00 = xx(ind(ii));
    x11 = xx(ind);
    
    % find all indices with positions closer than 400 bp
    ind_rem = find( abs(x00-x11)< 400 );
    
    % and consider only weaker peaks (i.e. larger than the current index)
    ind_rem = ind_rem(ind_rem>ii);
    
    % ero these indices out
    ind0(ind_rem) = 0;
end

% keep only the indices that haven't been zeroed out.
ind = ind0(ind0>0);
mll_min = mll_min(ind0>0);

% now the putative peaks have been pruned, make peak parameter vectors for
% each remaining peak position...

nind = numel(ind);

T_(nind) = T;
T0 = T;

for ii = 1:nind
    
    T = T0;
    xhat = xx(ind(ii));
    
    T.x0 = xhat;
    du = makedu( x__, T );
    
    [~,Ihat] = makeMLL( I, du, T.u0, T.s );
    
    T.I0 = T.I0*Ihat;
    T_(ii) = T;
    
    if 0
        figure(56);
        clf;
        
        plot( x__, I, '.' );
        hold on;
        plot( x__, Ihat*du+T.u0, 'LineWidth', 1 );
    end
end

end





%% makeP( lam, datap )
% this function makes the p value using the gumbel model
function [p] = makeP( lam, datap )


p = 1-exp( -exp( -(lam-datap.ull)/datap.sll ) );


end

%% makeWeight( data );
% this function estimates the total activity associated with each peak
% relative to the overall integrated activity.
function data = makeWeight( data )

L = data.T0.L0;

x0 = drill(data.Tvec, '.x0')';
lp = drill(data.Tvec, '.lpval')';

x0_real = x0(lp<log10(.05));

x_  = ones( size( x0))*data.x;
x0_ = x0*ones( size(data.x));

dx_ = abs(x_-x0_);
dx  = min( dx_,[],1);

no_peak = dx>16*L;

u_back = mean( data.D(no_peak));

N_back = u_back*numel( data.D );
N_tot  = sum( data.D );
p_back = u_back*numel( data.D );

data.p_back = p_back;
data.N_back = N_back;
data.N_tot  = N_tot;


for ii = 1:numel(data.Tvec)
        x0 = data.Tvec{ii}.x0;
        flag = abs(data.x-x0)<16*L;
        data.Tvec{ii}.N_peak = sum( data.D(flag)-u_back);
        data.Tvec{ii}.p_peak = data.Tvec{ii}.N_peak/N_tot;
end


end

