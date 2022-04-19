%% makeMLE_( I, x__, T)
% This function generates mll for a peak located in every position along
% the genome (fitting the amplitude only) in the interest of speed.
function [That, mll_min, mll, mll_0] = makeMLE_( I, x__, T)

% I are the conversion ratios at the genomic positions 

% x__ are the genomic positions

% T are default values for the peak parameters. In this function, only the
% position and amplituded are fit and the other "shape" parameters are left
% unchanged.

% make a copy of the genomic positions
xx = x__;

% minus-log-likelihood vector initialized with nans
mll = nan(size(xx));

% make a null hypothesis copy of the default parameter vector
T0 = T;
T1 = T;
% for which we set the peak amplitude to 0
T0.I0 = 0;

% make a model mean with peak amplitude 0
du0 = makedu( x__, T0 );

% calculate the MLL for this configuration
mll_0 = makeMLL( I, du0, T0.u0, T0.s );
% since the peak only has local support, it is most efficient to only 
% recalculate the likelihoods for positions close to the peak... we will 
% therefore use this as the "background" value and replace only parts of 
% this sum in the vicinity of the peak.    

% get the total number of positions considered
nn = numel(xx);

% loop through all genomic positions as putative positions for the peak
for ii = 1:nn
    
    % Get the neighborhood of values to replace
    ranger = ii+(-30:30);
    ranger(ranger<1) = ranger(ranger<1)+nn;
    ranger(ranger>nn) = ranger(ranger>nn)-nn;
    
    % Set the peak position to the current position
    T1.x0 = xx(ii);
    
    % make the model means for the null (0) and alternative (1) hypotheses
    du1 = makedu( x__(ranger), T1 );
    du0 = makedu( x__(ranger), T0 );
    
    % Use the background mll value and replace only the values in the
    % vicinity of the peak
    mll(ii) = mll_0 + ...
             +makeMLL( I(ranger), du1, T1.u0, T1.s )+...
             -makeMLL( I(ranger), du0, T0.u0, T0.s );
end

% Find out which genomic position considered minimized the
% minus-log-likelihood (mll)
[mll_min, ind] = min( mll );

% ind is the minimizing index... get that position and put it into T
xhat = xx(ind);
T.x0 = xhat;

% caculated the updated mean for this peak position.
du = makedu( x__, T );

% get the optimal peak amplitude for this peak position and load it into T
[~,Ihat] = makeMLL( I, du, T.u0, T.s );
T.I0 = T.I0*Ihat;

% That is the estimated peak parameters.
That = T;

end