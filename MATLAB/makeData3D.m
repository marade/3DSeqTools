function data = makeData3D( x, D, s_p_cut, datap )

if ~exist( 's_p_cut', 'var' ) || isempty( s_p_cut );
    s_p_cut = 4;
end

if ~exist( 'datap', 'var' ) || isempty( datap );
    datap = load( 'datap_emp.mat' );
end


% Peak shape parameters
x0 =  10;    % Position of peak
I0 =   0.1;
L0 = 400;    % Width of peak in bp
D0 =   1.5;  % Smoothness -- Shape param 1
a0 =   1;    % Tails -- Shape param 2
Lg = x(end);

% genome coarse-graining scale
% make an effective genome where every dL bases corresponds to an effective
% base.
dL = 250;

% T is the 'theta' model parameter vector for fitting peaks.
% load in the default values
T0    = [];
T0.x0 = x0;
T0.I0 = I0;
T0.L0 = L0;
T0.D0 = D0;
T0.a0 = a0;
T0.Lg = Lg;
T0.u0 = nan;
T0.s  = nan;

data    = [];

% copy in genomic positions
data.x  = x;
data.T0 = T0;
data.Lg = Lg;
data.dL = dL;

% use the salt and pepper filter?
data.salt_and_pep_flag = true;
data.s_p_cut = s_p_cut;

% copy in conversion fraction
data.D = D;

% create the coarse grained data.
[data] = decimateData3D( data );

data.max_peak_num = 10; 
data.datap = datap;
end