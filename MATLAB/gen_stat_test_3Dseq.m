%% Compupte P values from 3D-seq Data 
%
% The purpose of this script is to generate a statistical test in which p
% values, ie the probability that the fit values are more extreme than the
% observed values assume the null hypothesis. 
%
% In the current context, the distribution of the log likelihood ratio must 
% be computed explicitly since the peak model is statistically singular 
% (i.e. the Wilks Theorem does NOT apply.   
%

%% Where are we?

% This is where I am now
disp( pwd );

% What file is running:
disp( [mfilename('fullpath'),'.m'] );

%% Clear and close everything
clear all;
close all;
colordef white

%% Load in data from the xls file.
% We need to get the genomic positions into vector x and the conversion
% fraction at position x into D

data_ = load( '3Dseq_data.mat' );
D = data_.D;
x = data_.x;

%% Generate a null hypothesis from GacA data
Nbin = 100;

[y2,x2] = hist( D,Nbin );

figure(2);
clf;

semilogy( x2, y2/sum(y2), '.-' );

legend( seq_data.name(2)' );

ylabel( 'Probability' );
xlabel( 'Conversion fraction: R');

%% Create data structure
% We need to generate an emprical model for the null hypothesis. We will
% fit data from a control experiment
% We need to make decimated (reduced resolution) copy of the data. 

make_new_flag = true;

% Salt_and_Pepper filter will remove the isolated SNP's. 
cut_off_std = 1;
data_file_name = 'null_hyp_data.mat';

if make_new_flag   
    data = makeData3D( x, D, cut_off_std );
     
    % save this struture datafile so this step doesn't need to be repeated.
    save( data_file_name, '-struct', 'data' );
else
    % Load in data file if it isn't regenerated.
    data = load( data_file_name );
end


%% Set up ranger

make_new_flag = true;

num_dec = numel(data.Ddec);

if make_new_flag
    ranger = 1:num_dec;
else
    ranger = [1:2340,2360:16225,16235:num_dec];
end

%% Show the data we need to fit the distribution from
% This is the distribution we are attempting to model with an 
% empirical model.

figure(1);
clf;
plot(data.xdec,data.Ddec,'b.-');
hold on;
plot(data.xdec(ranger),data.Ddec(ranger),'r.-');

legend( {'All data', 'Included data' } );

xlim( [1,data.xdec(end)] );
xlabel( 'Genomic position' );
ylabel( 'Conversion fraction' );

% Generate histogram from background values

figure(2);
clf;

[yy_all,xx] = hist( data.Ddec(:),200 );
[yy_inc,xx] = hist( data.Ddec(ranger),xx );

semilogy( xx, yy_all/sum(yy_all), '.-' );
hold on;
semilogy( xx, yy_inc/sum(yy_all), '.-' );

ylabel( 'Probability' );
xlabel( 'Conversion Probability' );

legend( {'All data', 'Included data' } );


%% Run this to remove peaks 
% click on a peak you wish to remove on the genome plot
% and then run above section again.

figure(1);
x_rem = ginput( 1 );
x_rem = x_rem(1);

% remove indices within 50 kb of choice
ind = abs(x_rem-data.xdec(ranger))>5e4;

ranger = ranger( ind );


%% Generate histogram from background values

[yy,xx] = hist( data.Ddec(ranger),100 );
figure(2);
clf;

semilogy( xx, yy/sum(yy), '.-' );
hold on;

ylabel( 'Probability' );
xlabel( 'Conversion Probability' );


%% Make the empirical model
% Fit the empircal model and then compare it with the null hypothesis
% histogram.

null_model = makeEmpModel3D( data.Ddec(ranger) );

x = (0:.01:1)*1.0e-3;

mm = empProbFun( x, null_model );
mm(x>0) = mm(x>0)*(xx(2)-xx(1));

semilogy( x, mm, 'r.-' );

%% Generate simulated data.
% Compare the datad, simulated data, and model histograms
% As we can see, all the distributions match.

N = numel(data.Ddec(ranger));

FF = genF2(N, null_model );

[yyFF] = hist( FF,xx );

semilogy( xx, yyFF/sum(yyFF), 'g.-' );


legend( {'Data','Emprical model','Simulation'} ); 



%% Build the test statistic
% This will take a couple minutes

xdec = data.xdec;
T0   = data.T0;
T0.a0 = 1.5;
T0.D0 = 1;

% set the size of the coarse-grained genome
xx   = xdec;
Ngenome = numel( xx );

% set the number of simulated experiments
%Nloop = 10000; % 30 min
Nloop = 1000;  %  3 min

% 
Ddec_ = genF2(Ngenome, null_model );
[That, mll_min] = makeMLE_( Ddec_, xdec, T0);


T0.u0 = mean(Ddec_);
T0.s = std(Ddec_);

% test statistic vector
lam = nan( [1, Nloop] );
tic

%for ii = 1:Nloop
parfor ii = 1:Nloop
    
    disp( ['Loop fraction: ', num2str(ii/Nloop)] );
    
    Ddec_ = genF2(Ngenome, null_model );
    
    [That, mll_min,~,mll_0] = makeMLE_( Ddec_, xdec, T0);

    lam(ii) = mll_0-mll_min;

end
toc

%% Generate the test statistic CDF and plot it

lam = sort( lam );
p = 1-(0:(Nloop-1))/Nloop;

figure(11);
clf;
loglog( lam, p, '.-' );
hold on;

ylabel( 'P value' );
xlabel( 'Test statistic: \lambda' );

p_value_model = doGumbelFit2( lam );

u = p_value_model.u;
s = p_value_model.s;

loglog( lam, 1-exp(-exp( -(lam-u)/s)), 'r.-' );

legend( {['Sim. N = ', num2str(Nloop)], 'Gumbel model' } );
doPageFormat( [5,3] );
print -dpdf pval.pdf

% Save p value model

datap = [];

datap.lam = lam;
datap.p = p;
datap.ull = p_value_model.u;
datap.sll = p_value_model.s;

save( 'datap_emp.mat','-struct','datap');









