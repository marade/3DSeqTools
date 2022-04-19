%% Analyze 3D data 
% MATLAB script for modeling peaks in 3D-seq conversion frequency data
% Paul Wiggins (University of Washington)
%

%% Set things up.

% This is where I am now
disp( pwd );

% What file is running:
disp( [mfilename('fullpath'),'.m'] );

%%
clear all;
close all;
colordef white

%% Load in data from the xls file.
% We need to get the genomic positions into vector x and the conversion
% fraction at position x into D

data_ = load( '3Dseq_data.mat' );
D = data_.D;
x = data_.x;

%% Construct data 
% Get rid of SNP bands
% run data decimation to make a reduced resolution dataset.

data = makeData3D( x, D );


%% Fit peaks recursively  

data =  makeAnalysis3D( data );

%% Visual data fit
% Show plots of the entire genome
% as well as individual peaks, and 
% generating excel files of data and models

data = showAnalysis3D( data );


