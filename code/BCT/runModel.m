clear all
close all
clc
%% quick explanation
% the algorithm operates in multiple stages and works as follows:
%
% initialize with randomly-sampled parameters, generate networks, and
% calculate energies.
% for nlvls - 1
%   1. divide parameter space based on voronoi tesselation
%   2. assign energy to each parcel
%   3. sample new set of parameters, preferentially sampling from within
%      lower-energy parcels (the parameter pow controls the preference).
%   4. generate new networks, calculate energy.
% end
%
% in effect, the algorithm samples more densely from low-energy areas while
% avoiding high-energy areas.
%% load sample data
load('sampleData.mat');
%% set parameters
% set limits for eta (geometric) and gamma (topological) parameters
eta = [-15,1];
gam = [-2,2];
% parameters related to the optimization
pow = 2;        % severity
nlvls = 5;      % number of steps
nreps = 250;    % number of repetitions/samples per step
% specify model type and whether power-law or exponential (we really only
% used power-law).
modeltype = 'neighbors';
modelvar = [{'powerlaw'},{'powerlaw'}];
% number of edges
m = nnz(A)/2;
%% sample networks
[E,K,N,P] = fcn_sample_networks(A,Aseed,D,m,modeltype,modelvar,nreps,nlvls,eta,gam,pow);
%%
% E = energy
% K = KS stats
% N = list of edges (indices not in the seed matrix)
% P = parameters
scatter3(P(:,1),P(:,2),E,100,E,'filled');
xlabel('eta');
ylabel('gamma');
zlabel('energy');