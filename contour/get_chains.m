%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 04.04.2017-----%
% function: get chain1, chain2 based on 2 contours and 2 conjugate
% polygons

function [outChain1,inChain1,outChain2,inChain2] = get_chains(C,T,scale)

% debug = 1;

%% 1. get out-chain 
% n = length(T(1,:)); 
% m = [size(C{1},1) size(C{2},1)];
% to do in get_outChain.m: get_pieces.m 
direction = 'counter-clockwise';
[outChain1,~] = get_outChain(C{1},T(1,:),direction);

direction = 'clockwise';
[outChain2,~] = get_outChain(C{2},T(2,:),direction); 

%% 2. get in-chain
vs_T1 = C{1}(T(1,:),:); 
inChain1 = get_inChain(vs_T1,outChain2,scale);  
vs_T2 = C{2}(T(2,:),:);
inChain2 = get_inChain(vs_T2,outChain1,scale); 


