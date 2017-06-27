%-----copyRight(c) Shuhua Li<sue142857@gmail.com> 04.04.2017-----%
% function: 

function [outChain,pieces] = get_outChain(C,T,direction)

n = length(T);
m = size(C,1);
pieces = get_pieces(T,m,direction);
outChain = cell(n,1);
for i = 1:n
    outChain{i} = C(pieces{i},:);
end


