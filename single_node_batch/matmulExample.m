function matmulExample(n, m)
% function matmulExample(n, m)
% Purpose:
% Computes and collects timings data for a matrix-matrix multiply
% n  (input) row size
% m  (input) column size

tic
spmd
a = codistributed.rand(n,m); % distribute a random matrix by column
b = codistributed.rand(n,m); % distribute a random matrix by column
c = a*b; % distributed matrix multiply
end
Tp = toc;
disp(['Parallel job walltime is ' num2str(Tp)])
end
