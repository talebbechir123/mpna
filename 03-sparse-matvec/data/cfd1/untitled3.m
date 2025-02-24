

S = load('cfd1.mat');
A = S.Problem.A;  % A is now in your workspace
whos A
[V, D] = eigs(A, 1, 'largestabs');  % For real or complex matrices
% If A is symmetric, you could also do:
% [V, D] = eigs(A, 1, 'largestabs', 'IsSymmetric', true);

lambda = D(1,1);   % The eigenvalue
v = V(:,1);        % Corresponding eigenvector

fprintf('Largest eigenvalue by absolute value: %g\n', lambda);
disp('Eigenvector:');
disp(v);
