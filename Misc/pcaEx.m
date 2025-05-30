rng 'default'
close all

M = 7; % Number of observations
N = 2; % Number of variables observed

X = rand(M,N);

% De-mean
X = bsxfun(@minus,X,mean(X));

plot(X(:,1),X(:,2),'r*')

% Do the PCA
[coeff,score,latent] = pca(X);

hold on
plot(score(:,1),score(:,2),'bo')

% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix = cov(X);
[V,D] = eig(covarianceMatrix);

% "coeff" are the principal component vectors. These are the eigenvectors of the covariance matrix. Compare ...
coeff
V

% Multiply the original data by the principal component vectors to get the projections of the original data on the
% principal component vector space. This is also the output "score". Compare ...

dataInPrincipalComponentSpace = X*coeff
score

% The columns of X*coeff are orthogonal to each other. This is shown with ...

corrcoef(dataInPrincipalComponentSpace)

% The variances of these vectors are the eigenvalues of the covariance matrix, and are also the output "latent". Compare
% these three outputs

var(dataInPrincipalComponentSpace)'
latent
sort(diag(D),'descend')