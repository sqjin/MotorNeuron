function [numCluster, numCluster0, eigenvalues] = determineNumClusters(CM,fileName)
% Estimation of the number of clusters from the consensus matrix
%
% Input:
%   CM: consensus matrix constructed by running Seurat using multiple resolutions
%   fileName: a char giving the file name of the figure to save
% Output:
%   numCluster: Number of inferred clusters
%   numCluster0: the minimum number of inferred clusters based on the number of zero eigenvalues
%   eigenvalues: eigenvalues of the graph Laplacian
%
n = size(CM,1);
CM(1:n+1:end) = 1;

tol = 0.01;
numEigs = min(100,size(CM,1));

D = diag(CM*ones(n,1));
Prw = eye(n) - D^(-1/2)*CM*D^(-1/2);
all_eigs = real(eigs(Prw,numEigs,'sm'));
ZZ = sort(abs(real(all_eigs)));
numCluster0 = length(find(ZZ<=tol));
tau = 0.3;
CM(CM <= tau) = 0;
CM = (1/2)*(CM + CM');

D = diag(CM*ones(n,1));
Prw = eye(n) - D^(-1/2)*CM*D^(-1/2);
all_eigs = real(eigs(Prw,numEigs,'sm'));

zz = sort(abs(real(all_eigs)));

gap = zz(2:end) - zz(1:end-1);
[~,numCluster] = max(gap);


numCluster0 = length(find(zz<=tol));
display('Number of cluster based on zero eigenvalues & Largest gap ');
display([numCluster0 numCluster]);

eigenvalues = zz;

figure('position', [600, 200, 250, 200])
scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'k','filled');
hold on
scatter(numCluster,eigenvalues(numCluster),40,'r')
box off;
set(gca,'LineWidth',1);
set(gca,'FontSize',8,'FontName','Arial');
xlabel('Number of clusters','FontSize',10);
ylabel('Eigenvalue of graph Laplacian','FontSize',10);
title(['Inferred number of clusters: ', num2str(numCluster),'; Min number: ',num2str(numCluster0)],'FontSize',10)
saveas(gcf,['eigGap_',fileName,'.pdf'])


