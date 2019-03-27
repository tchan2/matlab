% Tammy Chan
% Bioinformatics
% April 25th, 2018

clear all;
close all;

load glioma_normalized_mirnaseq_data_.mat;
whos

% obtain all the singular values for dataset
X = data;
[U, S, V] = svd(X);
singular_values = diag(S);

% rank 2 approximation
rank_2_approximation = U(:, 1:2) * S(1:2, 1:2) * V(:, 1:2)';

% find top 20 ranked miRNA using SVD
u1 = U(:,1);
[sorted_u1, idx] = sort(-u1, 'descend');
gene_number = 20;

fprintf('\nThe gene contribution scores of top %d genes \n', gene_number);
selected_gene_contribution_scores = sorted_u1(1:gene_number);

disp(selected_gene_contribution_scores);

% plot contribution scores
figure ('name', 'Selected Gene Contribution Scores');
plot(selected_gene_contribution_scores);
grid on;
xlabel('Genes');
ylabel('Gene Contribution Scores');

% employ the two-sample t-test to find top 30 ranked microRNA
x = data(:, 1);
y = data(:, 2);
[h, p] = ttest2 (x, y, 'Alpha', 0.3)
