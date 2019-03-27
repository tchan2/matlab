% Tammy Chan
% Bioinformatics
% February 7th, 2018

clear all; % clear all memory
close all; % close all windows

% load data
load EColi.mat
seq = EColi;

% alphabet
ATCG = {'A', 'T', 'C', 'G'};

% use array p[]
for i = 1:4
    nt_count(i) = length(find(seq == ATCG{i}));
    p(i) = nt_count(i)/length(seq);
end

% calculate entropy
entropy = sum(-p.*log2(p));
    
% probabilities of ATCG
fprintf ('\n The entropy of this sequence is: \n');
disp (entropy);

% codon bias
cb = codonbias (EColi, 'PIE', true)




