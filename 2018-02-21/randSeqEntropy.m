% Tammy Chan
% Bioinformatics

clear all; % clear all memory
close all; % close all windows

% define variables
nt = 2001;
n = 200;

% generate 200 random sequences
for k = 1 : n
    rseq{k} = randseq(nt);
end

% alphabet
ATCG = {'A', 'T', 'C', 'G'};

% use vector p[] 
for k = 1 : n
    seq = rseq{k};
    for i = 1 : 4
        nt_count(i) = length(find(seq == ATCG{i}));
        p(i) = nt_count(i)/length(seq);
    end
    entropy(k) = sum(-p.*log2(p));
end
 
% display individual entropy values
fprintf ('\n All entropy values: \n', n);
disp (entropy);

% display average entropy value
avg_entropy = mean(entropy);
fprintf ('\n Average Entropy: %f\n', avg_entropy);
fprintf ('\n');



