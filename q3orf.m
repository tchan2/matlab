clear all; % clear all memory
close all; % close all windows

% random generation of RNA sequence
rna_seq = randseq (24, 'alphabet', 'rna');

% start and stop codons
start_codon = 'AUG';
end_codon = 'UAG';

% output
fprintf ('\nopen reading frame = \n');
fprintf ('\n');
disp(['     ', start_codon, rna_seq, end_codon]);
fprintf ('\n');

% save file
filename='q3orf.mat';
save (filename, 'rna_seq', 'start_codon', 'end_codon');
fprintf ('\n %s is saved!\n', filename);