% Tammy Chan
% Bioinformatics
% Central Dogma

clear all; % clear all memory
close all; % close all windows

% random generation of dna sequence
chr = randseq(77);
nt_seq = erase (chr, match);

fprintf ('\ndna sequence = \n');
fprintf ('\n');
disp(['     ', nt_seq]);
fprintf ('\n');

% remove stop codons in DNA sequence
match = 'TGA', 'TAG', 'TAA';
nt_seq = erase (chr, match);

% print new sequence
fprintf ('\nnew dna sequence = \n');
fprintf ('\n');
disp(['     ', rnt_seq]);
fprintf ('\n');

% change nucleotide seq to amino seq
%amino_seq = nt2aa(nt_seq);

% output
fprintf ('\nfinal amino sequence = \n');
fprintf ('\n');
disp(['     ', ]);
fprintf ('\n');
