% Tammy Chan
% Bioinformatics
% March 14th, 2018

clear all; % clear all
close all; % close all

Seq1 = 'TTATTCACCAAACGGGCAATTCTTTAAAA';
Seq2 = 'TTTTGCACTCGUCCCGGGGGGCCTGACAAAT';

[score, alignment] = nwalign (Seq1, Seq2, 'Alphabet', 'NT')

disp(score);
disp(alignment);
