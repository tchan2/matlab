% Tammy Chan
% Bioinformatics
% Due: 3/23/2018

clear all;
close all;

SeqX = 'CATTG';
SeqY = 'CACG'; 
x = length(SeqX) + 1;
y = length(SeqY) +1;

% create matrix filled with zeros (preallocated memory)
a = zeros(x, y); 

% set values 
a(1, :) = [ 0 -1 -2 -3 -4 ];
a(: , 1) = [0 -1 -2 -3 -4 -5];

for i = 2 : x
    for j  = 2 : y
        % when matched
        if SeqX(i-1) == SeqY(j-1)
            score = 2; 
        else
        % when mismatched
            score = 0;
        end
        
        % calculating scores to only retrieve the maximum score
        score = [a(i-1, j-1) + score, a(i-1, j) + (-1), a(i, j-1) + (-1)];
        a(i,j) = max(score);
    end
end

% flipping the matrix the elements are reversed
a = flipud(a)
a = fliplr(a)

NewSeqY = [];

% going through loop to check each element
    for j = 1 : y
        % checking the elements diagonal, below or next
        % if diagonal. then add the next letter of SeqY
      if x>1 && y>1 && a(x-1, y-1) >= a(x,y-1) && a(x-1,y)
           NewSeqY = strcat(SeqY(y-j), NewSeqY)
           % if the next largest element is in the same column, 
           % place a '-'
           if x>1 && y>1 && a(x-1, y-1) == a(x-1, y)
               NewSeqY = strcat('-', NewSeqY);
           end
      end
    end
      
      disp(SeqX);
      disp(NewSeqY);
