% function
function[orf] = createORF(length)
    if mod(length, 3) ~= 0
        error (' should be divisible by three!\n')
    end
    
    rna_seq = randseq (length-6, 'alphabet', 'rna');
    
    % start and stop codons
    start_codon = 'AUG';
    stop_codon = 'UAG'; 
    
    orf = [start_codon, rna_seq, stop_codon]
    fprintf ('ORF with length %d has been generated. \n', length);

end
