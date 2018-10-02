
clear all;
close all;

load('hexosaminidase.mat','humanHEXA')
sequence = humanHEXA.Sequence;

% Finding All Potential Forward Primers
N = length(sequence) % length of the target sequence
M = 20  % desired primer length
index = repmat((0:N-M)',1,M)+repmat(1:M,N-M+1,1);
fwdprimerlist = sequence(index);

% Finding All Potential Reverse Primers
comp_sequence = seqcomplement(sequence);
revprimerlist = seqreverse(comp_sequence(index));

for i = N-19:-1:1 % reverse order to pre-allocate structure
    fwdprimerprops(i) = oligoprop(fwdprimerlist(i,:));
end

comp_sequence = seqcomplement(sequence);
revprimerlist = seqreverse(comp_sequence(index));

for i = N-19:-1:1 % reverse order to preallocate structure
    revprimerprops(i) = oligoprop(revprimerlist(i,:));
end

% Convert to correct types for analyzing
fwdgc = [fwdprimerprops.GC]';
revgc = [revprimerprops.GC]';

fwdtm = cell2mat({fwdprimerprops.Tm}');
revtm = cell2mat({revprimerprops.Tm}');

fwddm  = cellfun('isempty',{fwdprimerprops.Dimers}');
revdm = cellfun('isempty',{revprimerprops.Dimers}');

fwdhp = cellfun('isempty',{fwdprimerprops.Hairpins}');
revhp = cellfun('isempty',{revprimerprops.Hairpins}');

fwd_clamp = (lower(fwdprimerlist(:,end)) == 'g') & (lower(fwdprimerlist(:,end)) == 'c');
rev_clamp = (lower(fwdprimerlist(:,end)) == 'g') & (lower(fwdprimerlist(:,end)) == 'c');

% Creating structures of all the properties
fwdprops = struct ('GC', fwdgc, 'MT', fwdtm, 'DM', fwddm, 'HP', fwdhp, 'CL', fwd_clamp);
revprops = struct ('GC', revgc, 'MT', revtm, 'DM', revdm, 'HP', revhp, 'CL', rev_clamp);

% Create a vector fwdPrimerData for which each element
% contains a structure with the following fields:
%     primerSeq = the sequence of the primer
%     primerProp = the properties of the primer
for i = 1:length(fwdprimerlist)
    fwdPrimerData(i) = struct('primerSeq', fwdprimerlist, 'primerProp', fwdprops);
end

% Create a vector revPrimerData for which each element
% contains a structure with the following fields:
%     primerSeq = the sequence of the primer
%     primerProp = the properties of the primer
for i = 1:length(revprimerlist)
    revPrimerData(i) = struct('primerSeq', revprimerlist, 'primerProp', revprops, 'deleteElement', false);
end

% Create a logical array to delete all elements that do not match
% conditions for GC content
fwd_gc = logical(((fwdPrimerData(i).primerProp.GC) > 45) & ((fwdPrimerData(i).primerProp.GC) < 55));
rev_gc = logical(((revPrimerData(i).primerProp.GC) > 45) & ((revPrimerData(i).primerProp.GC) < 55));

fwdPrimerData(i).primerProp.GC = fwdPrimerData(i).primerProp.GC(fwd_gc);
revPrimerData(i).primerProp.GC = revPrimerData(i).primerProp.GC(rev_gc);

% Check if data has been deleted
% disp (revPrimerData(i).primerProp(1))

% Create a logical array to delete all elements that do not match
% conditions for melting temperatures
fwd_mt = logical(((fwdPrimerData(i).primerProp.MT) > 50) & ((fwdPrimerData(i).primerProp.MT) < 60));
rev_mt = logical(((revPrimerData(i).primerProp.MT) > 50) & ((revPrimerData(i).primerProp.MT) < 60));

fwdPrimerData(i).primerProp(i).MT = fwdPrimerData(i).primerProp.MT(fwd_mt);
revPrimerData(i).primerProp(i).MT = revPrimerData(i).primerProp.MT(rev_mt);

% Check if data has been deleted
% disp (revPrimerData(i).primerProp(1))

disp (revPrimerData(i).primerProp.DM)

% Create a logical array to delete all elements that do not match
% conditions for self-dimerization and hairpin formation
fwd_dm = logical((fwdPrimerData(i).primerProp(i).DM) ~= 0)
rev_dm = logical((revPrimerData(i).primerProp(i).DM) ~= 0)

fwdPrimerData(i).primerProp(i).DM = fwdPrimerData(i).primerProp.DM(fwd_dm);
revPrimerData(i).primerProp(i).DM = revPrimerData(i).primerProp.DM(rev_dm);


%{
% Filtering Primers Based on GC Content
for i = 1: length(gc_primers)
    if gc_primers == 0
        gc_primers = [];
    end
end

disp (gc_primers)


% Filtering Primers Based on GC Content
for i = 1: length((revPrimerData(i).primerProp.GC))
    if ((revPrimerData(i).primerProp.GC) > 45)
        revPrimerData(i).deleteElement = true;
        revPrimerData(i) = [];
        i = i - 1;
    end
end
disp (revPrimerData(i).primerProp.GC);
revPrimerData = revPrimerData(revPrimerData.deleteElement(true));
%}
