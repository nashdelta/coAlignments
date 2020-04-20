

%This script finds all pairs of protein-protein interactions in the string
%database formatted as protein1 protein2 score. protein1 and protein2 are
%of the form TAXID.PROTEINID

%Be careful with textscan and MATLAB inferring scientific notation. For
%example, with a delimiter of '.' 123.E2 will not be read as 123, but 12300
%so instead of %d, the first field must be %s. Note that 123.D2 is the same
%not sure why 'D' returns exponent.

clear
clc

blockN = 10^6;
origID = fopen('STRING_p1_p2_score.txt', 'rt');
header = fgetl(origID);

c1 = {};
c2 = {};
s = [];

stp = 0;
count = 0;
tic
while ~stp  
    blockScan = textscan(origID,'%s %s %d',[blockN,3]);
    if ~isempty(blockScan{1})
        count = count + 1;
    
        writecell(blockScan{1},'tmp1.txt')
        writecell(blockScan{2},'tmp2.txt')
        
        fID1 = fopen('tmp1.txt');
        fID2 = fopen('tmp2.txt');
        
        blockScan1 = textscan(fID1,'%s %*[^\n]',[blockN,1],'Delimiter','.');
        blockScan2 = textscan(fID2,'%s %*[^\n]',[blockN,1],'Delimiter','.');
        
        fclose(fID1);
        fclose(fID2);
        
        bS1 = str2double(string(blockScan1{1}));
        bS2 = str2double(string(blockScan2{1}));
        
        nSame = ~(bS1==bS2);
        
        c1 = cat(1,c1,blockScan{1}(nSame));
        c2 = cat(1,c2,blockScan{2}(nSame));
        s = cat(1,s,blockScan{3}(nSame));
        
        toc
        disp(count*blockN)
    else    
        stp = 1;
    end
end

fclose(origID);

T1 = nan(length(c1),1);
T2 = nan(length(c1),1);
G1 = cell(length(c1),1);
G2 = cell(length(c1),1);
for i = 1:length(c1)
    f1 = find(c1{i}=='.',1,'first');
    f2 = find(c2{i}=='.',1,'first');
    T1(i) = str2double(c1{i}(1:f1-1));
    T2(i) = str2double(c2{i}(1:f2-1));
    G1{i} = c1{i}(f1+1:end);
    G2{i} = c2{i}(f2+1:end);
end
    
binaryInteract = cell(2,5);
binaryInteract{1,1} = 'taxID_1';
binaryInteract{1,2} = 'geneID_1';
binaryInteract{1,3} = 'taxID_2';
binaryInteract{1,4} = 'geneID_2';
binaryInteract{1,5} = 'score';
binaryInteract{2,1} = T1;
binaryInteract{2,2} = G1;
binaryInteract{2,3} = T2;
binaryInteract{2,4} = G2;
binaryInteract{2,5} = s;

allTaxID = unique([T1;T2]);
writematrix(allTaxID,'allTaxID.txt');

save('binaryInteract','binaryInteract')
