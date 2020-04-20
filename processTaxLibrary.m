

clear
clc

fid = fopen('taxLibrary.txt');
tS = textscan(fid,'%d %s','Delimiter','\t');

taxID = tS{1};
fullTax = tS{2};

lTax = nan(length(fullTax),1);
taxCell = cell(length(fullTax),1);
for i = 1:length(fullTax)
    taxCell{i} = strsplit(fullTax{i},';');
    lTax(i) = length(taxCell{i});
end

allNames = unique(cat(2,taxCell{:}));

taxMat = false(length(taxCell),length(allNames));
for i = 1:length(taxCell)
    taxMat(i,:) = ismember(allNames,taxCell{i});
end

speciesTax = {'taxID','allNames','taxMat';taxID,allNames,taxMat};

save('speciesTax','speciesTax')
