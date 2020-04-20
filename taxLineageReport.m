

clear
clc

lineageData = tdfread('tax_report.txt');
name = lineageData.name;
lin = lineageData.lineage;
nameCell = cell(length(name(:,1)),1);
linCell = nameCell;
for i = 1:length(nameCell)
    selector = ~(name(i,:)=='''');
    vect = name(i,selector);
    nameCell{i} = vect(1:find(~(vect==' '),1,'last'));
    linCell{i} = str2double(strsplit(lin(i,:),' '));
    linCell{i} = linCell{i}(~isnan(linCell{i}));
end

[nameCell,selector] = unique(nameCell);
linCell = linCell(selector);

lineageInfo = {nameCell,linCell};
save('lineageInfo','lineageInfo')




    