

clear
clc

vDwnld = readtable('ncbiVirusDownload.csv');

nCell = cell(1,width(vDwnld));
numMat = nan(size(vDwnld));
for i = 1:width(vDwnld(1,:))
    [nCell{1,i},~,numMat(:,i)] = unique(vDwnld(:,i));
    disp(width(vDwnld(1,:))-i)
end
[uNum,~,numIndex] = unique(numMat,'rows');

ncbiVirusDat = {nCell,uNum,numIndex};

save('ncbiVirusDat','ncbiVirusDat')

% clear
% clc
% 
% %%%Decompress
% 
% load('ncbiVirusDat.mat')
% 
% numMat = ncbiVirusDat{2}(ncbiVirusDat{3},:);
% datCell = cell(size(numMat));
% for i = 1:length(numMat(1,:))
%     tmp_i = table2cell(ncbiVirusDat{1}{i});
%     datCell(:,i) = tmp_i(numMat(:,i));
%     disp(length(numMat(1,:))-i)
% end
% 
% vDwnld2 = cell2table(datCell,'VariableNames',{'Species','Genus','Family','Protein','Host'});
% 
% vDwnld = readtable('ncbiVirusDownload.csv');
% 
% isequal(vDwnld,vDwnld2)


