

clear
clc

fID = fopen('virusAChost.txt','r');
vACH = textscan(fID,'%s %s','delimiter',',');
fclose(fID);

nCell = cell(1,2);
numMat = nan(length(vACH{1}),2);
for i = 1:2
    [nCell{i},~,numMat(:,i)] = unique(vACH{i});
end

ncbiVirusDat2 = {nCell,numMat};

save('ncbiVirusDat2','ncbiVirusDat2')

% clear
% clc
% 
% %%%Decompress
% 
% load('ncbiVirusDat2.mat')
% 
% numMat = ncbiVirusDat2{2};
% datCell = cell(size(numMat));
% for i = 1:length(numMat(1,:))
%     datCell(:,i) = ncbiVirusDat2{1}{i}(numMat(:,i));
% end
% 
% vACH2 = {datCell(:,1),datCell(:,2)};
% 
% fID = fopen('virusAChost.txt','r');
% vACH = textscan(fID,'%s %s','delimiter',',');
% fclose(fID);
% 
% isequal(vACH,vACH2)

