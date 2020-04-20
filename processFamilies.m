

clear
clc

% BLOCK 1

gName = 'plant';
gNs = 'Viridiplantae';

% gName = 'vertebrate';
% gNs = 'Vertebrata';

numFam = 1;%max number of families with which a "good" virus may interact
numInteract = 1;%number of interactions between host family and "good"
                %viruses must be greater than this. Could many host,
                %single virus, many virus single host.

fID = fopen([gName,'Keep.cls'],'r');
hostCLS = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);

clsCell = hostCLS{2};
for i = 1:length(clsCell)
    clsCell{i} = strsplit(clsCell{i}(1:end-1),' ')';%end-1 to remove trailing white space from .cls files
end

isE = 0;
if isempty(clsCell{end}{1})% if the "remain" field is empty
    clsCell = clsCell(1:end-1);
    isE = 1;
end

allHostCLS = vertcat(clsCell{:});

terry = clsCell(1:end-1);
terry = vertcat(terry{:});

load('binaryInteract')
load('speciesTax')

disp('loaded')

[~,isGroup] = ismember(gNs,speciesTax{2,2});
groupTax = speciesTax{2,1}(speciesTax{2,3}(:,isGroup));

[~,isVirus] = ismember('Viruses',speciesTax{2,2});
virusTax = speciesTax{2,1}(speciesTax{2,3}(:,isVirus));

hostVirus = logical(ismember(binaryInteract{2,1},groupTax).*ismember(binaryInteract{2,3},virusTax));
virusHost = logical(ismember(binaryInteract{2,1},virusTax).*ismember(binaryInteract{2,3},groupTax));

groupInteract = binaryInteract(:,1:4);
groupInteract{2,1} = cat(1,binaryInteract{2,1}(hostVirus),binaryInteract{2,3}(virusHost));
groupInteract{2,2} = cat(1,binaryInteract{2,2}(hostVirus),binaryInteract{2,4}(virusHost));
groupInteract{2,3} = cat(1,binaryInteract{2,3}(hostVirus),binaryInteract{2,1}(virusHost));
groupInteract{2,4} = cat(1,binaryInteract{2,4}(hostVirus),binaryInteract{2,2}(virusHost));
groupString = append(cellstr(num2str(groupInteract{2,1})),groupInteract{2,2},...
    cellstr(num2str(groupInteract{2,3})),groupInteract{2,4});
[~,uIndi] = unique(groupString);
for i = 1:4
    groupInteract{2,i} = groupInteract{2,i}(uIndi);
end
hostLabel = strtrim(append(cellstr(num2str(groupInteract{2,1})),'.',groupInteract{2,2}));%strtrim removes the leading whitespace from num2str
virusLabel = strtrim(append(cellstr(num2str(groupInteract{2,3})),'.',groupInteract{2,4}));

checkAll = min(sum(ismember(allHostCLS,hostLabel))/length(allHostCLS),...
    sum(ismember(hostLabel,allHostCLS))/length(hostLabel));
if lt(checkAll,1)
    disp('Something is Missing')
    keyboard
end

disp('Nothing Missing')

[hostLabel,~,hostID] = unique(hostLabel);
[virusLabel,~,virusID] = unique(virusLabel);
hostFamilyCell = cell(length(clsCell),1);
hostFID = nan(length(hostID),1);
disp('Processing Host Families')
L = length(hostFamilyCell);
mo1 = floor(L/10);
for i = 1:length(hostFamilyCell)
    [~,hostFamilyCell{i}] = ismember(clsCell{i},hostLabel);
    hostFID(ismember(hostID,hostFamilyCell{i})) = i;
    if isequal(mod(i,mo1),0)
        disp(L-i)
    end
end

pPI = zeros(length(hostFamilyCell),length(virusLabel));
for i = 1:length(pPI(:,1))
    pPI(i,virusID(ismember(hostID,hostFamilyCell{i}))) = pPI(i,virusID(ismember(hostID,hostFamilyCell{i})))+1;
end
if ~isE
    pPI = pPI(1:end-1,:);%ignore singletons for now
end
disp('Done')

vFamCount = sum(pPI);

vKeepIndex = find(le(vFamCount,numFam));
hKeepIndex = find(gt(sum(pPI(:,vKeepIndex),2),numInteract));

keepID = logical(ismember(hostFID,hKeepIndex).*ismember(virusID,vKeepIndex));
hostKeepID = hostID(keepID);
hostKFID = hostFID(keepID);
virusKeepID = virusID(keepID);
[uHK,uID] = unique(hostKeepID);
hCat = hostKFID(uID);
uVK = unique(virusKeepID);
sPPI = zeros(length(uHK),length(uVK));
for i = 1:length(sPPI(:,1))
    sPPI(i,ismember(uVK,virusKeepID(hostKeepID==uHK(i)))) = 1;
end

disp('plotting')

fEmpty = 0.5;
[~,catID] = ismember(hCat,unique(hCat));
virusName = 'virus';
hostCatNames = cellstr(num2str(unique(catID)));
iMat = sPPI;
hc = 3;
w1 = 1;
w2 = 1;
cThresh = 0;
maxR = 0.9;
rN = 0;
titleVar = gName;
numCol = 2;

plotBipartite1(fEmpty,catID,virusName,hostCatNames,iMat,hc,w1,w2,cThresh,maxR,rN,titleVar,numCol)

[~,gFID] = ismember(hostKFID,unique(hostKFID));
goodRecord = [gFID,hostKeepID,virusKeepID];
[~,ord] = sort(goodRecord(:,2),'descend');
goodRecord = goodRecord(ord,:);
[~,ord] = sort(goodRecord(:,3),'descend');
goodRecord = goodRecord(ord,:);
[~,ord] = sort(goodRecord(:,1),'ascend');
goodRecord = goodRecord(ord,:);

goodCell = cell(size(goodRecord));
goodCell(:,1) = cellstr(num2str(goodRecord(:,1)));
goodCell(:,2) = hostLabel(goodRecord(:,2));
goodCell(:,3) = virusLabel(goodRecord(:,3));

writecell(goodCell,[gName,'ReviewInteractions.txt'],'delimiter','tab')







