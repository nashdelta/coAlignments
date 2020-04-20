

clear
clc

% BLOCK 1

gName = 'plant';
gNs = 'Viridiplantae';

% gName = 'vertebrate';
% gNs = 'Vertebrata';

fID = fopen([gName,'Keep.cls'],'r');
hostCLS = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);

fID = fopen('virus3Keep.cls','r');
virusCLS = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);

clsCell = hostCLS{2};
for i = 1:length(clsCell)
    clsCell{i} = strsplit(clsCell{i}(1:end-1),' ')';%end-1 to remove trailing white space from .cls files
end

clsCellV = virusCLS{2};
for i = 1:length(clsCellV)
    clsCellV{i} = strsplit(clsCellV{i}(1:end-1),' ')';%end-1 to remove trailing white space from .cls files
end

isE = 0;
if isempty(clsCell{end}{1})% if the "remain" field is empty
    clsCell = clsCell(1:end-1);
    isE = 1;
end

isEV = 0;
if isempty(clsCellV{end}{1})% if the "remain" field is empty
    clsCellV = clsCellV(1:end-1);
    isEV = 1;
end

allHostCLS = vertcat(clsCell{:});

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

disp('Nothing Missing From Host')

[hostLabel,~,hostID] = unique(hostLabel);
[virusLabel,~,virusID] = unique(virusLabel);
hostFamilyCell = cell(length(clsCell),1);
virusFamilyCell = cell(length(clsCellV),1);
hostFID = nan(length(hostID),1);
virusFID = nan(length(virusID),1);
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

disp('Processing Virus Families')
L = length(virusFamilyCell);
mo1 = floor(L/10);
for i = 1:length(virusFamilyCell)
    [~,virusFamilyCell{i}] = ismember(clsCellV{i},virusLabel);
    virusFID(ismember(virusID,virusFamilyCell{i})) = i;
    if isequal(mod(i,mo1),0)
        disp(L-i)
    end
end

pairID = [hostID,virusID];
pPI = zeros(length(hostFamilyCell),length(virusFamilyCell));
for i = 1:length(pPI(:,1))
    currPair = pairID(ismember(pairID(:,1),hostFamilyCell{i}),:);
    for j = 1:length(pPI(1,:))
        cP2 = currPair(ismember(currPair(:,2),virusFamilyCell{j}),:);
        if ~isempty(cP2)
            [~,~,c1] = unique(cP2(:,1));
            [~,~,c2] = unique(cP2(:,2));
            m1 = max(c1);
            m2 = max(c2);
            cVect = zeros(m1*m2,1);
            cVect((c2-1)*m1+c1) = 1;
            pPI(i,j) = rank(reshape(cVect,m1,m2));
        end
    end
end

if ~isE
    pPI = pPI(1:end-1,:);%ignore singletons for now
end

if ~isEV
    pPI = pPI(:,1:end-1);%ignore singletons for now
end

disp('Done')

r = linspace(1,length(pPI(:,1)),length(pPI(:,1)));
c = linspace(1,length(pPI(1,:)),length(pPI(1,:)));
[gc,gr] = meshgrid(c,r);
kCC = cell(0,5);
for i = 1:max(pPI(:))
    colM = sum(ge(pPI,i-1),1);
    m2 = repmat(colM,length(pPI(:,1)),1);
    kPair = logical((pPI==i).*(m2==1));
    gP = pPI(:,gc(kPair));
    sP = nan(max(pPI(:)),length(gP(1,:)));
    if ~isempty(sP)
        for j = 1:length(sP(:,1))
            sP(j,:) = sum(gP==j,1);
        end
        vCell = clsCellV(gc(kPair));
        hCell = clsCell(gr(kPair));
        for j = 1:length(vCell)
            vKeep = ismember(virusID,find(ismember(virusLabel,vCell{j})));
            hKeep = ismember(hostID,find(ismember(hostLabel,hCell{j})));
            k = logical(vKeep.*hKeep);
            addK = cell(sum(k),5);
            addK(:,1) = cellstr(num2str(i*ones(length(addK(:,1)),1)));
            addK(:,2) = cellstr(num2str(j*ones(length(addK(:,1)),1)));
            addK(:,3) = repmat({strjoin(string(sP(:,j)),'|')},sum(k),1);
            addK(:,5) = virusLabel(virusID(k));
            addK(:,4) = hostLabel(hostID(k));
            kCC = cat(1,kCC,addK);
        end
    end
end

writecell(kCC,[gName,'ReviewFamilyPairs.txt'],'delimiter','tab')

