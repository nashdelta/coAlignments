

clear
clc

vDwnld = readtable('ncbiVirusDownload.csv');
load('lineageInfo.mat')
load('binaryInteract')
load('speciesTax')

%make binary interaction matrix unique (up to rows where order of columns
%doesn't matter).

bis = cell(length(binaryInteract{2,1}),1);
for i = 1:length(bis)
    s1 = [num2str(binaryInteract{2,1}(i)),'.',binaryInteract{2,2}{i}];
    s2 = [num2str(binaryInteract{2,3}(i)),'.',binaryInteract{2,4}{i}];
    s = sort({s1,s2});
    bis{i} = [s{1},'.',s{2},num2str(binaryInteract{2,5}(i))];
end
[~,sBI] = unique(bis);
for i = 1:length(binaryInteract(1,:))
    binaryInteract{2,i} = binaryInteract{2,i}(sBI); %#ok<SAGROW>
end

%Get binary interaction matrix with species on the axis and values the
%number of binary interactions identified between the pair of species
%(diagonal elements are zero per previous pruning).
interactMat = zeros(length(speciesTax{2,1}));
for i = 1:length(binaryInteract{2,1})
    [~,pos1] = ismember(binaryInteract{2,1}(i),speciesTax{2,1});
    [~,pos2] = ismember(binaryInteract{2,3}(i),speciesTax{2,1});
    interactMat(pos1,pos2) = interactMat(pos1,pos2)+1;
end

%Get number of interactions between viruses and the following groups,
%bacteria, plants, fungi, metazoa(ex vert.), vertebrates (ex. human),
%humans, else, and other viruses.

%Note "else" here is protists, Leishmania major 5664, and Trichomonas vaginalis 5722 are
%present in the database with 38 and 131 protein-protein interactions respectively

%note originally sought to look at all pairs, but turns out in this
%database besides self-pairs (within the same organism) only viral and
%other organism pairs are listed. Not sure if this is because I followed
%the link to the string database through the string/virus link but if
%that's true not sure why the self-pairs are there. In any case, this is
%what we want.

%no archaea in this data set as of 2/2020.
keyWordCol = nan(7,1);
[~,keyWordCol(1)] = ismember('Bacteria',speciesTax{2,2});
[~,keyWordCol(2)] = ismember('Viridiplantae',speciesTax{2,2});
[~,keyWordCol(3)] = ismember('Fungi',speciesTax{2,2});
[~,keyWordCol(4)] = ismember('Metazoa',speciesTax{2,2});
[~,keyWordCol(5)] = ismember('Vertebrata',speciesTax{2,2});
[~,keyWordCol(6)] = ismember('Homo sapiens',speciesTax{2,2});
[~,keyWordCol(7)] = ismember('Viruses',speciesTax{2,2});

%specify category names
hostCatNames = {'Bacteria','Protists','Plants','Fungi','Invertebrates','Vertebrates','Human'};

specKing = speciesTax{2,3}(:,keyWordCol);
specKing(:,4) = specKing(:,4)-specKing(:,5);
specKing(:,5) = specKing(:,5)-specKing(:,6);

%First get Virus-Virus data
vVMat = interactMat(specKing(:,7),:);
vVMat = vVMat(:,specKing(:,7));
vVMat = vVMat+vVMat';
for i = 2:length(vVMat(:,1))
    vVMat(i,1:i) = 0;
end

vSP = sum(logical(vVMat(:)));
vPP = sum(vVMat(:));
vN = sum(specKing(:,7));

%Now build hostVirus matrix
hostVirus = interactMat+interactMat';
hostVirus = hostVirus(~specKing(:,7),:);
hostVirus = hostVirus(:,specKing(:,7));
catID = sum(specKing.*repmat(linspace(1,length(specKing(1,:)),length(specKing(1,:))),length(specKing(:,1)),1),2);
catID = catID(lt(catID,7));

%correct for protist ordering
catID = catID+1;
pro = catID==1;
catID(catID==2) = 1;
catID(pro) = 2;

%get host-virus numbers
hN = nan(max(catID),1);
hSP = nan(max(catID),1);
hPP = nan(max(catID),1);
for i = 1:max(catID)
    currMat = hostVirus(catID==i,:);
    hSP(i) = sum(logical(currMat(:)));
    hPP(i) = sum(currMat(:));
    hN(i) = sum(catID==i);
end

%Output number of species in each group, number of host-virus species pairs
%in each group, and total number of binary interactions with viruses in
%each group.
saveMat1 = [[hN;vN]';[hSP;vSP]';[hPP;vPP]'];

%Construct bipartite graph with edges connecting host species and viral
%species

fEmpty = 0.7;
hc = 2;
w1 = 1;
w2 = 2;
cThresh = 10;
maxR = 0.9;
rN = 0;
titleVar = 'Number of Protein-Protein Interactions';
virusName = 'Viruses';
        
plotBipartite1(fEmpty,catID,virusName,hostCatNames,hostVirus,hc,w1,w2,cThresh,maxR,rN,titleVar)

rN = 1;%remove Human (accounts for about half of all protein-protein interactions)
plotBipartite1(fEmpty,catID,virusName,hostCatNames,hostVirus,hc,w1,w2,cThresh,maxR,rN,titleVar)

%Move on to process NCBI virus data:

%Host key taxa - note that viruses do not appear as "hosts" in this list.
keyID = [2,2157,33090,4751,33208,7742];%bacteria, archaea, plants, fungi, metazoa, vertebrates

[uViralGenus,~,gID] = unique(vDwnld.Genus);
viralHost = vDwnld.Host;
viralSpec = vDwnld.Species;

genusHostCount = nan(length(uViralGenus)-1,6);%ignore first blank genus
uSpecHost = nan(length(uViralGenus)-1,1);
uSpeciesHostPairs = nan(length(uViralGenus)-1,6);
numSpec = nan(length(uViralGenus)-1,6);
for i = 1:length(genusHostCount)
    gHost = viralHost(gID==(i+1));
    gSpec = viralSpec(gID==(i+1));
    [uHost,~,iHost] = unique(gHost);
    [uSpec,~,iSpec] = unique(gSpec);
    
    specHost = zeros(length(uSpec),length(uHost));
    for j = 1:length(gHost)
        specHost(iSpec(j),iHost(j)) = specHost(iSpec(j),iHost(j)) + 1;
    end
    specHost = logical(specHost);
    
    oneHost = specHost(sum(specHost,2)==1,:);
    oneSpec = oneHost(:,sum(oneHost,1)==1);
    
    uSpecHost(i) = sum(oneSpec(:));
    
    numSpec(i) = length(uSpec);
    
    uGenusHost = unique(viralHost(gID==(i+1)));
    for j = 1:length(uGenusHost)
        uGenusHost{j} = uGenusHost{j}(~(uGenusHost{j}==''''));
    end
    [~,gI] = ismember(uGenusHost,lineageInfo{1});
    uLin = lineageInfo{2}(gI);
    linMat = nan(length(uLin),6);
    for j = 1:length(linMat(:,1))
        linMat(j,:) = ismember(keyID,uLin{j});
    end
    genusHostCount(i,:) = sum(linMat,1);
    if isequal(mod(i,100),0)
        disp(i)
    end
end

genusHostCount(:,5) = genusHostCount(:,5)-genusHostCount(:,6);%remove vertebrates from metazoa
[~,ord] = sort(sum(genusHostCount,2),'descend');
genusHostCount = genusHostCount(ord,:);
uViralGenus = uViralGenus(2:end);
uViralGenus = uViralGenus(ord);
numSpec = numSpec(ord);
uSpecHost = uSpecHost(ord);

countBelow = nan(max(genusHostCount(:)),6);
countBelowSpec = nan(max(numSpec),1);
countBelowU = nan(max(uSpecHost),1);
for i = 1:max(genusHostCount(:))
    countBelow(i,:) = sum(ge(genusHostCount,i),1);
    countBelowSpec(i) = sum(ge(numSpec,i));
    countBelowU(i) = sum(ge(uSpecHost,i));
end

figure
loglog(countBelow)
xlabel('No. Unique Host Tax ID''s')
ylabel('No. Virus Genera w/ # or More Hosts')
legend({'Bacteria','Archaea','Plants','Fungi','Metazoa ex. Vert.','Vertebrates'})

figure
loglog(countBelowSpec,'k')
hold on
plot(countBelowU,'r')
xlabel('No. within Genus')
ylabel('No. Virus Genera w/ # or More')
legend({'No. Species','No. Unique Spec/Host Pairs'})



    