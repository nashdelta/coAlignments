

clear
clc

%BLOCK 1

% load('binaryInteract')
% 
% t1 = binaryInteract{2,1};
% t2 = binaryInteract{2,3};
% 
% g1 = binaryInteract{2,2};
% g2 = binaryInteract{2,4};
% 
% printCell = cell(length(t1),2);
% for i = 1:length(printCell)
%     printCell{i,1} = [num2str(t1(i)),'.',g1{i}];
%     printCell{i,2} = [num2str(t2(i)),'.',g2{i}];
% end
% 
% uPrintCell = unique(printCell(:));
% 
% writecell(uPrintCell,'subSTRINGnames.txt')

%BLOCK 2

% fID = fopen('subSTRINGnames.txt','r');
% fID2 = fopen('subSTRING.sr','r');
% 
% requested = textscan(fID,'%s');
% found = textscan(fID2,'%s %s');
% 
% fclose(fID);
% fclose(fID2);
% 
% remaining = requested{1}(~ismember(requested{1},found{1}));
% 
% writecell(remaining,'remainingTaxDotGene.txt');
% 
% fID = fopen('remainingTaxDotGene.txt','r');
% remaining = textscan(fID,'%s %s','delimiter','.');
% fclose(fID);
% 
% writecell(remaining{2},'remainingGene.txt');

%BLOCK 3

% clear
% clc
% 
% %first make sure browser is set so the uniprot column preferences are as
% %desired - needs "chain".
% 
% fID = fopen('uniprotUnmapped1.txt','r');
% searchIDs = textscan(fID,'%s');
% fclose(fID);
% 
% %Uniprot search url
% generalUrl = 'https://www.uniprot.org/uniprot/?query=';
% uniPrefix = 'id="checkbox_';
% 
% foundIDs = cell(length(searchIDs{1}),3);
% for i = 1:length(searchIDs{1})
%     queryUrl = [generalUrl,searchIDs{1}{i}];
%     commandString = strcat('curl',{' '},queryUrl);
%     [~,commandOutput] = system(commandString{1});
%     idRegEx = 'id="checkbox_[\S]+"';%[\S]+ means any number of non-white space characters
%     matchCell = regexp(commandOutput,idRegEx,'match');
%     if ~isempty(matchCell)
%         foundIDs{i,1} = matchCell{1}(length(uniPrefix)+1:end-1);
%         queryUrl2 = [generalUrl,foundIDs{i}];
%         
%         %open page referenced by found id to get chain information
%         options = weboptions;
%         pageContents = webread(queryUrl2,options);
%         proRegEx = [foundIDs{i},'[\S]+[\S]+',searchIDs{1}{i}];%at least two characters in between since foundIDs{i}#searchIDs{1}{i} appears also
%         matchCell = regexp(pageContents,proRegEx,'match');
%         openBracket = find(matchCell{1}=='[');
%         closedBracket = find(matchCell{1}==']');
%         sequenceRange = strsplit(matchCell{1}(openBracket+1:closedBracket-1),'-');
%         foundIDs{i,2} = sequenceRange{1};
%         foundIDs{i,3} = sequenceRange{2};
%     else
%         foundIDs{i,1} = 'NaN';
%         foundIDs{i,2} = 'NaN';
%         foundIDs{i,3} = 'NaN';
%     end
%     disp(length(searchIDs{1})-i)
% end
% 
% writecell(cat(2,searchIDs{1},foundIDs),'uniprotWebscrapResults.txt','Delimiter','tab');
% 
% uUniprot = unique(foundIDs(:,1));
% uUniprot = uUniprot(~ismember(uUniprot,'NaN'));
% writecell(uUniprot,'uUniprotWebscrap.txt');

% BLOCK 4

% clear
% clc
% 
% fID = fopen('uniprotPROsr.tab','r');
% uniprotPRO = textscan(fID,'%s %s');
% fclose(fID);
% 
% fID = fopen('uniprotWebscrapResults.txt','r');
% uniprotWebscrap = textscan(fID,'%s %s %d %d');
% fclose(fID);
% 
% uniprotPROnew = cell(length(uniprotWebscrap{1}),2);
% for i = 1:length(uniprotPROnew)
%     if gt(uniprotWebscrap{3}(i),0)
%         uniprotPROnew{i,1} = uniprotWebscrap{1}{i};
%         [~,pos] = ismember(uniprotWebscrap{2}{i},uniprotPRO{1});
%         uniprotPROnew{i,2} = uniprotPRO{2}{pos}(uniprotWebscrap{3}(i):uniprotWebscrap{4}(i));
%     else
%         uniprotPROnew{i,1} = 'NaN';
%         uniprotPROnew{i,2} = 'NaN';
%     end
% end
% 
% uniprotPROnew = uniprotPROnew(~ismember(uniprotPROnew(:,1),'NaN'),:);
% 
% fID = fopen('uniprotSSr.txt','r');
% uniprot = textscan(fID,'%s %*s %s');
% fclose(fID);
% 
% fID = fopen('subSTRING.sr','r');
% origSTRING = textscan(fID,'%s %s');
% fclose(fID);
% 
% fID = fopen('subSTRINGnames.txt','r');
% subSTRINGnames = textscan(fID,'%s %s %*[^\n]','delimiter','.');
% fclose(fID);
% 
% fID = fopen('subSTRINGnames.txt','r');
% subSTRINGnames2 = textscan(fID,'%s');
% fclose(fID);
% 
% [~,uPos] = unique(uniprot{1});
% uAdd0 = cat(2,uniprot{1}(uPos),uniprot{2}(uPos));
% 
% [~,uPos2] = unique(uniprotPROnew(:,1));%should always be unique anyway
% uAdd1 = uniprotPROnew(uPos2,:);
% 
% uAdd2 = cat(1,uAdd0,uAdd1);
% uAdd = {uAdd2(:,1),uAdd2(:,2)};
% 
% fullSr = cell(length(uAdd{1}),2);
% count0 = 0;
% countMore = 0;
% for i = 1:length(fullSr)
%     [~,pos] = ismember(subSTRINGnames{2},uAdd{1}{i});
%     pos = logical(pos);
%     if gt(sum(pos),1)
%         countMore = countMore + 1;
%     elseif isequal(sum(pos),0)
%         keyboard
%         count0 = count0 + 1;%Should not happen!
%     else
%         fullSr{i,1} = subSTRINGnames2{1}{pos};
%         fullSr{i,2} = uAdd{2}{i};
%     end
% end
% 
% combo = cat(2,cat(1,fullSr(:,1),origSTRING{1}),cat(1,fullSr(:,2),origSTRING{2}));
% writecell(combo,'blastSTRINGsr.txt','Delimiter','tab');
% writecell(combo(:,1),'tmp.txt');
% 
% fID = fopen('tmp.txt','r');
% cTaxa = textscan(fID,'%s %*[^\n]','delimiter','.');
% fclose(fID);
% 
% lostID = ~ismember(subSTRINGnames2{1},combo(:,1));
% lostNames = subSTRINGnames2{1}(lostID);
% lostGenes = unique(subSTRINGnames{2}(lostID));
% 
% lostTID = ~ismember(subSTRINGnames{1},cTaxa{1});
% lostTaxa = unique(subSTRINGnames{1}(lostTID));
% 
% writecell(lostTaxa,'lostTaxa.txt');
% writecell(lostNames,'lostNames.txt');
% writecell(lostGenes,'lostGenes.txt');

%BLOCK 5

% load('speciesTax')
% 
% fID = fopen('blastSTRINGsr.txt','r');
% allProteins = textscan(fID,'%s %s');
% fclose(fID);
% 
% writecell(allProteins{1},'tmp.txt');
% 
% fID = fopen('tmp.txt','r');
% taxID = textscan(fID,'%s %*[^\n]','delimiter','.');
% fclose(fID);
% 
% taxID = str2double(taxID{1});
% 
% [~,vPos] = ismember('Viruses',speciesTax{2,2});
% [~,metPos] = ismember('Metazoa',speciesTax{2,2});
% [~,vertPos] = ismember('Vertebrata',speciesTax{2,2});
% [~,bacPos] = ismember('Bacteria',speciesTax{2,2});
% [~,fPos] = ismember('Fungi',speciesTax{2,2});
% [~,pPos] = ismember('Viridiplantae',speciesTax{2,2});
% [~,protPos] = ismember('Eukaryota',speciesTax{2,2});
% [~,hPos] = ismember('Homo sapiens',speciesTax{2,2});
% 
% virusTax = speciesTax{2,1}(speciesTax{2,3}(:,vPos));
% hostTax = speciesTax{2,1}(~speciesTax{2,3}(:,vPos));
% 
% bTax = speciesTax{2,1}(speciesTax{2,3}(:,bacPos));
% protTax = speciesTax{2,1}(logical(speciesTax{2,3}(:,protPos)-speciesTax{2,3}(:,fPos)-speciesTax{2,3}(:,pPos)-speciesTax{2,3}(:,metPos)));
% pTax = speciesTax{2,1}(speciesTax{2,3}(:,pPos));
% fTax = speciesTax{2,1}(speciesTax{2,3}(:,fPos));
% mTax = speciesTax{2,1}(logical(speciesTax{2,3}(:,metPos)-speciesTax{2,3}(:,vertPos)));
% vTax = speciesTax{2,1}(speciesTax{2,3}(:,vertPos));
% hTax = speciesTax{2,1}(speciesTax{2,3}(:,hPos));
% 
% isHost = ismember(taxID,hostTax);
% isVirus = ismember(taxID,virusTax);
% isB = ismember(taxID,bTax);
% isProt = ismember(taxID,protTax);
% isP = ismember(taxID,pTax);
% isF = ismember(taxID,fTax);
% isM = ismember(taxID,mTax);
% isV = ismember(taxID,vTax);
% isH = ismember(taxID,hTax);
% 
% hostCell = cat(2,allProteins{1}(isHost),allProteins{2}(isHost));
% virusCell = cat(2,allProteins{1}(isVirus),allProteins{2}(isVirus));
% bCell = cat(2,allProteins{1}(isB),allProteins{2}(isB));
% protCell = cat(2,allProteins{1}(isProt),allProteins{2}(isProt));
% pCell = cat(2,allProteins{1}(isP),allProteins{2}(isP));
% fCell = cat(2,allProteins{1}(isF),allProteins{2}(isF));
% mCell = cat(2,allProteins{1}(isM),allProteins{2}(isM));
% vCell = cat(2,allProteins{1}(isV),allProteins{2}(isV));
% hCell = cat(2,allProteins{1}(isH),allProteins{2}(isH));
% 
% writecell(hostCell,'hostSTRINGsr.txt','delimiter','tab');
% writecell(virusCell,'virusSTRINGsr.txt','delimiter','tab');
% 
% writecell(bCell,'bacteriaSTRINGsr.txt','delimiter','tab');
% writecell(protCell,'protistSTRINGsr.txt','delimiter','tab');
% writecell(pCell,'plantSTRINGsr.txt','delimiter','tab');
% writecell(fCell,'fungiSTRINGsr.txt','delimiter','tab');
% writecell(mCell,'metazoaSTRINGsr.txt','delimiter','tab');
% writecell(vCell,'vertebrateSTRINGsr.txt','delimiter','tab');
% writecell(hCell,'humanSTRINGsr.txt','delimiter','tab');

%BLOCK 6

rector = dir('*MMC.txt');
names = cell(length(rector),1);
ls = nan(length(names),2);
for i = 1:length(names)
    names{i} = rector(i).name(1:end-7);
    fID = fopen(rector(i).name,'r');
    tmp = textscan(fID,'%s %s','delimiter','\t');
    fclose(fID);
    tmp2 = nan(length(tmp{1}),1);
    for j = 1:length(tmp{1})
        tmp2(j) = length(strsplit(tmp{2}{j},' '));
    end
    ls(i,1) = sum(tmp2);
    ls(i,2) = length(tmp2);
end

[ratio,ord] = sort(ls(:,1)./ls(:,2));
names = names(ord);
ls = ls(ord,:);

for i = 1:length(names)
    names{i} = [names{i},'[',num2str(round(ratio(i),2)),']'];
end

figure
hold on
bar(linspace(1,3*length(ls(:,1))-2,length(ls(:,1))),ls(:,1),0.33)
bar(linspace(2,3*length(ls(:,1))-1,length(ls(:,1))),ls(:,2),0.33)
set(gca,'yscale','log')
set(gca,'XTick',linspace(1,3*length(ls(:,1))-2,length(ls(:,1)))+0.5)
set(gca,'XTickLabel',names)
legend({'#proteins','#mmseqs clusters'})
ylim([1,10^5])
title('#proteins/#clusters')
ylabel('count')




















