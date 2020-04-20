

clear
clc

% %Block one compile accessions
% 
% rector = dir('*.csv');
% acCell = cell(length(rector),1);
% for i = 1:length(rector)
%     fID = fopen(rector(i).name,'r');
%     curr = textscan(fID,'%s %s %*[^\n]','delimiter',',');
%     fclose(fID);
%     acCell{i} = unique(curr{2});
% end
% acCell2 = unique(cat(1,acCell{:}));
% writecell(acCell2,'acc.txt');

%Block two match virus host and host taxa

fID = fopen('accTax.txt','r');
accTax = textscan(fID,'%s %s %s','delimiter',',');
fclose(fID);

fID = fopen('virusAChost.txt','r');
vACH = textscan(fID,'%s %s','delimiter',',');
fclose(fID);

fID = fopen('STRINGtargets.txt','r');
sTar = textscan(fID,'%s %s');
fclose(fID);

acCell = cell(length(sTar{1}),5);
for i = 1:length(acCell)
    fID = fopen([sTar{1}{i},'.csv'],'r');
    curr = textscan(fID,'%s %s %*[^\n]','delimiter',',');
    fclose(fID);
    writecell(curr{2},'tmp.txt');
    fID = fopen('tmp.txt','r');
    curr = textscan(fID,'%s %*[^\n]','delimiter','.');
    fclose(fID);
    hostAcc = curr{1};
    
    fID = fopen([sTar{2}{i},'.csv'],'r');
    curr = textscan(fID,'%s %s %*[^\n]','delimiter',',');
    fclose(fID);
    writecell(curr{2},'tmp.txt');
    fID = fopen('tmp.txt','r');
    curr = textscan(fID,'%s %*[^\n]','delimiter','.');
    fclose(fID);
    virusAcc = curr{1};
    
    %Since they may be called from blast at different times, remove the
    %accessions that don't appear in accTax. Also remove the virus
    %accessions missing from vACH. Sometimes virusAcc includes homologs
    %from cellular organisms which we want to exclude for the summary
    %statistics.
    hostAcc = hostAcc(ismember(hostAcc,accTax{1}));
    virusAcc = virusAcc(ismember(virusAcc,accTax{1}));
    virusAcc = virusAcc(ismember(virusAcc,vACH{1}));
    
    [~,virusHostID] = ismember(virusAcc,vACH{1});
    virusHost = vACH{2}(virusHostID);
    noHost = cellfun('isempty', virusHost);
    [~,hostNameID] = ismember(hostAcc,accTax{1});
    [~,virusNameID] = ismember(virusAcc,accTax{1});
    hostNames = accTax{3}(hostNameID);
    virusNames = accTax{3}(virusNameID);
    
    virusWHost = unique(virusNames(~noHost));
    
    hMatchAcc = cell(length(virusAcc),1);
    repVAcc = cell(length(virusAcc),1);
    for j = 1:length(virusAcc)
        hMatchAcc{j} = hostAcc(ismember(hostNames,virusHost{j}));
        repVAcc{j} = repmat(virusAcc(j),length(hMatchAcc{j}),1);
    end
    matchCell_1 = cat(1,hMatchAcc{:});
    matchCell_2 = cat(1,repVAcc{:});
    
    [~,fN1] = ismember(matchCell_1,accTax{1});
    [~,fN2] = ismember(matchCell_2,accTax{1});
    
    matchCell = cat(2,matchCell_1,matchCell_2);
    matchNames = join(cat(2,accTax{3}(fN1),accTax{3}(fN2)),'_');
    nTab = tabulate(matchNames);
    [~,ord] = sort(cell2mat(nTab(:,2)),'descend');
    nTab = nTab(ord,:);
    
    numHost = length(unique(accTax{3}(fN1)));
    numVirus = length(unique(accTax{3}(fN2)));
    numNoHostOrtho = length(virusWHost(~ismember(virusWHost,accTax{3}(fN2))));
    numNoHost = length(unique(virusNames(~ismember(virusNames,virusWHost))));
    
    labl = [num2str(i),'_',strrep(nTab{1,1},' ','-'),'.txt'];
    
    writecell(matchCell,labl,'delimiter','tab');
    acCell{i,1} = labl;
    acCell{i,2} = numHost;
    acCell{i,3} = numVirus;
    acCell{i,4} = numNoHostOrtho;
    acCell{i,5} = numNoHost;
    
end

recordCell = cat(2,sTar{1},sTar{2},acCell);
recordCell = cat(1,{'hostKey','virusKey','aliName','numHost','numVirus','vNoOrtho','vNoHost'},recordCell);
writecell(recordCell,'hostMatchIndex.txt','delimiter','tab')






