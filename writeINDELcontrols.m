

clear
clc

numRep = 1000;
dirName = 'littleMatchBatch';

%INDELible name, IQ-TREE name
modelTranslate = {'Poisson',        'Poisson';
                  'JTT',            'JTT';
                  'JTT-dcmut',      'JTTDCMut';
                  'Dayhoff',        'Dayhoff';
                  'Dayhoff-dcmut',  'DCMut';
                  'WAG',            'WAG';
                  'mtMAM',          'mtMAM';
                  'mtART',          'mtART';
                  'mtREV',          'mtREV';
                  'rtREV',          'rtREV';
                  'cpREV',          'cpREV';
                  'Vt',             'VT';
                  'Blosum',         'Blosum62';
                  'LG',             'LG';
                  'HIVb',           'HIVb';
                  'HIVw',           'HIVw'};

cd(dirName)
rector = dir('*.fa.iqtree');
fID = fopen('indel.tmp.txt','r');
inputLengths = textscan(fID,'%d %d');
fclose(fID);
inputLengths = cat(2,inputLengths{1},inputLengths{2});

strParams = cell(length(rector),8);
for i = 1:length(rector)
    mStr = 'Best-fit model according to BIC:';
    fStr = 'State frequencies: (empirical counts from alignment)';
    iStr = 'Proportion of invariable sites:';
    aStr = 'Gamma shape alpha:';
    
    tStr = 'Tree in newick format:';
    fID = fopen(rector(i).name,'r');
    lineBy = textscan(fID,'%s','delimiter','\n');
    fclose(fID);
    lineBy = lineBy{1};
    
    sName = rector(i).name;
    sName = sName(1:end-length('.fa.iqtree'));
    strParams{i,1} = sName;
    
    refCell = {'host','virus'};
    lenCol = 2-length(strfind(rector(i).name,'host'));
    lenRow = str2double(sName(length(refCell{lenCol})+2:end));
    seedL = num2str(inputLengths(lenRow,lenCol));
    strParams{i,2} = seedL;
    
    modelLine = lineBy{contains(lineBy,mStr)}(length(mStr)+2:end);
    modelCell = strsplit(modelLine,'+');
    subModel = modelCell{1};
    strParams{i,3} = modelTranslate{ismember(modelTranslate(:,2),subModel),1};
    
    isF = ismember('F',modelCell);
    if isF
       fLines = lineBy(find(contains(lineBy,fStr))+linspace(2,21,20));
       freqs = num2str(cellfun(@(x) str2double(x(9:end)),fLines)');
       strParams{i,4} = freqs;
    else
       strParams{i,4} = '';
    end
    
    isI = ismember('I',modelCell);
    if isI
        iVal = lineBy{contains(lineBy,iStr)}(length(iStr)+2:end);
        strParams{i,5} = iVal;
    else
        strParams{i,5} = '';
    end
    
    isG = isequal(sum(contains(modelCell,'G')),1+contains(subModel,'G'));
    if isG
        gNum = modelCell{find(contains(modelCell,'G'),1,'last')}(2:end);
        aVal = lineBy{contains(lineBy,aStr)}(length(aStr)+2:end);
        strParams{i,6} = gNum;
        strParams{i,7} = aVal;
    else
        strParams{i,6} = '';
        strParams{i,7} = '';
    end
    
    treeString = lineBy{find(contains(lineBy,tStr))+2}; 
    strParams{i,8} = treeString;
end

writeControl = {' ';' ';'// Anything on a line after two forward slashes is ignored.';...
    ' ';'[TYPE]  AMINOACID 1';' '};

for i = 1:length(strParams(:,1))
    writeControl = cat(1,writeControl,cat(2,'[MODEL] ',[strParams{i,1},'Model']));
    writeControl = cat(1,writeControl,cat(2,'[submodel] ',strParams{i,3}));
    writeControl = cat(1,writeControl,cat(2,'[rates] ',strParams{i,5},' ',...
       ' ',strParams{i,7},' ',strParams{i,6},' //pinv alpha ngamcat'));
    if ~isempty(strParams{i,4})
        writeControl = cat(1,writeControl,cat(2,'[statefreq] ',strParams{i,4},...
            ' //frequencies for A R N D C Q E G H I L K M F P S T W Y V, if model includes +F'));
    end
    writeControl = cat(1,writeControl,' ');
end

for i = 1:length(strParams(:,1))
    writeControl = cat(1,writeControl,cat(2,'[TREE] ',[strParams{i,1},'Tree '],strParams{i,8}));
end

writeControl = cat(1,writeControl,' ');

for i = 1:length(strParams(:,1))
    writeControl = cat(1,writeControl,cat(2,'[PARTITIONS] ',[strParams{i,1},'Partition //  [tree model root_length]']));
    writeControl = cat(1,writeControl,cat(2,'[',[strParams{i,1},'Tree '],...,
        [strParams{i,1},'Model '],strParams{i,2},']'));
    writeControl = cat(1,writeControl,' ');
end

writeControl = cat(1,writeControl,'[EVOLVE] // [partitionName #replicates outputName]');

for i = 1:length(strParams(:,1))
    writeControl = cat(1,writeControl,cat(2,[strParams{i,1},'Partition '],num2str(numRep),[' ',strParams{i,1},'INDEL']));
end

writeControl = cat(1,writeControl,' ');
writeControl = cat(1,writeControl,'// The true alignment will be output in a file named outputname_TRUE.phy');
writeControl = cat(1,writeControl,'// The unaligned sequences will be output in a file named outputname.fas');

cd ..

fID = fopen([dirName,'Control.dat'],'w');
for i = 1:length(writeControl)
    fprintf(fID,'%s\n',writeControl{i});
end
fclose(fID);

