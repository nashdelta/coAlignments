function plotBipartite1(fEmpty,catID,virusName,hostCatNames,iMat,hc,w1,w2,cThresh,maxR,rN,titleVar,numCol)
%Constructs bipartite graph, designed for hosts(rows) and viruses(columns)
%interaction. Hosts are placed along a circle with weighted spacing.
%Viruses are placed within that circle with distances from hosts defined by
%the number of protein-protein interactions identified between them. Hosts
%are divided up into pre-assigned categories and categories are separated
%by a uniform spacing.

%fEmpty is fraction of circle that's left as empty space between host
%categories

%catID is a vector indicating the category of the host must have no gaps
%i.e. if max(catID)==5, catID must have at least one element 1,2,3,4,5
%each.

%virusName is the name of the points connected to the categories in
%hotsCatNames, originally intended to be 'Viruses'

%hostCatNames is a cell array of names of categories must be matched with ID
%i.e. catName{1} is name of category with ID 1.

%iMat is matrix of interactions between hosts(rows) and viruses (columns),
%elements are integers

%hc specifies host colormap 1=parula, 2=jet, 3=hsv

%w1=[1 or 2] specifies the weighting rule for virus-host distances. Options are
%linear weights of number of interactions with host, log weights of number
%of interactions with host

%w2=[1,2,3 or 4] specifies the weighting rule for host distribution around the circle.
%Options are linear weights of total number of interactions including host,
%log weights total number of interactions, linear weights of total number
%of interacting viruses, log weights of total number of interacting
%viruses.

%cThresh=number Edges are colored light or dark gray. Edges corresponding to
%elements in iMat larger than cThresh are dark gray.

%maxR=+num<1 specifies the maximum ratio of a virus point from the origin
%to a host i.e. if == 1 viruses can overlap hosts in plot

%rN=integer removes top number of hosts specified with the highest host
%weighting as specified by w2. Note: if tied weights abritrarily picks out
%of tie. Removes columns with no nonzero entries after removal of specified
%rows.

%titleVar, title of plot

%numCol is the number of columns in the legend

%set w2
if isequal(w2,1)
    hW = sum(iMat,2);
elseif isequal(w2,2)
    hW = log(sum(iMat,2))+1;
elseif isequal(w2,3)
    hW = sum(gt(iMat,0),2);
else
    hW = log(sum(gt(iMat,0),2))+1;
end
hW = hW/sum(hW);

%remove rows if specified and reset w2
if gt(rN,0)
    [~,ord] = sort(hW,'descend');
    kIndex = true(size(hW));
    kIndex(ord(1:rN)) = false;
    
    catID = catID(kIndex);
    iMat = iMat(kIndex,:);
    
    hostCatNames = hostCatNames(unique(catID));
    
    [~,catID] = ismember(catID,unique(catID));
    
    iMat = iMat(:,gt(sum(iMat,1),0));
    
    if isequal(w2,1)
        hW = sum(iMat,2);
    elseif isequal(w2,2)
        hW = log(sum(iMat,2))+1;
    elseif isequal(w2,3)
        hW = sum(gt(iMat,0),2);
    else
        hW = log(sum(gt(iMat,0),2))+1;
    end
    hW = hW/sum(hW);
    
end

%set w1

if isequal(w1,1)
    vW = iMat;
else
    vW = zeros(size(iMat));
    vW(gt(iMat,0)) = log(iMat(gt(iMat,0)));
end
vW = vW./repmat(sum(vW,1),length(vW(:,1)),1);

%calculate host category weights

wC = nan(length(max(catID)),1);
for i = 1:max(catID)
    wC(i) = sum(hW(catID==i));
end

catIndex = cell(max(catID),1);
for i = 1:length(catIndex)
    catIndex{i} = find(catID==i);
end

%parametrize hosts along the circle
hT = nan(length(hW),1);
for i = 1:max(catID)
    hT(catIndex{i}(1)) = sum(wC(1:i-1))*(1-fEmpty)+(i-1)*fEmpty/max(catID)+1/2*hW(catIndex{i}(1))*(1-fEmpty);
    for j = 2:length(catIndex{i})
        hT(catIndex{i}(j)) = hT(catIndex{i}(j-1))+1/2*(hW(catIndex{i}(j-1))+hW(catIndex{i}(j)))*(1-fEmpty);
    end
end
hPos = [cos(2*pi*hT),sin(2*pi*hT)];

%calculate virus positions within circle
vPos(1,:) = sum(repmat(hPos(:,1),1,length(iMat(1,:))).*vW,1);
vPos(2,:) = sum(repmat(hPos(:,2),1,length(iMat(1,:))).*vW,1);
vPos = vPos*maxR;

%specify edge colors
grayMat = [0.9,0.9,0.9;0.2,0.2,0.2];

%specify host colors
if isequal(hc,1)
    hostColor = parula(max(catID));
elseif isequal(hc,2)
    hostColor = jet(max(catID));
else
    hostColor = hsv(max(catID));
end

%generate plot
figure
hold on
leL = plot([-1000,-1000],[-1000,-1000],'linewidth',0.1,'color',grayMat(1,:));
leD = plot([-1000,-1000],[-1000,-1000],'linewidth',0.1,'color',grayMat(2,:));
for i = 1:length(iMat(:,1))
    for j = 1:length(iMat(2,:))
        if gt(iMat(i,j),0)&&le(iMat(i,j),cThresh)
            plot([hPos(i,1),vPos(1,j)],[hPos(i,2),vPos(2,j)],'linewidth',0.1,'color',grayMat(1,:));
        end
    end
end
for i = 1:length(iMat(:,1))
    for j = 1:length(iMat(2,:))
        if gt(iMat(i,j),cThresh)
            plot([hPos(i,1),vPos(1,j)],[hPos(i,2),vPos(2,j)],'linewidth',0.1,'color',grayMat(2,:));
        end
    end
end 
vL = scatter(vPos(1,:),vPos(2,:),[],[0.5 0.5 0.5],'fill','markeredgecolor','k');

hL = nan(max(catID),1);
for i = 1:max(catID)
    hL(i) = scatter(hPos(catID==i,1),hPos(catID==i,2),[],hostColor(i,:),'fill','markeredgecolor','k');
end
numViruses = length(iMat(1,:));
for i = 1:length(hostCatNames)
    hostCatNames{i} = [hostCatNames{i},'[',num2str(sum(catID==i)),']'];
end
if gt(cThresh,0)
    legend([vL; hL; leL; leD],{[virusName,'[',num2str(numViruses),']'],hostCatNames{:},['<',num2str(cThresh)],['>',num2str(cThresh)]},'NumColumns',numCol) %#ok<CCAT>
else
    legend([vL; hL],{[virusName,'[',num2str(numViruses),']'],hostCatNames{:}},'NumColumns',numCol) %#ok<CCAT>
end
title(titleVar)

xlim([-1.1,1.1])
ylim([-1.1,1.1])
axis off

end

