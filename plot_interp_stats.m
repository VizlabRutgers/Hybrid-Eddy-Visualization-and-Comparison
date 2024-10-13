% look at eddy interpolation results
load interpolationEvaluation.mat
% variables: EddyIndex, start frame, total stride, 
%            interpolated frame, end frame, approximate depth
%            Euclidean distance, cosine distance (after scaled),
%            correlation, SSIM
%

eddyID=interpolationEvalutation.EddyHistoryIndex;
startF=interpolationEvalutation.('start frame');
endF=interpolationEvalutation.('end frame');
stride=interpolationEvalutation.('total stride');
thisF=interpolationEvalutation.('interpolated frame');
depth=interpolationEvalutation.('depth in ground truth');
Edist=interpolationEvalutation.('Euclidean distance');
Cosdist=interpolationEvalutation.('cosine distance (after scaled)');
Ecorr=interpolationEvalutation.('correlation');
SSIM=interpolationEvalutation.('SSIM');
DTW_distance=interpolationEvalutation.('DTW distance');

eddylist=unique(eddyID);
fprintf('%d eddies covered \n',length(eddylist))

colormap(colorcube)

figure(1)
swarmchart(stride,SSIM,depth,eddyID)
xlabel('Stride')
ylabel('SSIM')
ylim([-0.4,1])

figure(2)
thisEddy=13;
swarmchart(stride(eddyID==thisEddy),SSIM(eddyID==thisEddy),...
    depth(eddyID==thisEddy),eddyID(eddyID==thisEddy))
xlabel('Stride')
ylabel('SSIM')
title(['For Eddy ' num2str(thisEddy)])


meanEdist=zeros(length(eddylist),5);
meandepth=ones(length(eddylist),5);
minstart=zeros(length(eddylist),5);
maxend=zeros(length(eddylist),5);
eddydur=zeros(length(eddylist),5);
eddylab=zeros(length(eddylist),5);
eddystr=zeros(length(eddylist),5);
for i=1:length(eddylist)
    for j=2:6
        if sum(eddyID==eddylist(i) & stride==j)>0
    fprintf('eddy %d, stride %d: found %d entries \n',...
        eddylist(i),j,sum(eddyID==eddylist(i) & stride==j))
    meanEdist(i,j-1)=mean(Edist(eddyID==eddylist(i) & stride==j));
    meandepth(i,j-1)=mean(depth(eddyID==eddylist(i) & stride==j));
    minstart(i,j-1)=min(startF(eddyID==eddylist(i) & stride==j));
    maxend(i,j-1)=max(endF(eddyID==eddylist(i) & stride==j));
    eddydur(i,j-1)=maxend(i,j-1)-minstart(i,j-1);
    eddylab(i,j-1)=eddylist(i);
    eddystr(i,j-1)=j;
        else
    meanEdist(i,j-1)=NaN;
    meandepth(i,j-1)=NaN;
    minstart(i,j-1)=NaN;
    eddydur(i,j-1)=NaN;
    eddylab(i,j-1)=eddylist(i);
    eddystr(i,j-1)=j;
        end
    end
end

figure(3)
scatter(eddydur(:),meanEdist(:),20*eddystr(:),eddylab(:),'h')
xlabel('Duration (frames)')
ylabel('Euclidian Distance')

figure(4)
scatter(eddydur(:),meandepth(:),20*eddystr(:),eddylab(:),'h','filled')
xlabel('Duration (frames)')
ylabel('Mean Depth')


