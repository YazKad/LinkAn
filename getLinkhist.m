% getLinkHist
% Written By Yazan Kadkoy 5/2018
% This program takes _BonePart files and created a linkage histogram for 
% each possible linkage value. The program then calculated BVlink
% min,max,avg for each linkage value. place _link files into one folder and
% select that folder with prompted



clear all
clc
close all



for s_f=1:1

    selec_folder=s_f;
    
    if selec_folder == 2
    clearvars -except selec_folder bonepix zray_length percent
    end
    
% Load in directory and folder name
folder = uigetdir('C:\Users\Yazan\Documents\research');
dirListing = dir(folder);

for d = 3:length(dirListing)
    name{d-2}=dirListing(d).name;
end

% name=sort_nat(name);
end
fN = [];
for d = 1:length(name)
fileNames = fullfile(folder,char(name(d)));
fN = strvcat(fN,fileNames);
end
 
[num len]=size(fN);

for i=1:num;
 FN=char(fN(i,:));
 load (FN)


Bone(1).name=BonePart(1).name;

objCallus=BonePart(1).Callus.objconn;
objInner=BonePart(1).Inner.objconn;
objTotal=BonePart(1).Total.objconn;
objeffIn=BonePart(1).effInner.objconn;
objeffCal=BonePart(1).effCallus.objconn;
objeffTotal=BonePart(1).effTotal.objconn;


[Ex,Wy,Ze]=size(objCallus);
Object(1).obj=zeros(Ex,Wy,Ze);
Object(2).obj=zeros(Ex,Wy,Ze);
Object(3).obj=zeros(Ex,Wy,Ze);
Object(4).obj=zeros(Ex,Wy,Ze);
Object(5).obj=zeros(Ex,Wy,Ze);
Object(6).obj=zeros(Ex,Wy,Ze);

Object(1).obj=objCallus;
Object(2).obj=objInner;
Object(3).obj=objTotal;
Object(4).obj=objeffCal;
Object(5).obj=objeffIn;
Object(6).obj=objeffTotal;


for I=1:numel(Object)

for z=1:Ze
   im=Object(I).obj(:,:,z);
   
   bound=bwboundaries(im);
   sslice=cellfun('length',bound);
  

% hold on
% imshow(im);

   stats = regionprops('table',im,'Centroid',...
    'MajorAxisLength','MinorAxisLength','Area');
    centers = stats.Centroid;
      

bound=bwboundaries(im);

    if(isempty(centers))

     
     tv=[];

    else
count=1;
 bound2=bwboundaries(bwareaopen(im,25)); 
       for k=1:numel(bound2)
    boundaryx2= bound2{k}(:,2);
    boundaryy2= bound2{k}(:,1);
        
        for q=1:numel(boundaryx2)
            bx(count)=boundaryx2(q);
            by(count)=boundaryy2(q);
                count=count+1;
        end
      end
  cx=BonePart(1).Centroid(1);
  cy=BonePart(1).Centroid(2);
  %calculate all distaces and associated angles
       allDistances = sqrt((bx-cx).^2 + (by-cy).^2);
       allAngles = round(atan2d((bx-cx),(by-cy))); 
        
       if(~isnan(allAngles))
uangles=unique(allAngles);
%find largest distance for every angle
  for w=1:numel(uangles)
      angle=uangles(w);
      ind=find(allAngles==angle);
      [mxl mxi]=max(allDistances(ind));
      distin2(w)=find(allDistances==mxl,1);
      [ml mi]=min(allDistances(ind));
      distin(w)=find(allDistances==ml,1);
  end
end
  for e=1:numel(distin2)
      outerx(e)=bx(distin2(e));
      outery(e)=by(distin2(e));
      XY(e,1)=bx(distin2(e));
      XY(e,2)=by(distin2(e));
  end
%%code based off ellipse fit (2009, Nikolai Chernov)
%ellipse approximated to determine TV for the slice
      centroid=mean(XY);
    D1=[(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
      (XY(:,2)-centroid(2)).^2];
  D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
    S1 = D1'*D1;
    S2 = D1'*D2;
    S3 = D2'*D2;
    T = -inv(S3)*S2';
    M = S1 + S2*T;
    M = [M(3,:)./2; -M(2,:); M(1,:)./2];
    [evec,eval] = eig(M);
    cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
    A1 = evec(:,find(cond>0));
    A = [A1; T*A1];
    A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
    A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
    A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
         A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
    A(4) = A4;  A(5) = A5;  A(6) = A6;
    A = A/norm(A);
    syms x  y
    B=A(2);
    C=A(3);
    D=A(4);
    E=A(5);
    F=A(6);
    A=A(1);
     %% Roger Stafford 19 Jul 2014
     e = 4*A*C-B^2; if e<=0, error('This conic is not an ellipse.'), end
   x0 = (B*E-2*C*D)/e; y0 = (B*D-2*A*E)/e;   % Ellipse center
   F0 = -2*(A*x0^2+B*x0*y0+C*y0^2+D*x0+E*y0+F);
   g = sqrt((A-C)^2+B^2); a = F0/(A+C+g); b = F0/(A+C-g);
   if (a<=0)|(b<=0), error('This is a degenerate ellipse.'), end
   a = sqrt(a);  b = sqrt(b); % Major & minor axes
   t = 1/2*atan2(B,A-C); ct = cos(t); st = sin(t);   % Rotation angle
   p = linspace(0,2*pi,100); cp = cos(p); sp = sin(p);   % Variable parameter
   x = x0+a*ct*cp-b*st*sp; y = y0+a*st*cp+b*ct*sp;   % Generate points on ellipse


   %Get TV 
    Bone(I).TV(z)=pi*a*b;
    tv(z)=pi*a*b;
    %get histogram each bin 0-26
    for g=1:27
        bin=0;
        rcount=0;
        rcount=sum(im==g);
        bin=sum(rcount);
        hist(g)=bin;
    end
    %get summation of histogram for every slize
    for h=1:27
        
       sumhist(z,h)=sum(hist(h:end)); 
        sumhist2(z,h)=sum(hist(1:h));
    end
    [~,ii]=max(hist);
    
  if I==1
Bone(z).peak1=ii+1;
Bone(z).hist1=hist;
  elseif I==2
 Bone(z).peak2=ii+1;
Bone(z).hist2=hist;     
  elseif I==3
 Bone(z).peak3=ii+1;
Bone(z).hist3=hist;     
  elseif I==4
Bone(z).peak4=ii+1;
Bone(z).hist4=hist;      
  elseif I==5
 Bone(z).peak5=ii+1;
Bone(z).hist5=hist;     
  elseif I==6
 Bone(z).peak6=ii+1;
Bone(z).hist6=hist;                 
  end
end
end
Bone(1).ROI="Callus";
Bone(2).ROI="Inner";
Bone(3).ROI="Total";
Bone(4).ROI="Effective Callus";
Bone(5).ROI="Effective Inner";
Bone(6).ROI="Effective Total";
 %get max min and ave for every possible linkage value
 for ii=1:27
     holdsumhist=[];
     holdsumhist2=[];
     
     holdsumhist=sumhist(1:end,ii);
     holdsumhist2=sumhist2(1:end,ii);
     
    
  
    maxhist(ii)=max(holdsumhist);
    minhist(ii)=min(holdsumhist(holdsumhist>0));
    avehist(ii)=nanmean(holdsumhist);
    
    maxhist2(ii)=max(holdsumhist2);
    minhist2(ii)=min(holdsumhist2(holdsumhist>0));
    avehist2(ii)=nanmean(holdsumhist2);
    
    
    
    
 end

Bone(I).TVmax=max(tv);
Bone(I).TVmin=min(tv(tv>0));
Bone(I).TVave=nanmean(tv);
Bone(I).histmax=maxhist;
Bone(I).histmin=minhist;
Bone(I).avehist=avehist;
Bone(I).histmax2=maxhist2;
Bone(I).histmin2=minhist2;
Bone(I).avehist2=avehist2;

end
savename=strcat(Bone(1).name,'_ConnHist.mat')
save(savename,'Bone')
clearvars -except fn name folder fN
       
              

end

