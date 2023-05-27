% segBone 
% Written by Yazan Kadkoy 3/2018
% This program segments linkage matracies representing cortical and callus of a fracture callus
%All desired _link files should placed into a folder. This folder should be
%selected when prompted


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

Table=["name";"max";"remainder";"percent bridged"];
for i=1:num
    FN=char(fN(i,:));
    load (FN)
count=1;
[Ro Co]=size(rotatedbone(2).vox);
totmask=zeros(Ro,Co);
%Centroid is determined as the average centroid of the top and bottom 10
%slices. 
for z=2:numel(rotatedbone)-1
    
    im=rotatedbone(z).vox;
    
    stat2=regionprops(true(size(im)), imfill(im,'holes'),  'WeightedCentroid');
   
   
    
    
    if((z<=10)||(z>=numel(rotatedbone)-11))
        
     cx1(count)=stat2.WeightedCentroid(1);
    cy1(count)=stat2.WeightedCentroid(2); 
    count=count+1;
    end
    
end
    cx=nanmean(cx1);
    cy=nanmean(cy1);
    
for z=2:numel(rotatedbone)-1
    distin2=[];
  im=rotatedbone(z).vox;
  [Ro Co]=size(im);

  count=1;
%   calculate distances and angles of every positive voxel to the centroid
    [bx by]=find(rotatedbone(z).vox);
       allDistances = sqrt((bx-cx).^2 + (by-cy).^2);
       allAngles = round(atan2d((bx-cx),(by-cy))); 
      if(~isnan(allAngles))
uangles=unique(allAngles);


  for w=1:numel(uangles)
      angle=uangles(w);
      ind=find(allAngles==angle);
  % determine the largest point for every angle
    [mx mxi]=max(allDistances(ind));

      distin2(w)=ind(mxi);
      
  
      
      De(w)=mx;
      x1=bx(ind(mxi));
      y1=by(ind(mxi));
      
   %create a slightly larger outline  
       XY(w,1)=x1+(mx*.3*(x1-cx)/mx);
      XY(w,2)=y1+(mx*.3*(y1-cy)/mx); 
           
      
  end  
       
      end   

    %reorder points  
       ang = atan2( (XY(:,2)-cy),(XY(:,1)-cx));
    [sortedAngles, sortIndices] = sort(ang);
    xnew=XY(:,1);
    ynew=XY(:,2);
     rx = xnew(sortIndices);
    ry = ynew(sortIndices);
    
   
   %create mask
    
    mask1=poly2mask(rx,ry,Ro,Co);
    mask=imfill(mask1,'holes');
    leftover=rotatedbone(z).vox.*~mask;
    r=[];
    c=[];
    %cover any leftover points
    [r,c]=find(leftover);
    if(~isempty(r))
    for i=1:numel(r)
        mask(r(i),c(i))=1;
    end
    end
    
    
    totmask=totmask+mask;
end
for z=2:numel(rotatedbone)
    objact(:,:,z)=rotatedbone(z).vox;
    objmask(:,:,z)=totmask;
%display mask for quality control
    subplot(1,3,1)
    imshow(rotatedbone(z).vox)
    pause(0.1)
    hold on
    subplot(1,3,2)
    imshow(mask)
    subplot(1,3,3)
    imshow(rotatedbone(z).vox.*~mask)
end
      for xx=1:Ro
          for yy=1:Co
             
    
                  
                  
              zfind=find(objact(xx,yy,:));
              %find sum every linkage value in a given row and column
              zrayact(xx,yy)=sum(objact(xx,yy,zfind));
              %caculate the maximum possible linkage value summation
              zrayest(xx,yy)=numel(find(objmask(xx,yy,:)))*27;
              
              if(zrayest(xx,yy)==0)
                zraynorm(xx,yy)=0;  
              else
               
              zraynorm(xx,yy)=zrayact(xx,yy)/zrayest(xx,yy);
          end
             
          end
      end

        [r,c,v]=find(zraynorm);
     allDistances2 = sqrt((c-cx).^2 + (r-cy).^2);
       allAngles2 = round(atan2d((c-cx),(r-cy)));
       
       
       
        if(~isnan(allAngles))
uangles2=unique(allAngles2);

  for w=1:numel(uangles2)
      angle=uangles2(w);
      ind=find(allAngles2==angle);
%find largest percentage for every angle
[mx mxi]=max(v(ind));

      
      xmaxconn(w)=c(ind(mxi));
      ymaxconn(w)=r(ind(mxi));
      
      maxconndis(w)=allDistances2(ind(mxi));
  
      XY2(w,1)=xmaxconn(w);
      XY2(w,2)=ymaxconn(w); 
      
  end  
    
        end 
%%code based off ellipse fit (2009, Nikolai Chernov)
%create ellipse approximating locations of max linkage value percentage for
%every angle from centroid
       centroid=mean(XY2);
    D1=[(XY2(:,1)-centroid(1)).^2, (XY2(:,1)-centroid(1)).*(XY2(:,2)-centroid(2)),...
      (XY2(:,2)-centroid(2)).^2];
  D2 = [XY2(:,1)-centroid(1), XY2(:,2)-centroid(2), ones(size(XY2,1),1)];
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
    syms xx  yy
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
   xx = x0+a*ct*cp-b*st*sp; yy = y0+a*st*cp+b*ct*sp;   % Generate points on ellipse   
   
  
    
   %create mask to seperate cortical bone from callus
    
    
    mask1=poly2mask(xx,yy,Ro,Co);

    maskout=imfill(mask1,'holes'); 
    maskin=~maskout;
   
    for zzz=2:numel(rotatedbone)-1   
        objinner(:,:,zzz)=rotatedbone(zzz).vox.*maskout;
        objcallus(:,:,zzz)=rotatedbone(zzz).vox.*maskin;
        objTotal(:,:,zzz)=rotatedbone(zzz).vox;
    end   
     for xx=1:Ro
          for yy=1:Co
             
    
                  
                  
              zfind2=find(objcallus(xx,yy,:));
              zrayact2(xx,yy)=sum(objcallus(xx,yy,zfind2));
              
             
              
              if(zrayest(xx,yy)==0)
                zraynorm2(xx,yy)=0;  
              else
              zraynorm2(xx,yy)=zrayact2(xx,yy)/zrayest(xx,yy);
          end
             
          end
      end

        [r2,c2,v2]=find(zraynorm2);
     allDistances3 = sqrt((c2-cx).^2 + (r2-cy).^2);
       allAngles3 = round(atan2d((c2-cx),(r2-cy)));
       
       
       
        if(~isnan(allAngles3))
uangles3=unique(allAngles3);
%record location of max linkage value for each angle
  for w=1:numel(uangles3)
      angle=uangles3(w);
      ind=find(allAngles3==angle);

[mx mxi]=max(v(ind));
[md mdi]=max(allDistances3(ind));


      XY3(w,1)=x1+(mx*.3*(x1-cx)/mx);
      XY3(w,2)=y1+(mx*.3*(y1-cy)/mx);

      
      xmaxconn2(w)=c2(ind(mxi));
      ymaxconn2(w)=r2(ind(mxi));

      maxconndis2(w)=allDistances3(ind(mxi));
  
      
      
  end  
    
        end 
%%code based off ellipse fit (2009, Nikolai Chernov)        
     centroid=mean(XY3);
    D1=[(XY3(:,1)-centroid(1)).^2, (XY3(:,1)-centroid(1)).*(XY3(:,2)-centroid(2)),...
      (XY3(:,2)-centroid(2)).^2];
  D2 = [XY3(:,1)-centroid(1), XY3(:,2)-centroid(2), ones(size(XY3,1),1)];
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
    syms xx  yy
    B=A(2);
    C=A(3);
    D=A(4);
    E=A(5);
    F=A(6);
    A=A(1);
 %appoximate ellipse to determine TV 
  %% Roger Stafford 19 Jul 2014
     e = 4*A*C-B^2; if e<=0, error('This conic is not an ellipse.'), end
   x0 = (B*E-2*C*D)/e; y0 = (B*D-2*A*E)/e;   % Ellipse center
   F0 = -2*(A*x0^2+B*x0*y0+C*y0^2+D*x0+E*y0+F);
   g = sqrt((A-C)^2+B^2); a = F0/(A+C+g); b = F0/(A+C-g);
   if (a<=0)|(b<=0), error('This is a degenerate ellipse.'), end
   a = sqrt(a);  b = sqrt(b); % Major & minor axes
    
    
    objcallusbw=objcallus>=1;
    objinnerbw=objinner>=1;
    objTotalbw=objTotal>=1;
    statscal=regionprops3(objcallusbw,'Volume','Image','Centroid','Orientation','VoxelList');
    statsinner=regionprops3(objinnerbw,'Volume','Image','Centroid','Orientation','VoxelList');
    statstotal=regionprops3(objTotalbw,'Volume','Image','Centroid','Orientation','VoxelList');
    
    
    


str=strsplit(rotatedbone(1).name,'_');
str=string(str(1));
BonePart(1).name=str;
BonePart(1).TV=pi*a*b*numel(rotatedbone);
BonePart(1).zrayact=zrayact;
BonePart(1).zrayest=zrayest;
BonePart(1).zraynorm=zraynorm;
BonePart(1).zraynorm2=zraynorm2;
BonePart(1).maxconndist=maxconndis;
BonePart(1).xmaxconn=xmaxconn;
BonePart(1).ymaxconn=ymaxconn;
BonePart(1).maxconndistcal=maxconndis2;
BonePart(1).xmaxconncal=xmaxconn2;
BonePart(1).ymaxconncal=ymaxconn2;
BonePart(1).Callus.obj=objcallusbw;
BonePart(1).Callus.objconn=objcallus;
BonePart(1).Callus.stat=statscal;
BonePart(1).Inner.obj=objinnerbw;
BonePart(1).Inner.objconn=objinner;
BonePart(1).Inner.stat=statsinner;
BonePart(1).Centroid=[cx,cy];
BonePart(1).Total.obj=objTotalbw;
BonePart(1).Total.objconn=objTotal;
BonePart(1).Total.stat=statstotal;
% determinatino of break type
    Calvol=BonePart(1).Callus.stat;
    calmax=max(Calvol.Volume);
    indcal=find(Calvol.Volume>=.15*calmax);
    
    calsum=sum(Calvol.Volume(indcal(:)));
    
    Innervol=BonePart(1).Inner.stat;
    Inmax=max(Innervol.Volume);
    indin=find(Innervol.Volume>=.15*Inmax);
    tv=BonePart(1).TV;
    Insum=sum(Innervol.Volume(indin(:)));
   if(Inmax/Insum>=.95)
       BT="Union";
   else
       BT="non-union";
   end
   
    
       siz=size(BonePart(1).Inner.obj);
       obj=zeros(siz);
       obj(:,:,1:5)=1;
       obj(:,:,end-5:end)=1;
       objintemp=zeros(siz);

countin=0;
countout=0;
%create effective objects

       for iii=1:numel(indin)
           objtemp=obj;
        vox=(Innervol.VoxelList{indin(iii),1});
       inds=sub2ind(siz,vox(:,2),vox(:,1),vox(:,3));
       objtemp(inds)=1;
       objintemp(inds)=1;
       objtemp=objtemp>=1;
       objintemp=objintemp>=1;
       
       
       statobj=regionprops3(objtemp,'Volume','VoxelList');
       
       numo=3;
       numt=numel(statobj.Volume);
       if(numo-numt>=2)
          countin=countin+1;
          effin(countin)=iii;
           
           
       end
       end
       statobj=regionprops3(objintemp,'Volume','VoxelList');
      
       
       if(countin==0)
           BT="Partial union";
        numo=numel(statobj.Volume);
        for ii=1:numel(Calvol.VoxelList)
                
           objtemp=objintemp;
            
            vox2=Calvol.VoxelList{ii,1};
             inds2=sub2ind(siz,vox2(:,2),vox2(:,1),vox2(:,3));
            objtemp(inds2)=1;
            objtemp=objtemp>=1;
            
            stattemp=regionprops3(objtemp,'Volume');
            numt=numel(stattemp.Volume);
            
            if(numo-numt>=2)
                countout=countount+1;
                effout(countout)=ii;
            end
        end
        
       else
           numo=3;
            for ii=1:numel(Calvol.VoxelList)
                
           objtemp=obj;
            
            vox2=Calvol.VoxelList{ii,1};
             inds2=sub2ind(siz,vox2(:,2),vox2(:,1),vox2(:,3));
            objtemp(inds2)=1;
            objtemp=objtemp>=1;
            
            stattemp=regionprops3(objtemp,'Volume');
            numt=numel(stattemp.Volume);
            
            if(numo-numt>=2)
                countout=countout+1;
                effout(countout)=ii;
            end
        
        
       end
        
       end
       effinobj=zeros(siz);
       effoutobj=zeros(siz);
       effobj=zeros(siz);
       
       if(countin>=1)
        for xx=1:countin
            effvox=Innervol.VoxelList{effin(xx),1};
            effinds=sub2ind(siz,effvox(:,2),effvox(:,1),effvox(:,3));
         effinobj(effinds)=1;
         effobj(effinds)=1;
         effinobj=effinobj>=1;
         effobj=effobj>=1;
        end
       end
       if(countout>=1)
           for yy=1:countout
               effvox=Calvol.VoxelList{effout(yy),1};
               effinds=sub2ind(siz,effvox(:,2),effvox(:,1),effvox(:,3));
            effoutobj(effinds)=1;
            effobj(effinds)=1;
            effoutobj=effoutobj>=1;
            effobj=effobj>=1;
           end
       end
    statseffcal=regionprops3(effoutobj,'Volume','Image','Centroid','Orientation','VoxelList');
    statseffinner=regionprops3(effinobj,'Volume','Image','Centroid','Orientation','VoxelList');
    statseff=regionprops3(effobj,'Volume','Image','Centroid','Orientation','VoxelList');
    
    
    
    BonePart(1).effTotal.obj=effobj;
    BonePart(1).effTotal.stat=statseff;
    BonePart(1).effTotal.objconn=objTotal.*effobj;
    
    BonePart(1).effInner.obj=effinobj;
    BonePart(1).effInner.objconn=effinobj.*objinner;
    BonePart(1).effInner.stat=statseffinner;
    
    BonePart(1).effCallus.obj=effoutobj;
    BonePart(1).effCallus.objconn=effoutobj.*objcallus;
    BonePart(1).effCallus.stat=statseffcal;
    
    BonePart(1).BreakType=BT;   
       
       
savename=strcat(str, '_BonePart')

save(savename,'BonePart','-v7.3')

clearvars -except fn name folder Table fN
end