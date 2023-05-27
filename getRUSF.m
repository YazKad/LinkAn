% getRUSF
% Written by Yazan Kadkoy 5/2018
% This program takes _BonePart files and scores them according to mRUST scoring system
%Place all _BonePat files you wish to be scored into one folder and select
%that folder when prompted


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

Total=["Name";"BreakType";"cort 1";"Cort 2";"Cort 3";"Cort 4"; "view 1";"view 2";"total";"force"];
                    
for i=1:num
    FN=char(fN(i,:));
    load (FN)
 
    name=BonePart(1).name;
    bt=BonePart(1).BreakType;
    
%     [~,~,zee]=size(obj);
    cx=round(BonePart(1).Centroid(1));
    offsetx=0;
    cy=round(BonePart(1).Centroid(2));
    offsety=0;
    os2=25;
    % creates 4 objects representing the 4 corticies of bone. This is done
    % for the Inner object, Inner object>26, Callus object and Callus
    % object >26
    Bridge(1).Inner(1).obj=BonePart(1).Inner.obj(cy-os2:cy+os2,1:cx-offsetx,5:end-5);
    Bridge(1).Inner(1).stat=regionprops3(Bridge(1).Inner(1).obj,'Volume','VoxelList');
    Bridge(1).Inner(1).size=size(Bridge(1).Inner(1).obj);
    Bridge(1).Inner(2).obj=BonePart(1).Inner.obj(cy-os2:cy+os2,cx+offsetx:end,5:end-5);
    Bridge(1).Inner(2).stat=regionprops3(Bridge(1).Inner(2).obj,'Volume','VoxelList');
    Bridge(1).Inner(2).size=size(Bridge(1).Inner(2).obj);
    Bridge(1).Inner(3).obj=BonePart(1).Inner.obj(1:cy-offsety,cx-os2:cx+os2,5:end-5);
    Bridge(1).Inner(3).stat=regionprops3(Bridge(1).Inner(3).obj,'Volume','VoxelList');
    Bridge(1).Inner(3).size=size(Bridge(1).Inner(3).obj);
    Bridge(1).Inner(4).obj=BonePart(1).Inner.obj(cy+offsety:end,cx-os2:cx+os2,5:end-5);
    Bridge(1).Inner(4).stat=regionprops3(Bridge(1).Inner(4).obj,'Volume','VoxelList');
    Bridge(1).Inner(4).size=size(Bridge(1).Inner(4).obj);
     

    Bridge(1).Inner22(1).obj=BonePart(1).Inner.objconn(cy-os2:cy+os2,1:cx-offsetx,5:end-5)>=26;
    Bridge(1).Inner22(1).stat=regionprops3(Bridge(1).Inner22(1).obj,'Volume','VoxelList');
    Bridge(1).Inner22(1).size=size(Bridge(1).Inner22(1).obj);
    Bridge(1).Inner22(2).obj=BonePart(1).Inner.objconn(cy-os2:cy+os2,cx+offsetx:end,5:end-5)>=26;
    Bridge(1).Inner22(2).stat=regionprops3(Bridge(1).Inner22(2).obj,'Volume','VoxelList');
    Bridge(1).Inner22(2).size=size(Bridge(1).Inner22(2).obj);
    Bridge(1).Inner22(3).obj=BonePart(1).Inner.objconn(1:cy-offsety,cx-os2:cx+os2,5:end-5)>=26;
    Bridge(1).Inner22(3).stat=regionprops3(Bridge(1).Inner22(3).obj,'Volume','VoxelList');
    Bridge(1).Inner22(3).size=size(Bridge(1).Inner22(3).obj);
    Bridge(1).Inner22(4).obj=BonePart(1).Inner.objconn(cy+offsety:end,cx-os2:cx+os2,5:end-5)>=26;
    Bridge(1).Inner22(4).stat=regionprops3(Bridge(1).Inner22(4).obj,'Volume','VoxelList');
    Bridge(1).Inner22(4).size=size(Bridge(1).Inner22(4).obj);
    
    
    Bridge(1).Callus(1).obj=BonePart(1).Callus.obj(cy-os2:cy+os2,1:cx-offsetx,5:end-5);
    Bridge(1).Callus(1).stat=regionprops3(Bridge(1).Callus(1).obj,'Volume','VoxelList');
    Bridge(1).Callus(1).size=size(Bridge(1).Callus(1).obj);
    Bridge(1).Callus(2).obj=BonePart(1).Callus.obj(cy-os2:cy+os2,cx+offsetx:end,5:end-5);
    Bridge(1).Callus(2).stat=regionprops3(Bridge(1).Callus(2).obj,'Volume','VoxelList');
    Bridge(1).Callus(2).size=size(Bridge(1).Callus(2).obj);
    Bridge(1).Callus(3).obj=BonePart(1).Callus.obj(1:cy-offsety,cx-os2:cx+os2,5:end-5);
    Bridge(1).Callus(3).stat=regionprops3(Bridge(1).Callus(3).obj,'Volume','VoxelList');
    Bridge(1).Callus(3).size=size(Bridge(1).Callus(3).obj);
    Bridge(1).Callus(4).obj=BonePart(1).Callus.obj(cy+offsety:end,cx-os2:cx+os2,5:end-5);
    Bridge(1).Callus(4).stat=regionprops3(Bridge(1).Callus(4).obj,'Volume','VoxelList');
    Bridge(1).Callus(4).size=size(Bridge(1).Callus(4).obj);
    Bridge(3).size=size(Bridge(1).Callus(1).obj);
    
    Bridge(1).Callus22(1).obj=BonePart(1).Callus.objconn(cy-os2:cy+os2,1:cx-offsetx,5:end-5)>=26;
    Bridge(1).Callus22(1).stat=regionprops3(Bridge(1).Callus22(1).obj,'Volume','VoxelList');
    Bridge(1).Callus22(1).size=size(Bridge(1).Callus22(1).obj);
    Bridge(1).Callus22(2).obj=BonePart(1).Callus.objconn(cy-os2:cy+os2,cx+offsetx:end,5:end-5)>=26;
    Bridge(1).Callus22(2).stat=regionprops3(Bridge(1).Callus22(2).obj,'Volume','VoxelList');
    Bridge(1).Callus22(2).size=size(Bridge(1).Callus22(2).obj);
    Bridge(1).Callus22(3).obj=BonePart(1).Callus.objconn(1:cy-offsety,cx-os2:cx+os2,5:end-5)>=26;
    Bridge(1).Callus22(3).stat=regionprops3(Bridge(1).Callus22(3).obj,'Volume','VoxelList');
    Bridge(1).Callus22(3).size=size(Bridge(1).Callus22(3).obj);
    Bridge(1).Callus22(4).obj=BonePart(1).Callus.objconn(cy+offsety:end,cx-os2:cx+os2,5:end-5)>=26;
    Bridge(1).Callus22(4).stat=regionprops3(Bridge(1).Callus22(4).obj,'Volume','VoxelList');
    Bridge(1).Callus22(4).size=size(Bridge(1).Callus22(4).obj);
        
    
    line=zeros(7,1);
    count=1;
    for x=1:4
        %creates imaginary plate on either end of break site
        plate1=zeros(Bridge(1).Inner(x).size);
    plate1(:,:,1)=1;
    plate1(:,:,end-5:end)=1;
    
    plate2=zeros(Bridge(1).Inner22(x).size);
    plate2(:,:,1)=1;
    plate2(:,:,end-5:end)=1;
    plate3=zeros(Bridge(1).Callus(x).size);
    plate3(:,:,1)=1;
    plate3(:,:,end-5:end)=1;
    plate4=zeros(Bridge(1).Callus22(x).size);
    plate4(:,:,1)=1; 
    plate3(:,:,end-5:end)=1;    
        
        in=0;
        in22=0;
        cal=0;
        cal22=0;
        % Add in poritons of inner and callus to determine if they span
        % fracture site. for each portion added. Volum fraction is checked.
        % must be greater than 85% to be considered in calculation
       [~,lc1]=max(Bridge(1).Inner(x).stat.Volume);
       perin=max(Bridge(1).Inner(x).stat.Volume)/sum(Bridge(1).Inner(x).stat.Volume);
        objtemp1=plate1;
            
            vox2=Bridge(1).Inner(x).stat.VoxelList{lc1,1};
             inds2=sub2ind(Bridge(1).Inner(x).size,vox2(:,2),vox2(:,1),vox2(:,3));
            objtemp1(inds2)=1;
            objtemp1=objtemp1>=1;
            
            stattemp=regionprops3(objtemp1,'Volume');
            numt=numel(stattemp.Volume);
           if(perin>.85) 
            if(numt==1)
             in=1;
            end
           end
           
         [~,lc2]=max(Bridge(1).Inner22(x).stat.Volume);
         perin22=max(Bridge(1).Inner22(x).stat.Volume)/sum(Bridge(1).Inner22(x).stat.Volume);
        objtemp2=plate2;
            
            vox2=Bridge(1).Inner22(x).stat.VoxelList{lc2,1};
             inds2=sub2ind(Bridge(1).Inner22(x).size,vox2(:,2),vox2(:,1),vox2(:,3));
            objtemp2(inds2)=1;
            objtemp2=objtemp2>=1;
            
            stattemp=regionprops3(objtemp2,'Volume');
            numt=numel(stattemp.Volume);
           if(perin22>.85) 
            if(numt==1)
             in22=1;
            end
           end
         [~,lc3]=max(Bridge(1).Callus(x).stat.Volume);
         percal=max(Bridge(1).Callus(x).stat.Volume)/sum(Bridge(1).Callus(x).stat.Volume);
        objtemp3=plate3;
            
            vox2=Bridge(1).Callus(x).stat.VoxelList{lc3,1};
             inds2=sub2ind(Bridge(1).Callus(x).size,vox2(:,2),vox2(:,1),vox2(:,3));
            objtemp3(inds2)=1;
            objtemp3=objtemp3>=1;
            
            stattemp=regionprops3(objtemp3,'Volume');
            numt=numel(stattemp.Volume);
            if(percal>.85)
            if(numt==1)
             cal=1;
            end
            end
      [~,lc4]=max(Bridge(1).Callus22(x).stat.Volume);
      percal22=max(Bridge(1).Callus22(x).stat.Volume)/sum(Bridge(1).Callus22(x).stat.Volume);
        objtemp4=plate4;
            
            vox2=Bridge(1).Callus22(x).stat.VoxelList{lc4,1};
             inds2=sub2ind(Bridge(1).Callus22(x).size,vox2(:,2),vox2(:,1),vox2(:,3));
            objtemp4(inds2)=1;
            objtemp4=objtemp4>=1;
            
            stattemp=regionprops3(objtemp4,'Volume');
            numt=numel(stattemp.Volume);
            if(percal22>.85)
            if(numt==1)
             cal22=1;
            end
            end
            
            if(in22&& cal22)
                line(x,1)=4;
            elseif(in && cal)
                line(x,1)=3;
            elseif(~in&& cal)
                line(x,1)=2;
            elseif(~in&&~cal)
                line(x,1)=1;
            else
                line(x,1)=0;
            end
            
            
            
   
    end
   line(5,1)=line(1,1)+line(2,1);
   line(6,1)=line(3,1)+line(4,1);
   line(7,1)=line(5,1)+line(6,1);
    
   look=find(forces(1,:)==name);
   next=forces(2,look);
    
    col=[name;bt;line;next];
    
    Total=[Total col];
    
 
    
end
   %perform regression
   
    two=str2double(Total(end,2:end));
    ncol=["R^2";"--"];
for q=3:numel(Total(:,1))
    one=str2double(Total(q,2:end));
    
   r=regression(one,two);
    
    ncol=[ncol;r];
     
 
 
end
%plot points
Total=[Total ncol];



a=str2double(Total(end-1,2:end-1));
b=str2double(Total(end,2:end-1));


ex=str2double(Total(end-1,2:end-1));
ey=str2double(Total(end,2:end-1));
coefficients = polyfit(ex, ey, 1);
xFit = linspace(min(ex), max(ex), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;
scatter(a,b,30,'k','filled')
title('RUSF Score vs Force')
ylabel('Peak Torque (N*mm)')
xlabel('RUSF Score')
