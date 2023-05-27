% % getLink
% % Written by Yazan Kadkoy 9/2017
% % 
% % This program takes _vox files and converts binary matrices into linkage matracies.
% % A linkage matrix consits of linkage values, defined as the number of postive neighbors within the 3X3X3 
% % space around any postive voxel. 




clear all
clc
close all
T=[];


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
iter=num2str(i);
[~,var,~]=fileparts(FN);

file(i).name(:)=char(var);

total(i).bone=bone;
end
for ii=1:num
   bone=total(ii).bone;
[R C]=size(bone(1).vox);
rotatedbone=repmat(struct('vox',zeros(R,C)),numel(bone),1);

for z=2:(numel(bone)-1)
    con=0;
    tot=0;
    
    for x=1:R
        for y=1:C
            %check next vox for 1 and 3X3 array
            %exceptions for going out of bounds
            %spit out least connectend and highest C/J if D>thresh
            connect=0;
                  if(bone(z).vox(x,y)==1)
                      
            xx=x;
            yy=y;
            xc=1;
            yc=1;
            for q=xx-1:xx+1
             if((q>0)&&(q<=R))
                 xs(xc)=q;
                 xc=xc+1;
             end
            end
             for p=yy-1:yy+1
             if((p>0)&&(p<=C))
                 ys(yc)=p;
                 yc=yc+1;
             end
            end
          for zz=z-1:z+1
              found=find(bone(zz).vox(xs(:),ys(:)));

              if(~isempty(found))
              connect=connect+numel(found);
              end
                  
          end
              
              
              
                  end 
         rotatedbone(z).vox(x,y)=connect;
         
        end
    end
    
 
end 
rotatedbone(1).name=bone(1).name;
savename=[bone(1).name '_Link.mat']
save(savename,'rotatedbone')

clearvars -except fN FN fN1 ii total

end
    
    
    
    