%% Clear all data
clear all
clc
close all

uiopen
% load('percent.mat');

%%
f=percent.full;
[xf yf]=size(f);
full=f;
full(full~=0)=255;

[B,L,N,A] = bwboundaries(full,8);

figure(1)
imshow(f)
hold on

for ii = 1:length(B)
    sB(ii,1:2)=size(B{ii});
end

[bsize bindex]=max(sB(:,1));
outer=B{bindex};

BW_out = poly2mask(outer(:,2),outer(:,1),xf,yf);

plot(outer(:,2),outer(:,1), 'r', 'LineWidth', 2)
hold on
plot(mean(outer(:,2)),mean(outer(:,1)),'.')

% inner ROI
cent1=mean(outer(:,1));
cent2=mean(outer(:,2));
inner(:,1)=outer(:,1)-cent1;
inner(:,2)=outer(:,2)-cent2;

width=input('Percent of ROI width (0 - 100%):  ')
width=(100-width)/100;
inner=inner*width;

inner(:,1)=inner(:,1)+cent1;
inner(:,2)=inner(:,2)+cent2;

BW_in = poly2mask(inner(:,2),inner(:,1),xf,yf);

BW_ROI=BW_out-BW_in;

plot(inner(:,2),inner(:,1), 'g', 'LineWidth', 2)

figure(2)
subplot(1,3,1)
imshow(f)

subplot(1,3,2)
imshow(f)
hold on
plot(outer(:,2),outer(:,1), 'r', 'LineWidth', 2)
plot(mean(outer(:,2)),mean(outer(:,1)),'.')
plot(inner(:,2),inner(:,1), 'g', 'LineWidth', 2)

subplot(1,3,3)
imshow(BW_ROI)


%%

BW_ROI = abs(double(BW_ROI));

area=regionprops(BW_ROI,'area');
area=area.Area;

ppp=BW_ROI.*percent.defect;
ppp=ppp*100;

[s1 s2]=size(ppp);

for i=1:s2
    a=ppp(:,i);
    a=round(a*10^2)/(10^2);
    ppp(:,i)=a;
end



iter=1;
x=[];
val=[];
pr=[];

for i=iter:iter:100
    num = sum(sum(ppp>=i-(iter/2) & ppp<i+(iter/2)));
    per=(num/area)*100;
    x=[x; i];
    val=[val; num];
    pr=[pr; per];
end

figure 
plot(x,val,'.', 'MarkerSize', 10)
set(gca, 'FontSize', 18)
title('Number of Z-Rays vs. Percent Filled', 'FontSize', 24)
xlabel('Percent Z-Ray Fill', 'FontSize', 22)
ylabel('Number of Z-Rays', 'FontSize', 22)
axis([0 100 0 max(val)+100])
grid on

figure
plot(x,pr,'.', 'MarkerSize', 10)
set(gca, 'FontSize', 18)
title('Percent of Z-Rays vs. Percent Filled', 'FontSize', 24)
xlabel('Percent Z-Ray Fill', 'FontSize', 22)
ylabel('Number of Z-Rays', 'FontSize', 22)
xlabel('Percent Z-Ray Fill')
ylabel('Percent of Total Z-Rays')
axis([0 100 0 max(pr)+5])
grid on

area(1:length(x),:)=area;


for i=1:100
zzz=ppp;
zzz(ppp<=i-(iter/2))=0;
zzz(ppp>i+(iter/2))=0;
BI2=zzz;
BI2(zzz>0.001)=1;
BI2=BI2.*zray_length.defect;
total_t_pix(i)=sum(sum(BI2));
[row_t,col_t,v_t] = find(BI2);
avgt(i)=mean(v_t);
st_devt(i)=std(v_t);

yyy=(zzz/100).*zray_length.defect;
yyy=round(yyy);
total_b_pix(i)=sum(sum(yyy));
[row_b,col_b,v_b] = find(yyy);
avgb(i)=mean(v_b);
st_devb(i)=std(v_b);
end

data=[x val area pr total_b_pix' avgb' st_devb' total_t_pix' avgt' st_devt'];
data_cells=num2cell(data);

col_header={'Percent Iteration','Number of Z-Rays','Total Selected Z-Rays','Percent of Region',...
    'Total Bone Pixels','Mean Bone','STD Bone','Total Tissue Pixels','Mean Tissue','STD Tissue'};

matrix_out=[col_header; data_cells];

d=date;
c=clock;
c3=round(c(6));
c1=num2str(c(4));
c2=num2str(c(5));
c3=num2str(c3);
dd=[date '_' c1 '_' c2 '_' c3]

xlswrite(dd,matrix_out)
