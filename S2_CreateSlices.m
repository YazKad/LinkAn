[s p]=size(rotated2.slice(1).pic);
[m n]=size(rotated2.slice);

for j=1:1:s
for i=1:n
cut(:,i)=(rotated2.slice(i).pic(j,:));
end

imshow(cut)
% pause(0.1)

filename=num2str(j);
filename=[filename '.bmp'];

imwrite(cut,filename)


end
