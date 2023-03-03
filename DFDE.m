%% Flood Depth Estimation in Agricultural Lands from L and C-Band SAR Images and Digital Elevation Model%%
%  DATE OF CREATION OF FILE: 23/03/2022 
%  AUTHOR: SAMVEDYA SURAMPUDI @Microwave lab, VIT University %%

%########################DESCRIPTION##########################################################
% INPUTS:              C-Band data (Sentine-1A)
%                      L-Band data (ALOS-2 PALSAR-2)
%                      Digital Elevation Model (SRTM, 1arc Sec)
%      Derived data:   L-Band derived NDFI flood boundary (Binary raster)
                    
%#############################################################################################
clc
clear all
% Data sets
% Reading Sentinel-1A and ALOS-2 data 
sent=imread('Sentinel.tif'); % Read C-band data
%sent=sent(1:674,1:1037); (Clip to same dimensions)
alos=imread('ALOS2.tif'); % Read L-band data
alos=10*log10(alos.^2)-83; % L-band Calibration
dem=imread('DEM.tif'); % Read DEM
dem10=dem(:,:,2);

% Flood boundary extraction
% Wavelet image fusion and segmentation
alos1=alos;
sent1=sent;
figure(1); imshow(alos1,[]); impixelinfo
figure(2);imshow(sent1,[]); impixelinfo
xfus2 = wfusimg(sent1,alos1,'haar',5,'mean','mean');
figure(3);imshow(xfus2,[]); colormap bone; impixelinfo; 
thresh = multithresh(double(xfus2),5);
seg_I = imquantize(double(xfus2),thresh);
figure(4);imshow(seg_I,[]); impixelinfo;

depth1=(seg_I==1); depth1=double(depth1);
depth2=(seg_I==2); depth2=double(depth2);
depth3=(seg_I==3); depth3=double(depth3);
depth4=(seg_I==4); depth4=double(depth4);
depth5=(seg_I==5); depth5=double(depth5);
% Finding all elevation (dem10) that corresponds to extent in each boundary
dem10_lev1=dem10(find(depth1==1));
dem10_lev2=dem10(find(depth2==1));
dem10_lev3=dem10(find(depth3==1));
dem10_lev4=dem10(find(depth4==1));
dem10_lev5=dem10(find(depth5==1));
dem10=dem10(:);
% clipping DEM10 in flood extent
d_lev1=zeros(length(dem10),1);
d_lev1(find(depth1==1))=dem10_lev1;
d_lev1=reshape(d_lev1,[734,1037]);
d_lev1(d_lev1==0)=NaN;
figure(1)
imshow(d_lev1,[]);colormap jet; impixelinfo;

d_lev2=zeros(length(dem10),1);
d_lev2(find(depth2==1))=dem10_lev2;
d_lev2=reshape(d_lev2,[734,1037]);
imshow(d_lev2,[]);colormap jet; impixelinfo;

d_lev3=zeros(length(dem10),1);
d_lev3(find(depth3==1))=dem10_lev3;
d_lev3=reshape(d_lev3,[734,1037]);
imshow(d_lev3,[]);colormap jet; impixelinfo;

d_lev4=zeros(length(dem10),1);
d_lev4(find(depth4==1))=dem10_lev4;
d_lev4=reshape(d_lev4,[734,1037]);
imshow(d_lev4,[]);colormap jet; impixelinfo;

d_lev5=zeros(length(dem10),1);
d_lev5(find(depth5==1))=dem10_lev5;
d_lev5=reshape(d_lev5,[734,1037]);
imshow(d_lev5,[]);colormap jet; impixelinfo;
%% 
% Section -2
%###############################################################################################

[r1,c1] = ind2sub(size(d_lev1),find(d_lev1>0));
un=unique(c1);
bounds=cell(length(un),1);
for i=1:length(un)
  bounds{i,1}=r1(find(c1==un(i)));
 bounds{i,2}=c1(find(c1==un(i)));
end

for i=1:length(bounds)

   block{i} = d_lev1(bounds{i,1},bounds{i,2});  
   block{i}=block{i}(:);
  
end
%########################################

for i=1:length(block)
re{i}=reshape((block{1,i}(1:length(block{1,i}))),size(bounds{i},1),size(bounds{i},1));
end

re2=d_lev1;
for i=1:length(block)
re2(bounds{i,1},bounds{i,2})=re{i};
end

%########################################

for i=1:length(block)

p1=100:-(100/length(block{i})):0; 
cs1=block{i};
pertil11 = prctile(cs1,p1);

% Finding out flood depth. 
for j=length(p1)-6:-1:1
    if (pertil11(length(p1)-j)-pertil11(length(p1)-5-j))<30
        depth11(j)=(pertil11(length(p1)-j)+ pertil11(length(p1)-5-j))/2;
    end
end

depth11(length(depth11)+1:length(pertil11))=0;
flooddepth1{i}=pertil11-depth11(1:length(pertil11));% Y vector
flooddepth1{i}=flooddepth1{i}(:);

end


% Reshaping flood vectors
for i=1:length(flooddepth1)
new_fd{i}=reshape((flooddepth1{1,i}(1:length(flooddepth1{1,i})-1)),size(bounds{i},1),size(bounds{i},1));
end

%writing values into d_lev1 matrix
fd_fin1=zeros(size(d_lev1));
for i=1:length(bounds)

   fd_fin1(bounds{i,1},bounds{i,2})=new_fd{i};
     
end
figure(2)
imshow(fd_fin1,[]);colormap jet; impixelinfo;

% Creating of flood depth values

for i=1:length(flooddepth1)
    
    x{i}=flooddepth1{1,i};
    xx{i}=0:(flooddepth1{1,i}(1))/(length(flooddepth1{1,i})-1):(flooddepth1{1,i}(1));
    pd_cs1{i} = fitdist(double(xx{i})','Normal');
    pd{i} = makedist('Normal','mu',pd_cs1{i}.mu,'sigma',pd_cs1{i}.sigma);
    y{i}=pdf(pd{i},xx{i});
    inv_dat{i} = icdf('Normal',y{i},pd_cs1{i}.mu,pd_cs1{i}.sigma);
    inv_dat{i}=inv_dat{i}';
     
end


% Reshaping flood vectors
for i=1:length(flooddepth1)
fdepth_fin1{i}=reshape((inv_dat{1,i}(1:length(inv_dat{1,i})-1)),size(bounds{i},1),size(bounds{i},1));
end

%writing values into Layer_lev1 matrix
Layer_lev1=zeros(size(d_lev1));
for i=1:length(bounds)

   Layer_lev1(bounds{i,1},bounds{i,2})=-(fdepth_fin1{i});
     
end
figure(3)
imshow(Layer_lev1);colormap jet; impixelinfo;

Lay_L1_write=-Layer_lev1';
% Interpolation

depths1=Lay_L1_write(find(p1>50&p1<60));
depths2=Lay_L2_write(find(p2>50&p2<54));
depths3=Lay_L3_write(find(p3>50&p3<54));
depths4=Lay_L4_write(find(p4>50&p4<54));
depths5=Lay_L5_write(find(p5>50&p5<54));

points1=54:-((54-50)/length(dem1_req)):50;
points2=54:-((54-50)/length(dem2_req)):50;
points3=54:-((54-50)/length(dem3_req)):50;
points4=54:-((54-50)/length(dem4_req)):50;
points5=54:-((54-50)/length(dem5_req)):50;

y1=interp1(Lay_L1_write(find(p1>50&p1<60)),depths1,points1,'linear');
y2=interp1(Lay_L2_write(find(p2>50&p2<54)),depths2,points2,'linear');
y3=interp1(Lay_L3_write(find(p3>50&p3<54)),depths3,points3,'linear');
y4=interp1(Lay_L4_write(find(p4>50&p4<54)),depths4,points4,'linear');
y5=interp1(Lay_L5_write(find(p5>50&p5<54)),depths5,points5,'linear');

% Writing the depth values
level=dem10(:);
level(dem1_ind1)=y1(1:length(dem1_ind1));
level(dem2_ind1)=y2(1:length(dem2_ind1));
level(dem3_ind1)=y3(1:length(dem3_ind1));
level(dem4_ind1)=y4(1:length(dem4_ind1)); 
level(dem5_ind1)=y5(1:length(dem5_ind1));

shapewrite(fld_shp, 'level.shp')
shapewrite(fld_shp, 'level2.shp')
shapewrite(fld_shp, 'level3.shp')
shapewrite(fld_shp, 'level4.shp')
shapewrite(fld_shp, 'level5.shp')

% Note:
% Repeat section- 2 for all levels from level-1(L1) to level-5 L5 to get
% Lay_L()_write
Lay_fin_write=Lay_L1_write+Lay_L2_write+Lay_L3_write+Lay_L4_write+Lay_L5_write; % Final layer
Lay_fin=Lay_fin_write(:);
% Visulaization for single sample 
figure(3)
surf(new_fd{1})
figure(4)
surf(-fdepth_fin1{50})
shading flat
figure(5)
surf(d_lev1(bounds{1,1},bounds{1,2}))
figure(6)
surf((Layer_lev1(bounds{50,1},bounds{50,2})))
qqplot(fdepth_fin1{50})


