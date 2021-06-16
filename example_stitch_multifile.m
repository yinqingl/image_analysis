% automatic stitching of FISH images
%
% use DAPI as stitching reference
% based on SIFT descriptor for matching
% also keep track of the dapi intensity as a mask for later processing
%
% Yinqing Li
% yinqingl@csail.mit.edu
% 2016

%%
%initialization

%analyze fish image
addpath('imread')
addpath('helper_function')
addpath('IAT_v0.9.2')
addpath('TSP')
addpath('bfmatlab')

%{
ATTENTION: run the following code once
iat_setup('forever')
%}
iat_setup('forever')

%%
%list nd2 files
path = '..\Slide1\'; %system path to the images
path = strrep(path,'\',filesep);
datname = [path,'20161021_Ecel1T4_Spp1T1_Slide1_Section1_TRNL.nd2'];

%list all nd2 files
datname_s = {};
path_sets = {'..\Slide7\'}; %system path to the images
for path_i = 1:length(path_sets),
    path = path_sets{path_i};
    path = strrep(path,'\',filesep);
    listing = dir(path);
    for i = 1:length(listing),
        filename = listing(i).name;
        if length(filename)>3 && strcmp(filename(end-3:end),'.nd2'),
            datname_s{end+1} = [path,filename];
        end
    end
end

%print
for i_img = 1:length(datname_s)
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname))
end


%%
%configuration

%set reference channel
ref_channel = 1;

%filter setup
h_stdfilt = ones(5); %std filter size

is_z_halfdepth_filter_used = 0; %use dapi to filter z stacks, set to 0 for weak dapi
z_halfdepth_around_ref_pl = 3; %the half depth of z stacks containing signal
med2filter_n = 5;   %number of x grids for z stack filter
med2filter_m = 5;   %number of y grids for z stack filter

%background adjustment
% is_d_bk_sel = {'all'};
% is_d_bk_sel = {'on',[1]}; %only adjust background for the selected channels
is_d_bk_sel = {'off',[1]}; %do not adjust background for the selected channels


%%
%process all image files
%series are splitted in different images files

%read meta file to determine series coordinates
i_img_s = [34,35];
pos_xy = [];
i_img_pos_xy = [];
s_pos_xy = [];
for i_img = i_img_s,
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname));

% datname = [path,'2016-10-3_Slide1_Section2_TRN_3Color.nd2'];
% meta = imreadBFmeta_nd2(datname);
[T,meta] = evalc('imreadBFmeta_nd2(datname);');

meta

width = meta.width;
height = meta.height;
zsize = meta.zsize;
nframes = meta.nframes;
nseries = meta.nseries;
nchannels = meta.channels;

pos_xy = [pos_xy; double(cell2mat(meta.result(:,2)))];
i_img_pos_xy = [i_img_pos_xy; ones(nseries,1)*i_img];
s_pos_xy = [s_pos_xy, [1:nseries]];
end

pos_xy = pos_xy + rand(size(pos_xy))*min(pos_xy(:))/1e4;
pos_xy = pos_xy + rand(size(pos_xy))*min(pos_xy(:))/1e4;
%determine the order of series for alignment based on the x,y coordinates
pos_xy_dist = dist(pos_xy', 'euclidean');
UG = sparse(pos_xy_dist);
[~,s_tour_pred] = graphminspantree(UG);
%BFS on series tour
s_tour = [];
queue = [find(s_tour_pred==0)];
while(~isempty(queue)),
    s_tour = [s_tour, queue(1)];
    queue(1) = [];
    queue = [queue, find(s_tour_pred == s_tour(end))];
end

% figure; scatter(pos_xy(:,1),pos_xy(:,2),20,[1:size(pos_xy,1)],'filled')
% view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))
% view(biograph(s_tour,[],'ShowArrows','off','ShowWeights','on'))

%%
%stitch images 
if ~is_z_halfdepth_filter_used,
    z_halfdepth_around_ref_pl = zsize; %the half depth of z stacks containing signal
    med2filter_n = 1;   %number of x grids for z stack filter
    med2filter_m = 1;   %number of y grids for z stack filter    
end

img_mask_s = {};

% for s = s_tour,
for s = s_tour(6:end),
    disp(sprintf('loading frame: %d',s));
    
    i_img = pos_xy_i_img(s);
    s_img = s_pos_xy(s);
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname));    
    
    %for each channel
    img = zeros(width, height, nchannels);
    for c = 1:nchannels,
        %read in all z stacks
%         [vol]=imreadBF_nd2(datname,[1:zsize],t,c);
        [T,vol] = evalc('imreadBF_nd2(datname,[1:zsize],s_img,c);');
%         if meta.nd2,
%             [T,vol] = evalc('imreadBF_nd2(datname,[1:zsize],t,c);');
%         else
%             [T,vol] = evalc('imreadBF(datname,[1:zsize],t,c);');
%         end

        % figure;
        % imagesc(vol);
        % colormap(gray);

        %max projection
%         I_zmax = max(vol,[],3);
        %max projection
        if c == ref_channel,
            [I_zmax,idx_zmax] = max(vol,[],3);
            if min(idx_zmax(:)) == max(idx_zmax(:)),
                z_ref_pl = ones(med2filter_n,med2filter_m);
                z_ref_pl(:,:) = min(idx_zmax(:));
            else,
                stdfilt_J = stdfilt(idx_zmax,h_stdfilt); %find nuclei by looking at projected z, because the max z at no nuclei position should be random
                [kmeans_idx,kmeans_c] = kmeans(stdfilt_J(:),2);
                kmeans_idx = reshape(kmeans_idx,height,width);
        %         figure;imagesc(kmeans_idx);
                [~,kmeans_idx_min] = min(kmeans_c);
                z_ref_mask = (kmeans_idx == kmeans_idx_min);
                z_ref_pl = zeros(med2filter_n,med2filter_m);
                n_range = round(linspace(1,height,med2filter_n+1));
                m_range = round(linspace(1,width,med2filter_m+1));
                idx_zmax_msk = idx_zmax.*z_ref_mask;
                for n = 1:med2filter_n,
                    for m = 1:med2filter_m,
                        pl = idx_zmax_msk(n_range(n):n_range(n+1),m_range(m):m_range(m+1));
                        z_ref_pl(n,m) = median(pl(pl>0));
                    end
                end
                %clear stdfilt_J kmeans_idx idx_zmax_msk;
            end
        else
            %for other channels, take max projection within a range of z
            %the range is determined based on the reference channel
            I_zmax = zeros(height,width);
            n_range = round(linspace(1,height,med2filter_n+1));
            m_range = round(linspace(1,width,med2filter_m+1));
            for n = 1:med2filter_n,
                for m = 1:med2filter_m,
                    z_range = round([max(z_ref_pl(n,m)-z_halfdepth_around_ref_pl,1):1:...
                        min(z_ref_pl(n,m)+z_halfdepth_around_ref_pl,zsize)]);
                    vol_nm = vol(n_range(n):n_range(n+1),m_range(m):m_range(m+1),z_range);
                    [I_zmax_nm,idx_zmax] = max(vol_nm,[],3);
                    I_zmax(n_range(n):n_range(n+1),m_range(m):m_range(m+1)) = I_zmax_nm;
                end
            end
        end

        %substract background
        I = squeeze(I_zmax);
        background = imopen(I,strel('disk',25));
        I = I - background;
        I_bk_adj = I;
        img(:,:,c) = I_bk_adj;
    end
    %clear vol;
    %clear I_zmax I_bk_adj;
    
    

    %for the first position in the series, skip image alignment
    if s==s_tour(1),
        img_ref = img;
        img_ref_mask = ones(size(img_ref(:,:,ref_channel)));
        z_ref_mask_ref = z_ref_mask;
        
        img_mask_s{s} = img_ref_mask;
        continue
    end
    
    s_pred = s_tour_pred(s);
    img_ref_mask = img_mask_s{s_pred};
    
    %register two images
    try
        [t,tform] = get_sift_affine_t(img(:,:,ref_channel),img_ref_mask.*img_ref(:,:,ref_channel));
    catch
        display(sprintf('cannot find matching between series %d and %d', s, s_pred));
        %update s_pred
        s_tour_pred(s_tour_pred == s) = s_pred;
        s_tour
        s_tour_pred(s_tour)
        continue
    end

    %construct a canvas for the reference and the image
    [~,xdata,ydata] = imtransform(img(:,:,ref_channel), tform,...
        'UData',[1,width], 'VData',[1,height]);
    xdata = ceil(abs(xdata)).*sign(xdata);
    ydata = ceil(abs(ydata)).*sign(ydata);
    t_ref = [1 0 0; 0 1 0; -min(0,xdata(1)), -min(0,ydata(1)), 1];
    t.T(3,1) = t.T(3,1) + t_ref(3,1);
    t.T(3,2) = t.T(3,2) + t_ref(3,2);
    xdata_max = max(size(img_ref,2),xdata(2)) + t_ref(3,1);
    ydata_max = max(size(img_ref,1),ydata(2)) + t_ref(3,2);

    Rcb = imref2d([ydata_max,xdata_max]);
    t_ref = affine2d(t_ref);

    %apply transformation to all channels and construct composite image
    img_w = imwarp(img,t,'outputview',Rcb);
    img_ref_w = imwarp(img_ref,t_ref,'outputview',Rcb);
    
    %update series mask
    for ss = find(cellfun(@(x) ~isempty(x), img_mask_s)),
        img_mask_s{ss} = imwarp(img_mask_s{ss},t_ref,'outputview',Rcb);
    end
    img_mask_s{s} = img_w(:,:,ref_channel) ~= 0;
    
    % figure, imshowpair(img_w(:,:,1), img_ref_w(:,:,1))
    
    %adjust background based on the difference in the overlapping area
    %the background area is found by excluding the dapi high area
    img_mask_w = imwarp(~z_ref_mask,t,'outputview',Rcb);
    img_mask_ref_w = imwarp(~z_ref_mask_ref,t_ref,'outputview',Rcb);
    img_overlap_mask = (img_mask_w + img_mask_ref_w)==2;
%     figure; imagesc(img_w(:,:,2).*img_overlap_mask)
%     figure; imagesc(img_ref_w(:,:,2).*img_overlap_mask)
    
    for c = 1:nchannels,
        img_w_masked = img_w(:,:,c).*img_overlap_mask;
        img_ref_w_masked = img_ref_w(:,:,c).*img_overlap_mask;
        
%         med_bk = median(img_w_masked(img_w_masked>0));
%         med_bk_ref = median(img_ref_w_masked(img_ref_w_masked>0));        
        
        img_w_bk = img_w_masked(img_w_masked>0);
        img_ref_w_bk = img_ref_w_masked(img_ref_w_masked>0);
        
        %norm using kde on linear scale hist local maximum
        I_bk = img_w_bk;
        I_th_5 = quantile(I_bk,0.95);
        [~,n,x_bk,~] = kde((I_bk(I_bk<I_th_5)),256); n_bk = n/sum(n);
        [~,x_bk_i] = max(n_bk);
        med_bk = (x_bk(x_bk_i));
        
        I_bk = img_ref_w_bk;
        I_th_5 = quantile(I_bk,0.95);
        [~,n,x_bk,~] = kde((I_bk(I_bk<I_th_5)),256); n_bk = n/sum(n);
        [~,x_bk_i] = max(n_bk);
        med_bk_ref = (x_bk(x_bk_i));
        
%         figure; plot(x_bk, n_bk); hold on; plot(x_bk,n_bk,'r'); 
        
%         %norm on log scale hist local maximum
%         [n,x] = hist(log(img_w_bk),256);
%         n_bk = histc(log(img_w_bk),x); n_bk = n_bk/sum(n_bk);
%         n_bk_ref = histc(log(img_ref_w_bk),x); n_bk_ref = n_bk_ref/sum(n_bk_ref);
% %         figure; plot(x,n_bk); hold on; plot(x,n_bk_ref,'r');
%         [~,x_bk_i] = max(n_bk);
%         med_bk = exp(x(x_bk_i));
%         [~,x_bk_i] = max(n_bk_ref);
%         med_bk_ref = exp(x(x_bk_i));
        
        d_bk = med_bk - med_bk_ref;
        
        if ~strcmp(is_d_bk_sel{1},'all'),
            if strcmp(is_d_bk_sel{1},'on'),
                if ~ismember(c,is_d_bk_sel{2})
                    d_bk = 0;
                end
            elseif strcmp(is_d_bk_sel{1},'off'),
                if ismember(c,is_d_bk_sel{2})
                    d_bk = 0;
                end
            end
        end
        
%         d_bk = (img_w_masked./img_ref_w_masked);
%         d_bk_mask = img_w_masked==0 | img_ref_w_masked==0;
%         d_bk = d_bk(~d_bk_mask);
% %         figure; hist(log(d_bk),256);
%         [n,x] = hist(log(d_bk),256);
%         [~,i_x] = max(n);
%         d_bk = x(i_x);
        
        %reduce background
        if d_bk > 0,
            img_d_bk = zeros(size(img_w_masked));
            img_d_bk = d_bk;
%             img_d_bk(img_mask_w) = d_bk;
            img_w(:,:,c) = max(img_w(:,:,c) - img_d_bk,0);
%             img_w(:,:,c) = max(img_w(:,:,c)/exp(img_d_bk),1);
%             img_ref_w(:,:,c) = max(img_ref_w(:,:,c)*exp(img_d_bk),1);
        else
            img_d_bk = zeros(size(img_ref_w_masked));
            img_d_bk = -d_bk;
%             img_d_bk(img_mask_ref_w) = -d_bk;            
            img_ref_w(:,:,c) = max(img_ref_w(:,:,c) - img_d_bk,0);
%             img_ref_w(:,:,c) = max(img_ref_w(:,:,c)/exp(img_d_bk),1);
%             img_w(:,:,c) = max(img_w(:,:,c)*exp(img_d_bk),1);
        end
    end
    
    %stitch two images
    img_c = max(img_w, img_ref_w);
    % figure, imagesc(img_c(:,:,2));
    
    %update reference mask
    img_mask_w = imwarp(z_ref_mask,t,'outputview',Rcb);
    img_mask_ref_w = imwarp(z_ref_mask_ref,t_ref,'outputview',Rcb);
    z_ref_mask_ref = max(img_mask_w, img_mask_ref_w);
    img_ref = img_c; 
    %clear img_c img_w img_ref_w;
end
% figure; imagesc(img_ref(:,:,1)); daspect([1 1 1]);
% figure; imagesc(z_ref_mask_ref(:,:)); daspect([1 1 1]);

% %adjust contrast
% I_int_max = double(intmax('uint16'));
% for c = 1:nchannels,
%     I = double(img_ref(:,:,c));
%     I = I/max(I(:))*I_int_max;
%     img_ref(:,:,c) = I;
% end

%write to tiff
datname = datname_s{i_img_s(1)};
datname_output = [datname '.tif'];
% data = uint16(round(img_c));
data = uint16(round(img_ref));
while size(data,3)<3,   %imagej cannot open tif with less than 3 channels
    data(:,:,end+1) = zeros(size(data(:,:,1)));
end
t = Tiff(datname_output, 'w');
tagstruct.ImageLength = size(data, 1);
tagstruct.ImageWidth = size(data, 2);
tagstruct.Compression = Tiff.Compression.None;
% tagstruct.Compression = Tiff.Compression.LZW;        % compressed
tagstruct.SampleFormat = Tiff.SampleFormat.UInt; 
% tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;                        % float data
tagstruct.SamplesPerPixel = size(data,3);
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
t.setTag(tagstruct);
t.write(data);
t.close();

%write reference channel mask
datname_output = [datname_output '.mat'];
save(datname_output,'ref_channel','z_ref_mask_ref','-v7.3');

