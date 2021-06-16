function count_retrobeads()

% addpath('imread')
% addpath('helper_function')
% addpath('IAT_v0.9.2')
% addpath('TSP')
% addpath('bfmatlab')
% addpath('adaptivethreshold')
% addpath('bradley')

s = datestr(now,'yyyy-mm-dd-HH-MM-SS');
str_name = sprintf('quant_c_%s.mat',s);

fprintf('\nthe results will be saved in file: %s\n', str_name)

quant_s = {};

while(1)
    prompt = 'Take a screen shot from slide scanner viewer, enter y or empty to finish? ';
    x = input(prompt,'s');
    if isempty(x),
        break
    end
    q = process_image();
    quant_s(end+1) = {q};
    save(str_name,'quant_s','-v7.3')
    q{1}
    prompt = 'data saved';
end

fprintf('\nthe results are saved in file: %s\n', str_name)


end

function q = process_image()
I = imclipboard('paste');
I = I(:,:,1);
fh_I = figure; imshow(I);
background = imopen(I,strel('disk',25));

% figure; imshow(background)
% figure; imshow(I - background)
% mask = adaptivethreshold(I,[150 150],0,0);
% mask = adaptivethresh(I);
% mask = bradley(I, [125 125], 10);
% mask = adaptthresh(I, 0.4);
% figure; imshow(mask)
% BW = imbinarize(I,mask);
% figure; imshow(BW);

I_bk = I - background;
I_bk_u = unique(I_bk(:));
[n,bin]=histc(I_bk(:),I_bk_u);
% figure; plot(I_bk_u,a,'.');
f_hist = figure; 
plot(I_bk_u,[0,diff(n)'],'.');
f_bk = figure; 
figure(f_hist);
while(1)
    prompt = 'What is the background threshold, enter a number or empty to finish? ';
    x = input(prompt);
    if isempty(x),
        break
    end
    figure(f_bk);
    I_c = I_bk - I_bk_u(x) > 1;
    imshow(I_c==0);
end
xy_s = get_roi(f_bk, I_c==0);

bw_s = {};
[m,n] = size(I_bk);
for i = 1:size(xy_s,1)
    x = xy_s{i,1};
    y = xy_s{i,2};
    bw = poly2mask(x,y,m,n);
    bw_s{i} = bw;
end

q_s = [];
for i = 1:size(xy_s,1)
    bw = bw_s{i};
    c = I_c.*bw;
    q_s(i) = sum(c(:));
end

q = {q_s, bw_s, I_c, I};

close(fh_I)
close(f_hist)
close(f_bk)
end

function xy_s = get_roi(fh, I_c)
display('select roi, left click: set a polygon for roi, d: delete, n:next roi, q:finish')
xy_s = {};

c_line_de = {'b','r','g'}';
c_line  = c_line_de;

xs = [];
ys = [];
figure(fh)
title(sprintf('select roi #%i, select roi, left click: set a polygon for roi, d: delete, n:next roi, q:finish',size(xy_s,1)+1))
while 1,
    [x,y, key] = ginput(1);
    switch key
        case 1  %selecting roi
            xs(end+1) = x;
            ys(end+1) = y;
            c = c_line{1};
            hold on;
            plot(xs,ys,c);
        case 110    %n
            xy_s(end+1,:) = {xs,ys};
            c_line  = c_line_de;
            for i=1:size(xy_s,1)
                c = c_line{1};
                xs_p = xy_s{i,1};xs_p(end+1) = xs_p(1);
                ys_p = xy_s{i,2};ys_p(end+1) = ys_p(1);
                hold on;
                plot(xs_p,ys_p,c)
                c_line = circshift(c_line,1);
            end     
            xs = [];
            ys = [];
            figure(fh)
            title(sprintf('select roi #%i, select roi, left click: set a polygon for roi, d: delete, n:next roi, q:finish',size(xy_s,1)+1))
        case 100    %d
            figure(fh)
            imshow(I_c)
            xy_s = {};
            c_line  = c_line_de;
            
            xs = [];
            ys = [];
            figure(fh)
            title(sprintf('select roi #%i, select roi, left click: set a polygon for roi, d: delete, n:next roi, q:finish',size(xy_s,1)+1))
        case 113    %q
            xy_s(end+1,:) = {xs,ys};
            break
        otherwise
            1;
    end
end

end