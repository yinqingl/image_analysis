%{
take tapestation electrophoresis images
calculate the molar concentration for each sample
%}

img = imread('hsd1k_yl_170204-HSD1000_HS D1000.png');
img = squeeze(img(:,:,1));

addpath('Y:\Shared\Projects\YL_Shared\Nucseq_analysis\Matlab_analysis\common_code')

fh = custom_analysis_functions;

fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

%identify lanes
l = sum(img < 255,1);
l = l>max(l)/3;

il = find(l);
l = kmeans(il,2);

%sort lanes from left to right
il_m = zeros(1,length(unique(l)));
for i = hv(unique(l)),
    il_m(i) = mean(il((l==i)));
end
[~,l_s] = sort(il_m);
l_sorted = l;
for i = hv(unique(l)),
    l_sorted(l==i) = l_s(i);
end
l = l_sorted;

%identify ladder
i = 1;
pos_c = round(quantile(il(l==i),[0.3,0.7]));
e = img(:,pos_c(1):pos_c(2));
e = fliplr(hv(255-mean(e,2)));
mask = zeros(size(e));
mask(100:650) = 1;
e = e.*mask;

[MAXTAB, MINTAB] = peakdet(e, 9);  %50 is hard coded
%{
figure; plot(e);
hold on;
plot(MAXTAB(:,1), MAXTAB(:,2), 'r*');
%}
size_l = [25,50,100,200,300,400,500,700,1000,1500];

%interpolate ladder
Xq = [MAXTAB(1,1):MAXTAB(end,1)];
Vq = interp1(MAXTAB(:,1)',size_l,Xq,'spline');
dVq = diff(Vq);
%{
figure; plot(Vq, e(Xq));
figure; plot(Xq,Vq,'*-');
%}

%iterate through samples
q = zeros(size(unique(l)));
for i = hv(unique(l)),
    if i == 1,
        continue
    end
    
pos_c = round(quantile(il(l==i),[0.3,0.7]));
e = img(:,pos_c(1):pos_c(2));
e = fliplr(hv(255-mean(e,2)));

%the gel electrophoresis shows the total amount
%convert total amount to density
e = e(Xq);
e = e./[dVq(1) dVq];

%find a dilution to 4nM library
%1nM is 1pmol dna in 1ml of solution
%1pmol dna of 1kbp is 653.69ng
%calculate the amount of sample for 10ul of 4nM library, which is 40fmol

%select a range of sample to calculate effective molecular weight
figure; plot(Vq, e);
[x,~] = ginput(2); 
ylim = get(gca,'ylim');
hold on;
plot([x(1),x(1)],ylim,'r')
plot([x(2),x(2)],ylim,'r')

f = @(x) interp1(Vq,e,x,'linear');
g = @(x) f(x)./(x/1000*0.65369);
q(i) = integral(g,x(1),x(2))/integral(f,x(1),x(2)); %q fmol/ng
title(sprintf('sample %d, molar concentration: %f fmol/ng',i,q(i)));
end

%update concentration
c = [0,26.9];

%show dilution
v = zeros(size(q));
display('to prepare 10ul of 4nm sample, the total sample needed is 40fmol')
for i = hv(unique(l)),
    if i == 1,
        continue
    end
    v(i,1)=40/(c(i)*q(i));
    v(i,2)=10-v(i,1);
end
v
