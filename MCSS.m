%@Author:   Michael Desch
%@Date:     1/6/2021
%@Project:  Sensitivity Analysis on S Monte Carlo Simulation

%Housekeeping
clc;
clear;
close all;
tic;


%%Part 1%%


%S-N test data from Fisher et.al. 1969
%Tests on full scale A-36 cover plate beams, ends welded
%Col 1: Stress range (ksi)
%Col 2: Minimum stress (ksi)
%Col 3: Cycles to Failure (kilo-cycles)
FisherData = [16,-6,392.500000000000;16,-6,393.300000000000;16,-6,336.700000000000;20,-6,192.200000000000;20,-6,168.100000000000;20,-6,288.200000000000;20,-6,176.100000000000;24,-6,114.400000000000;24,-6,93.7000000000000;24,-6,85;12,2,797.700000000000;12,2,654.500000000000;12,2,724.300000000000;16,2,276.900000000000;16,2,316.500000000000;16,2,328.600000000000;16,2,325;20,2,197.700000000000;20,2,159;20,2,147.800000000000;8,10,2227.40000000000;8,10,2693.10000000000;8,10,2453.20000000000;12,10,675.600000000000;12,10,777.600000000000;12,10,657.800000000000;12,10,738.600000000000;16,10,300.700000000000;16,10,344.100000000000;16,10,297.200000000000;20,10,107.700000000000;20,10,180.300000000000;20,10,172;20,10,166];
%Ultimate Strength 
Su = 61*6.89476;

%Extracting Data to variable vectors
S = FisherData(:,1)*6.89476;
Smin = FisherData(:,2)*6.89476;
N = 1000*FisherData(:,3);
n = length(FisherData);

%Adjust mean stress to zero using Gerber Eq
Smean = Smin + S./2;
SadjGerber = S.*(1./(1-((Smean/Su).^2)));

%Transform to log-log space
logSGerber = log(SadjGerber);
logN = log(N);

%Mean values
XbarGerber = sum(logSGerber)/n;
Ybar = sum(logN)/n;

%Intermediate Sums
SSxGerber = sum((logSGerber-XbarGerber).^2);
SxyGerber = sum((logSGerber-XbarGerber).*(logN-Ybar));

%Calculating Regression Coefficients
b1Gerber = SxyGerber/SSxGerber;
b0Gerber = Ybar - b1Gerber*XbarGerber;

%Calculating Residuals
eGerber = logN - (b0Gerber + b1Gerber*logSGerber);

%Calculate Determinant parameters C & m
CGerberIP = exp(b0Gerber);
mGerberIP = -1*b1Gerber;
disp("C & m (Gerber) = " + CGerberIP + " & " + mGerberIP);

%Calculate interval C bounds
ClGerberIP = exp(b0Gerber+min(eGerber));
CuGerberIP = exp(b0Gerber+max(eGerber));
disp("C Interval (Gerber) = [" + ClGerberIP + " ," + CuGerberIP + "]");

%Stress time history Field Data (Imperial)
sth = [9.98664781200000,10.3284979700000,9.98664781200000,10.5871373500000,10.3557577000000,10.6885113400000,10.4992739200000,10.9589314500000,10.8289854500000,11.0319123800000,10.7924800700000,11.0319123800000,10.8289854500000,11,13.2459453900000,12.7000647900000,14.9972191100000,13.2459453900000,15.8932709400000,10.2204432600000,10.9793911700000,9.98664781200000]*6.89476;
%Cycle Size threshold
cst = 0.1*6.89476;
stressrangedataIP = A3(A2(A1(sth), cst));

SminDet = stressrangedataIP(:,2);
SmaxDet = stressrangedataIP(:,3);
SDet = SmaxDet-SminDet;
SmDet = (SmaxDet+SminDet)./2;
SaDet = (SmaxDet-SminDet)./2;

SapDetGerber = SaDet./(1-((SmDet./Su).^2));
SpDetGerber = 2*SapDetGerber;
ADTT = 200;
r = 0.01;
td = 1/365;
ta = 20;
ts = 18;

trGerbermin = zeros(26,1);
trGerbermax = zeros(26,1);
for k = 1:26
    frac = (k-1)/100;
    for i = 1:10000000

        %Random C in interval bounds
        CGerber(1,i) = rand(1,1)*(CuGerberIP-ClGerberIP)+ClGerberIP;

        for j = 1:10
            Sminf(j,i) = rand(1,1)*((SminDet(j)+frac*SaDet(j))-(SminDet(j)-frac*SaDet(j)))+(SminDet(j)-frac*SaDet(j));
            Smaxf(j,i) = rand(1,1)*((SmaxDet(j)+frac*SaDet(j))-(SmaxDet(j)-frac*SaDet(j)))+(SmaxDet(j)-frac*SaDet(j));
            Sf(j,i) = Smaxf(j,i)-Sminf(j,i);
            Smf = (Smaxf(j,i)+Sminf(j,i))/2;
            Saf = (Smaxf(j,i)-Sminf(j,i))/2;
            SpfGerber(j,i) = 2*(Saf/(1-((Smf/Su)^2)));
        end

        DdGerber(1,i) = ADTT*sum(SpfGerber(:,i).^mGerberIP)/CGerber(1,i);
        DmGerber(1,i) = DdGerber(1,i)/td;
        DeGerber(1,i) = DmGerber(1,i)*(((1+r)^(ta+1-ts)-(1+r)^(1-ts))/r);
        tlGerber(1,i) = log(((r/DmGerber(1,i))*(1+r)^(ts-1)+1))/log(1+r);
        trGerber(1,i) = tlGerber(1,i)-ta;

    end

    Sminmin = min(Sminf,[],2);
    Sminmax = max(Sminf,[],2);
    Smaxmin = min(Smaxf,[],2);
    Smaxmax = max(Smaxf,[],2);
    Sfmin = min(Sf,[],2);
    Sfmax = max(Sf,[],2);
    SpfGerbermin = min(SpfGerber,[],2);
    SpfGerbermax = max(SpfGerber,[],2);

    DdGerbermin = min(DdGerber);
    DdGerbermax = max(DdGerber);
    %disp("MCS Dd Interval (Gerber) = [" + DdGerbermin + " ," + DdGerbermax + "]");

    DmGerbermin = min(DmGerber);
    DmGerbermax = max(DmGerber);
    %disp("MCS Dm Interval (Gerber) = [" + DmGerbermin + " ," + DmGerbermax + "]");

    DeGerbermin = min(DeGerber);
    DeGerbermax = max(DeGerber);
    %disp("MCS De Interval (Gerber) = [" + DeGerbermin + " ," + DeGerbermax + "]");

    tlGerbermin = min(tlGerber);
    tlGerbermax = max(tlGerber);
    %disp("MCS tl Interval (Gerber) = [" + tlGerbermin + " ," + tlGerbermax + "]");

    trGerbermin(k,1) = min(trGerber);
    trGerbermax(k,1) = max(trGerber);
    %disp("MCS tr Interval (Gerber) = [" + trGerbermin + " ," + trGerbermax + "]");

    %disp("COV = " + std(trGerber)/mean(trGerber));
end
toc;

%Stress time history analysis algorithms

%Algorithm 1: Primary Data compression
%Return local minima and maxima of stress time history
function ret = A1(sth)
    tic;
    keep = sth(1);
    buffer = sth(2);
    for i = 3:length(sth)
        next = sth(i);
        if(buffer == keep(length(keep)))
            buffer = next;
        elseif(buffer > keep(length(keep)) && next >= buffer)
            buffer = next;
        elseif(buffer > keep(length(keep)) && next < buffer)
            keep(length(keep)+1) = buffer;
            buffer = next;
        elseif(buffer < keep(length(keep)) && next > buffer)
            keep(length(keep)+1) = buffer;
            buffer = next;
        elseif(buffer < keep(length(keep)) && next <= buffer)
            buffer = next;
        end
    end
    keep(length(keep)+1) = buffer;
    ret = keep;
    toc;
end

%Algorithm 2: Secondary Data Compression
%Remove all cycles smaller than threshold from stress time history
function ret = A2(sth, cst)
    tic;
    keep = sth(1);
    for i = 2:length(sth)
        buffer = max([sth(i-1) sth(i)])-min([sth(i-1) sth(i)]);
        if(buffer >= cst)
            keep(length(keep)+1) = sth(i);
        elseif(buffer < cst && i > 2 && i < length(sth))
            keep = keep(1:(length(keep)-1));
        end
    end
    ret = keep;
    toc;
end

%Algorithm 3: Stress Range Identification
%Will return a [(n-1)/2 x 4] matrix including stress ranges, minimum 
%stresses, and maximum stresses & full or half cycle from a compressed 
%input stress time history. Runs in O(n^2) time, possibly able to rewrite 
%to run in O(n) time
function ret = A3(sth)
    tic;
    S = [];
    asth = sth;
    cycle = 1;
    for j = 1:((length(sth))/2)
        ranges = zeros(1,length(asth)-1);
        for i = 2:length(asth)
           ranges(1,i-1) = max([asth(i-1) asth(i)])-min([asth(i-1) asth(i)]);
        end
        [mv,minid] = min(ranges);
        if(length(asth) == 2)
            cycle = 0.5;
        end
        S(j,[1 2 3 4]) = [mv min([asth(minid) asth(minid+1)]) max([asth(minid) asth(minid+1)]) cycle];
        asth([minid minid+1])=[];
    end
    ret = S;
    toc;
end