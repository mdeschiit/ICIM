%@Author:   Michael Desch
%@Date:     10/18/2021
%@Project:  Sensitivity Analysis of S using Optimization

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
Su = 61;

%Extracting Data to variable vectors
S = FisherData(:,1);
Smin = FisherData(:,2);
N = 1000*FisherData(:,3);
n = length(FisherData);

%Adjust mean stress to zero using Gerber Eq
Smean = Smin + S./2;
Sadj = S.*(1./(1-((Smean/Su).^2)));

%Transform to log-log space
logS = log(Sadj);
logN = log(N);

%Mean values
Xbar = sum(logS)/n;
Ybar = sum(logN)/n;

%Intermediate Sums
SSx = sum((logS-Xbar).^2);
Sxy = sum((logS-Xbar).*(logN-Ybar));

%Calculating Regression Coefficients
b1 = Sxy/SSx;
b0 = Ybar - b1*Xbar;

%Calculating Residuals
e = logN - (b0 + b1*logS);

%Calculate Determinant parameters C & m
CIP = exp(b0);
mIP = -1*b1;
disp("C & m (IP) = " + CIP + " & " + mIP);

%Calculate interval C bounds
ClIP = exp(b0+min(e));
CuIP = exp(b0+max(e));
disp("C Interval (IP) = [" + ClIP + " ," + CuIP + "]");
 

%%Part 2%%


%Stress time history Field Data
sth = [9.98664781200000,10.3284979700000,9.98664781200000,10.5871373500000,10.3557577000000,10.6885113400000,10.4992739200000,10.9589314500000,10.8289854500000,11.0319123800000,10.7924800700000,11.0319123800000,10.8289854500000,11,13.2459453900000,12.7000647900000,14.9972191100000,13.2459453900000,15.8932709400000,10.2204432600000,10.9793911700000,9.98664781200000];
%Cycle Size threshold
cst = 0.1;
stressrangedataIP = A3(A2(A1(sth), cst));

SminDet = stressrangedataIP(:,2);
SmaxDet = stressrangedataIP(:,3);
SaDet = (SmaxDet-SminDet)./2;
    
%Percent variation in stress range intervals +1
for k = 1:26

frac = (k-1)/100;
SminL = SminDet-frac.*SaDet;
SminU = SminDet+frac.*SaDet;
SmaxL = SmaxDet-frac.*SaDet;
SmaxU = SmaxDet+frac.*SaDet;
disp("Interval Min Stress:");
disp("[ "+SminL+", "+SminU+" ]");
disp("Interval Max Stress:");
disp("[ "+SmaxL+", "+SmaxU+" ]");

nn = 200;
m = 2.9004;
Su = 61;
x0 = zeros(21,1);
LB = zeros(21,1);
LB(1,1) = ClIP;
LB(2:11,1) = SminL;
LB(12:21,1) = SmaxL;
UB = zeros(21,1);
UB(1,1) = CuIP;
UB(2:11,1) = SminU;
UB(12:21,1) = SmaxU;

fun1 = @(x)sum(nn.*((x(12:21)-x(2:11))./(1-((x(12:21)+x(2:11))./(2*Su)).^2)).^m)./x(1);
fun2 = @(x)1-((sum(nn.*((x(12:21)-x(2:11))./(1-((x(12:21)+x(2:11))./(2*Su)).^2)).^m)./x(1)));

options = optimoptions('fmincon','Algorithm','interior-point','OptimalityTolerance',1.0e-20,'Display','iter','StepTolerance',1.0e-20,'ConstraintTolerance',1.0e-20,'MaxFunctionEvaluations',50000);

DdLx = fmincon(fun1,x0,[],[],[],[],LB,UB,[],options);
DdUx = fmincon(fun2,x0,[],[],[],[],LB,UB,[],options);

DdL = sum(nn.*((DdLx(12:21)-DdLx(2:11))./(1-((DdLx(12:21)+DdLx(2:11))./(2*Su)).^2)).^m)/DdLx(1);
DdU = sum(nn.*((DdUx(12:21)-DdUx(2:11))./(1-((DdUx(12:21)+DdUx(2:11))./(2*Su)).^2)).^m)/DdUx(1);

%Problem parameters
r = 0.01;       %Cycle growth rate / year
td = 1/365;     %Length of sample period in years
ta = 20;        %Age of detail at time of evaluation
ts = 18;        %Age of detail at 

%Mean Damage
DmIntGerberLIP = DdL/td;
DmIntGerberUIP = DdU/td;
disp("Interval Mean Damage (IP): ["+DmIntGerberLIP+" ,"+DmIntGerberUIP+"]");

%Existing Damage
DeIntGerberLIP = DmIntGerberLIP*(((1+r)^(ta+1-ts)-(1+r)^(1-ts))/r);
DeIntGerberUIP = DmIntGerberUIP*(((1+r)^(ta+1-ts)-(1+r)^(1-ts))/r);
disp("Interval Existing Damage (IP): ["+DeIntGerberLIP+" ,"+DeIntGerberUIP+"]");

%Total Life
tlIntGerberLIP = log(((r/DmIntGerberUIP)*(1+r)^(ts-1)+1))/log(1+r);
tlIntGerberUIP = log(((r/DmIntGerberLIP)*(1+r)^(ts-1)+1))/log(1+r);
disp("Interval Total Life (IP): ["+tlIntGerberLIP+", "+tlIntGerberUIP+"]");

%Remaining Life
rlIntGerberLIP = tlIntGerberLIP - ta;
rlIntGerberUIP = tlIntGerberUIP - ta;
disp("Interval Remaining Life (IP): ["+rlIntGerberLIP+", "+rlIntGerberUIP+"]");
LBS(k)=rlIntGerberLIP;
UBS(k)=rlIntGerberUIP;
end
toc;

%Stress time history analysis algorithms

%Algorithm 1: Primary Data compression
%Return local minima and maxima of stress time history
function ret = A1(sth)

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

end

%Algorithm 2: Secondary Data Compression
%Remove all cycles smaller than threshold from stress time history
function ret = A2(sth, cst)

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

end

%Algorithm 3: Stress Range Identification
%Will return a [(n-1)/2 x 4] matrix including stress ranges, minimum 
%stresses, and maximum stresses & full or half cycle from a compressed 
%input stress time history. Runs in O(n^2) time, possibly able to rewrite 
%to run in O(n) time
function ret = A3(sth)

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

end