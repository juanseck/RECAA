%  Reversible Cellular Automata Algortihm (RCAA)  
%
%  Developed in MATLAB R2015a(7.08)                                                                   
%                                                                                                     
%  Author and programmer: Juan Carlos Seck Tuoh Mora
%                                                                                                     
%       email:   jseck@uaeh.edu.mx
%                juanseck@gmail.com                                                             
%                                                                                                     
%       Homepage: 
%                                                                                                     
%  Main paper:                                                                                        
%  A new algorithm inspired on reversible elementary cellular automata for global optimization
%  IEEE Access, DOI: http://
%
% This function containts full information and implementations of the benchmark 
% functions from the literature 
% Codes of every test function can be consulted in:
%
% http://benchmarkfcns.xyz/fcns
% https://www.sfu.ca/~ssurjano/index.html
% https://github.com/luclaurent/optiGTest
% http://infinity77.net/global_optimization/test_functions_nd_C.html
%
% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = benchmark_functions(F)


switch F
    case 'F1'
        % Brown Function
        fobj = @F1;
        lb=-1;
        ub=4;
        dim=30;
        
    case 'F2'
        % Chung Reynolds Function Function
        fobj = @F2;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F3'
        % Exponential Function
        % N dimensional
        % Unimodal
        % Minimo en f(0,0,...,0) = -1
        fobj = @F3;
        lb=-1;
        ub=1;
        dim=30;
        
    case 'F4'
        % Powell Singular 1 Function
        fobj = @F4;
        lb=-4;
        ub=5;
        dim=30;
        
    case 'F5'
        % Powell Singular 2 Function
        fobj = @F5;
        lb=-4;
        ub=5;
        dim=30;
        
    case 'F6'
        % Powell Sum Function
        fobj = @F6;
        lb=-1;
        ub=1;
        dim=30;
        
    case 'F7'
        % Quartic
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=30;

    case 'F8'
        % Ridge Function
        % Parametros  d=1; alpha=0.5;
        % N dimensional
        % Unimodal
        % Minimo en f(lb,0,...,0) = lb
        % http://benchmarkfcns.xyz/benchmarkfcns/ridgefcn.html
        fobj = @F8;
        lb=-5;
        ub=5;
        dim=30;
        
    case 'F9'
        %Schwefel 1.2 Function
        fobj = @F9;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F10'
        %Schwefel 2.20 Function
        % N dimensional
        % Unimodal
        % Minimo en f(0,0,...,0) = 0
        % http://benchmarkfcns.xyz/benchmarkfcns/schwefel220fcn.html
        fobj = @F10;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F11'
        %Schwefel 2.21 Function
        fobj = @F11;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F12'
        %Schwefel 2.22 Function
        fobj = @F12;
        lb=-100;
        ub=100;
        dim=30;
   
    case 'F13'
        %Schwefel 2.23 Function
        % N dimensional
        % Unimodal
        % Minimo en f(0,0,...,0) = 0
        % http://benchmarkfcns.xyz/benchmarkfcns/schwefel223fcn.html
        fobj = @F13;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F14'
        % Sphere Function
        fobj = @F14;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F15'
        % Step 2 Function
        fobj = @F15;
        lb=-100;
        ub=100;
        dim=30;
        
     case 'F16'
         % Sum Squares Function
        fobj = @F16;
        lb=-100;
        ub=100;
        dim=30;    
        
        
     case 'F17'
        % Ackley 1 Function
        fobj = @F17;
        lb=-35;
        ub=35;
        dim=30; 
        
    case 'F18'
        % Alpine 1 Function 
        fobj = @F18;
        lb=-10;
        ub=10;
        dim=30;
        
    case 'F19'
        % Alpine 2 Function
        fobj = @F19;
        lb=0;
        ub=10;
        dim=30;
        
    case 'F20'
        % High Conditioned Elliptic Function
        fobj = @F20;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F21'
        % Cosine Mixture Function
        fobj = @F21;
        lb=-1;
        ub=1;
        dim=30;
        
    case 'F22'
        % Griewank Function
        fobj = @F22;
        lb=-100;
        ub=100;
        dim=30;

    case 'F23'
        % Levy
        fobj = @F23;
        lb=-10;
        ub=10;
        dim=30;
        
    case 'F24'
        % Pathological Function
        fobj = @F24;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F25'
        % Penalty 1
        fobj = @F25;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F26'
        % Penalty 2
        fobj = @F26;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F27'
        % Rastrigin Function
        fobj = @F27;
        lb=-5.12;
        ub=5.12;
        dim=30;
        
    case 'F28'
        % Salomon Function
        fobj = @F28;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F29'
        % Schwefel 2.26 Function
        fobj = @F26;
        lb=-500;
        ub=500;
        dim=30;
        
    case 'F30'
        % Shubert 3 Function
        % N dimensional
        % Multimodal
        % Minimo en f(x) = −29.6733337
        % http://benchmarkfcns.xyz/benchmarkfcns/shubert3fcn.html
        fobj = @F30;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F31'
        % Streched V Sine Wave
        fobj = @F31;
        lb=-10;
        ub=10;
        dim=30;
        
    case 'F32'
        % Wavy Function
        fobj = @F32;
        lb=-pi;
        ub=pi;
        dim=30;
        
    case 'F33'
        % Xin-She Yang Function
        % N dimensional
        % Multimodal
        % Minimo en f(0*) = 0
        % http://benchmarkfcns.xyz/benchmarkfcns/xinsheyangn1fcn.html
        fobj = @F33;
        lb=-5;
        ub=5;
        dim=30;

    case 'F34'
        % Beale
        fobj = @F34;
        lb=-4.5;
        ub=4.5;
        dim=2;
        
    case 'F35'
        % Corana
        fobj = @F35;
        lb=-500; 
        ub=500;  
        dim=4;
        
    case 'F36'
        % Cross-in-Tray Function
        % 2 dimensional
        % Multimodal
        % Minimo en f(±1.349406685353340, ±1.349406608602084) = −2.06261218
        % http://benchmarkfcns.xyz/benchmarkfcns/crossintrayfcn.html
        fobj = @F36;
        lb=-100;
        ub=100;
        dim=2;
        
    case 'F37'
        % Drop-Wave Function
        % 2 dimensional
        % Multimodal
        % Minimo en f(0,0) = −1
        % http://benchmarkfcns.xyz/benchmarkfcns/dropwavefcn.html
        fobj = @F37;
        lb=-100;
        ub=100;
        dim=2;

    case 'F38'
        % Foxholes
        fobj = @F38;
        lb=-65.536;
        ub=65.536;
        dim=2;
        
    case 'F39'
        % Hartman1
        fobj = @F39;
        lb=0;   % 1
        ub=1;  % 3
        dim=3;
        
    case 'F40'
        % Hartman2
        fobj = @F40;
        lb=0;
        ub=1;
        dim=6;
        
    case 'F41'
        % Helical Valley
        fobj = @F41;
        lb=-10;
        ub=10;
        dim=3;
        
    case 'F42'
        % Hump
        fobj = @F42;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F43'
        % Kowalik
        fobj = @F43;
        lb=-5;
        ub=5;
        dim=4;
        
    case 'F44'
        %Matyas Function
        %http://benchmarkfcns.xyz/benchmarkfcns/matyasfcn.html
        fobj = @F44;
        lb=-10;
        ub=10;
        dim=2;
        
    case 'F45'
        % Miele-Cantrell's Function
        fobj = @F45;
        lb=-1.1; 
        ub=1.1;  
        dim=4;
        
    case 'F46'
        % SHEKEL 5
        fobj = @F46;
        lb=0;
        ub=10;
        dim=4;
        
    case 'F47'
        % SHEKEL 7
        fobj = @F47;
        lb=0;
        ub=10;
        dim=4; 
    
    case 'F48'
        % SHEKEL 10
        fobj = @F48;
        lb=0;
        ub=10;
        dim=4;   
    
    case 'F49'
        % Watson
        fobj = @F49;
        lb=-10;
        ub=10;
        dim=6;
        
    case 'F50'
        % Wolfe
        fobj = @F50;
        lb=0;
        ub=2;
        dim=3;       
        
end

end

% F1
% Brown Function
function o = F1(x)
    n = size(x, 2);  
    o = 0;
    x = x .^ 2;
    for i = 1:(n-1)
        o = o + x(:, i) .^ (x(:, i+1) + 1) + x(:, i+1).^(x(:, i) + 1);
    end
end

% F2
% Chung Reynolds Function
function o = F2(x)
o=(sum(x.^2))^2;
end


% F3
% Exponential
function o = F3(x)
   x2 = x .^2;
   o = -exp(-0.5 * sum(x2, 2));
end 

% F4
% Powell Singular 1 Function
function o = F4(x)
d = length(x);
sum = 0;
for ii = 1:(d/4)
	term1 = (x(4*ii-3) + 10*x(4*ii-2))^2;
	term2 = 5 * (x(4*ii-1) - x(4*ii))^2;
	term3 = (x(4*ii-2) - 2*x(4*ii-1))^4;
	term4 = 10 * (x(4*ii-3) - x(4*ii))^4;
	sum = sum + term1 + term2 + term3 + term4;
end
o = sum;
end

% F5
% Powell Singular 2 Function
function o = F5(x)
d = length(x);
sum = 0;
for ii = 2:(d-2)
	term1 = (x(ii-1) + 10*x(ii))^2;
	term2 = 5 * (x(ii+1) - x(ii+2))^2;
	term3 = (x(ii) - 2*x(ii+1))^4;
	term4 = 10 * (x(ii-1) - x(ii+2))^4;
	sum = sum + term1 + term2 + term3 + term4;
end
o = sum;
end

% F6
% Powell Sum Function
function o = F6(x)
    n = size(x, 2);
    absx = abs(x);
    o = 0;
    for i = 1:n
        o = o + (absx(:, i) .^ (i + 1));
    end
end

% F7
% Quartic
function o = F7(x)
dim=size(x,2);
a=1:dim;
o=sum(a.*(x.^4));
end

% F8
% Ridge
function o = F8(x)
    d=1;
    alpha=0.5;
    x1 = x(:, 1);
    o = x1 + d * (sum(x(:, 2:end).^2, 2) ^ alpha);
end 

% F9
% Schwefel 1.2 Function
function o = F9(x)
dim=size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F10
% Schwefel 2.20 Function
function o = F10(x)
    o = sum(abs(x), 2);
end

% F11
% Schwefel 2.21 Function
function o = F11(x)
o=max(abs(x));
end

% F12
% Schwefel 2.22 Function
function o = F12(x)
o=sum(abs(x))+prod(abs(x));
end

% F13
% Schwefel 2.23 Function
function o = F13(x)
    o = sum(x .^10, 2);
end

% F14
% Sphere Function
function o = F14(x)
o=sum(x.^2);
end

% F15
% Step 2 Function
function o = F15(x)
o=sum(floor((x+.5)).^2);
end

% F16
% Sum Squares Function
function o = F16(x)
i=1:length(x);
o=sum((x.^2).*(i));
end

% F17
% Ackley 1 Function
function o = F17(x)
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F18
% Alpine 1 Function
function o = F18(x)
o=sum(abs(x.*sin(x)+0.1*(x)));
end

% F19
% Alpine 2 Function
function o = F19(x)
o = -prod(sqrt(x) .* sin(x), 2);
end

% F20
% Conditioned Elliptic Function
function o = F20(x)
n=length(x);
o=0;
for i=1:n
    o=o+(((10^6)^((i-1)/(n-1)))*(x(i)^2));
end
end

% F21
% Cosine Mixture Function
function o = F21(x)
o=sum(0.1*length(x)) - ((0.1*sum(cos(5*pi*x)))-(sum(x.^2)));
end

% F22
% Griewank Function
function o = F22(x)
dim=size(x,2);
o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F23
% Levy
function o = F23(x)
d = length(x);
for ii = 1:d
	w(ii) = 1 + (x(ii) - 1)/4;
end
term1 = (sin(pi*w(1)))^2;
term3 = (w(d)-1)^2 * (1+(sin(2*pi*w(d)))^2);
sum = 0;
for ii = 1:(d-1)
	wi = w(ii);
        new = (wi-1)^2 * (1+10*(sin(pi*wi+1))^2);
	sum = sum + new;
end
o = term1 + sum + term3;
end


% F24
% Pathological Function
function o = F24(x)
%constants
a=100;
b=1/2;
c=1e-3;
d=ones(1,size(x,2)-1)*0.5;
%evaluation and derivatives
xi=x(:,1:end-1);
xii=x(:,2:end);
%
xA=a*xii.^2+xi.^2;
sqA=sqrt(xA);
sA=sin(sqA);
pA=sA.^2-b;
%
xB=xi.^2-2*xi.*xii+xii.^2;
xB2=xB.^2;
%
pB=c*xB2+1;
%
pt=pA./pB+d;
%
o=sum(pt,2);
end

% F25
% Penalty 1
function o = F25(x)
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F26
% Penalty 2
function o = F26(x)
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F27
% Rastrigin Function
function o = F27(x)
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F28
% Salomon Function
function o = F28(x)
o=1-cos((2*pi)*sqrt(sum(x.^2)))+(0.1*sqrt(sum(x.^2)));
end

% F29
% Schwefel 2.26 Function
function o = F29(x)
o=sum(-x.*sin(sqrt(abs(x))));
end

% F30
% Shubert 3 Function
function o = F30(x)
    n = size(x, 2);
    o = 0;
    for i = 1:n
        for j = 1:5
            o = o + j * sin(((j + 1) * x(:, i)) + j);
        end
    end
end

% F31
% Streched V Sine Wave
function o = F31(x)
o=0;
for i=1:length(x)-1
    o=o+(((x(i)^2+2*x(i+1)^2)^.25) * (1+(sin((50*(x(i)^2+x(i+1)^2))^0.1)^2)));
end
end

% F32
% Wavy Function
function o = F32(x)
%constants
a=1;
b=1/2;
k=10;
nbvar=size(x,2);
%evaluation and derivatives
cx=cos(k*x);
ex=exp(-b*x.^2);
pcx=cx.*ex;
%
o=a-1/nbvar.*sum(pcx,2);
end

% F33
% Xin-She Yang Function
function o = F33(x)    
    n = size(x, 2);
    o = 0;
    for i = 1:n
        o = o + rand * (abs(x(:, i)) .^ i);
    end
end

% F34
% Beale
function o = F34(x)
    n = size(x, 2);
    assert(n == 2, 'Beale''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);
    o = (1.5 - X + (X .* Y)).^2 + ...
             (2.25 - X + (X .* (Y.^2))).^2 + ...
             (2.625 - X + (X .* (Y.^3))).^2;
end

% F35
% Corana
function o = F35(x)
%constants
a=0.05;
b=0.2;
c=0.15;
d=[1,1e3,10,1e2];
e=0.49999;
%evaluation and derivatives
zz=b*floor(abs(x/b)+e).*sign(x);
vv=abs(x-zz);
maskV=(vv<a);
%
pzz=(zz+(-a*sign(zz)));
pz=(c*pzz.^2).*d;
%
pd=d.*x.^2;
%
pS=zeros(1,size(x,2));
pS(maskV)=pz(maskV);
pS(~maskV)=pd(~maskV);
%
o=sum(pS,2);
end

% F36
% Cross-in-Tray Function
function o = F36(x)
    n = size(x, 2);
    assert(n == 2, 'The Cross-in-tray function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);
    expcomponent = abs(100 - (sqrt(X .^2 + Y .^2) / pi));
    o = -0.0001 * ((abs(sin(X) .* sin(Y) .* exp(expcomponent)) + 1) .^ 0.1);
end

%F37
% Drop-Wave Function 
function o = F37(x)
    n = size(x, 2);
    assert(n == 2, 'Drop-Wave function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);
    numeratorcomp = 1 + cos(12 * sqrt(X .^ 2 + Y .^ 2));
    denumeratorcom = (0.5 * (X .^ 2 + Y .^ 2)) + 2;
    o = - numeratorcomp ./ denumeratorcom;
end

% F38
% Foxholes
function o = F38(x)
aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
for j=1:25
    bS(j)=sum((x'-aS(:,j)).^6);
end
o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F39
% Hartman1
function o = F39(x)
aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F40
% Hartman2
function o = F40(x)
aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
cH=[1 1.2 3 3.2];
pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F41
% Helical Valley
function o = F41(x)
%constants
a=100;
b=10;
c=1;
d=1/(2*pi);
%variables
xxx=x(1);
yyy=x(2);
zzz=x(3);
%evaluation and derivatives
th=atan2(yyy,xxx)*d;
pa=yyy-b*th;
pb=sqrt(xxx.^2+yyy.^2)-c;
%
o=a*(pa.^2+pb.^2)+zzz.^2;
end

% F42
% Hump
function o = F42(x)
o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F43
% Kowalik
function o = F43(x)
aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

%F44
%Matyas Function
%http://benchmarkfcns.xyz/benchmarkfcns/matyasfcn.html
function [o] = F44(x)
    X = x(:, 1);
    Y = x(:, 2);
    o = 0.26 * (X .^ 2 + Y.^2) - 0.48 * X .* Y;
end

% F45
% Miele-Cantrell's Function
function o = F45(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
o = ((exp(x1)-x2)^4)+(100*(x2-x3)^6)+((tan(x3-x4))^4)+(x1^8);
end

% F46
% SHEKEL 5
function o = F46(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
o=0;
for i=1:5
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F47
% SHEKEL 7
function o = F47(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
o=0;
for i=1:7
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F48
% SHEKEL 10
function o = F48(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
o=0;
for i=1:10
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F49
% Watson
function o = F49(x)
%constants
o=x(1)^2;
for i=0:29
    a=i/29;
    aux2=0;
    for j=0:4
        aux2=aux2+((j+1)*(a^j)*x(j+2));
    end
    aux3=0;
    for j=0:5
        aux3=aux3+((a^j)*x(j+1));
    end
    aux4=(aux2-aux3^2-1);
    o=o+(aux4)^2;
end
end

% F50
% Wolfe
function o = F50(x)
    n = size(x, 2);
    assert(n == 3, 'The Wolfe function is defined only on the 3-D space.')
    X = x(:, 1);
    Y = x(:, 2);
    Z = x(:, 3);
    o = (4/3)*(((X .^ 2 + Y .^ 2) - (X .* Y)).^(0.75)) + Z;
end 

function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end



















