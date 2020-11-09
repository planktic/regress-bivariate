%% Regress_Bivariate
%% 1. About
%
% *Regress_Bivariate*: 
% 
% This code produces a weighted least squares fit of a straight line to
% a set of points with error in both coordinates. It can handle bivariate
% regression where the errors in both coordinates are correlated and is
% capable of performing force-fit regression.
%
% For details on the algorithm used in the code, see _Unified Equations 
% for the slope, intercept, and standard errors of the best straight line_ 
% (<http://is.gd/jVA5fE York et al., 2004>)
%
% Written by Kaustubh Thirumalai (kau@ig.utexas.edu) at the University of
% Texas Institute for Geophysics in July 2009. Latest Update: 08/28/2014
% Kaustubh is currently at Brown University (kaustubh_thirumalai@brown.edu)
%
% Citation: Thirumalai, K., A. Singh, and R. Ramesh (2011), _A MATLAB code 
% to perform weighted linear regression with (correlated or uncorrelated) 
% errors in bivariate data_, Journal of the Geological Society of India, 
% 77(4), 377–380, doi: <http://is.gd/sk1hMu 10.1007/s12594-011-0044-1>.
%
%=====================================
% INPUT
%=====================================
% Excel file/array with the following format:
% 
% Column 1 | Column 2 | Column 3 | Column 4 | Column 5    |
% X data   | sigX     | Y data   | sigY     | r^2(sigmas) |
% 
%=====================================
% OUTPUT
%=====================================
% 
% •	Graph of Y vs. X data (with errorbars)
% •	Comparison of Weighted Linear Regression vs. Simple Linear Regression
% •	Errors on slope and intercept of the line of best-fit
% 
% NOTE: Correlation coefficient between errors in X(i) and Y(i) has been
% taken as zero in this program. If there is any correlation between the
% errors, then you can include it directly in the worksheet in column 5
% and you can de-comment it in the code below. By default it will be taken
% as zero. For force-fit regression, place a sufficiently low error (i.e.
% ~infinite weight) on the coordinates of the point at which the line is 
% to be forced.

%% 2. Input
% A .xlsx (.xls) file with the following format:
% 
% * Column 1 - |X Data|
% * Column 2 - |X Error|
% * Column 3 - |Y Data|
% * Column 4 - |Y Error|
% * Column 5 - |Corr.Coef. between errors (if needed)|
%
% Otherwise, similar input for an array in |mat|

tic

clc;
clear;

tol = 1e-8;    % Tolerance

%---- Excel Data Extraction ----

mat = xlsread('Pearson.xls');

% Note that data can be entered directly as an array in |mat|

X = mat(:,1);
sigX = mat(:,2);
Y = mat(:,3);
sigY = mat(:,4);

%% 3. Simple Linear Regression
% Performs simple linear regression (SLR) for comparison.

R = corrcoef(X,Y);
n = length(X);

ri = 0;
X1 = [n,sum(X);sum(X),sum(X.*X)];             % LSE Value of 'b'
Y1 = [sum(Y);sum(X.*Y)];                      % Polyfit can also be used
Z1 = X1\Y1;                                     % Matrix division
a1 = Z1(1);
b1 = Z1(2);
sigres = sum((Y - a1 - b1.*X).^2)/(n-2);        % Sigma Residual
delta = det(X1);                                % Determinant
varx = var(X);
sigb1sq = (n^2)*(n-1)*varx*sigres/(delta^2);    % Sigma(b) without weights
sigb1 = sqrt(sigb1sq);
Xbar = mean(X);
siga1sq = (sigres/(varx*n))*(varx + Xbar^2);
siga1 = sqrt(siga1sq);
A_SLR = [a1 siga1];                             % SLR Intercept
B_SLR = [b1 sigb1];                             % SLR Slope

%---- SLR p-value ----

B1 = 0;
t1=(b1-B1)/sigb1;
Pval1=2*(1-tcdf(abs(t1),1));

%% 4. Maximum Likelihood Method (Weighted Linear Regression)
% Implements the _York et al._ [2004] algorithm to perform weighted linear
% regression using the maximum likelihood method (MLM).

%---- Weighting Errors ----

wX = abs(1./(sigX.^2));                     
wY = abs(1./(sigY.^2));
alpha = sqrt(wX.*wY);
%ri = mat(:,5);         % CorrCoef between sig(X) and sig(Y)

b = b1;                             
d = tol;                            
i = 0; 

%---- York et al. (2004) Algorithm ----

while (d > tol || d == tol)         % Tolerance check loop
    i = i+1;
    b2 = b;
    W = wX.*wY./((wX) + ((b^2).*wY) - (2.*b.*alpha.*ri));
    meanX = sum(W.*X)/sum(W);
    meanY =  sum(W.*Y)/sum(W);
    U = X(:) - meanX;
    V = Y(:) - meanY;
    Beta = W.*((U./wY)+((b.*V)./wX) - (b.*U + V).*(ri./alpha));
    meanBeta = sum(W.*Beta)/sum(W);
    b = sum(W.*Beta.*V)/sum(W.*Beta.*U);
    dif = b - b2;
    d = abs(dif);
end

U2 = U.^2;
V2 = V.^2;
a = meanY - b.*meanX;
x = meanX + Beta;
meanx = sum(W.*x)/sum(W);
u = x - meanx;
sigbsq = 1./(sum(W.*(u.*u)));
sigb = sqrt(sigbsq);
sigasq = 1./(sum(W)) + meanx^2.*(sigbsq);
siga = sqrt(sigasq);
S = sum(W.*((Y - b.*X - a)).^2);

%---- MLE p-value ----
B = 0;
t=(b-B)/sigb;
Pval=2*(1-tcdf(abs(t),length(X)-2));

%% 5. Plot

figure(1)
clf;
hold on;
errorbar_x(X,Y,sigX,'o');
set(gca,'FontSize',12,'FontName','Myriad Pro');
% See http://www.mathworks.com/matlabcentral/fileexchange/12751)
h1 = errorbar(X,Y,sigY,'o','markersize',...
    6,'linewidth',1); grid on;   
p = min(X)-(max(X)-min(X))/10:(max(X)-min(X))/10:max(X)+(max(X)-min(X))/10;
q = a + (b.*p);
h2 = plot(p,q,'linewidth',2);
p = min(X)-(max(X)-min(X))/10:(max(X)-min(X))/10:max(X)+(max(X)-min(X))/10;
q = a1 + (b1.*p);
h3 = plot(p,q,'--');
pl = legend([h1,h2,h3],'Data','Weighted Line (MLM)',...
    'Simple Line (SLR)','Location','NorthEast');
htext=findobj(get(pl,'children'),'type','text');
set(htext,'FontName','Avenir','Fontsize',14);
xl = xlabel('X');
set(xl,'FontName','Avenir','FontSize',14);
yl = ylabel('Y');
set(yl,'FontName','Avenir','FontSize',14);

%% 6. Output

display(' ');
display('Regress_Bivariate:');
display(' ');

D = [X sigX Y sigY];
display('       X       sigX       Y        sigY');
disp(D);

A = [a siga];
B = [b sigb];
disp('--------------------------------------');
display('Weighted Linear Regression:');disp(' ');
wr = sum(U.*V)./sqrt((sum(U2).*sum(V2)));
display('Y = a + bX');display(' ');
disp('       a      siga');disp(A);
disp('       b      sigb');disp(B);
disp(['r: ',num2str(wr)]);
disp(['r^2: ',num2str(wr^2)]);
disp(['p-value: ',num2str(Pval)]);

disp('--------------------------------------');
display('Simple Linear Regression:');disp(' ');
display('Y = a + bX');display(' ');
disp('       a      siga');disp(A_SLR);
disp('       b      sigb');disp(B_SLR);
disp(['r: ',num2str(R(1,2))]);
disp(['r^2: ',num2str(R(1,2)^2)]);
disp(['p-value: ',num2str(Pval1)]);
disp('--------------------------------------');
disp(['Total no. of iterations: ',num2str(i)]);
toc
