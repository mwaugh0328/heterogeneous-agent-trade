function [tssdmat, rtausd, ssd, err_var]=waughtgt_sim(trdshrs,distance,border,talkblock,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program uses the approach from Waugh (2009) on the data set from EK
% (2002). Using the same exact data it estimates S's and trade costs and
% then computes the implied price indicies and plots them versus wages. 
% 
% To run this program, load tgtdata.mat and replicatedata.mat in that order
% into the work space, then run the following:
%
% waughtgt(trdshrs,tradex,distance,border,talkblock,w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fist creates the dummy variable matrix for S's to be recovered, this
% matrix is called dummy
theta = theta;

N = 19;

xrow = -ones(1,N-1);
eee=eye(N-2,N-2);

dummy = [-ones(N-2,1) eee  ];
xrow(1,1)=-2;
dummy = [dummy;xrow];

size(dummy);
for i = 1:N-2
    xrow = -ones(1,N-1);
    dum = [eee(:,1:i), -ones(N-2,1), eee(:,i+1:end)];
    xrow(1,i+1) = -2 ;
    dum = [dum;xrow];
    size(dum);
    dummy = [dummy;dum];
end

dummy = [dummy; ones(N-1,N-1) + eye(N-1,N-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This then creates the dummy variable matrix for the exporter fixed effects
% to be recvoerd, this matrix is called asym

asym = eye(N-1,N-1);
for i = 1:N-1
    ass = eye(N-1,N-1);
    ass(i,:) = [];
    ass = [-ones(1,N-1);ass];
    asym = [asym;ass];
end
asym = rot90(asym,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organize the data, distance and other varibles in the appropriate way and
% then run least squares

trdata = (trdshrs);
hh = (exp(trdata(:))~=0 & exp(trdata(:))~=1);
qq = (exp(trdata(:))~=0);
ff = (exp(trdata(:))~=1);

d_nz = distance(ff,:);
b_nz = border(ff);

trdata = trdata(ff);
gg = (exp(trdata(:))~=0);

trdata = trdata(gg);
dummy = dummy(gg,:);
asym = asym(gg,:);

d = distance(hh,:);
b = border(hh);
tb = talkblock;


% trdata = trdshrs;
% qq = (trdata(:)~=0);
% trdata = trdata(qq);
% d = distance(qq,:);
% b = border(qq);
% tb = talkblock;

[bsd,bintsd,rsd,rintsd,statssd] = regress((trdata),[dummy asym d b tb]);
disp(bsd)
% disp('R-squared')
% disp(statssd(1))

err_var = statssd(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct Fitted Values
yfit =exp([dummy asym d b ]*bsd);

msdmat = zeros(N^2,1);
msdmat(hh) = yfit;
msdmat(~ff) = 1;
msdmat = reshape(msdmat,N,N);
sdc = sum(msdmat);
sdc = repmat(sdc,N,1);
tssdmat = msdmat./sdc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This just computes bootstrap standard errors, these will be different
% than EK because they use an alternative approach.

% disp('Running Bootstrap')
% 
% resid = trdata - yfit;
% straps = 1000;
% 
% [btg]=bootstrp(straps, @(bootr) regress(yfit+bootr,[dummy asym d b tb ]),resid);
% 
% seSr = .01.*round(100*std([(-sum(btg(1:straps,1:N-1)'))',btg(1:straps,1:N-1)]));
% 
% seDer = .01.*round(100*std([(-sum(btg(1:straps,N:end-10)'))',btg(1:straps,N:end-10)]));
% seBr = .01.*round(100*std([btg(1:straps,end-3)]));
% seTr = .01.*round(100*std([btg(1:straps,end-2:end)]));
% seDr = .01.*round(100*std([btg(1:straps,end-9:end-4)]));
% 
% disp('Boot Strap Complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reports Results --- compare to Table III in EK 2002 
% notice that except for the S_is all the esimates are exactly the same. 
% But the trade cost estimates have a different meaning.

ssd = [bsd(1:N-1);-sum(bsd(1:N-1))];

bsasym = bsd(N:end-7);

bsasym = bsd(N:end-7);
bsasym = [bsd(N:end-7);-sum(bsasym)];
beasym = exp(-theta.*bsasym)-1;
arr = 100.*beasym;

dspd = exp(-theta.*bsd(end-6:end-1))-1;

befsd = exp(-theta.*bsd(end))-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('Theta')
% disp(theta)
% 
% disp('Distance Effects')
% dpc = 100.*(exp(-theta.*Dr)-1);
% disp(.01.*round(100*dpc))
% 
% disp('Exporter Fixed Effect')
% depc = 100.*(exp(-theta.*Der)-1);
% disp(.01.*round(100*depc))
% 
% disp('Border Effects')
% bpc = 100.*(exp(-theta.*Br)-1);
% disp(.01.*round(100*bpc))
% 
% disp('Talk-Trade Effects')
% tpc = 100.*(exp(-theta.*Tr)-1);
% disp(.01.*round(100*tpc))
% 
% vvv = 032878;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrange estimates to recover the trade cost matrix 

beasym=(beasym(:,:)*ones(1,N));
beasym = beasym(ff);

tau = zeros(N^2,1);
tau(~ff) = 1;
tau(ff) = (1+d_nz*dspd).*(1+b_nz.*befsd).*(1+beasym);
tau = reshape(tau,N,N);
rtausd = tau;

% disp('log((Xni/Xn)/(Xii/Xi))-log((Xin/Xi)/(Xnn/Xn))')
% disp('US and Denmark')
% disp(tradex(5,19)-tradex(19,5))
% disp('Difference in Arrival Effects')
% disp(-best(19) + best(5))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To see the observable implications for my approach relative to EK 2002,
% here are esimates of the price indicies plotted versus wages in each
% country.

% expshat = (exp(ssd));
% pmat = repmat(expshat, 1,N);
% ipmat = pmat.*tau.^(-1./theta);
% price = ipmat.^-theta;
% price = sum(price);
% 
% axes1 = axes('Parent',figure,...
% 'FontWeight','bold','FontSize',16);
% 
% hold('all');
% plot(log2(exp(wages)), log2(price./price(end)),'bo','LineWidth',3)
% xlabel('Log 2 Relative Mfg. Wages, U.S. = 0','FontWeight','bold','FontSize',16);
% ylabel('Log 2 Model Implied Price of Tradeables,  U.S. = 0','FontWeight','bold','FontSize',16);






