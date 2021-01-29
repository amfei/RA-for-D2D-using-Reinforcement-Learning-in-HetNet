

function model = modelHetNet(I,C,D,K,CR,crp,dd,pathloss_parameter,FREQ,BW,BOLTZ)

%% users location
%cell R and BS
xBS=[0 300 -150 -150];
yBS=[0  0   260 -260];


%circle
tf = linspace(0, 2*pi, 100);
xr = CR*cos(tf) + xBS(1);
yr = CR*sin(tf) + yBS(1);


xr1 = crp*cos(tf) + xBS(2);
yr1 = crp*sin(tf) + yBS(2);

xr2 = crp*cos(tf) + xBS(3);
yr2 = crp*sin(tf) + yBS(3);

xr3 = crp*cos(tf) + xBS(4);
yr3 = crp*sin(tf) + yBS(4);

% Eav in cell
E=2*pi*rand(1,1);
rE  = CR*rand(1,1);
xE =rE.*cos(E) ;
yE = rE.*sin(E);
%CUs in circle cell
centerC = [0 ,0];
anglec = 2*pi*rand(1,C);
rc = CR*sqrt(rand(1,C));
xc = rc.*cos(anglec)+ centerC(1) ;
yc = rc.*sin(anglec)+ centerC(2) ;
%DU transmiter in R
centerD = [0 ,0];
angledt = 2*pi*rand(1,D);
rdt = 0.8*CR*sqrt(rand(1,D));
xdt = rdt.*cos(angledt)+ centerD(1);
ydt = rdt.*sin(angledt)+ centerD(2);
%DU reciver in distance dd from transmitter
for d=1:D
    centerDr=[xdt(d) ydt(d)];
    angledr= 2*pi*rand; %returns a single uniformly distributed random number in the interval (0,1).
    xdr(d) = dd.*cos(angledr)+ centerDr(1);
    ydr(d) = dd.*sin(angledr)+ centerDr(2);
    
end


%% CUs  to BSi
d_mi=zeros();
for bs=1:I
    for b=1:C
        d_mi(b,bs)=sqrt((xc(b)-xBS(bs)).^2 + (yc(b)-yBS(bs)).^2);
    end
end

%%  c * RB * BSi
PL_iB=10*pathloss_parameter*log10(d_mi)+ 22.7 + 26*log10(FREQ);
shadowing_iB = 4*ones(size(d_mi));
shadowing_iB = shadowing_iB.*randn(size(shadowing_iB));
sub_power_array_iB = 10.^((shadowing_iB-PL_iB)./10);% convert to linear
h_iB = sqrt(0.5)*randn(C,K,I) +sqrt(0.5)* 1j*randn(C,K,I);

for ee=1:I
    for cc=1:C
        g_mBi(cc,:,ee)=sub_power_array_iB(cc,ee).*(abs(h_iB(cc,:,ee))).^2;
    end
end

%% pathloss c *RB * e
d_mE= zeros(1,C);
for v=1:C
    d_mE(v)=  sqrt((xc(v)-xE).^2 + (yc(v)-yE).^2); % 1xN
end

PL_ie=10*pathloss_parameter*log10(d_mE)+ 22.7 + 26*log10(FREQ);
shadowing_ie = 4*ones(size(d_mE));
shadowing_ie = shadowing_ie.*randn(size(shadowing_ie));
sub_power_array_ie = 10.^((shadowing_ie-PL_ie)./10);% convert to linear
h_ie = sqrt(0.5)*randn(C,K) +sqrt(0.5)* 1j*randn(C,K);
for ci=1:C
    g_mE(ci,:)=sub_power_array_ie(ci)*(abs(h_ie(ci,:))).^2;
end

%% pathloss d*RB*Bi
d_di=zeros();
for aa=1:I
    for t=1:D
        d_di(t,aa)=sqrt(  (xdt(t)-xBS(aa))^2+(ydt(t)^2 -yBS(aa)));
    end
end
PL_jtB=10*pathloss_parameter*log10(d_di)+ 22.7 + 26*log10(FREQ);
shadowing_jt_B = 4*ones(size(d_di));
shadowing_jt_B = shadowing_jt_B.*randn(size(shadowing_jt_B));
sub_power_array_jt_B = 10.^((shadowing_jt_B-PL_jtB)./10);% convert to linear
h_jtB = sqrt(0.5)*randn(D,K,I) +sqrt(0.5)* 1j*randn(D,K,I);
for fk=1:I
    for dtt=1:1:D
        g_dtBi(dtt,:,fk)=sub_power_array_jt_B(dtt,fk)*(abs(h_jtB(dtt,:,fk))).^2;
    end
end

%% pathloss d * RB * e
d_dE= zeros();
for f=1:D
    d_dE(f)=  sqrt(   (xdt(f)-xE).^2 + (ydt(f)-yE).^2); % 1xN
end
PL_jt_e=10*pathloss_parameter*log10(d_dE)+ 22.7 + 26*log10(FREQ);
shadowingjt_e = 4*ones(size(d_dE));
shadowingjt_e = shadowingjt_e.*randn(size(shadowingjt_e));
sub_power_array_jt_e = 10.^((shadowingjt_e-PL_jt_e)./10);% convert to linear
h_jt_e = sqrt(0.5)*randn(D,K) +sqrt(0.5)* 1j*randn(D,K);
for dttt=1:1:D
    g_dtE(dttt,:)=sub_power_array_jt_e(dttt)*(abs(h_jt_e(dttt,:))).^2;
end

%% pathloss c * RB * d
d_md=zeros();
for s=1:C
    for x=1:D
        d_md(s,x)=sqrt(   (xc(s)-xdr(x)).^2 + (yc(s)-ydr(x)).^2);
    end
end

PL_dijr=10*pathloss_parameter*log10(d_md)+ 22.7 + 26*log10(FREQ);
shadowing_ijr = 4*ones(size(d_md));
shadowing_ijr = shadowing_ijr.*randn(size(shadowing_ijr));
sub_power_array_ijr = 10.^((shadowing_ijr-PL_dijr)./10);% convert to linear
h_ij = sqrt(0.5)*randn(C,K,D) +sqrt(0.5)* 1j*randn(C,K,D);
for ci=1:C
    for drr=1:D
        g_mdr(ci,:,drr)=sub_power_array_ijr(ci,drr).*(abs(h_ij(ci,:,drr))).^2 ; %C*K*D
    end
end

%% pathloss d *RB* d
d_d=zeros();
for dtt=1:D
    for drr=1:D
        d_d(dtt,drr)=  sqrt( (xdt(dtt)-xdr(drr)).^2 +(ydt(dtt)-ydr(drr)).^2 );
    end
end
PL_djJr=10*pathloss_parameter*log10(d_d)+ 22.7 + 26*log10(FREQ);
shadowing_djJr = 4*ones(size(d_d));
shadowing_djJr = shadowing_djJr.*randn(size(shadowing_djJr));
sub_power_array_djJr = 10.^((shadowing_djJr-PL_djJr)./10);% convert to linear
h_jj = sqrt(0.5)*randn(D,K,D) +sqrt(0.5)* 1j*randn(D,K,D);

for dti=1:D
    for drj=1:D
           g_dtdr(dti,:,drj)=sub_power_array_djJr(dti,drj).*(abs(h_jj(dti,:,drj))).^2 ; %D*D*N_RB
    end
end

 %% pathloss d * RB
PL_D2D=10*pathloss_parameter*log10(dd*ones(1,D))+ 22.7 + 26*log10(FREQ);
shadowing_D2D = 4*ones(size(dd*ones(1,D)));
shadowing_D2D = shadowing_D2D.*randn(size(shadowing_D2D));
sub_power_array_dd = 10.^((shadowing_D2D-PL_D2D)./10);% convert to linear
h_D2D = sqrt(0.5)*randn(D,K) +sqrt(0.5)* 1j*randn(D,K);
for dtttt=1:1:D
    g_J(dtttt,:)=sub_power_array_dd(dtttt)*(abs(h_D2D(dtttt,:))).^2;
end


%%

model.xdt=xdt;
model.ydt=ydt;

model.d_md=d_md;
model.g_mdr=g_mdr;
model.g_mBi=g_mBi;
model.d_mi=d_mi;
model.g_dtdr=g_dtdr;
model.g_J=g_J;
model.d_d=d_d;
model.g_dtE=g_dtE;
model.d_dE=d_dE;
model.g_dtBi=g_dtBi;
model.g_mE =g_mE;
model.d_mE =d_mE;
model.C=C;
model.D=D;
model.K=K;
model.CR=CR;
model.dd=dd;
model.pathloss_parameter=pathloss_parameter;
model.FREQ=FREQ;
model.BW=BW;
model.BOLTZ=BOLTZ;
end
