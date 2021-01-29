
function  d2d_throughput=D2D_throughput(CU_CB,DUs_CB,k,d,pds,Pm,g_dtdr,g_mdr,BOLTZ,BW,K)

DUs_CoC=find (DUs_CB(:,2)==k);
for du=1:size(DUs_CoC,1)
    I_dpd(du)=pds(DUs_CoC(du)).*(g_dtdr(DUs_CoC(du),d,k));
end
I_dpd(d)=0;


CUs_CoC=find (CU_CB(:,2)==k);
for cu=1:size(CUs_CoC,1)
    I_md(cu)=Pm.*g_mdr(CUs_CoC(cu),d,k);
end

SINR_D=(pds(d).*g_dtdr(d,d,k))./(BOLTZ*293*BW/K+sum(I_dpd)+sum(I_md));
d2d_throughput = log2(1+SINR_D);
