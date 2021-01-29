function CU_CB=CUsAssignments(g_mBi,I,C)
assignment=zeros(I,C);
CU_CB=zeros(C,4);
channel=zeros(1,C);
gmiBest=zeros(C,I);
for s=1:I
    assignment(s,:)=munkres(1./g_mBi(:,:,s)); % find strongest channel of all user for all BS
end

Assignment=assignment';
for i=1:I
    for q=1:C
        gmiBest(q,i)=g_mBi(q,Assignment(q,i),i);
    end
end

[maxGains,I_ms]=max(gmiBest,[],2);
for f=1:C
    channel(f)=Assignment(f,I_ms(f));
end
%     %  [maxGains,I_ms]=max(g_mi(:,:,s),[],2); % strongest channel of all user to a BS
%     All_BS_maxGains(:,s)=maxGains;% strongest channel to ALL BSs
%  Chs_BSi(:,s)=I_ms; % stongest channels of users to all BS
%
%
%
% [BestGain,BSi]=max(All_BS_maxGains,[],2); % strongest channel  to best BSs
% for c=1:C
%
%     CU_ch(c)=Chs_BSi(c,BSi(c));
% end
% CU_CB=[BestGain  CU_ch' BSi];
for c=1:C
    CU_CB(c,:)=[c maxGains(c) channel(c) I_ms(c)];
end

end