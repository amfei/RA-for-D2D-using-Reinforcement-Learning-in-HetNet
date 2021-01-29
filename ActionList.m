function Act_list=ActionList(Pd,C,CU_CB)
action_pd_ch_I=zeros(C,3);
kk=0;
actionT{1,numel(Pd)}=[];
for te=1:numel(Pd)
    for j=1:C
        action_pd_ch_I(j,:)=[Pd(te) CU_CB(j,3) CU_CB(j,4)];
    end
    kk = kk+1;
    actionT{kk}=action_pd_ch_I;
end
Act_list=cell2mat(actionT'); % each D2D has a selection on this set
