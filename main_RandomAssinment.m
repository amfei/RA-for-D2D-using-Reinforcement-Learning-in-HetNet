
clc
clear
close all
%Pm=0.200
 alpha = 0.95;    % learning rate    % TODO : we need exploration rate schedule

   Pm= 0.2% [0.0316 0.0630 0.125 0.25 ]  % [15 18 21 25]
    Qdu= 3
    Pd=[0.1 0.02 0.08];
    beta = 0.9;    % discount factor  % TODO : we need learning rate schedule
        epsilon = 0.9;  % exploration probability (1-epsilon = exploit / epsilon = explore)
    %pathloss_parameter = 3.5;
    Rc=500; % level of distance between Cu and d2d for assignment
    Qcu=3; % min sinr guaranteeing QoS cus
    
    T=500;% time loop with same model
    tau=0.01;
    rd=150;
    C=50 ;%15 20 25 30 35 40 45];
    D=50;%10 15 20 25 30 35 40];
    K=50;%15 20 25 30 35 40 45]; % K: number of RB
    I=4;
    dd=20;
    pathloss_parameter = 3.5;
    FREQ = 2.15;
    BW = 10e6;
    BOLTZ = 1.3806488e-23;
    CR = 500;   %CellRadius
    crp=200; %radious of pico
    
    loop=10;
    % Mont carlo loop
    for l=1:loop
        
        %model=SelectModel();
        model = modelHetNet(I,C,D,K,CR,crp,dd,pathloss_parameter,FREQ,BW,BOLTZ);
        g_dtdr=model.g_dtdr;              %Dt*Dr*RB
        g_mdr=model.g_mdr;                % C*K*D
        g_mBi=model.g_mBi;
        g_dtBi=model.g_dtBi;
        g_dtE=model.g_dtE;
        g_mE=model.g_mE;
        C=model.C;
        D=model.D;
        K=model.K;
        pathloss_parameter=model.pathloss_parameter;
        FREQ=model.FREQ;
        BW=model.BW;
        BOLTZ=model.BOLTZ;
        R=zeros(C,D);
        d_md=model.d_md;
        xdt=model.xdt; %location of d2d transmiter
        ydt=model.ydt;
        dt=[xdt;ydt];
        % Initalization
        next_rewarD=0;
        FF=0;
        TT=0;
        states=[0,1];
        state_idx=1;
        % CU Channel and BS Assignment
        Bi=randi(4,1,D); 
        Kk=randperm(D);
        for ccc=1:C
            G_mBi(ccc)=g_mBi(ccc,Kk(ccc),Bi(ccc));
        end
        CU_CB=[G_mBi' Kk' Bi' ];
        
        %CU_CB=CUsRandomAssignments(g_mBi,I,D); %[CU_gains k BSi] find strongest channel of all cu to a BS using munkler
        %Create D2D action list from D2D CU_CB
        Act_list=ActionList(Pd,C,CU_CB);  % 3 level of power for all d2d
        %Binnay co-channels users
       % BiCoCH=BinaryCoChannel(CU_CB(:,2));
       
        % creat  D Q-table for initialixzation
        Q = zeros(numel(states),size(Act_list,1),D);
        %% Q-learning loop
        for d=1:D
           for t=1:T
                 for D2D=1:D
                [current_action,act_idx]= RandomAction(Act_list);
               
               
                    D2D_RAM(D2D,:)=[D2D current_action act_idx]; %
                end
                %Action of current d2d d
                k=D2D_RAM(d,3);  %  d2d channel d at t
                BS=D2D_RAM(d,4);%  d2d  BSi  at t
                pd=D2D_RAM(d,2);  %  d2d power d at t
                % Acton of All d2d pairs
                ks=D2D_RAM(:,3);
                pds=D2D_RAM(:,2);
                dI=D2D_RAM(:,4);%  d2ds use resource od CUs that is belong to BSi
                DUs_CB=D2D_RAM(:,[1 3:4]) ; % [d2d k BSi]
                % current D2D Throughput
                d2d_throughput=D2D_throughput(CU_CB,DUs_CB,k,d,pds,Pm,g_dtdr,g_mdr,BOLTZ,BW,K);
                % Chek D2D QoS condition
                if  d2d_throughput > Qdu
                    % Update next state
                    next_state=1;
                    % Eav Throughput
                    Eav_throughput=EAV_throughput(DUs_CB,k,d,pds,g_dtE,CU_CB,Pm,g_mE,BOLTZ,BW,K);
                    % D2D Secracy capacity
                    SC_DUS=max(d2d_throughput-Eav_throughput,0);
                    % Reward
                    reward_t=SC_DUS;
                else
                    next_state=0;
                    reward_t=0;
                end
                %% Find Next State Reward
                next_rewarD=next_rewarD+reward_t;
                %% Find Next State Index
                next_state_idx = find(states==next_state);  % id of the next state
                %%  Update Q
                % Independent Q
                Q(state_idx,act_idx ,d)=(1-alpha)*Q(state_idx,act_idx,d)+...
                    alpha*(next_rewarD + beta* (max(Q(next_state_idx,:,d))));
                %S;
                %% Update Netx State Index
                state_idx = next_state_idx;
                [Sec,II]=max(Q(:,:,d),[],2)    ;     % finding the max value of each row of dth matrix
                G(t)=Sec(1);
                
            end
            TT=TT+G;
        end
        AA(l,:)=TT;
    end
    BB=round(sum(AA)/loop)
    Ave_sec_capacity_loops=max(BB) 
   BBB=num2str(BB)