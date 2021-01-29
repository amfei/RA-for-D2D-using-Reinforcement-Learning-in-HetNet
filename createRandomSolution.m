%% Create Random solution that meet resource allocation constraints
function f=createRandomSolution(C,D,K)


x_ik=zeros(C,K);
x_jk=zeros(D,K);
CK=randperm(C);
DK=randperm(D);

for i=1:C
x_ik(CK(i), i)=1;
end
for j=1:D
x_jk(DK(j), j)=1;
end

f1= x_ik(:,1:K);
f2=x_jk(:,1:K);
f=[f1;f2];