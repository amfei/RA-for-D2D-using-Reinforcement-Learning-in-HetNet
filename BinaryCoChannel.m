function  B=BinaryCoChannel(A)
for col=1:size(A,1)
    for i=1:size(A,1)
        if A(col)==A(i)
            B(col,i)=1 ;
            B(col,col)=1;
        else
            B(col,i)=0;
        end
    end
end