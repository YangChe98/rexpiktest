clc
clear
N=2;

B=[1   2   17;
3   4   18;
5   6   19;
7   8   20;
9   10  21;
11  12  22;
13  14  23;
15  16  24;
25  26  27];
A(1,:,:)=B(1:3,:);
A(2,:,:)=B(4:6,:);
A(3,:,:)=B(7:9,:);
for i=1:3
    for j=1:3
        for k=1:3
            X=[i,j,k,A(i,j,k)];
            display(X)
            
        end
    end
end
a_real=fopen('/home/yangche/vscodecode/rexpiktest/r_real10.txt', 'wt');

for i=1:3
    for j=1:3
        for k=1:3
            fprintf(a_real, '%e\t', A(i,j,k));
        end
         fprintf(a_real, '\n');
    end

end
%a_real=fopen('/home/yangche/matlab_code/numerical singular/r_real10.txt', 'wt');
for k1=1:N
    for k2=1:N
        for k3=1:N
            for j1=1:N
                for j2=1:N
                    for j3=1:N
                        kk11=k1+j1-N;
                        kk22=k2+j2-N;
                        kk33=k3+j3-N;
                        
                        kk1=abs(k1+j1-N-2)+1;
                       kk2=abs(k2+j2-N-2)+1;
                       kk3=abs(k3+j3-N-2)+1;
                       X=[k1,k2,k3,j1,j2,j3,kk1,kk2,kk3]
                        % fprintf(a_real, '%e\t', A(kk1,kk2,kk3));
                        A_re((k1-1)*N^2+(k2-1)*N+k3,(j1-1)*N^2+(j2-1)*N+j3)=A(kk1,kk2,kk3);
                       % if(kk11==0)
                      %  A_im((k1-1)*N^2+(k2-1)*N+k3,(j1-1)*N^2+(j2-1)*N+j3)
                    end
                end
            end
                %fprintf(a_real, '\n');
        end
    end
end