function  ret=sdr1cvx_unsensor(Q0,B,C,a,sensor_sel_number,sensor_number,dim)
%input Q0,B,C,a,sensor_sel_number
%凸优化
tic;
cvx_clear
cvx_begin sdp
cvx_solver sdpt3
cvx_precision best
cvx_quiet(1)
variable Select(sensor_number,1);
variable Z(dim,dim);
variable W(sensor_number,sensor_number);
minimize(trace(Z));
subject to
[W,Select;Select',1]>=0; %减少（sensor_number-1+...+2）个变量
trace(W)==sensor_sel_number;
diag(W)==Select;
[inv(Q0)-B/(C)*B'+(1/a)*(diag(Select)-(1/sensor_sel_number)*W),B/(C);(C)\B',Z]>=0; %约束中包含的变量有Z，Select,W;而TOA方法中变量只有Z和Select
cvx_end
toc;
%传感器选择向量（高斯随机化）
Select_copy=Select;
N=100;
Dw=W-Select*Select';
Dw=(Dw+Dw')./2;
%判断rank1
landa=eig(Dw);
landa=sort(landa,'descend');
lambda=eig(W);
lambda=sort(lambda,'descend');
if (lambda(1)/lambda(2))>100000
    %传感器选择向量01优化
    Select_copy=Select;
    Select_sort=sort(Select,'descend');
    select_max_ssn=Select_sort(sensor_sel_number);
    Select(Select<select_max_ssn)=0;
    Select(Select>0)=1;
    %-------------------------------------
    %选择任意选中的传感器作为参照传感器
    ref1=unidrnd(sensor_sel_number);%随机生成最大为sensor_sel_number的正整数
    ref2=find(Select==1);
    ref=ref2(ref1);%随机选择参照传感器的序号
    Select_normal=Select;
    Select_ref=zeros(sensor_number,1);
    Select_ref(ref)=1;
    Select_normal(ref)=0;
    Fisher_information_matrix_inverse_01n=inv(C-B'/(inv(Q0)+1/(a)*(diag(Select_ref)+diag(Select_normal)-(1/sensor_sel_number)*(Select_ref+Select_normal)*(Select_ref+Select_normal)'))*B);
    crlb_rank1=trace(Fisher_information_matrix_inverse_01n);
    crlb_01=crlb_rank1;
    ret=[crlb_01;crlb_rank1;cvx_cputime];
else
    if landa>=0
        for i=1:N
        SELECT_gaussian_ramdomization(:,i)=mvnrnd(Select_copy,Dw);
        Select_gaussian_ramdomization_sort=sort(SELECT_gaussian_ramdomization(:,i),'descend');
        select_gaussian_ramdomization_max_ssn=Select_gaussian_ramdomization_sort(sensor_sel_number);
        T=SELECT_gaussian_ramdomization(:,i).*0;
        T(find(SELECT_gaussian_ramdomization(:,i)>=select_gaussian_ramdomization_max_ssn))=1;
        SELECT_gaussian_ramdomization(:,i)=T;
        %tdoa模型的FIM
        %选择任意参照传感器
        ref1_g=unidrnd(sensor_sel_number);%随机生成最大为sensor_sel_number的正整数
        ref2_g=find(SELECT_gaussian_ramdomization(:,i)==1);
        ref_g=ref2_g(ref1_g);%随机选择参照传感器的序号
        Select_normal_g=SELECT_gaussian_ramdomization(:,i);
        Select_ref_g=zeros(sensor_number,1);
        Select_ref_g(ref_g)=1;
        Select_normal_g(ref_g)=0;
        Fisher_information_matrix_inverse_rn=inv(C-B'/(inv(Q0)+1/(a)*(diag(Select_ref_g)+diag(Select_normal_g)-(1/sensor_sel_number)*(Select_ref_g+Select_normal_g)*(Select_ref_g+Select_normal_g)'))*B);
        crlb_01_TDOArn(i)=trace(Fisher_information_matrix_inverse_rn);
        
        end
    else
        [V,D]=eig(Dw);
        D(D<0)=0;
        Dw=V*D*V';
        for i=1:N
                SELECT_gaussian_ramdomization(:,i)=mvnrnd(Select_copy,Dw);
                Select_gaussian_ramdomization_sort=sort(SELECT_gaussian_ramdomization(:,i),'descend');
                select_gaussian_ramdomization_max_ssn=Select_gaussian_ramdomization_sort(sensor_sel_number);
                T=SELECT_gaussian_ramdomization(:,i).*0;
                T(find(SELECT_gaussian_ramdomization(:,i)>=select_gaussian_ramdomization_max_ssn))=1;
                SELECT_gaussian_ramdomization(:,i)=T;
                %tdoa模型的FIM
                %选择任意参照传感器
                ref1_g=unidrnd(sensor_sel_number);%随机生成最大为sensor_sel_number的正整数
                ref2_g=find(SELECT_gaussian_ramdomization(:,i)==1);
                ref_g=ref2_g(ref1_g);%随机选择参照传感器的序号
                Select_normal_g=SELECT_gaussian_ramdomization(:,i);
                Select_ref_g=zeros(sensor_number,1);
                Select_ref_g(ref_g)=1;
                Select_normal_g(ref_g)=0;
                Fisher_information_matrix_inverse_rn=inv(C-B'/(inv(Q0)+1/(a)*(diag(Select_ref_g)+diag(Select_normal_g)-(1/sensor_sel_number)*(Select_ref_g+Select_normal_g)*(Select_ref_g+Select_normal_g)'))*B);
                crlb_01_TDOArn(i)=trace(Fisher_information_matrix_inverse_rn);
        end
    end
    crlb_radom=min(crlb_01_TDOArn);
    %未高斯随机化的CRLB
    Select_sort=sort(Select,'descend');
    select_max_ssn=Select_sort(sensor_sel_number);
    Select(Select<select_max_ssn)=0;
    Select(Select>0)=1;
    
    %选择任意选中的传感器作为参照传感器
    ref1=unidrnd(sensor_sel_number);%随机生成最大为sensor_sel_number的正整数
    ref2=find(Select==1);
    ref=ref2(ref1);%随机选择参照传感器的序号
    Select_normal=Select;
    Select_ref=zeros(sensor_number,1);
    Select_ref(ref)=1;
    Select_normal(ref)=0;
    Fisher_information_matrix_inverse_01n=inv(C-B'/(inv(Q0)+1/(a)*(diag(Select_ref)+diag(Select_normal)-(1/sensor_sel_number)*(Select_ref+Select_normal)*(Select_ref+Select_normal)'))*B);
    crlb_01=trace(Fisher_information_matrix_inverse_01n);
    %返回值 第一个为w/o GR
    ret=[crlb_01;crlb_radom];
        
end
