function  ret=sdr2cvx_unsensor(range,noise_cmatrix,sensor_sel_number,sensor_number,dim)
%凸优化
cvx_clear
cvx_begin sdp
cvx_solver sdpt3
cvx_precision best
cvx_quiet(1)
variable Select(sensor_number,1);
variable Z(dim+1,dim+1);
variable W(sensor_number,sensor_number);
minimize(trace(Z(1:dim,1:dim)))
subject to
[Z,eye(dim+1);eye(dim+1),range'*(diag(Select)*inv(noise_cmatrix)*range)]>=0;%增加了12个变量;恒为15个约束
[W,Select;Select',1]>=0;%减少(sensor_number-1+...+2)个变量
trace(W)==sensor_sel_number;
diag(W)==Select;
cvx_end
%传感器选择向量（高斯随机化）
Select_copy=Select;
N=100;
Dw=W-Select*Select';
%判断rank1
landa=eig(Dw);
landa=sort(landa,'descend');
if (landa(1)/landa(2))>100000
    Select_sort=sort(Select,'descend');
    select_max_ssn=Select_sort(sensor_sel_number);
    Select(Select<select_max_ssn)=0;
    Select(Select>0)=1;
    
    Select_r=Select;
    Fisher_information_matrix_inverse_rank1=inv(range'*diag(Select_r)/(noise_cmatrix)*diag(Select_r)*range);
    crlb_rank1=trace(Fisher_information_matrix_inverse_rank1(1:dim,1:dim));
    crlb_01=crlb_rank1;
    ret=[crlb_01;crlb_rank1;cvx_cputime];
else
    if landa>0
    for i=1:N
        SELECT_gaussian_ramdomization(:,i)=mvnrnd(Select_copy,Dw);
        Select_gaussian_ramdomization_sort=sort(SELECT_gaussian_ramdomization(:,i),'descend');
        select_gaussian_ramdomization_max_ssn=Select_gaussian_ramdomization_sort(sensor_sel_number);
        T=SELECT_gaussian_ramdomization(:,i).*0;
        T(find(SELECT_gaussian_ramdomization(:,i)>=select_gaussian_ramdomization_max_ssn))=1;
        SELECT_gaussian_ramdomization(:,i)=T;

        Fisher_information_matrix_inverse_radoms=inv(range'*diag(SELECT_gaussian_ramdomization(:,i))*inv(noise_cmatrix)*diag(SELECT_gaussian_ramdomization(:,i))*range);
        Crlb_radom(i)=trace(Fisher_information_matrix_inverse_radoms(1:dim,1:dim));
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
            %toa模型的FIM
            Fisher_information_matrix_inverse_radoms=inv(range'*diag(SELECT_gaussian_ramdomization(:,i))*inv(noise_cmatrix)*diag(SELECT_gaussian_ramdomization(:,i))*range);
            Crlb_radom(i)=trace(Fisher_information_matrix_inverse_radoms(1:dim,1:dim));
        end
    end
        crlb_radom=min(Crlb_radom);
        number_min=find(Crlb_radom==crlb_radom);
        Select_r=SELECT_gaussian_ramdomization(:,number_min(1));
        %未高斯随机化的CRLB
        Select_sort=sort(Select,'descend');
        select_max_ssn=Select_sort(sensor_sel_number);
        Select(Select<select_max_ssn)=0;
        Select(Select>0)=1;
    
        Select_r=Select;
        Fisher_information_matrix_inverse_rank1=inv(range'*diag(Select_r)/(noise_cmatrix)*diag(Select_r)*range);
        crlb_01=trace(Fisher_information_matrix_inverse_rank1(1:dim,1:dim));
        %返回值 第一个为w/o GR
        ret=[crlb_01;crlb_radom;cvx_cputime];
        
end
