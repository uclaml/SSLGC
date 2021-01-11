function [classifier, error, run_time, mistakes, SVs, TMs]=  GPA(X,Y,options,id_list)

[d n] = size(X);

eta = options.eta;
w = zeros(d,1);
T = length(id_list);
SV = [];
mistakes = zeros(1,T);
error = zeros(1,T);
SVs = zeros(1,T);
TMs=zeros(1,T);
err_count = 0;
tic;
for t = 1:T
    i = id_list(t);
    x_t = X(:,i);
    y_t = Y(i);
    f_t=w'*x_t;
    hat_y_t = sign(f_t);
    if (hat_y_t==0)
        hat_y_t=1;
    end
        
    %% Error Counting and Updating
    if (hat_y_t~=y_t)
        w = w + eta*x_t*y_t;
        err_count = err_count + 1;     
        SV = [SV i];
    end
    %if Y(i)*w'*X(:,i) <= 0
        
    
        %err_count = err_count+1;
        
    %end
    
    run_time = toc;
    %if (mod(t,t_tick)==0)
    error(t) = err_count;%[error ];
    mistakes(t) = err_count/t;%[mistakes err_count/t];
    SVs(t) = length(SV);%[SVs length(SV)];
    TMs(t)= run_time;
end

classifier.SV = SV;
classifier.w = w;
fprintf(1,'GPA: The number of mistakes = %d\n', err_count);
run_time = toc;