function [classifier, error, run_time, mistakes, SVs, TMs] = OLLGC(X, Y, options, id_list)

%% initialize parameters
ID = id_list;
d = size(X,1);
T = length(ID);
SV = [];
%ID = id_list;
err_count = 0;
%t_tick = options.t_tick;
mistakes = zeros(1,T);
error = zeros(1,T);
SVs = zeros(1,T);
TMs=zeros(1,T);
v=zeros(d,1);
a = options.a;
%A = 1/a*eye(d);
invA = 1/a*eye(d);

%% loop
tic;
for t = 1:T,
    id = ID(t);
    y_t = Y(id);
    x_t = X(:,id);
    
    %% Prediction
    tmp=x_t'*invA*x_t;
    %beta_t=1/(tmp+1);
    %f_t= v'*(invA-invA*x_t*x_t'*invA/(1+tmp))*x_t;
    f_t= v'*invA*x_t;
    hat_y_t = sign(f_t);
    if (hat_y_t==0)
        hat_y_t=1;
    end
    
        
    %% Error Counting and Updating
    if (hat_y_t~=y_t),
        err_count = err_count + 1;    
        v = v + y_t*x_t;
        invA=invA-invA*x_t*x_t'*invA/(1+tmp);
        SV = [SV id];
    end
        
    %end
    
    run_time = toc;
    %if (mod(t,t_tick)==0)
    error(t) = err_count;%[error err_count];
        mistakes(t) = err_count/t;%[mistakes err_count/t];
        SVs(t) = length(SV);
        %TMs=[TMs run_time];
        TMs(t) = run_time;
    %end
end

classifier.SV = SV;
classifier.v = v;
fprintf(1,'OLLGC: The number of mistakes = %d\n', err_count);
run_time = toc;


