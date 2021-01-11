clear
close all
clc

%===========================================

path = ['.\PubMed\'];
name = 'PubMed';
load([path name]);
%[n d] = size(fea);
%times = 1;

n = size(A,1);


c = length(unique(gnd));

k = 100;


X = M';



trial = 20;
err = zeros(trial,c,n);
mistakes = zeros(trial,c,n);
SVs = zeros(trial,c,n);
TMs = zeros(trial,c,n);
time = zeros(trial,c,n);
Query = zeros(trial,c,n);

method = 'GPA'
options.eta = 0.1;
for i = 1:trial
    id_list = ID_ALL(i,:);
    for j=1:c
        Y = -ones(n,1);
        index_pos = find(gnd==j);
        %index_neg = setdiff(1:n,index_pos);       
        Y(index_pos) = 1;
        [classifier, err(i,j,:),run_time, mistakes(i,j,:), SVs(i,j,:), TMs(i,j,:)] = GPA(X,Y,options,id_list);
    end
end
meanErr2 = mean(err,2);
meanErr2 = squeeze(meanErr2);
meanErr = mean(meanErr2,1);
stdErr = std(meanErr2,0,1);
meanMis2 = mean(mistakes,2);
meanMis2 = squeeze(meanMis2);
meanMis = mean(meanMis2,1);
stdMis = std(meanMis2,0,1);
meanSV2 = mean(SVs,2);
meanSV2 = squeeze(meanSV2);
meanSV = mean(meanSV2,1);
stdSV = std(meanSV2,0,1);
meanTM2 = mean(TMs,2);
meanTM2 = squeeze(meanTM2);
meanTM = mean(meanTM2,1);
stdTM = std(meanTM2,0,1);
save([path method],'options','meanErr','stdErr','meanMis','stdMis','meanSV','stdSV','meanTM','stdTM');
fid = fopen([path method '.txt'],'a+');
fprintf(fid, '%d %.4f %.4f %.4f %.4f %.4f %.4f\n',k,meanErr(end),stdErr(end),meanMis(end),stdMis(end),meanTM(end),stdTM(end));
fclose(fid);


method = 'OLLGC'  
options.a = 0.001;
for i = 1:trial
    id_list = ID_ALL(i,:);
    for j=1:c
        Y = -ones(n,1);
        index_pos = find(gnd==j);
        %index_neg = setdiff(1:n,index_pos);       
        Y(index_pos) = 1;
        [classifier, err(i,j,:),run_time, mistakes(i,j,:), SVs(i,j,:), TMs(i,j,:)] = OLLGC(X,Y,options,id_list);
    end
end
meanErr2 = mean(err,2);
meanErr2 = squeeze(meanErr2);
meanErr = mean(meanErr2,1);
stdErr = std(meanErr2,0,1);
meanMis2 = mean(mistakes,2);
meanMis2 = squeeze(meanMis2);
meanMis = mean(meanMis2,1);
stdMis = std(meanMis2,0,1);
meanSV2 = mean(SVs,2);
meanSV2 = squeeze(meanSV2);
meanSV = mean(meanSV2,1);
stdSV = std(meanSV2,0,1);
meanTM2 = mean(TMs,2);
meanTM2 = squeeze(meanTM2);
meanTM = mean(meanTM2,1);
stdTM = std(meanTM2,0,1);
save([path method],'options','meanErr','stdErr','meanMis','stdMis','meanSV','stdSV','meanTM','stdTM');
fid = fopen([path method '.txt'],'a+');
fprintf(fid, '%d %.4f %.4f %.4f %.4f %.4f %.4f\n',k,meanErr(end),stdErr(end),meanMis(end),stdMis(end),meanTM(end),stdTM(end));
fclose(fid);


method = 'SSLGC'  
options.a = 0.001;
options.k = 0.4;
for i = 1:trial
    id_list = ID_ALL(i,:);
    for j=1:c
        Y = -ones(n,1);
        index_pos = find(gnd==j);
        %index_neg = setdiff(1:n,index_pos);       
        Y(index_pos) = 1;
        [classifier, err(i,j,:),run_time, mistakes(i,j,:), SVs(i,j,:), TMs(i,j,:),Query(i,j,:)] = SSLGC(X,Y,options,id_list);
    end
end
meanErr2 = mean(err,2);
meanErr2 = squeeze(meanErr2);
meanErr = mean(meanErr2,1);
stdErr = std(meanErr2,0,1);
meanMis2 = mean(mistakes,2);
meanMis2 = squeeze(meanMis2);
meanMis = mean(meanMis2,1);
stdMis = std(meanMis2,0,1);
meanSV2 = mean(SVs,2);
meanSV2 = squeeze(meanSV2);
meanSV = mean(meanSV2,1);
stdSV = std(meanSV2,0,1);
meanTM2 = mean(TMs,2);
meanTM2 = squeeze(meanTM2);
meanTM = mean(meanTM2,1);
stdTM = std(meanTM2,0,1);
meanQU2 = mean(Query,2);
meanQU2 = squeeze(meanQU2);
meanQU = mean(meanQU2,1);
stdQU = std(meanQU2,0,1);
save([path method],'options','meanErr','stdErr','meanMis','stdMis','meanSV','stdSV','meanTM','stdTM','meanQU','stdQU');
fid = fopen([path method '.txt'],'a+');
fprintf(fid, '%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n',k,meanErr(end),stdErr(end),meanMis(end),stdMis(end),meanTM(end),stdTM(end),meanQU(end),stdQU(end));
fclose(fid);


figure
method = 'GPA'

load([path method]);
plot(meanMis,'r','LineWidth',2);
hold on

method = 'OLLGC'  

load([path method]);
plot(meanMis,'b--','LineWidth',2);

method = 'SSLGC'  

load([path method]);
plot(meanMis,'k-.','LineWidth',2);

legend('GPA','OLLGC','SSLGC',1)
xlabel('#rounds')
ylabel('cumulative error rate')

axis([1 length(meanMis) 0.15 0.3])

figure

method = 'GPA'

load([path method]);
n = length(meanErr);
plot(1:n,'r','LineWidth',2);
hold on

method = 'OLLGC'  

load([path method]);
plot(1:n,'b--','LineWidth',2);

method = 'SSLGC'  

load([path method]);
plot(meanQU,'k-.','LineWidth',2);

legend('GPA','OLLGC','SSLGC',2)
xlabel('#rounds')
ylabel('#cumulative queried nodes')

figure

method = 'GPA'

load([path method]);
plot(meanTM,'r','LineWidth',2);
hold on

method = 'OLLGC'  

load([path method]);
plot(meanTM,'b--','LineWidth',2);

method = 'SSLGC'  

load([path method]);
plot(meanTM,'k-.','LineWidth',2);

legend('GPA','OLLGC','SSLGC',2)
xlabel('#rounds')
ylabel('cumulative time (in second)')
