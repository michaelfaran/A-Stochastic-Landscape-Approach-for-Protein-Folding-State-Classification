% A code for the paper- "A Stochastic Landscape Approach for Protein
% Folding State Classification", by Michael Faran, Dhiman Ray, Shubhadeep Nag, Umberto Raucci,
% Michele Parrinello, and Gili Bisker.
% This code was written by Michael Faran, 28/3/2024
% For any questions or inquiries, please send an email to: faranmic@mail.tau.ac.il
% 
% User Guide:
% 
% 1. Activate this script from the its origin folder. Please Note that
% this code is working for MATLAB version 23.2 (product number 1), and
% requires also the 'Statistics and Machine Learning Toolbox'(product
% number 19). It might work for previous versions, but it is not certain.
% 
% 2.Download beforehand to this main folder the following .m files:
%     a.nmi.m from:
%         https://www.mathworks.com/matlabcentral/fileexchange/29047-normalized-mutual-information
%     b.Knee_pt.m from:
%         https://www.mathworks.com/matlabcentral/fileexchange/35094-knee-point
%     c. munkres.m from:
%         https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
%     d.clustering_comparison.m from:
%         https://www.mathworks.com/matlabcentral/fileexchange/45222-hierarchical-cluster-comparison
%     e.randindex.m from:
%         https://www.mathworks.com/matlabcentral/fileexchange/130779-rand-and-adjusted-rand-index-calculator-for-cluster-analysis
%     f. Install to this folder the BEAST algorithm .m files, saving them in a
%         subfolder of this folder named "Beast" under this folder. Follow the instructions here:
%         https://www.mathworks.com/matlabcentral/fileexchange/72515-bayesian-changepoint-detection-time-series-decomposition
%         Be sure to have in the sub folder 'Beast/Matlab' the following functions:
%         beast.m,plotbeast.m, extractbeast.m and installbeast.m 
% 
% 3. The user inputs are assumed to be in the 'Main' folder, where this code
%     is activated from. They are listed as:
%     a. CV_Matrix_name-A N*M matrix, where N is the CV trajectory length and M is the
%         number of CVs to segment
%     b. M-Insert explictly M as the number of CVs of interest
%     c. DSI-The downsample index DSI, which the CV needs to be under the length
%         of 10000 samples (for the BEAST algorithm to converge in a reasonable
%         time).
%     d. The ground truth CV-vector- A 1*N vector of the correct state vs.
%          sample of the protein thoughout its dynamics. It should be a
%          vector where each sample attain a number correspnding to a
%          certain state. For example- "Folded" is 1, "Unfolded" is 2,
%          "Misfolded" is 3. The name of the vector inside the loaded
%          structure should be "ground_truth_vec". We will use "0" as
%          the label for unclassified samples.
%     e. Protein_name- a string containing the protein name. As
%          examples or test cases, there are two proteins-the Chingolin and the Trp
%          Cage.
%     f. N_states- The number of protein states in the groun truth vector.
% 
% 4. The code outputs are for each CV seperately (i index runs from i=1...M for each CV), assuming for
%     simplictly the protein name is "Protein_name":
% 
%     a. "BEAST_RAW_Protein_name_New_i"- The raw CV data figure, with the
%         trend changepoints marks detected by the BEAST algorithm.
%     b. "SLM_Clustering_Protein_namen_New_No_color_i"- The yielded stochastic landscape, after 
%         PCA and normalization, with the no cluster color marks on the
%         segment data (The segments are points in this space). 
%     c."Stochastic_Landscape_Unnormalized_Matrixi"- The raw segment data
%         of the segments stochastic coordinates for CV i
%     d. "SLM_Clustering_Protein_name_New_i"- The yielded stochastic landscape, after 
%         PCA and normalization, with the cluster color marks on the
%         segment data (The segments are points in this space).
%     e. "Stochastic_Landscape_PCA_Matrixi"- The PCA result on the
%         normalized segment data, which is projected later into the stochastic
%         landscape in the .png file of
%         "SLM_Clustering_Protein_name_New_No_color_i"
%     f. "BEAST_Summary_Protein_name_New_i"- The CV data clustered according to
%         the DBSCAN algorithm, with vertical lines indicating the trend
%         changepoints times along the CV.
%     g. "CV_vs_time_Protein_name_New_i"-  Protein state labeling, determined by the Kuhn-Munkres 
%         algorithm applied to the DBSCAN clustering results, projected onto
%         the original CV trajectory. This is the main output of interest
%         in a "black box" manner.
%     h."BEAST_Summary_Protein_name_New_i_workspace"- The MATLAB workspace in
%         the end of the run for each CV i
%     i. "Truth_mat"- The coed output of all the classification matrices 
%         (NMI  RI  ARI  Dice  FM  Jaccard), compared with the ground truth,
%         The rows correspond to each CV by their order in the input matrix
%         ("CV_Matrix_name").
% 
% Please refer to the paper for additional details, and examaine the repository:
% "https://github.com/luigibonati/deep-learning-slow-modes"
% 
% References:
% 
% This code uses the BEAST algorithmas a baseline-
% Zhao, Kaiguang, et al. "Detecting change-point, trend, and seasonality in satellite time series data
% to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm." Remote sensing of Environment 232 (2019): 111181.
% in conjuction with the Stochastic Landscape, first published in:
% Faran, Michael, and Gili Bisker. "Nonequilibrium Self-Assembly Time Forecasting by the Stochastic Landscape Method." The Journal of Physical Chemistry B 127.27 (2023): 6113-6124.
% On protein folding MD collective variables data, published in:
% Ray, D.; Trizio, E.; Parrinello, M. Deep learning collective variables from transition
% path ensemble. The Journal of Chemical Physics 2023, 158, 204102
% Bonati, Luigi, Giovanni,Maria Piccini, and Michele Parrinello. "Deep learning the slow modes for rare events sampling."
% Proceedings of the National Academy of Sciences 118.44 (2021): e2113533118.

%User inputs (put the CV_Matrix_name.mat in the 'Main' folder):
CV_Matrix_name='CVs_Trp_Cage.mat';%CV_Matrix_name-A N*M matrix, change to CVs_Chignolin.mat for the second shared example.
ground_truth_vec_name='GT_Tr_cage.mat'; %Change to GT_Chingolin.mat for the second shared example.
M=9; %Number of CVs
DSI=1;%The downsample index
Protein_name='Trp_Cage'; %A string containing the protein name. Change to Chignolin for the second shared example.
N_states=3; %Number of states in the ground truth

load(ground_truth_vec_name); %The ground truth CV-vector
CV_content=load(CV_Matrix_name);
field_name = fieldnames(CV_content);
CV_MAT = CV_content.(field_name{1});
folder_name = horzcat(pwd,'\',Protein_name,'_Results');

if ~exist(folder_name, 'dir')
    mkdir(folder_name);
    disp(['Folder "', folder_name, '" created successfully.']);
else
    disp(['Folder "', folder_name, '" already exists. Skipping creation.']);
end

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
save_here=folder_name;

% DSI=1;%was 8
% load('GT_Chingolin_New.mat'); %this yields ground_truth_vec
% Truth_mat=zeros(4,6); %The size of the 

for ok=1:1:M

    E_vec=CV_MAT(:,ok);
    Ipo=find(isnan((E_vec)));
    E_vec(Ipo)=E_vec(Ipo+1);
    CV_data=E_vec(1:DSI:end);
    o=beast(CV_data,'start', 0,'season','none','torder.minmax', [1,1], 'tcp.minmax',    [0,10000]); %this is the activated BEAST algorithm on the CV data
        
    %this command plots the BEAST algorithm on the
    %data
    yyy=figure;
    xxx=plotbeast(o);
    h = findobj(xxx(1),'Type','line');
    Y=h(end-1,:).YData; 
    A={};
    cp=sort(o.trend.cp(1:o.trend.ncp_median));    
    if sum(isnan(cp))~=0
       cp=cp(1:(find (isnan(cp),1)-1));     
    end   
    cp=[1 ;cp; length(CV_data)];
    if ~isempty(find (cp==0))
    cp=cp(find (cp~=0)); 
    end     
    %Here, we extract relveant vectors out of the CVs segmented data
    %The skewness option is left for future possible research on its
    %contribution
    [mu_vec,std_vec,Skewness_vec,trend_vec,times_vec]=give_vecs(CV_data,Y,cp);
    cc={ mu_vec std_vec  trend_vec times_vec };
    AA=[ mu_vec ;std_vec;  trend_vec ;times_vec]';

    A=cc;
    %This is a structure that save the unormailizied CV segments data
    save(horzcat(save_here,'\',horzcat('Stochastic_Landscape_Unnormalized_Matrix',regexprep(num2str(ok, '%5.0f'),'\.','_'))),'A');
    M_reduced=AA(:,1:3);
    Szz=size(AA,1);
    Mm_reduced=(M_reduced-mean(M_reduced,1))./std(M_reduced,1);
    [coeff,score,latent] = pca(Mm_reduced);
    %This is a structure that save the noramlizied+ PCA CV segments data
    save(horzcat(save_here,'\',horzcat('Stochastic_Landscape_PCA_Matrix',regexprep(num2str(ok, '%5.0f'),'\.','_'))),'coeff','score','latent','Mm_reduced');
    
    %The figure below plots the stochastic landscape of the segmented data 
    figgg=figure; 
    bbb=get(figgg,'Position');
    new_width=8.7;
    set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
    X=score;
    %This uses the elbow method to decide on DBscan heuristic parameters,
    %using the knee properity of the NN heuristics, displayed here:
    
    %Sander, J.; Ester, M.; Kriegel, H.-P.; Xu, X. Density-based clustering in spatial
    %databases: The algorithm gdbscan and its applications. Data mining and knowledge
    %discovery 1998, 2, 169–194

    %Kaplan, D. Knee Point. https://www.mathworks.com/matlabcentral/
    %fileexchange/35094-knee-point, 2024; [Online; accessed January 22, 2024]
    
    elbow_search=zeros(1,5);
    kk=6;
    Idx = knnsearch(X,X,'K',(kk));
    distancee=zeros(kk-1,size(Idx,1));
    distance=zeros(1,size(Idx,1));
    for jj=2:1:(kk)
        for ii=1:1:size(Idx,1)  
        distancee(jj-1,ii)=norm(X(ii,:)-X(Idx(ii,jj),:));   
        end
    end
    for gg=1:1:size(distancee,2)
        distance(gg)=mean(distancee(:,gg));
    end
    distance=sort(distance);
    ee=knee_pt(1:1:length(distance),distance);
    CV=CV_data;
    T1 = dbscan(X, ee,5);
    scatter3(X(:,1),X(:,2),X(:,3),10,T1,'filled');
    xlabel('y1');
    ylabel('y2');
    zlabel('y3');
    camorbit(90, 0);
    set(gca,'FontSize',6)    
    Namexs=horzcat('SLM_Clustering_',Protein_name,'_New_',num2str(ok));
    print (horzcat(save_here,'\',Namexs),'-dpng','-r300');
    
    %This figure plots the yielded stochastic landscape without the
    %clustering 
    figgg=figure; 
    bbb=get(figgg,'Position');
    new_width=8.7;
    set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
    scatter3(X(:,1),X(:,2),X(:,3),10, gca().Colormap(1,:),'filled');   
    %title('Result of Clustering');
    xlabel('y1');
    ylabel('y2');
    zlabel('y3');
    camorbit(90, 0);
    set(gca,'FontSize',6)
    Namexs=horzcat('SLM_Clustering_',Protein_name,'_New_No_Color',num2str(ok));
    print (horzcat(save_here,'\',Namexs),'-dpng','-r300');
    
    %The code below plots the CV with the post Munkres algorithm
    %clustering, as described by:
    
    %Kuhn, H. W. The Hungarian method for the assignment problem. Naval research logis-
    %tics quarterly 1955, 2, 83–97.
    %Chen, Yewang, et al. "Decentralized clustering by finding loose and distributed density cores." Information Sciences 433 (2018): 510-526.

    PostPCA_coor=[score AA(:,4) ];
    timec=[];  

    for ii=1:1:Szz
    pivot=T1(ii).*ones(1,AA(ii,4));

    timec=[timec pivot];
    end
    CV=CV_data;
    c = [timec NaN]; % Must end in NaN to filling a solid
    if length(c)~=length(ground_truth_vec)
        c=c(1:end-1);
        c(end)=NaN;
        CV=CV(1:end-1);
        CV(end)=NaN;
        CV_data=CV_data(1:end-1);
        CV_data(end)=NaN;
    end

    c_updated=-ones(1,length(c));
    c_updated(end)=1;
    num_of_states_minkres=length(unique(c))-2;
    M=zeros(num_of_states_minkres,N_states);
    c=c';
      
    for ii=1:1:num_of_states_minkres
           M(ii,1)=length(find (ground_truth_vec==1 & c'==ii));
           M(ii,2)=length(find (ground_truth_vec==2 & c'==ii));
    end
    
    Khun=munkres(-M);
    new_labels_vec=zeros(1, length(N_states));
    for jj=1:1:N_states
       if ~isempty(find(Khun(:,jj)==1))
       new_labels_vec(jj)=find(Khun(:,jj)==1);
       else
       new_labels_vec(jj)=nan;
       end
    end
    new_labels_vec=new_labels_vec(~isnan(new_labels_vec));
    for ii=1:1:length(c_updated)
        if c(ii)~=-1 && ~isnan(c(ii)) && c(ii)<N_states+1
            c_updated(ii)= new_labels_vec(c(ii));
        end    
    end
    
    %Here, we calculate the classification metric values comparred with the
    %ground-truth
    NMI=nmi(c_updated,ground_truth_vec);
    [RI,ARI,Dice,JD]=randindex(c_updated, ground_truth_vec);
    FM=clustering_comparison(c_updated, ground_truth_vec );
    Truth_mat(ok,:)=[NMI,RI,ARI,Dice,FM, JD];
    CV_data=CV;
    %This figure plots the CV with the clustered new colors
    figgg=figure; 
    bbb=get(figgg,'Position');
    new_width=8.7;
    set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
    c_updated(c_updated==-1)=0;
    scatter((1:1:length(CV_data))*DSI,CV_data,1,c_updated);
    ylabel('CV');
    xlabel('Time');
    set(gca,'FontSize',6);
    Namexs=horzcat('CV_vs_time_',Protein_name,'_New_',num2str(ok));
    % Get the handle of the scatter plot (assuming one scatter plot)
    scatter_handle = findobj(gcf, 'Type', 'scatter');
    % Extract unique colors from the scatter plot
    A=colororder;
    CCOLOR=A(1:length(unique(c_updated)),:);
    % % Create a colormap with segments for each unique color
    colormap(CCOLOR);
    unique_colors=size(CCOLOR,1);
    % Create the colorbar
    cb = colorbar();
    % Set tick marks at the boundaries between colors
    tick_positions = linspace(0.5, numel(unique_colors)+0.5, numel(unique_colors) + 1);
    set(cb, 'Ticks', tick_positions, 'TickLabels', []);  % Remove default tick labels
    % Add custom tick labels below the colorbar
    tick_labels = 0:(unique_colors-1);  % Convert colors to strings
    yticks = get(cb, 'YTick');
    range_me=cb.Limits(2)-cb.Limits(1);
    tick_positions=(1:1:unique_colors).*range_me/unique_colors-0.5.*range_me/unique_colors;
    set(cb, 'Ticks', tick_positions, 'TickLabels', tick_labels);
    print (horzcat(save_here,'\',Namexs),'-dpng','-r300');
    CV=CV_data;
    given_vec=E_vec;
    save(horzcat(save_here,'\',horzcat(Namexs,'.mat')),'CV','DSI','c_updated','given_vec')

    %This figure plots the CV with ticks on the trend-changepoints found by
    %the BEAST algorithm

    figgg=figure; 
    bbb=get(figgg,'Position');
    new_width=8.7;
    set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
    plot(1:1:length(CV_data),CV_data,'.', 'LineWidth', 0.5);
    ylimw=ylim;
    ylimits=[min(CV_data)-1 ylimw(2)];
    hold on;
    for i = 1 : length(cp)
        s=scatter(cp(i),ylimits(1),[],[0.8500 0.3250 0.0980],"|",'filled','MarkerEdgeColor',[0.8500 0.3250 0.0980]);
        s.SizeData = 20;
        ylabel('CV')
        xx0=xlabel('MC Steps');
        ylim([min(CV_data)-2 ylimits(2)])
    end 

    xlabel('Time','FontSize',6)  ;    
    ylabel('CV','FontSize',6);
    set(gca,'FontSize',6)     
    Namexs=horzcat('BEAST_RAW_',Protein_name,'_New_',num2str(ok));
    print (horzcat(save_here,'\',Namexs),'-dpng','-r300');
    
    %This figure plots the segmented CV data with the cluster correspnding colors

    figgg=figure; 
    bbb=get(figgg,'Position');
    h_factor=bbb(3)/bbb(4);
    new_width=8.7;
    set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 0.5*new_width]);
    patch((1:1:length(CV_data))*DSI,CV_data,c, 'EdgeColor', 'interp', 'LineWidth', 1);
    title('Segment changes along y axis');
    hold on;
    for i = 1 : length(cp)
        plot( [cp(i),cp(i)]*DSI, get(gca,'Ylim'), 'color', 'k');
    end   
    ylabel('CV');
    xlabel('Time');
    set(gca,'FontSize',6);
    Namexs=horzcat('BEAST_Summary_',Protein_name,'_New_',num2str(ok));
    print (horzcat(save_here,'\',Namexs),'-dpng','-r300');

    %This part of the code saves the run information\ MATLAB workspace
    save(horzcat(save_here,'\',Namexs,'_workspace')); 
    clearvars -except Truth_mat save_here DSI ground_truth_vec ground_truth_vec_name N_states M CV_Matrix_name CV_MAT Protein_name; 
    load(ground_truth_vec_name);
    load(CV_Matrix_name);
    close all

end

%this prints the classifcation metric values, comparing between the
%ground-truth to the CV data

truth_mat=Truth_mat;
%print the header
fprintf('CV  NMI  RI  ARI  Dice  FM  Jaccard\n');
% print the results
for r=1:size(Truth_mat,1)
    fprintf('%.2f %.2f %.2f %.2f %.2f %.2f\n', truth_mat(r,1), truth_mat(r,2), truth_mat(r,3), truth_mat(r,4), truth_mat(r,5),truth_mat(r,6));
end
