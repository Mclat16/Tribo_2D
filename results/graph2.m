colors = {"#0072BD", "#D95319", "#EDB120","#7E2F8E"};
% Define force values
force = [2,5,10,20,30,50,70,85,100,130];
temproom = 300;
materials = {'b-AlSb','b-antimonene','b-arsenene','b-bismuthene'};
angle = 0:4:88;
% Initialize empty arrays to store data

% Loop through each force value
for t= 1:4
    m = materials{t};
    disp(m)
    cf = zeros(length(force),1);
    for i = 1:length(force)
        % Build filename based on force value
        filename = [m '/afm/100x_100y/sub_aSi_tip_Si_r25/K300/results/fc_ave_slide_' num2str(force(i)) 'nN_0angle_2ms_l1' ];
        %filename = [m '/afm/100x_100y/sub_aSi_tip_Si_r25/K300/results/fc_ave_slide_30nN_' num2str(angle(i)) 'angle_2ms_l1' ];

        % Open the file for reading
        fid = fopen(filename, 'r');
        formatspec = '%f %f %f %f %f %f %f %f %f %f'; % Skip columns other than 4 and 5
        % Read the data, skipping the first line (header)
        data = textscan(fid, formatspec, 'HeaderLines', 2);
        fid = fclose(fid);
        time = data{1};% first column
        nf = data{2};
        lfx = data{3};
        lfy = data{4};
        comz = data{7};
        tipx = data{8};
        tipy = data{9};
        lf = sqrt(lfx.^2 + lfy.^2);


        ave_nf = abs(mean(nf));
        ave_lf = abs(mean(comz));
        

        %cf(i) = ave_lf/ave_nf;
        cf(i) = ave_lf;
    end


    % Plot cf vs force on the same graph
    plot(force, cf,'-o','MarkerEdgeColor',colors{t},'MarkerSize', 10, 'LineWidth', 2, 'DisplayName',m);
    hold on

    title('Displacement (A) vs. Load (nN) at 30nN NF');
    legend;


end

hold off
set(gcf,'units','points','position',[10,10,830,420])
set(gca,'FontSize',26)
grid minor
xlabel('Normal force (nN)');
ylabel('COF');


