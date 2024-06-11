%%% LAB 08 - SCIENTIFIC DATA VISUALIZATION
clear all
close all
clc

load Data_Lab03.mat

% divide male and female subjects
CT_male = CT(find(gender == 1),:);
CT_female = CT(find(gender == 2),:);

% divide left and right hemisphere
CT_left = CT(:,1:31);
CT_right = CT(:,32:end);

% both
CT_male_left = CT_male(:,1:31);
CT_female_left = CT_female(:,1:31);
CT_male_right = CT_male(:,32:end);
CT_female_right = CT_female(:,32:end);

%% Comparison left and right CT for ROI n. 1 with boxplot
ROI = ROIs(1);
figure(1)
boxplot([CT_left(:,1), CT_right(:,1)], 'Labels', {'left', 'right'})
ylabel('CT [mm]')
xlabel('hemisphere')
title('Comparison CT: left and right hemisphere for ' + ROI)

%% Comparison left and right CT for every ROI with boxplot
for i = 1:3
    ROI = ROIs(i);
    figure(i)
    boxplot([CT_left(:,i), CT_right(:,i)], 'Labels', {'left', 'right'})
    ylabel('CT [mm]')
    xlabel('hemisphere')
    title('Comparison CT: left and right hemisphere for ' + ROI)
end

saveas(figure(1), 'caudal anterior cingulate', 'jpg')
saveas(figure(2), 'caudal middle frontal', 'jpg')
saveas(figure(3), 'cuneus', 'jpg')


%% Comparison left and right CT for every ROI with plot
avg_left = zeros(1,31);
avg_right = zeros(1,31);
for i = 1:31
    ROI = ROIs(i);
    avg_left(i) = mean(CT_left(:,i));
    avg_right(i) = mean(CT_right(:,i));
end

figure
plot(1:31, avg_left, 1:31, avg_right, 'LineWidth', 2)
title('Comparison average CT: left and right hemispheres')
xlabel('ROIs')
ylabel('CT [mm]')
xlim([1, 31])
legend('average CT left', 'average CT right')

%% Difference left and right CT for every ROI with plot/scatter
difference = avg_left-avg_right;
x_left = [0; 0; 31; 31];
y_left = [0; 0.3; 0.3; 0];
x_right = [0; 0; 31; 31];
y_right = [-0.3; 0; 0; -0.3];

idx_left = find(difference > 0);
idx_right = find(difference < 0);

f_diff = figure('Position', [50 100 1000 400]); hold on
stem(1:31, difference, 'k', 'filled')
yline(0)
patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
title('Comparison averaged CT: left and right hemispheres')
xlabel('ROIs')
ylabel('CT difference left-right [mm]')
xlim([1, 31])
ylim([-0.25, 0.25])
legend('CT_{left} - CT_{right}', 'threshold', 'CT_{left} bigger', 'CT_{right} bigger', ...
    'Location', 'eastoutside')
ax = gca;
xticks(1:31)
xticklabels(ROIs)
for l = 1:length(idx_left)
    ax.XTickLabel{idx_left(l)} = ['\color{blue}' ax.XTickLabel{idx_left(l)}];
end
for r = 1:length(idx_right)
    ax.XTickLabel{idx_right(r)} = ['\color{red}' ax.XTickLabel{idx_right(r)}];
end

saveas(f_diff, 'difference', 'jpg')

%% Comparison CT of left hemisphere between male and female subjects
left_hist = figure('Position', [50 100 500 400]);
hold on
histogram(CT_male_left(:,1))
histogram(CT_female_left(:,1))
legend('male', 'female')
title('CT_{left} for ' + ROI(1))
xlabel('CT [mm]')
ylabel('number of subjects')

right_hist = figure('Position', [50 100 500 400]);
hold on
histogram(CT_male_right(:,1))
histogram(CT_female_right(:,1))
legend('male', 'female')
title('CT_{right} for ' + ROI(1))
xlabel('CT [mm]')
ylabel('number of subjects')

saveas(left_hist, 'left_gender', 'jpg')
saveas(right_hist, 'right_gender', 'jpg')

%% Comparison CT of both hemispheres between male and female subjects

mat_diff_male = CT_male_left - CT_male_right;
mat_diff_female = CT_female_left - CT_female_right;

difference_male = mean(mat_diff_male);
difference_female = mean(mat_diff_female);

std_diff_male = std(mat_diff_male);
std_diff_female = std(mat_diff_female);

[sorted_diff_male, idx_male] = sort(difference_male);
[sorted_diff_female, idx_female] = sort(difference_female);

ROIs_male = ROIs(idx_male);
ROIs_female = ROIs(idx_female);

mat_diff_male = mat_diff_male(:, idx_male);
mat_diff_female = mat_diff_female(:, idx_female);

boxchart_diff = figure('Position', [50 70 700 700]), hold on
% male
boxchart(mat_diff_male, 'MarkerStyle', 'none', 'BoxFaceColor', 'b')
ylim([-0.8, 0.8])
% female
mat_diff_female_visual = mat_diff_female(:, idx_male);
boxchart(mat_diff_female_visual, 'MarkerStyle', 'none', 'BoxFaceColor', 'r')
% visualization
view(-90,90)
xlabel('ROIs')
ylabel('CT difference [mm]')
ylim([-0.8, 0.8])
title('Difference CT_{left} - CT_{right} for male and female')
xticklabels(ROIs_male)
yline(0)
legend('male', 'female', 'zero difference')

saveas(boxchart_diff, 'difference boxchart', 'jpg')

plot_diff = figure('Position', [50 70 700 700]); hold on
title('Difference CT_{left} - CT_{right} for male and female')
xlabel('CT difference [mm]')
yyaxis left
plot(mean(mat_diff_male), 1:31, 'b')
ylabel('ROIs')
yticks(1:31)
yticklabels(ROIs_male)
ylim([1,31])
ax = gca;
ax.YColor = 'b';
yyaxis right
plot(mean(mat_diff_female), 1:31, 'r')
ylabel('ROIs')
yticks(1:31)
yticklabels(ROIs_female)
ylim([1,31])
ax = gca;
ax.YColor = 'r';
legend('male', 'female', 'Location', 'southeast')

saveas(plot_diff, 'difference plot', 'jpg')



% mat_map = [mean(CT_male_right); mean(CT_female_right)];
% figure
% heatmap(mat_map)
