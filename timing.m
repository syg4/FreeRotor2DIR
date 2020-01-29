states = [1 2 3 4 5 6 7 8 9 10 12 15 17 20];
times = zeros(1, length(states));

for ii = 1:length(states)
    tic
    [transDipoleMatrixZ, transDipoleMatrixX, transDipoleMatrixY, Cz, Az, Cx, Ax, Cy, Ay] = getCO2RoVibTransDipoleMatrix(2, states(ii));
    times(ii) = toc
end

%%
fig1 = figure(1);
clf
plot(((states+1).^2.*3).^2, times, 'LineWidth', 2);
% fig1.Position = [435 126 1120 840];
box off
xlabel('Total Number of M States');
ylabel('Run Time (s)');

% ax1 = gca; % current axes
% 
% 
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position', ax1_pos, 'XAxisLocation','top','YAxisLocation','right','Color', 'none', 'XColor', 'b');


% line(states,times,'Parent',ax2,'Color','b')
% xlabel('Total Number of J States');
% ylabel('Run Time (min)')

%%
states_vec = [1 2 3 4 5 6 7 8 9 10 12 15 17 20 25 30 35 40 45];
times_vec = zeros(1, length(states_vec));

for ii = 1:length(states_vec)
    tic
    [transDipoleMatrixZ, transDipoleMatrixX, transDipoleMatrixY, Cz, Az, Cx, Ax, Cy, Ay] = getCO2RoVibTransDipoleMatrix_Vectorized(2, states_vec(ii));
    times_vec(ii) = toc
end

%%
% fig2 = figure(2);
hold on
plot(((states_vec+1).^2.*3).^2, times_vec, 'LineWidth', 2);
% fig2.Position = [435 126 1120 840];

xlabel('Total Number of M States');
ylabel('Run Time (s)');
legend({'For loop', 'Array Function/Vectorization'}, 'Location', 'northwest')
% ax1 = gca; % current axes
% 
% 
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position', ax1_pos, 'XAxisLocation','top','YAxisLocation','right','Color', 'none', 'XColor', 'b');


% line(states_vec,times_vec./60,'Parent',ax2,'Color','b')
% xlabel('Total Number of J States');
% ylabel('Run Time (min)')
