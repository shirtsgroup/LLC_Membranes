function plot_states(S, mat, traj_no, linked)

data = load(mat);  % load original data
if nargin == 4
    traj = data.linked_data;
else
    traj = data.traj;  % (npoints x ntraj)
end

nT = size(S.stateSeq(traj_no).z, 2);  % nT < len(traj)
t = linspace(1, nT, nT); 

%first_order_difference = traj(2:nT, traj_no, 3) - traj(1:nT-1, traj_no, 3);
%plot(first_order_difference)
plot(traj(:, traj_no, 3));  % plot original data
hold on;

estimated_states = S.stateSeq(traj_no).z;  % estimated state sequence
states = unique(estimated_states);  % unique estimated state labels

% if nargin == 4
%     fprintf('\nFound %d states\n\n', size(states, 2))
% 
%     true_states = data.labels(:, traj_no);
%     true_state_labels = unique(true_states);
% 
%     mismatches = [];
%     p = perms(states);
%     for i=1:size(p, 1)
%         M = containers.Map(true_state_labels, p(i, :));
%         wrong = 0;  % number of states misidentified with this combo of state labels
%         for s=1:size(estimated_states, 2)
%             if estimated_states(s) ~= M(true_states(s + 2))
%                 wrong = wrong + 1;
%             end
%         end
%         mismatches = [mismatches, wrong];
%     end
% 
%     M = containers.Map(true_state_labels, p(argmin(mismatches), :));
% 
%     wrong_label = [];
%     for s=1:size(estimated_states, 2)
%         if estimated_states(s) ~= M(true_states(s + 2))
%             wrong_label = [wrong_label , s + 2];
%         end
%     end
%     
%     fprintf('%.1f percent of states identified correctly\n\n', 100*(nT - size(wrong_label, 2)) / nT)
%     scatter(t(wrong_label), traj(wrong_label, traj_no, 3), 100, 'r', 'x')
% end

colors = rand(3, max(states));
colorSeq = zeros(3, nT);

for i=1:nT
    colorSeq(:, i) = colors(:, S.stateSeq(traj_no).z(i));
end

% transition matrix
T = S.dist_struct(traj_no).pi_z(states, states);
disp('Estimated Transition Matrix:')
disp(T)

scatter(t, traj(1:nT, traj_no, 3), 4, colorSeq')
%scatter(t(1:nT-1), first_order_difference, 4, colorSeq(:, 1:nT-1)')
