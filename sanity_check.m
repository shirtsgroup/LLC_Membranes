function sanity_check(S, mat, traj_no, linked)

data = load(mat);  % load original data

if nargin == 4
    traj = data.linked_data;
else
    traj = data.traj;  % (npoints x ntraj)
end

nT = size(S.stateSeq(traj_no).z, 2);  % nT < len(traj)
t = linspace(1, nT, nT); 

plot(traj(:, traj_no, 3));  % plot original data
hold on;

estimated_states = S.stateSeq(traj_no).z;  % estimated state sequence
found_states = unique(estimated_states);  % unique estimated state labels

fprintf('\nFound %d states\n\n', size(found_states, 2))

if nargin  == 4
    true_states = data.linked_state_labels(:, traj_no);
else
    true_states = data.labels(:, traj_no);
end

true_state_labels = unique(true_states);

shift = size(true_states, 1) - size(estimated_states, 2);

mismatches = [];

subsets = nchoosek(found_states, size(true_state_labels, 1));

% handles case where less states are found than were used to generate data
dummy_states = [];
if size(subsets, 1) == 0
    diff = size(true_state_labels, 1) - size(found_states, 2);
    for i=1:diff
        dummy = sum(found_states);
        found_states = [found_states, dummy];  % generate unique fake states
        dummy_states = [dummy_states, dummy];
    end
    subsets = zeros(1, size(found_states, 2));
    subsets(1, :) = found_states(1, :);
end

nsubsets = size(subsets, 1);
nperm = size(perms(subsets(1, :)), 1);
mismatches = zeros(nsubsets, nperm);  %k, i

for k=1:nsubsets
    p = perms(subsets(k, :));
    for i=1:size(p, 1)
        M = containers.Map(true_state_labels, p(i, :));
        wrong = 0;  % number of states misidentified with this combo of state labels
        for s=1:size(estimated_states, 2)
            if estimated_states(s) ~= M(true_states(s + shift))
                wrong = wrong + 1;
            end
        end
        mismatches(k, i) = wrong;
    end
end

[subset, perm] = find(mismatches==min(mismatches(:)));

[subset, perm] = find(mismatches==min(mismatches(:)));
p = perms(subsets(subset, :));
states = p(perm, :);
    
M = containers.Map(true_state_labels, states);
reverseM = containers.Map(states, true_state_labels);

fake_states = found_states(find(0==ismember(found_states, states)));

wrong_label = [];
for s=1:size(estimated_states, 2)
    if estimated_states(s) ~= M(true_states(s + shift))
        wrong_label = [wrong_label , s + shift];
    end
end

fprintf('%.1f percent of states identified correctly\n\n', 100*(nT - size(wrong_label, 2)) / nT)

scatter(t(wrong_label-shift), traj(wrong_label, traj_no, 3), 100, 'r', 'x')

colors = rand(3, max(found_states));
colorSeq = zeros(3, nT);

for i=1:nT
    colorSeq(:, i) = colors(:, S.stateSeq(traj_no).z(i));
end

if nargin == 4
    phantom_state = M(M.Count - 1); % identify phantom state
    states(states==phantom_state) = [];  % remove phantom state
end

for i=find(1==ismember(states, dummy_states))
    states(states==states(i)) = [];
end

% transition matrix
T = S.dist_struct(traj_no).pi_z(states, states);
% make rows of transition matrix sum to 1
for i=1:size(states)
    T(i, :) = T(i, :) / sum(T(i, :));
end
disp('Estimated Transition Matrix:')
disp(T)

disp('Actual Transition Matrix:')
disp(data.T)

% calculate error from actual matrix
if size(T, 1) == size(data.T, 1)
    err = (T - data.T).^2;
    mean(mean(err))
    disp(norm(err))
else
    disp('Too few states were found, so no error is reported')
end

scatter(t, traj(1:nT, traj_no, 3), 25, colorSeq', 'filled')

figure
legend_labels = {};
for i=1:size(found_states, 2)
    plot([0, 1], [i, i], 'color', colors(:, found_states(i)), 'LineWidth', 5);
    hold on;
    if nargin == 4 & found_states(i) == phantom_state
        legend_labels{i} = 'Phantom State';
    elseif ismember(found_states(i), fake_states)
        legend_labels{i} = 'Fake State';
    else
        legend_labels{i} = sprintf('State %d', reverseM(found_states(i)));
    end
end
ylim([0, size(found_states, 2) + 1])
legend(legend_labels)
