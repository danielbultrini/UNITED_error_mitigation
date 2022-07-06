tag = '_LinGrow';
%{
Q=4;  p=4;
for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end

Q=6;  p=6;

for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end

Q=8;  p=8;

for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end
%}


Q=10;  p=10;

for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end



%
Q=4;  p=64;
for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end

%
Q=6;  p=96;

for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end

%
Q=8;  p=128;

for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end


%{
Q=10;  p=160;

for seed = 0:29
    for nlsp = 1:3
        for budget_no = 1:4
            check_shot_noise_budget(budget_no, Q, p, seed, nlsp, tag)
        end
    end
end
%}


        