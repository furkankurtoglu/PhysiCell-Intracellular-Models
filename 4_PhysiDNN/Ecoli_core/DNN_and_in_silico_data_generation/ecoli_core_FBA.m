load("ecoli_core_model.mat")



num_of_met_run = 200;
glucose_bounds = -1 * linspace(0, 10, num_of_met_run);
oxygen_bounds = -1 * linspace(0, 1000, num_of_met_run);

sol_vector = zeros(1,5);

save_data = 'Y';


for i = 1:num_of_met_run
    for j = 1:num_of_met_run
        oxy_lb = oxygen_bounds(j);
        glc_lb = glucose_bounds(i);
        model.lb(28) = glc_lb;
        model.lb(36) = oxy_lb;
        sol = optimizeCbModel(model);
        if ~(strcmp(sol.origStat,'INFEASIBLE'))
            sol_vector = [-1*oxy_lb, -1*glc_lb, -1*sol.x(36), -1*sol.x(28), sol.f];
            if save_data == "Y"
                writematrix(sol_vector,'ecoli_in_silico_data.csv',Delimiter=',',WriteMode='append');
            end
        end
    end
end

fclose('all');