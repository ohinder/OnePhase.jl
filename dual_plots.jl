
####################
###################

function summarize_by_iteration(hist::Array{alg_history,1})
    sum_hist = Array{alg_history,1}()
    for i = 2:length(hist)
      if hist[i].t > hist[i-1].t
        push!(sum_hist, hist[i-1])
      end
    end
    push!(sum_hist, hist[end])

    return sum_hist
end

sh = summarize_by_iteration(hist)
data_us = Dict()
field_list = Dict("max dual variable" => :y_norm, "primal feas" => :primal_residual, "dual feas" => :norm_grad_lag, "complementarity" => :dot_sy);

for (label, field) in field_list
    data_us[field] = zeros(length(sh));
end

for i = 1:length(sh)
    for (label, field) in field_list
      data_us[field][i] = getfield(sh[i],field);
    end
end

using PyPlot

subplot(121)
title("A well-behaved IPM")
ylabel("L-infinity norm")
xlabel("iterations")

for (label, field) in field_list
  semilogy(data_us[field],label=label)
end
ax = gca()
ylims = [1e-9, 1e10]
ax[:set_ylim](ylims)

legend()

subplot(122)
title("IPOPT with perturbed constraints")
xlabel("iterations")

for (label, field) in field_list
  semilogy(ds_ipopt[field],label=label)
end
ax = gca()
ax[:set_ylim](ylims)
