
using Revise

using Tiger
using JuMP

# =================== TESTING ===================

cobra = read_cobra("/Users/jensen/Dropbox/repos/COBRA_models/iSMU.mat", "iSMU")
cobra.lb[cobra.lb .> 0.0] .= 0.0  # remove NGAM
model = build_base_model(cobra)
optimize!(model)

set_media_bounds!(model, "./test/CDM.toml")
optimize!(model)

add_gprs_cnf!(model, cobra)
optimize!(model)

#set_silent(model)
#@time single_deletions(model, variable_by_name.(model, cobra.genes))

ext = extend_cobra_cnf(cobra)
m_ext = build_base_model(ext)
set_media_bounds!(m_ext, "./test/CDM.toml")
optimize!(m_ext)

