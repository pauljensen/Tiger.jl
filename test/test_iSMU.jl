
using Tiger

# =================== TESTING ===================

b = parse_boolean("a & (b | c) & d")

print(
    dnf(b)
)

print(
    cnf(b)
)

cobra = read_cobra("/Users/jensen/Dropbox/repos/COBRA_models/iSMU.mat", "iSMU")
model = build_base_model(cobra)
optimize!(model)

print(atoms(b))
print(names(b))
