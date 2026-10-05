const β = 0.92
const r = 0.017 
const θ = 0.282 
const y_const = -0.5
const y_curv = 0.53

const ρy = 0.945 
const σy = 0.025
const stds_out = 3.0 
const nquad = 15

const b_L = 0.0 
const b_H = 0.6
const q_L = 0.0
const q_H = 1/(1+r)

const nb = 201 
const ny = 31

const σy_uncond = σy/sqrt(1-ρy^2)
const y_L = exp(-stds_out*σy_uncond)
const y_H = exp(stds_out*σy_uncond)

const bgrid = range(b_L, b_H, length = nb);
const ygrid = range(y_L, y_H, length = ny);

const gh_nodes, gh_weights = gausshermite(nquad)
const y_nodes = √2 * σy * gh_nodes
const y_weights = gh_weights / √π
