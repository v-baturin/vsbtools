def binding_energy_1comp(n_vs_e_formatted, start=0, stop=None):
    all_nats = n_vs_e_formatted[:, 0]
    all_e_s = n_vs_e_formatted[:, 1]
    e_1 = all_e_s[all_nats == 1]
    n_at_s = all_nats[slice(start, stop)]
    e_s = all_e_s[slice(start, stop)]
    return (e_s - n_at_s * e_1) / n_at_s, n_at_s