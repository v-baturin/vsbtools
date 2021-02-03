def stoich2formula(stoich, el_symbols):
    formula = ''
    for n_at,symbol in zip(stoich, el_symbols):
        if n_at == 1:
            formula += symbol
        elif n_at > 1:
            formula += symbol + str(int(n_at))
    return formula