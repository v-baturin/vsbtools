from ase.formula import Formula

def stoich2formula(stoich, el_symbols):
    formula = ''
    for n_at, symbol in zip(stoich, el_symbols):
        if n_at == 1:
            formula += symbol
        elif n_at > 1:
            formula += symbol + str(int(n_at))
    return formula



def get_chemical_formula_custom_order(atoms, ordered_elements):
    fla_dict = Formula(atoms.get_chemical_formula()).count()
    result = ''
    for el in ordered_elements:
        if el in fla_dict:
            result += el
            cnt = fla_dict[el]
            result += str(int(cnt)) if cnt !=1 else ''
    return result
