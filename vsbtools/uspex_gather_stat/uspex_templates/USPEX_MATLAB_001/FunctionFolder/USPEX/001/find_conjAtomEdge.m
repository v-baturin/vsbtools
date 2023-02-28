function [Atom_conj, Edge_conj] = find_conjAtomEdge(Structure, Atom, Edge)

global AR_VARIANTS

Atom_conj = AR_VARIANTS.Structure(Structure).Atom(Atom).Edge(Edge).Atom;

Edge_conj = -1;
for k = 1 : length(AR_VARIANTS.Structure(Structure).Atom(Atom_conj).Edge)
    if AR_VARIANTS.Structure(Structure).Atom(Atom_conj).Edge(k).Atom == Atom
        Edge_conj = k;
        break;
    end
end