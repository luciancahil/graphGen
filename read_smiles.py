from rdkit import Chem
from rdkit.Chem import Draw
import networkx as nx

f = open("SMILES.txt", "r")

s = [line.strip() for line in f.readlines()]
print(s)

def rdkit_to_networkx(mol):
    """Convert an RDKit molecule to a NetworkX graph."""
    G = nx.Graph()
    
    # Add nodes with atom properties
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), atomic_num=atom.GetAtomicNum(), symbol=atom.GetSymbol())
    
    # Add edges with bond properties
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType())
    
    return G

# Example usage
smiles = s[0]
mol = Chem.MolFromSmiles(smiles)
G = rdkit_to_networkx(mol)

# Display the nodes and edges
print("Nodes:", G.nodes(data=True))
print("Edges:", G.edges(data=True))