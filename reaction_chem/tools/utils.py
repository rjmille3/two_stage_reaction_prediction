from openeye import oechem
from rpCHEM.Common.Util import molBySmiles, smi_to_unique_smi_fast, mol_to_unique_smi_fast, exact_mass

def calculate_molecular_weight(smi):
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smi)
    return exact_mass(mol)