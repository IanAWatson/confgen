# Conformer generator copied from
# https://iwatobipen.wordpress.com/2021/01/31/generate-conformers-script-with-rdkit-rdkit-chemoinformatics/
# Input option is now a smiles file, and all smiles in the file are processed.

import click

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdmolops

def make_conformers(mol, writer, prunermsthresh, numconf, add_ref):
    mol = Chem.AddHs(mol, addCoords=True)
    refmol = mol
    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = prunermsthresh
    cids = rdDistGeom.EmbedMultipleConfs(mol, numconf, param)
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, mmffVariant='MMFF94s')
    if add_ref:
        refmol.SetProp('CID', '-1')
        refmol.SetProp('Energy', '')
        writer.write(refmol)
    res = []
 
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))
    sorted_res = sorted(res, key=lambda x:x[1])
    rdMolAlign.AlignMolConformers(mol)
    for cid, e in sorted_res:
        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))
        writer.write(mol, confId=cid)
 
@click.command()
@click.option('--input', '-i', help='inputfile smiles file', required=True)
@click.option('--output', '-o', help='output sdf file path', default='gen_confs.sdf')
@click.option('--prunermsthresh', '-t', default=0.1, type=float, help='RMS threshold for conformer retention')
@click.option('--numconf', default=50, type=int, help='Number of conformations to produce (def 50)')
@click.option('--add_ref', '-r', default=False, type=bool, help='Also write the starting molecule')
@click.option('--rmchiral/--normchiral', default=False, help='Remove chirality before processing')
@click.option('--verbose', '-v', count=True, help='Verbosity')
def confgen(input, output, prunermsthresh, numconf, add_ref, rmchiral, verbose):
  """Performs 2D-3D conversion and conformer generation on smiles input"""
  writer = Chem.SDWriter(output)
  with open(input, 'r') as input:
    for line in input:
      m = Chem.MolFromSmiles(line)
      if m is None:
        print(f'Cannot parse {line}, ignored')
        continue
      
      if verbose > 0:
        print(f'Processing {m.GetProp("_Name")}')

      if rmchiral:
        rdmolops.RemoveStereochemistry(m)

      make_conformers(m, writer, prunermsthresh, numconf, add_ref)

  writer.close()

if __name__=='__main__':
    confgen()
