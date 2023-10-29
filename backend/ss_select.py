from rdkit.Chem import AllChem as Chem
from backend import smiles_syntax
apply_glue = smiles_syntax.convert_labeled_wildcards_to_closures
import re
import copy

# Every query needs a "starting material" and "product", either user-defined, or None if not user-defined
# The content from the sketcher can come in as a rxnfile or molfile
def resolve_sketched_content(sketched_content: str):

    if '$RXN' in sketched_content:
        components = sketched_content.split('$MOL\n')
        header = components[0]
        # Rxnfile keeps track of whether or not starting material, product, and reagents exist using the below pattern searched for as a regex
        # For example, '  1  0' means starting material but no product, and '  0  1  1' means product and reactants but no starting material
        # (The third number only shows up as 1 if a structure is drawn above the arrow, and we ignore it; should we send back a message telling the user it was ignored?)
        components_inventory = re.search('(\n  [0-1]  [0-1]\n|\n  [0-1]  [0-1]  1\n)', header).group(0)[1:-1]
        starting_material, product = int(components_inventory[2]), int(components_inventory[5])
        assert starting_material or product

        if starting_material and product:
            return (components[1], components[2])
        elif starting_material:
            return (components[1], None)
        elif product:
            return (None, components[1])
    else:
        # If the user doesn't draw a reaction arrow, treat the drawn structure as a starting material
        return (sketched_content, None)

def prepare_mol(molfile):
    mol = Chem.MolFromMolBlock(molfile, sanitize=False) if molfile is not None else None
    if mol is not None:
        # Sanitization seems to be necessary for the mol.GetRingInfo().AtomRings() function to find the rings
        # However, sanitizing at the Chem.MolFromMolBlock function stage seems to remove explicit H
        # But calling Chem.SanitizeMol on the mol seems to NOT remove explicit H, while still enabling ring finding
        Chem.SanitizeMol(mol)
    return mol

def prepare_reaction_mols(sketched_content):
    starting_material, product = resolve_sketched_content(sketched_content)
    starting_material = Chem.MolFromMolBlock(starting_material, sanitize=False) if starting_material is not None else None
    product = Chem.MolFromMolBlock(product, sanitize=False) if product is not None else None
    for mol in (starting_material, product):
        if mol is not None:
            # Sanitization seems to be necessary for the mol.GetRingInfo().AtomRings() function to find the rings
            # However, sanitizing at the Chem.MolFromMolBlock function stage seems to remove explicit H
            # But calling Chem.SanitizeMol on the mol seems to NOT remove explicit H, while still enabling ring finding
            Chem.SanitizeMol(mol)
    return starting_material, product

# Often the user may select atoms that will not lead to any results, based on how the MMP data is organized. This isn't the user's fault, we need to help them out
# Therefore, below, we try and catch any selection that will never lead to results, and modify the selection to the closest selection that will generate a fruitful search
# E.g., If any ring contains the selected atoms, include all atom indices from all such rings together with the original selected atoms (because we didn't split open rings in our MMP definitions)
# If user has drawn a reaction, and has previous selections on one side of the arrow, and is making a new selection on the opposite side, then preserve the previous selection
def validate_selection(selection_type: str, sketched_content: dict, mol1_selected_atoms: str, mol1_selected_bonds: str, mol2_selected_atoms: str, mol2_selected_bonds: str):

    mol1 = Chem.MolFromMolBlock(sketched_content.mol1_molfile, sanitize=False) if sketched_content.mol1_molfile != '' else None
    mol2 = Chem.MolFromMolBlock(sketched_content.mol2_molfile, sanitize=False) if sketched_content.mol2_molfile != '' else None
    for mol in (mol1, mol2):
        if mol is not None:
            # Sanitization seems to be necessary for the mol.GetRingInfo().AtomRings() function to find the rings
            # However, sanitizing at the Chem.MolFromMolBlock function stage seems to remove explicit H
            # But calling Chem.SanitizeMol on the mol seems to NOT remove explicit H, while still enabling ring finding
            Chem.SanitizeMol(mol)

    extracted_sketched_content = {}
    for key, value in sketched_content:
        if key not in ('mol1_molfile', 'mol2_molfile'):
            split = value.split(',')
            extracted_sketched_content[key] = [int(x) for x in split] if split != [''] else []
        else:
            extracted_sketched_content[key] = value
    sketched_content = extracted_sketched_content

    mol1_selected_atoms = set([int(x) for x in mol1_selected_atoms.split(',')]) if mol1_selected_atoms != '' else set()
    mol1_selected_bonds = set([int(x) for x in mol1_selected_bonds.split(',')]) if mol1_selected_bonds != '' else set()
    mol2_selected_atoms = set([int(x) for x in mol2_selected_atoms.split(',')]) if mol2_selected_atoms != '' else set()
    mol2_selected_bonds = set([int(x) for x in mol2_selected_bonds.split(',')]) if mol2_selected_bonds != '' else set()

    selection_dict = {'mol1': {}, 'mol2': {}}

    for name, mol, selected_atoms, selected_bonds in (('mol1', mol1, mol1_selected_atoms, mol1_selected_bonds), ('mol2', mol2, mol2_selected_atoms, mol2_selected_bonds)):

        if mol is None:
            selection_dict[name]['entire_molecule_selected'] = "False"
            selection_dict[name]['autoenv_atoms'] = []
            selection_dict[name]['selected_atoms'] = []
            selection_dict[name]['expanded_bond_group'] = []
            continue

        ### Disabling this feature with ketcher, because it feels more accurate when lassoing to not include atoms touched by bonds
        # Sometimes the user may select bonds without selecting all attached atoms
        # Find all atoms touched by selected bonds, and make sure all these atoms are included in selected_atoms
        """
        atoms_touched_by_bonds = set()
        for bond_index in selected_bonds:
            bond = mol.GetBondWithIdx(bond_index)
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                atoms_touched_by_bonds.add(atom.GetIdx())
        selected_atoms = selected_atoms.union(atoms_touched_by_bonds)
        """

        if selection_type == 'variable':

            rings = mol.GetRingInfo().AtomRings()
            start_length = len(selected_atoms)
            end_length = 0
            skipped_indices = ['this is here so length evaluates to > 0 on first iteration of the below while loop']

            # By default we don't cut rings during MMP definition, so automatically expand the user's selection to include entire rings if they select only part of a ring
            # Expand the selected atoms with every loop iteration, to include an entire ring if even one atom of that ring is in the selection
            # There can be a domino effect where expanding selection to an entire ring then causes another ring to be included in the selection
            while (end_length != start_length) and len(skipped_indices) > 0:
                skipped_indices = []
                start_length = len(selected_atoms)
                for i, ring in enumerate(rings):
                    ring = set(ring)
                    if selected_atoms == ring:
                        continue
                    elif len(selected_atoms.intersection(ring)) > 0:
                        selected_atoms = selected_atoms.union(ring)
                    else:
                        skipped_indices.append(i)
                rings = [rings[i] for i in skipped_indices]
                end_length = len(selected_atoms)
    
            # If explicit H are adjacent to a variable atom, include that H
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    if len(set([a.GetIdx() for a in atom.GetNeighbors()]).intersection(selected_atoms)) > 0:
                        selected_atoms = selected_atoms.union(set([atom.GetIdx()]))

            # Find single atoms that border the variable selection
            autoenv_atoms = []
            for atom in mol.GetAtoms():
                if atom.GetIdx() not in selected_atoms:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in selected_atoms:
                            autoenv_atoms.append(atom.GetIdx())
                            break
            selection_dict[name]['autoenv_atoms'] = [i for i in autoenv_atoms]

        elif selection_type == 'environment':
            selection_dict[name]['autoenv_atoms'] = []
            
        # Find all bonds that are completely enclosed by the expanded set of selected_atoms
        expanded_bond_group = []
        for bond in mol.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            if begin_atom_idx in selected_atoms and end_atom_idx in selected_atoms:
                expanded_bond_group.append(bond.GetIdx())

        selection_dict[name]['selected_atoms'] = [i for i in selected_atoms]
        selection_dict[name]['expanded_bond_group'] = expanded_bond_group

        #if len(selected_atoms) == len(mol.GetAtoms()):
        selection_dict[name]['entire_molecule_selected'] = "True" if len(selected_atoms) == len(mol.GetAtoms()) else "False"

    mol1_variable_atoms, mol1_variable_bonds, mol1_environment_atoms, mol1_environment_bonds, mol1_entire_molecule_selected = [], [], [], [], selection_dict['mol1']['entire_molecule_selected']
    mol2_variable_atoms, mol2_variable_bonds, mol2_environment_atoms, mol2_environment_bonds, mol2_entire_molecule_selected = [], [], [], [], selection_dict['mol2']['entire_molecule_selected']

    for mol, variable_atoms, variable_bonds, environment_atoms, environment_bonds, sketched_variable_atoms, sketched_variable_bonds, sketched_environment_atoms, sketched_environment_bonds in (
        ('mol1', mol1_variable_atoms, mol1_variable_bonds, mol1_environment_atoms, mol1_environment_bonds, sketched_content['mol1_variable_atoms'], sketched_content['mol1_variable_bonds'], sketched_content['mol1_environment_atoms'], sketched_content['mol1_environment_bonds']),
        ('mol2', mol2_variable_atoms, mol2_variable_bonds, mol2_environment_atoms, mol2_environment_bonds, sketched_content['mol2_variable_atoms'], sketched_content['mol2_variable_bonds'], sketched_content['mol2_environment_atoms'], sketched_content['mol2_environment_bonds'])
    ):
        # Include the previous variable/environment selections, if user did not make any new selection on the same side of the reaction arrow
        if selection_dict[mol]:
            if not selection_dict[mol]['selected_atoms']:
                variable_atoms += sketched_variable_atoms
                variable_bonds += sketched_variable_bonds
                environment_atoms += sketched_environment_atoms
                environment_bonds += sketched_environment_bonds
            elif selection_type == 'environment':
                if not set(selection_dict[mol]['selected_atoms']).intersection(set(sketched_variable_atoms)):
                    variable_atoms += sketched_variable_atoms
                    variable_bonds += sketched_variable_bonds

            if selection_type == 'variable':
                variable_atoms += selection_dict[mol]['selected_atoms']
                variable_bonds += selection_dict[mol]['expanded_bond_group']
                environment_atoms += selection_dict[mol]['autoenv_atoms']
            elif selection_type == 'environment':
                environment_atoms += selection_dict[mol]['selected_atoms']
                environment_bonds += selection_dict[mol]['expanded_bond_group']

    return (mol1_variable_atoms, mol1_variable_bonds, mol1_environment_atoms, mol1_environment_bonds, mol1_entire_molecule_selected,
            mol2_variable_atoms, mol2_variable_bonds, mol2_environment_atoms, mol2_environment_bonds, mol2_entire_molecule_selected)

# This function is designed for rxnfile input, with MarvinJS
def validate_rxnfile_selection(sketched_content: str, selection_type: str, selected_atoms: str, selected_bonds: str, variable_atoms: str, variable_bonds: str, environment_atoms: str, environment_bonds: str):

    starting_material, product = prepare_reaction_mols(sketched_content)

    # Need to decrement atom indices by 1, because MarvinJS starts from 1, but RDKit starts from 0
    selected_atoms = set(map(lambda x: x-1, [int(x) for x in selected_atoms.split(',')])) if selected_atoms != '' else set()
    variable_atoms = set(map(lambda x: x-1, [int(x) for x in variable_atoms.split(',')])) if variable_atoms != '' else set()
    environment_atoms = set(map(lambda x: x-1, [int(x) for x in environment_atoms.split(',')])) if environment_atoms != '' else set()
    # Marvin bond notation takes the form of atom_index-atom_index,atom_index-atom_index,...
    # We need to unpack these bonds to lists of integers that we can work with, and decrement atom indices as above

    variable_bonds = list(map(lambda x: x.split('-'), variable_bonds.split(','))) if variable_bonds else []
    environment_bonds = list(map(lambda x: x.split('-'), environment_bonds.split(','))) if environment_bonds else []
    for bonds in variable_bonds, environment_bonds:
        for i, (x,y) in enumerate(bonds):
            bonds[i] = [int(x)-1, int(y)-1]

    # Sometimes the user may select bonds without selecting all attached atoms
    # Find all atoms touched by selected bonds, and make sure all these atoms are included in selected_atoms
    atoms_touched_by_bonds = []
    # Marvin bond notation takes the form of atom_index-atom_index,atom_index-atom_index,...
    for bond in selected_bonds.split(','):
        atoms_touched_by_bonds += bond.split('-')
    atoms_touched_by_bonds = set(map(lambda x: x-1, [int(x) for x in atoms_touched_by_bonds])) if selected_bonds != '' else set()
    selected_atoms = selected_atoms.union(atoms_touched_by_bonds)

    # We use sm_index_end to figure out whether atom indices are in the starting_material or product
    if starting_material is not None:
        sm_index_end = len(starting_material.GetAtoms()) if starting_material is not None else 0
        sm_atoms = set(i for i in selected_atoms if i < sm_index_end)
    else:
        sm_index_end = 0
        sm_atoms = None

    if product is not None:
        product_index_end = sm_index_end + len(product.GetAtoms())
        product_atoms = set(i - sm_index_end for i in selected_atoms if i >= sm_index_end and i < product_index_end)
    else:
        product_atoms = None

    selection_dict = {'starting_material': {}, 'product': {}}
    entire_molecule_selected = "False"

    for name, mol, selected_atoms, index_start in (('starting_material', starting_material, sm_atoms, 0), ('product', product, product_atoms, sm_index_end)):

        if mol is None:
            continue
        elif selection_type == 'variable':

            rings = mol.GetRingInfo().AtomRings()
            start_length = len(selected_atoms)
            end_length = 0
            skipped_indices = ['this is here so length evaluates to > 0 on first iteration of the below while loop']

            # By default we don't cut rings during MMP definition, so automatically expand the user's selection to include entire rings if they select only part of a ring
            # Expand the selected atoms with every loop iteration, to include an entire ring if even one atom of that ring is in the selection
            # There can be a domino effect where expanding selection to an entire ring then causes another ring to be included in the selection
            while (end_length != start_length) and len(skipped_indices) > 0:
                skipped_indices = []
                start_length = len(selected_atoms)
                for i, ring in enumerate(rings):
                    ring = set(ring)
                    if selected_atoms == ring:
                        continue
                    elif len(selected_atoms.intersection(ring)) > 0:
                        selected_atoms = selected_atoms.union(ring)
                    else:
                        skipped_indices.append(i)
                rings = [rings[i] for i in skipped_indices]
                end_length = len(selected_atoms)
    
            # If explicit H are adjacent to a variable atom, include that H
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    if len(set([a.GetIdx() for a in atom.GetNeighbors()]).intersection(selected_atoms)) > 0:
                        selected_atoms = selected_atoms.union(set([atom.GetIdx()]))

            # Find single atoms that border the variable selection
            autoenv_atoms = []
            for atom in mol.GetAtoms():
                if atom.GetIdx() not in selected_atoms:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in selected_atoms:
                            autoenv_atoms.append(atom.GetIdx())
                            break
            selection_dict[name]['autoenv_atoms'] = [i + index_start for i in autoenv_atoms]

        elif selection_type == 'environment':
            selection_dict[name]['autoenv_atoms'] = []
            
        # Find all bonds that are completely enclosed by the expanded set of selected_atoms
        expanded_bond_group = []
        for bond in mol.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            if begin_atom_idx in selected_atoms and end_atom_idx in selected_atoms:
                # Need to increment atom indices by 1, because we're passing the atom indices back to MarvinJS
                # We need to adjust the indices by adding '+ index_start' so that they match the molfile/rxnfile atom indices exported from the sketcher
                expanded_bond_group.append(str(min(begin_atom_idx + index_start + 1, end_atom_idx + index_start + 1)) + '-' + str(max(begin_atom_idx + index_start + 1, end_atom_idx + index_start + 1)))

        selection_dict[name]['selected_atoms'] = [i + index_start for i in selected_atoms]
        selection_dict[name]['expanded_bond_group'] = expanded_bond_group

        if len(selected_atoms) == len(mol.GetAtoms()):
            entire_molecule_selected = "True"

    # Include the previous variable/environment selections, if user did not make any new selection on the same side of the reaction arrow
    validated_variable_atoms, validated_variable_bonds, validated_environment_atoms, validated_environment_bonds = [], [], [], []


    if selection_dict['starting_material']:
        if not selection_dict['starting_material']['selected_atoms']:
            validated_variable_atoms += [i for i in variable_atoms if i < sm_index_end]
            validated_variable_bonds += list(map(lambda x: '-'.join([str(x[0] + 1), str(x[1] + 1)]), [index_pair for index_pair in variable_bonds if index_pair[0] < sm_index_end]))
            validated_environment_atoms += [i for i in environment_atoms if i < sm_index_end]
            validated_environment_bonds += list(map(lambda x: '-'.join([str(x[0] + 1), str(x[1] + 1)]), [index_pair for index_pair in environment_bonds if index_pair[0] < sm_index_end]))
        elif selection_type == 'environment':
            if not set(selection_dict['starting_material']['selected_atoms']).intersection(set(variable_atoms)):
                validated_variable_atoms += [i for i in variable_atoms if i < sm_index_end]
                validated_variable_bonds += list(map(lambda x: '-'.join([str(x[0] + 1), str(x[1] + 1)]), [index_pair for index_pair in variable_bonds if index_pair[0] < sm_index_end]))

    if selection_dict['product']:
        if not selection_dict['product']['selected_atoms']:
            validated_variable_atoms += [i for i in variable_atoms if i >= sm_index_end]
            validated_variable_bonds += list(map(lambda x: '-'.join([str(x[0] + 1), str(x[1] + 1)]), [index_pair for index_pair in variable_bonds if index_pair[0] >= sm_index_end]))
            validated_environment_atoms += [i for i in environment_atoms if i >= sm_index_end]
            validated_environment_bonds += list(map(lambda x: '-'.join([str(x[0] + 1), str(x[1] + 1)]), [index_pair for index_pair in environment_bonds if index_pair[0] >= sm_index_end]))
        elif selection_type == 'environment':
            if not set(selection_dict['product']['selected_atoms']).intersection(set(variable_atoms)):
                validated_variable_atoms += [i for i in variable_atoms if i >= sm_index_end]
                validated_variable_bonds += list(map(lambda x: '-'.join([str(x[0] + 1), str(x[1] + 1)]), [index_pair for index_pair in variable_bonds if index_pair[0] >= sm_index_end]))

    for name in selection_dict:
        if selection_dict[name]:
            if selection_type == 'variable':
                validated_variable_atoms += selection_dict[name]['selected_atoms']
                validated_variable_bonds += selection_dict[name]['expanded_bond_group']
                validated_environment_atoms += selection_dict[name]['autoenv_atoms']
            elif selection_type == 'environment':
                validated_environment_atoms += selection_dict[name]['selected_atoms']
                validated_environment_bonds += selection_dict[name]['expanded_bond_group']

    # Need to increment atom indices by 1, because we're passing the atom indices back to MarvinJS
    validated_variable_atoms = map(lambda x: x+1, validated_variable_atoms)    
    validated_variable_atoms = ','.join([str(a) for a in validated_variable_atoms])
    validated_variable_bonds = ','.join(validated_variable_bonds)

    validated_environment_atoms = map(lambda x: x+1, validated_environment_atoms)
    validated_environment_atoms = ','.join([str(a) for a in validated_environment_atoms])
    validated_environment_bonds = ','.join(validated_environment_bonds)

    return validated_variable_atoms, validated_variable_bonds, validated_environment_atoms, validated_environment_bonds, entire_molecule_selected

# label variable atoms with map number = 1
# label constant atoms with map number = 2
def attach_maps(molfile, variable_maps, constant_maps):

    # RDKit cartridge will throw error if we use A as wildcard, so we need to use *
    molfile = molfile.replace(' A ', ' * ')

    # The maps are strings of atom indices selected by the user, e.g. 1,2,7
    variable_maps = variable_maps.split(',')
    for a,b in enumerate(variable_maps):
        variable_maps[a] = int(b)

    if constant_maps != '':
        constant_maps = constant_maps.split(',')
        for a,b in enumerate(constant_maps):
            constant_maps[a] = int(b)
    else:
        constant_maps = []

    num_patt = '([0-9][0-9][0-9]| [0-9][0-9]|  [0-9])'
    # the last number in the num_patt*10 is the one we need to change to the desired map
    atom_line_patt = '([A-Z]|\*)([a-z]| )' + num_patt*10

    for i, match in enumerate(re.finditer(atom_line_patt, molfile)):
        middle = '  0'

        if i in variable_maps:
            middle = '  1'
        elif i in constant_maps:
            middle = '  2'
        else:
            continue

        beginning = molfile[0:match.end() - 3]
        end = molfile[match.end():]
        molfile = beginning + middle + end

    return molfile


# functions find_border_bonds through extract_variable_constant have initially been developed
# to extract desired variable and constant regions, based on the assumption that variable atoms
# will be labeled with atom map = 1, and constant atoms will be labeled with atom map = 2

def find_border_bonds(mol):
    
    varAtoms = []
    constAtoms = []

    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 1:
            varAtoms.append(atom)
        elif atom.GetAtomMapNum() == 2:
            constAtoms.append(atom)

    var_const_contact_bonds = []
    var_unlabeled_contact_bonds = []

    # Find bonds where the variable atoms (marked with atom map number 1) makes contact with the rest of the molecule
    # The bond indexes are loaded to above lists
    for atom in varAtoms:
        for neighbor in atom.GetNeighbors():
            # Don't look at dummy atoms (isotope > 0) that are there because of a previous fragmentation event
            if neighbor.GetIsotope() == 0:
                # Constant atoms are marked with atom map number 2
                if neighbor.GetAtomMapNum() == 2:
                    var_const_contact_bonds.append(mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetIdx())
                # By default, atom map number is 0 if there is no atom map number assigned
                elif neighbor.GetAtomMapNum() == 0:
                    var_unlabeled_contact_bonds.append(mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetIdx())

    const_unlabeled_contact_bonds = []

    # Find bonds where the constant atoms (marked with atom map number 2) makes contact with the rest of the molecule
    for atom in constAtoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIsotope() == 0:
                if neighbor.GetAtomMapNum() == 0:
                    const_unlabeled_contact_bonds.append(mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetIdx())
                    
    return {
        'variable_constant' : var_const_contact_bonds,
        'variable_unlabeled' : var_unlabeled_contact_bonds,
        'constant_unlabeled' : const_unlabeled_contact_bonds
    }

def cleave_bonds(mol, bond_list, dummyLabels=None, customLabelSet=None):
    num_bonds = len(bond_list)
    if num_bonds > 0:
        if dummyLabels:
            mol = Chem.FragmentOnBonds(mol, bond_list, dummyLabels = [dummyLabels]*num_bonds)
        elif customLabelSet:
            mol = Chem.FragmentOnBonds(mol, bond_list, dummyLabels = customLabelSet[0:num_bonds])
    return mol
        
def find_next_isotope_atom(RWFrag, isotope):
    for RWAtom in RWFrag.GetAtoms():
        if RWAtom.GetIsotope() == isotope:
            return RWAtom.GetIdx()

def find_lone_aromatic_atom_indices(mol, only_first = False):
    lone_aromatic_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() == True:
            atom_is_lone_aromatic = True
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic() == True:
                    atom_is_lone_aromatic = False
                    break
            if atom_is_lone_aromatic == True:
                lone_aromatic_atoms.append(atom.GetIdx())
                if only_first == True:
                    return lone_aromatic_atoms
    return lone_aromatic_atoms

def find_lone_generic_atom_index(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and len(atom.GetNeighbors()) == 0:
            return atom.GetIdx()


def extract_variable_constant(molfile1, molfile2):

    mol1, mol2 = [prepare_mol(molfile) for molfile in (molfile1, molfile2)]
    observations = {
        #"all_atoms_selected_as_variable" : "False",
        #"variable_selection_over_limit" : "False",
        "max_variable_atoms" : 20 
    }

    # Check if the whole molecule is selected as variable
    # Also check if number of selected variable atoms is greater than max_variable_atoms
    # Current DB settings limit variable fragments to 20 heavy atoms at most for >=2 cuts, and 15 atoms for single cut
    max_variable_atoms = observations["max_variable_atoms"]

    for mol in (mol1, mol2):
        if mol is None:
            continue
        total_atoms = len(mol.GetAtoms())
        num_selected_variable_atoms = 0
        num_selected_heavy_atoms = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                num_selected_variable_atoms += 1
                if atom.GetAtomicNum() != 1:
                    num_selected_heavy_atoms += 1
        if num_selected_variable_atoms == total_atoms:
            observations["all_atoms_selected_as_variable"] = "True"
        if num_selected_heavy_atoms > max_variable_atoms:
            observations["variable_selection_over_limit"] = "True"

    # Abort if the query will never give results
    if "all_atoms_selected_as_variable" in observations or "variable_selection_over_limit" in observations:
        return (None, None, None, observations)

    extracted_content = {'mol1': {}, 'mol2': {}}

    for name, mol in (('mol1', mol1), ('mol2', mol2)):
        if mol is None:
            continue

        mol = cleave_bonds(mol, find_border_bonds(mol)['variable_constant'], customLabelSet = ((1,1), (2,2), (3,3)))
        mol = cleave_bonds(mol, find_border_bonds(mol)['variable_unlabeled'], dummyLabels = (6, 6))
        mol = cleave_bonds(mol, find_border_bonds(mol)['constant_unlabeled'], dummyLabels = (7, 7))
    
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                if atom.GetIsotope() in [1,2,3,5,6]:
                    # currently using Rb (atomic number 37) as dummy atom in database for connection point
                    atom.SetAtomicNum(37)
                    if atom.GetIsotope() not in [1,2,3]:
                        # We need to 'remember' if the cut occurred on a variable-constant border
                        # In such a case (where isotope is in [1,2,3]), we keep track of the matching atoms between variable and constant, based on the isotope
                        atom.SetIsotope(0)
                #elif atom.GetIsotope() == 7:
                #atom.SetAtomicNum(1)
                    #atom.SetIsotope(0)
                
        variable_frags, constant_frags = [], []

        for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
            for atom in frag.GetAtoms():
                if atom.GetAtomMapNum() == 1:
                    variable_frags.append(frag)
                    break
                elif atom.GetAtomMapNum() == 2:
                    # Check if any attachments to non-variable atoms (isotope number 7), if so remove them
                    RWFrag = Chem.RWMol(frag)
                    while True:
                        delete_index = find_next_isotope_atom(RWFrag, 7)
                        if delete_index is not None:
                            RWFrag.RemoveAtom(delete_index)
                        else:
                            break
                    constant_frags.append(RWFrag)
                    break
                elif atom.GetAtomicNum() != 37:
                    break

        for mol_frags in [variable_frags, constant_frags]:
            while len(mol_frags) > 1:
                grouped = Chem.CombineMols(mol_frags[0], mol_frags[1])
                del mol_frags[1]
                mol_frags[0] = grouped
    
        # This block is designed specifically for when there are lone aromatic atoms in the environment
        # The code will ensure that the output constant fragment molblock encodes the atom as aromatic
        # Aromaticity is encoded in the bond block in molfiles, not the atom block. Therefore we need to add a second aromatic atom and build an aromatic bond
        if len(constant_frags) > 0:
            constant_mol = Chem.RWMol(constant_frags[0])
            # For each lone aromatic atom, we need to add a generic aromatic neighbor, for recognition by the query
            num_lone_aromatics = len(find_lone_aromatic_atom_indices(constant_mol))
            for i in range(num_lone_aromatics):
                constant_mol = Chem.CombineMols(constant_mol, Chem.MolFromSmarts('[c,n,o,s]'))
            constant_mol = Chem.RWMol(constant_mol)

            # Now start fusing the lone aromatic atoms with the newly added generic aromatic atoms
            for i in range(num_lone_aromatics):
                lone_aromatic_index = find_lone_aromatic_atom_indices(constant_mol, only_first=True)[0]
                lone_generic_index = find_lone_generic_atom_index(constant_mol)
                constant_mol.AddBond(lone_aromatic_index, lone_generic_index)
                constant_mol.GetBondBetweenAtoms(lone_aromatic_index, lone_generic_index).SetIsAromatic(True)
                constant_mol.GetBondBetweenAtoms(lone_aromatic_index, lone_generic_index).SetBondType(Chem.BondType.AROMATIC)

            constant_smiles_allBondsExplicit = Chem.MolToSmiles(constant_mol, allBondsExplicit=True)
        else:
            constant_mol = None

        wildcard_frag_smiles = ['*[Rb]', '*[Rb].*[Rb]', '*[Rb].*[Rb].*[Rb]']
        one_atom_frag_smiles = ['C[Rb]', '*c[Rb]', 'O[Rb]', 'N[Rb]', 'F[Rb]', 'S[Rb]', 'Cl[Rb]']

        ### TODO: Fix bug where drawing A -> B, and highlighting variable atoms in A (but nothing in B) leads to an exception being thrown on below line:
        for atom in variable_frags[0].GetAtoms():
            atom.SetAtomMapNum(0)
        variable_mol = variable_frags[0]
        variable_pkl = variable_frags[0].ToBinary()
        variable_smiles = Chem.MolToSmiles(variable_frags[0])
        #print(variable_smiles)
        # Validation for whether or not the query is too general -- these queries could return many gigabytes
        if variable_smiles in wildcard_frag_smiles:
            observations["has_wildcard_variable"] = "True"
        elif variable_smiles in one_atom_frag_smiles:
            observations["has_one_atom_variable"] = "True"

        #if len(constant_frags) > 0:
        if constant_mol is not None:
            for atom in constant_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            constant_smiles = Chem.MolToSmiles(constant_mol)

            # Using allBondsExplicit=True so that cases such as a single aromatic environment atom, instead of getting converted from mol to *c[Rb] which loses its aromatic bond information upon later mol generation,
            # will instead be converted to *:c-[Rb], which is properly recognized by the Chem.MolFromSmiles(mol, sanitize=False) function
            #constant_smiles_allBondsExplicit = Chem.MolToSmiles(constant_mol, allBondsExplicit=True)
            #print('constant_smiles allBondsExplicit=True')
            #print(constant_smiles_allBondsExplicit)

            constant_pkl = constant_mol.ToBinary()
            # Seems necessary to convert to Smarts then back to mol, in order for the wildcard to be recognized - should investigate this further, to avoid a conversion
            # The below line also forces recognition of aromatic vs. aliphatic
            # The below line also causes SO2 recognition in constant fragment to fail
            #constant_pkl = Chem.MolFromSmarts(constant_smiles_allBondsExplicit).ToBinary()

            if constant_smiles in wildcard_frag_smiles:
                observations["has_wildcard_constant"] = "True"
            elif constant_smiles in one_atom_frag_smiles:
                observations["has_one_atom_constant"] = "True"
        else:
            constant_smiles = None
            constant_smiles_allBondsExplicit = None
            constant_pkl = None

        extracted_content[name]['variable'] = variable_mol
        extracted_content[name]['constant'] = constant_mol

    # We don't enforce constaints on what the user selects for the constant environment, except that it can't overlap with the variable atoms
    # The user might select different constant environments in mol1 vs. mol2; in this case, just use the larger selection (or the starting material selection if tied); should we change this?

    variable1, variable2 = extracted_content['mol1'].get('variable'), extracted_content['mol2'].get('variable')
    sm_constant, p_constant = extracted_content['mol1'].get('constant'), extracted_content['mol2'].get('constant')

    if sm_constant is not None and p_constant is not None:
        constant = sm_constant if len(sm_constant.GetAtoms()) >= len(p_constant.GetAtoms()) else p_constant
    # If user draws a reaction, and has an environment selection on one side (but not the other), then tell the user this isn't allowed
    # If any environment is selected on one side of the arrow, an environment must be selected on the other side of the arrow
    # TODO We can probably remove this restriction if we modify extract_variable_constant to apply the right tags to atoms in the variable that border non-environment atoms
    # TODO But this could introduce more bugs and is not a priority right now
    elif (sm_constant is not None or p_constant is not None) and (variable1 is not None and variable2 is not None):
        observations["missing_constant"] = "True"
        constant = None
    else:
        constant = sm_constant if sm_constant is not None else p_constant

    return (variable1, variable2, constant, observations)

# Works for rxnfiles
def extract_rxnfile_variable_constant(sketched_content):

    starting_material, product = prepare_reaction_mols(sketched_content)
    observations = {
        #"all_atoms_selected_as_variable" : "False",
        #"variable_selection_over_limit" : "False",
        "max_variable_atoms" : 20 
    }

    # Check if the whole molecule is selected as variable
    # Also check if number of selected variable atoms is greater than max_variable_atoms
    # Current DB settings limit variable fragments to 20 heavy atoms at most for >=2 cuts, and 15 atoms for single cut
    max_variable_atoms = observations["max_variable_atoms"]

    for mol in (starting_material, product):
        if mol is None:
            continue
        total_atoms = len(mol.GetAtoms())
        num_selected_variable_atoms = 0
        num_selected_heavy_atoms = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                num_selected_variable_atoms += 1
                if atom.GetAtomicNum() != 1:
                    num_selected_heavy_atoms += 1
        if num_selected_variable_atoms == total_atoms:
            observations["all_atoms_selected_as_variable"] = "True"
        if num_selected_heavy_atoms > max_variable_atoms:
            observations["variable_selection_over_limit"] = "True"

    # Abort if the query will never give results
    if "all_atoms_selected_as_variable" in observations or "variable_selection_over_limit" in observations:
        return (None, None, None, observations)

    extracted_content = {'starting_material': {}, 'product': {}}

    for name, mol in (('starting_material', starting_material), ('product', product)):
        if mol is None:
            continue

        mol = cleave_bonds(mol, find_border_bonds(mol)['variable_constant'], customLabelSet = ((1,1), (2,2), (3,3)))
        mol = cleave_bonds(mol, find_border_bonds(mol)['variable_unlabeled'], dummyLabels = (6, 6))
        mol = cleave_bonds(mol, find_border_bonds(mol)['constant_unlabeled'], dummyLabels = (7, 7))
    
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                if atom.GetIsotope() in [1,2,3,5,6]:
                    # currently using Rb (atomic number 37) as dummy atom in database for connection point
                    atom.SetAtomicNum(37)
                    if atom.GetIsotope() not in [1,2,3]:
                        # We need to 'remember' if the cut occurred on a variable-constant border
                        # In such a case (where isotope is in [1,2,3]), we keep track of the matching atoms between variable and constant, based on the isotope
                        atom.SetIsotope(0)
                #elif atom.GetIsotope() == 7:
                #atom.SetAtomicNum(1)
                    #atom.SetIsotope(0)
                
        variable_frags, constant_frags = [], []

        for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
            for atom in frag.GetAtoms():
                if atom.GetAtomMapNum() == 1:
                    variable_frags.append(frag)
                    break
                elif atom.GetAtomMapNum() == 2:
                    # Check if any attachments to non-variable atoms (isotope number 7), if so remove them
                    RWFrag = Chem.RWMol(frag)
                    while True:
                        delete_index = find_next_isotope_atom(RWFrag, 7)
                        if delete_index is not None:
                            RWFrag.RemoveAtom(delete_index)
                        else:
                            break
                    constant_frags.append(RWFrag)
                    break
                elif atom.GetAtomicNum() != 37:
                    break

        for mol_frags in [variable_frags, constant_frags]:
            while len(mol_frags) > 1:
                grouped = Chem.CombineMols(mol_frags[0], mol_frags[1])
                del mol_frags[1]
                mol_frags[0] = grouped
    
        # This block is designed specifically for when there are lone aromatic atoms in the environment
        # The code will ensure that the output constant fragment molblock encodes the atom as aromatic
        # Aromaticity is encoded in the bond block in molfiles, not the atom block. Therefore we need to add a second aromatic atom and build an aromatic bond
        if len(constant_frags) > 0:
            constant_mol = Chem.RWMol(constant_frags[0])
            # For each lone aromatic atom, we need to add a generic aromatic neighbor, for recognition by the query
            num_lone_aromatics = len(find_lone_aromatic_atom_indices(constant_mol))
            for i in range(num_lone_aromatics):
                constant_mol = Chem.CombineMols(constant_mol, Chem.MolFromSmarts('[c,n,o,s]'))
            constant_mol = Chem.RWMol(constant_mol)

            # Now start fusing the lone aromatic atoms with the newly added generic aromatic atoms
            for i in range(num_lone_aromatics):
                lone_aromatic_index = find_lone_aromatic_atom_indices(constant_mol, only_first=True)[0]
                lone_generic_index = find_lone_generic_atom_index(constant_mol)
                constant_mol.AddBond(lone_aromatic_index, lone_generic_index)
                constant_mol.GetBondBetweenAtoms(lone_aromatic_index, lone_generic_index).SetIsAromatic(True)
                constant_mol.GetBondBetweenAtoms(lone_aromatic_index, lone_generic_index).SetBondType(Chem.BondType.AROMATIC)

            constant_smiles_allBondsExplicit = Chem.MolToSmiles(constant_mol, allBondsExplicit=True)
        else:
            constant_mol = None

        wildcard_frag_smiles = ['*[Rb]', '*[Rb].*[Rb]', '*[Rb].*[Rb].*[Rb]']
        one_atom_frag_smiles = ['C[Rb]', '*c[Rb]', 'O[Rb]', 'N[Rb]', 'F[Rb]', 'S[Rb]', 'Cl[Rb]']

        ### TODO: Fix bug where drawing A -> B, and highlighting variable atoms in A (but nothing in B) leads to an exception being thrown on below line:
        for atom in variable_frags[0].GetAtoms():
            atom.SetAtomMapNum(0)
        variable_mol = variable_frags[0]
        variable_pkl = variable_frags[0].ToBinary()
        variable_smiles = Chem.MolToSmiles(variable_frags[0])
        #print(variable_smiles)
        # Validation for whether or not the query is too general -- these queries could return many gigabytes
        if variable_smiles in wildcard_frag_smiles:
            observations["has_wildcard_variable"] = "True"
        elif variable_smiles in one_atom_frag_smiles:
            observations["has_one_atom_variable"] = "True"

        #if len(constant_frags) > 0:
        if constant_mol is not None:
            for atom in constant_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            constant_smiles = Chem.MolToSmiles(constant_mol)

            # Using allBondsExplicit=True so that cases such as a single aromatic environment atom, instead of getting converted from mol to *c[Rb] which loses its aromatic bond information upon later mol generation,
            # will instead be converted to *:c-[Rb], which is properly recognized by the Chem.MolFromSmiles(mol, sanitize=False) function
            #constant_smiles_allBondsExplicit = Chem.MolToSmiles(constant_mol, allBondsExplicit=True)
            #print('constant_smiles allBondsExplicit=True')
            #print(constant_smiles_allBondsExplicit)

            constant_pkl = constant_mol.ToBinary()
            # Seems necessary to convert to Smarts then back to mol, in order for the wildcard to be recognized - should investigate this further, to avoid a conversion
            # The below line also forces recognition of aromatic vs. aliphatic
            # The below line also causes SO2 recognition in constant fragment to fail
            #constant_pkl = Chem.MolFromSmarts(constant_smiles_allBondsExplicit).ToBinary()

            if constant_smiles in wildcard_frag_smiles:
                observations["has_wildcard_constant"] = "True"
            elif constant_smiles in one_atom_frag_smiles:
                observations["has_one_atom_constant"] = "True"
        else:
            constant_smiles = None
            constant_smiles_allBondsExplicit = None
            constant_pkl = None

        extracted_content[name]['variable'] = variable_mol
        extracted_content[name]['constant'] = constant_mol

    # We don't enforce constaints on what the user selects for the constant environment, except that it can't overlap with the variable atoms
    # The user might select different constant environments in starting_material vs. product; in this case, just use the larger selection (or the starting material selection if tied); should we change this?

    variable1, variable2 = extracted_content['starting_material'].get('variable'), extracted_content['product'].get('variable')
    sm_constant, p_constant = extracted_content['starting_material'].get('constant'), extracted_content['product'].get('constant')

    if sm_constant is not None and p_constant is not None:
        constant = sm_constant if len(sm_constant.GetAtoms()) >= len(p_constant.GetAtoms()) else p_constant
    # If user draws a reaction, and has an environment selection on one side (but not the other), then tell the user this isn't allowed
    # If any environment is selected on one side of the arrow, an environment must be selected on the other side of the arrow
    # TODO We can probably remove this restriction if we modify extract_variable_constant to apply the right tags to atoms in the variable that border non-environment atoms
    # TODO But this could introduce more bugs and is not a priority right now
    elif (sm_constant is not None or p_constant is not None) and (variable1 is not None and variable2 is not None):
        observations["missing_constant"] = "True"
        constant = None
    else:
        constant = sm_constant if sm_constant is not None else p_constant

    return (variable1, variable2, constant, observations)

# Permute both variable1 and variable2 simultaneously to ensure same mappings from each to the same set of constant permutations
def permute_variables_constant(variable1=None, variable2=None, constant=None, has_labeled_wildcards=False):

    # prevent modification in place
    variable1, variable2, constant = tuple(map(copy.deepcopy, (variable1, variable2, constant)))

    # Example input for has_labeled_wildcards=False: mol whose structure is [Rb]CO[Rb]
    # Example input for has_labeled_wildcards=True: mol whose structure is [*:1]CO[*:2], which we convert to [Rb]CO[Rb]
    if has_labeled_wildcards:
        for wildcard_mol in (variable1, variable2, constant):
            if wildcard_mol is not None:
                for atom in wildcard_mol.GetAtoms():
                    if atom.GetAtomicNum() == 0 and atom.GetIsotope() in (1,2,3):
                        #atom.SetAtomMapNum(0)
                        atom.SetAtomicNum(37)

    assert variable1 is not None or variable2 is not None
    
    variable1_permutations = []
    variable1_permutable_indices = []
    variable2_permutations = []
    variable2_permutable_indices = []
    constant_permutations = []
    constant_permutable_indices = []
    num_cuts = 0
    
    if variable1:
        for atom in variable1.GetAtoms():
            if atom.GetAtomicNum() == 37:
                variable1_permutable_indices.append(atom.GetIdx())
        num_cuts = len(variable1_permutable_indices)
            
    if variable2:
        for atom in variable2.GetAtoms():
            if atom.GetAtomicNum() == 37:
                variable2_permutable_indices.append(atom.GetIdx())
        num_cuts = len(variable2_permutable_indices)
            
    if constant:
        for atom in constant.GetAtoms():
            if atom.GetAtomicNum() == 37:
                constant_permutable_indices.append(atom.GetIdx())
            
    assert num_cuts > 0 and num_cuts < 4
    # Atomic numbers: He = 2, Ne = 10, Ar = 18
    if num_cuts == 1:
        patterns = ((2,),)
    if num_cuts == 2:
        patterns = ((2, 10), (10, 2))
    if num_cuts == 3:
        patterns = ((2, 10, 18), (2, 18, 10), (10, 2, 18), (10, 18, 2), (18, 2, 10), (18, 10, 2))
        
    for i in range(len(patterns)):
        variable1_permutations.append(copy.deepcopy(variable1))
        variable2_permutations.append(copy.deepcopy(variable2))
        constant_permutations.append(copy.deepcopy(constant))
        
    # We need to order things correctly for the for loop and return statement; we can't start the for loop with an empty variable
    if variable1:
        variable_permutations, variable_permutable_indices = variable1_permutations, variable1_permutable_indices
        other_variable_permutations, other_variable_permutable_indices = variable2_permutations, variable2_permutable_indices
    elif variable2 and not variable1:
        variable_permutations, variable_permutable_indices = variable2_permutations, variable2_permutable_indices
        other_variable_permutations, other_variable_permutable_indices = variable1_permutations, variable1_permutable_indices
    
    for variable_copy, constant_copy, other_variable_copy, pattern in zip(variable_permutations, constant_permutations, other_variable_permutations, patterns):
        vc_map = {}
        for idx, new_atomic_num in zip(variable_permutable_indices, pattern):
            atom = variable_copy.GetAtomWithIdx(idx)
            atom.SetAtomicNum(new_atomic_num)
            # Atoms in variable touching the constant, and vice versa, have an isotope in (1,2,3)
            # The cv_map for each variable_copy will let us apply matching He, Ne, Ar to atoms with matching isotope in the corresponding constant_copy (if the user selected a constant portion)
            isotope = atom.GetIsotope()
            if isotope in (1,2,3):
                vc_map[isotope] = new_atomic_num
            atom.SetIsotope(0)
                
        if vc_map:
            for idx in constant_permutable_indices:
                atom = constant_copy.GetAtomWithIdx(idx)
                isotope = atom.GetIsotope()
                atom.SetAtomicNum(vc_map[isotope])
                atom.SetIsotope(0)
        elif constant is None:
            # Relevant when constant is completely detached from variable
            # In such a case, vc_map will be empty because no variable atoms were touching constant, but there could still be a constant that we want to search with, so we check if it's None
            constant_permutations = None
            
        if variable1 and variable2 and vc_map:
            # Above we mapped from variable1 to constant, here we map from constant to variable2, to ensure we can use the same single ordering of constant permutations for both variable1/variable2
            atomic_nums_not_in_vc_map = [a for a in pattern if a not in set(vc_map.values())]
            counter = 0
            for idx in other_variable_permutable_indices:
                atom = other_variable_copy.GetAtomWithIdx(idx)
                isotope = atom.GetIsotope()
                if isotope in (1,2,3):
                    atom.SetAtomicNum(vc_map[isotope])
                else:
                    atom.SetAtomicNum(atomic_nums_not_in_vc_map[counter])
                    counter += 1
                atom.SetIsotope(0)

        elif variable1 and variable2:
            for idx, new_atomic_num in zip(other_variable_permutable_indices, pattern):
                atom = other_variable_copy.GetAtomWithIdx(idx)
                atom.SetAtomicNum(new_atomic_num)
                atom.SetIsotope(0)
        else:
            other_variable_permutations = None
            
    if variable1:
        return (variable_permutations, other_variable_permutations, constant_permutations)
    elif variable2 and not variable1:
        return (other_variable_permutations, variable_permutations, constant_permutations)
    
def map_unique_permutations(mols):
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    smiles_ids = [i for i in range(len(smiles))]
    for i in range(len(mols)):
        for j in range(i+1, len(mols)):
            # This means we already found these were equal in separate comparisons
            if smiles_ids[i] == smiles_ids[j]:
                continue
            if smiles[i] == smiles[j]:
                smiles_ids[j] = smiles_ids[i]
    return smiles_ids

def extract_variable_and_enumeration_core(molfile):

    mol = Chem.MolFromMolBlock(molfile, sanitize=False)
    Chem.SanitizeMol(mol)

    bonds_touching_variable = find_border_bonds(mol)['variable_constant'] + find_border_bonds(mol)['variable_unlabeled']
    mol = cleave_bonds(mol, bonds_touching_variable, customLabelSet = ((1,1), (2,2), (3,3)))

    variable_frags = []
    non_variable_frags = []
    # Separate the variable fragments from the rest of the fragments
    for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
        found_variable_frag = False
        for atom in frag.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                found_variable_frag = True
                variable_frags.append(frag)
                break
        if found_variable_frag == False:
            non_variable_frags.append(frag)

    for frag_list in (variable_frags, non_variable_frags):
        # Combine the fragments into one molecule object
        while len(frag_list) > 1:
            grouped = Chem.CombineMols(frag_list[0], frag_list[1])
            del frag_list[1]
            frag_list[0] = grouped

    variable = variable_frags[0]
    enumeration_core = non_variable_frags[0]

    return variable, enumeration_core

def get_num_cuts(mol):
    num_cuts = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() in (1,2,3):
            num_cuts += 1
    return num_cuts

def wildcards_in_core(molfile):
    mol = Chem.MolFromMolBlock(molfile)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == 2:
            return True
    return False

def mapped_wildcards_to_nobles(mol):
    mol = copy.deepcopy(mol)
    mapping = {1:2, 2:10, 3:18} # Convert a wildcard labeled with map in (1,2,3) to He, Ne, or Ar
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            map_num = atom.GetAtomMapNum()
            if map_num in (1,2,3):
                atom.SetAtomMapNum(0)
                atom.SetAtomicNum(mapping[map_num])
    return mol

def isotopes_to_mapped_wildcards(mol):
    mol = copy.deepcopy(mol)
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        if isotope in (1,2,3) and atom.GetAtomicNum() == 0:
            atom.SetIsotope(0)
            atom.SetAtomMapNum(isotope)
    return mol

def nobles_to_mapped_wildcards(mol):
    mol = copy.deepcopy(mol)
    mapping = {2:1, 10:2, 18:3} # Convert a He, Ne, or Ar atom to wildcard with map in (1,2,3)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num in (2, 10, 18):
            atom.SetAtomicNum(0)
            atom.SetAtomMapNum(mapping[atomic_num])
    return mol

def clear_atom_maps(mol, non_wildcards_only=False):
    mol = copy.deepcopy(mol)
    for atom in mol.GetAtoms():
        if non_wildcards_only:
            if atom.GetAtomicNum() == 0:
                continue
        atom.SetAtomMapNum(0)
    return mol

def get_gluable_smiles(mol):
    converted_smiles = ""
    for root_atom_index in range(-1, mol.GetNumAtoms()):
        # Some smiles variations will not successfully undergo closure identification with smiles_syntax; try various smiles representations until one works
        # Wildcards from query, *, should not be mistaken for attachment point atoms, therefore convert query wildcards to something else
        # Also convert the He/Ne/Ar dummy atoms to labeled wildcards for closure generation; we needed He/Ne/Ar for RDKit cartridge differentiation of multicut dummy atoms
        candidate_smiles = Chem.MolToSmiles(mol, rootedAtAtom=root_atom_index).replace('He', '*:1').replace('Ne', '*:2').replace('Ar', '*:3').replace('[*:7]', '[Rn:7]')
        try:
            converted_smiles = smiles_syntax.convert_labeled_wildcards_to_closures(candidate_smiles)
            break
        except Exception as e:
            pass
    if not converted_smiles:
        raise Exception(f"No way to represent this smiles with closures using smiles_syntax dependency: {candidate_smiles}")

    return candidate_smiles

def get_multifragment_gluable_smiles(mol):

    canonical_smiles = Chem.MolToSmiles(mol, rootedAtAtom=-1)
    # The below code will fail in some cases with multifrag molecules, but will succeed on individual fragments, therefore work with fragments
    frag_mols = [Chem.MolFromSmiles(frag, sanitize=False) for frag in canonical_smiles.split('.')]
    frag_smiles_list = []
    for frag in frag_mols:
        frag_smiles_list.append(get_gluable_smiles(frag))

    return '.'.join(frag_smiles_list)

def apply_transforms_to_input(variable, core, transform_data, aggregation_type, num_cuts):

    #core = isotopes_to_mapped_wildcards(clear_atom_maps(core, non_wildcards_only=True))
    enumerations = {}
    assert aggregation_type in ('individual_transforms', 'group_by_fragment') and num_cuts > 0 and num_cuts < 4
    if num_cuts == 1:
        core = isotopes_to_mapped_wildcards(clear_atom_maps(core))

        # The presence of a mapped wildcard in the middle of a SMILES string, for a single-cut fragment, sometimes
        # causes smiles_syntax.convert_labeled_wildcards_to_closures to raise NotImplementedError("intermediate groups not supported" 
        core = Chem.MolToSmiles(clear_atom_maps(core))
        core_to_weld = smiles_syntax.convert_wildcards_to_closures(core)

    for row in transform_data:
        if num_cuts == 1:
            # Trivial because the core and to_smiles can only be glued together in one orientation

            # Important to clear atom maps due to smiles_syntax, see above where we clear atom maps for core
            variable_to_weld = Chem.MolToSmiles(clear_atom_maps(Chem.MolFromSmiles(row['to_smiles'])))
            variable_to_weld = smiles_syntax.convert_wildcards_to_closures(variable_to_weld)
            welded = Chem.MolFromSmiles(f"{core_to_weld}.{variable_to_weld}")
            products = [welded]
        else:
            # Multiple products are possible depending on symmetry and environment
            products = []
            # Multiple identical products will be generated if the from_smiles and to_smiles share the same symmetry, for example [He]c1ccc([Ne])cc1 >> [Ne]c1cc([He])ccc1
            unique_product_smiles = set()
            # For multicut queries, to make sure stuff is glued in the right orientation, we need to flip around the attachment points, considering all possibilities of how to match variable to core
            # First we generate all attachment point permutations of user-defined variable and core, using permute_variables_constant
            # Attachment points in input variable/core are atomic mass=0, isotope in (1,2,3). Below function converts these attachment points to He/Ne/Ar for substructure matching
            variable_perms, is_None, core_perms = permute_variables_constant(variable1=variable, variable2=None, constant=core, has_labeled_wildcards=True)

            for variable_perm, core_perm in zip(variable_perms, core_perms):
                if aggregation_type == 'individual_transforms':
                    # We need to know that the orientation of the from_smiles is proper before we glue on the corresponding to_smiles (otherwise we could get the wrong isomer)
                    from_mol = mapped_wildcards_to_nobles(Chem.MolFromSmiles(row['from_smiles']))
                    # We use a substruct match in case the user put wildcards in their input. In any case, all we really care about here are the attachment points matching
                    if not from_mol.HasSubstructMatch(variable_perm):
                        continue
                elif aggregation_type == 'group_by_fragment':
                    to_mol = mapped_wildcards_to_nobles(Chem.MolFromSmiles(row['to_smiles']))
                    # The user-defined variable must be a substructure of the to_smiles
                    # We convert from labeled wildcard smiles, e.g. having [*:1], to mols with noble gas dummy atoms to match variable_perms above, for substructure matching
                    if not to_mol.HasSubstructMatch(variable_perm):
                        continue

                # The environment_smarts (if exists) must be a substructure of the core
                if row.get('environment_smarts') is not None:
                    environment = Chem.MolFromSmarts(row['environment_smarts'])
                    if not core_perm.HasSubstructMatch(environment):
                        continue

                #core_perm_smiles = Chem.MolToSmiles(nobles_to_mapped_wildcards(core_perm))

                core_perm_smiles = get_multifragment_gluable_smiles(nobles_to_mapped_wildcards(core_perm))
                product = Chem.MolFromSmiles(f"{apply_glue(row['to_smiles'])}.{apply_glue(core_perm_smiles)}")
                #product = clear_atom_maps(product)
                product_smiles = Chem.MolToSmiles(product)
                if product_smiles not in unique_product_smiles:
                    unique_product_smiles.add(product_smiles)
                    products.append(product)

        if len(products) == 0:
            # Should we put some kind of error molecule output here? How do we communicate to the user that gluing fragments failed in a case?
            # continue
            products = ['error']

        enumerations[row['transform_id']] = {'products': products}

    return enumerations
