"""Test validate_selection function."""

from backend.ss_select import validate_selection
from backend.models import SketchedContent

def test_validate_selection():

    mol1_molfile = """
  Ketcher 103023 8282D 1   1.00000     0.00000     0

  2  1  0  0  0  0            999 V2000
    4.9167   -5.8333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7827   -5.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0     0  0
M  END
"""
    selection_type = 'variable'
    sketched_content = SketchedContent.parse_obj({'mol1_molfile': mol1_molfile})
    mol1_selected_atoms = '1'
    mol1_selected_bonds = ''
    mol2_selected_atoms = ''
    mol2_selected_bonds = ''

    result = validate_selection(selection_type, sketched_content, mol1_selected_atoms, 
        mol1_selected_bonds, mol2_selected_atoms, mol2_selected_bonds)

    assert result == ([1], [], [0], [], 'False', [], [], [], [], 'False')
