from pydantic import BaseModel, Field, conlist
from typing import Optional, Callable, List, Literal, Union, Dict
import copy

sketched_content_examples = {
                "mol1_molfile": '\n  Ketcher  8112214332D 1   1.00000     0.00000     0\n\n  7  7  0  0  0  0            999 V2000\n    6.1042   -4.2292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.9702   -4.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.9702   -5.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.1042   -6.2292    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    5.2382   -5.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.2382   -4.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.1042   -3.2292    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0     0  0\n  2  3  2  0     0  0\n  3  4  1  0     0  0\n  4  5  2  0     0  0\n  5  6  1  0     0  0\n  6  1  2  0     0  0\n  1  7  1  0     0  0\nM  END\n',
                "mol1_variable_atoms": "6",
                "mol1_variable_bonds": "",
                "mol1_environment_atoms": "0",
                "mol1_environment_bonds": "",
                "mol2_molfile": '\n  Ketcher  8112214352D 1   1.00000     0.00000     0\n\n  7  7  0  0  0  0            999 V2000\n    6.4583   -4.4792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3243   -4.9792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3243   -5.9792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4583   -6.4792    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5923   -5.9792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5923   -4.9792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4583   -3.4792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0     0  0\n  2  3  2  0     0  0\n  3  4  1  0     0  0\n  4  5  2  0     0  0\n  5  6  1  0     0  0\n  6  1  2  0     0  0\n  1  7  1  0     0  0\nM  END\n',
                "mol2_variable_atoms": "6",
                "mol2_variable_bonds": "",
                "mol2_environment_atoms": "0",
                "mol2_environment_bonds": "",
}

class SketchedContent(BaseModel):
    mol1_molfile: Optional[str] = ''
    mol1_variable_atoms: Optional[str] = ''
    mol1_variable_bonds: Optional[str] = ''
    mol1_environment_atoms: Optional[str] = ''
    mol1_environment_bonds: Optional[str] = ''
    mol2_molfile: Optional[str] = ''
    mol2_variable_atoms: Optional[str] = ''
    mol2_variable_bonds: Optional[str] = ''
    mol2_environment_atoms: Optional[str] = ''
    mol2_environment_bonds: Optional[str] = ''

no_selection_sketched_content_examples = copy.deepcopy(sketched_content_examples)
no_selection_sketched_content_examples['mol1_variable_atoms'] = ""
no_selection_sketched_content_examples['mol1_environment_atoms'] = ""
no_selection_sketched_content_examples['mol2_variable_atoms'] = ""
no_selection_sketched_content_examples['mol2_environment_atoms'] = ""

class ValidateSelectionInput(BaseModel):
    selection_type: Literal['variable', 'environment']
    sketched_content: SketchedContent
    mol1_selected_atoms: Optional[str] = ''
    mol1_selected_bonds: Optional[str] = ''
    mol2_selected_atoms: Optional[str] = ''
    mol2_selected_bonds: Optional[str] = ''

    class Config:
        schema_extra = {
            "example": {
                "selection_type": "variable",
                "sketched_content": no_selection_sketched_content_examples,
                "mol1_selected_atoms": "6",
                "mol1_selected_bonds": "",
                "mol2_selected_atoms": "6",
                "mol2_selected_bonds": "",
            }
        }

advanced_options_examples = {
                "variable_min_heavies": 0,
                "variable_max_heavies": 0,
                "compound_min_heavies": 0,
                "compound_max_heavies": 0,
                "aggregation_type": "individual_transforms",
}

class AdvancedOptions(BaseModel):
    variable_min_heavies: Optional[int] = 0
    variable_max_heavies: Optional[int] = 0
    compound_min_heavies: Optional[int] = 0
    compound_max_heavies: Optional[int] = 0
    aggregation_type: Optional[Literal['individual_transforms', 'group_by_fragment']] = 'individual_transforms'

    class Config:
        schema_extra = {
            "examples": advanced_options_examples,
        }

class QueryInput(BaseModel):
    snapquery_id: Optional[str] = ''
    query_id: Optional[int] = -1
    query_type: Literal['exact', 'substructure']
    transform_order: Literal['first_order', 'second_order']
    sketched_content: SketchedContent
    REQUIRED_properties: Optional[str] = ''
    OPTIONAL_properties: Optional[str] = ''
    advanced_options: Optional[AdvancedOptions] = AdvancedOptions.parse_obj(advanced_options_examples)
    snapfilter_id: Optional[str] = ''
    snapfilter_string: Optional[str] = ''

    class Config:
        schema_extra = {
            "example": {
                "snapquery_id": "",
                "query_id": -1,
                "query_type": "exact",
                "transform_order": "first_order",
                "sketched_content": sketched_content_examples, 
                "REQUIRED_properties": "",
                "OPTIONAL_properties": "hERG_pIC50",
                "snapfilter_id": "",
                "snapfilter_string": "",
            }
        }


class repQueryInput(BaseModel):
    rule_environment_statistics_id: int
    #rep_property: str = Field(..., min_length=1)

class ExampleQuery(QueryInput):
    # For use when we pre-load example query inputs to the database
    # No query has actually been run, therefore there is no query_id that has been assigned
    # If the rest of the code in backend/frontend works correctly, a query_id should be assigned when the query is first run
    query_id: Optional[Literal['NULL']] = 'NULL'


# Filter MMPs based on the unaggregated data
# Example: A Papp < 15, meaning the MMP's starting compound's Papp is less than 15
class RangeFilter(BaseModel):
    # A refers to starting compound (LHS) in the MMP, B is ending compound (RHS)
    # Example for methanol to ethanol, CO >> CCO: CO is A, CCO is B
    compound: Literal['A', 'B']
    property_name: str = Field(..., min_length=1)
    operator: Literal['<', '>', '<=', '>=', '=', '!=']
    base: Optional[Literal['raw', 'log', 'negative_log']] = 'raw'
    units: Optional[Literal['uM', 'M']] = 'uM'
    value: float = Field(...)


# Desired statistics to calculate
class Statistic(BaseModel):
    # can add more statistics here, if desired
    statistic: Literal['median']
    property_name: str = Field(..., min_length=1)
    change_type: Literal['delta', 'fold_change']
    base: Optional[Literal['raw', 'log']] = 'raw'
    units: Optional[Literal['uM', 'M']] = None


class CreateSnapshotByRule(BaseModel):
    rule_environment_id: int
    query_type: Literal['exact', 'substructure'] = 'exact'
    transform_order: Literal['first_order', 'second_order']
    advanced_options: AdvancedOptions = AdvancedOptions(
        variable_min_heavies=0,
        variable_max_heavies=0,
        compound_max_heavies=0,
        aggregation_type='group_by_fragment'
    )
    REQUIRED_properties: Optional[str] = ''
    OPTIONAL_properties: Optional[str] = ''


class AggregationParameters(BaseModel):
    query_id: int = Field(...)
    aggregation_type: Optional[Literal['individual_transforms', 'group_by_fragment']] = 'individual_transforms'
    # Require the MMPs to have data for these properties
    range_filters: Optional[List[RangeFilter]] = []
    statistics: List[Statistic]

    class Config:
        schema_extra = {
            "example": {
                "query_id": 1,
                "aggregation_type": "individual_transforms",
                "statistics": [
                    {
                        "statistic": "median",
                        "property_name": "hERG_pIC50",
                        "change_type": "fold_change",
                        "base": "raw",
                        "units": "uM",
                    }
                ]
            }
        }


class PlotDataDimension(BaseModel):
    property_name: str = Field(..., min_length=1)
    change_type: Literal['A', 'B', 'average', 'fold_change', 'delta']
    base: Optional[Literal['raw', 'log']] = 'raw'
    units: Optional[Literal['M', 'uM']] = None


class PlotParameters(BaseModel):
    query_id: int = Field(...)
    aggregation_type: Optional[Literal['individual_transforms', 'group_by_fragment']] = 'individual_transforms'
    grouped_by_environment: Optional[bool] = False
    ids: Union[List[int], List[str]]
    range_filters: Optional[List[RangeFilter]] = []
    x_data: PlotDataDimension = Field(...)
    y_data: PlotDataDimension = Field(...)


class PairPropChange(PlotDataDimension):
    change_type: Literal['fold_change', 'delta']


class PairsCondensed(BaseModel):
    prop_changes: List[PairPropChange]
    construct_id_pairs: List[conlist(int, min_items=2, max_items=2)]


class GetAllRawData(BaseModel):
    query_id: int
    query: QueryInput


class EnumerationData(BaseModel):
    query: QueryInput
    query_id: int = Field(...)
    row_data: List[Dict]
    link_address: Optional[str] = ''
