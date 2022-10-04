var backend_root = "{{ external_backend_root }}";
var frontend_root = "{{ external_frontend_root }}";
var dash_path = frontend_root + '/dash/';

// When the page initially loads, the Dash app submit button will be automatically fired, which is necessary for automatic query submission when a snapshot is loaded in via template
// However, we want to use this same html file even without loading a snapshot in via template, in which case the Dash submit button firing would trigger error messages, without the below global variable
var suppress_initial_errors = true;
var snapquery_id = "{{ snapquery_id }}";
var snapfilter_id = "{{ snapfilter_id }}";
var query_id = parseInt("{{ query_id }}");
query_id = (!isNaN(query_id) ? query_id : -1);
var snapfilter_string = "{{ snapfilter_string }}";

var schema = "{{ schema }}";
var ketcher1 = undefined;
var ketcher2 = undefined;

// Highlights are included in the editor's history stack, meaning the user can undo a highlight with Ctrl+z (or Cmd+z)
// This is a problem, because the undo function is outside the scope of our control, and when a user undos highlights, the sketcher highlights will not be synchronized with our highlights variable in this file
// To solve this problem, we use the below two variables to 'understand' when we did (or did not) cause a highlight change with code in this file, vs. the editor undo function
// When a highlight change occurs in ketcher1, but ketcher1_legit_highlight_change is false, then we know that highlight change was caused by the undo button, and we can take appropriate action
var ketcher1_legit_highlight_change = false;
var ketcher2_legit_highlight_change = false;

function createHighlight (ketcher_id, input) {
    if (ketcher_id === 'ketcher1') {
        ketcher1_legit_highlight_change = true;
        ketcher1.editor.highlights.create(input);
        ketcher1_legit_highlight_change = false;
    } else if (ketcher_id === 'ketcher2') {
        ketcher2_legit_highlight_change = true;
        ketcher2.editor.highlights.create(input);
        ketcher2_legit_highlight_change = false;
    }
}

function parseHighlights (input) {
    let output = input.split(',');
    if (JSON.stringify(output) === JSON.stringify([""])) {return [];}
    output = output.map(x => parseInt(x))
    return output;
}

const empty_highlights = {
    variable: {
        atoms: [], bonds: []
    },
    environment: {
        atoms: [], bonds: []
    },
    maps: {
        atoms: {}, bonds: {}
    },
    inverse_maps: {
        atoms: {}, bonds: {}
    }
};

var highlights = {
    mol1: JSON.parse(JSON.stringify(empty_highlights)),
    mol2: JSON.parse(JSON.stringify(empty_highlights))
};

function startup() {
    get_prop_options();
    bind_ketcher_editors();
    bind_event_listeners();
    initialize_content();
}

function initialize_content() {

    if (snapquery_id === "") {
        // null is the default value, but might be annoying to deal with across JS/python/SQL, in all use cases, therefore set '' as default value
        document.getElementById("ketcher1_content").setAttribute("value", '');
        document.getElementById("ketcher2_content").setAttribute("value", '');
        document.getElementById("outputPlot").setAttribute("src", dash_path);
        return;
    }

    let optPropBoxes = document.getElementsByClassName("optPropBox");

    // Only populate the sketchers after they have loaded, same with property boxes
    if (ketcher1 !== undefined && ketcher2 !== undefined && optPropBoxes.length > 0) {
        // The goal with this template is to execute a query with the below input state, upon page load
        // Initialize radio buttons
        var initial_query_type = "{{ query_type }}";
        document.getElementById(initial_query_type).click();
        var initial_transform_order = "{{ transform_order }}";
        document.getElementById(initial_transform_order).click();

        if (`{{ mol2_molfile }}` !== "") {
            document.getElementById("starting_molecule_and_product").click();
        }

        // Current version of Dash does not support promises in clientside callbacks, therefore we need the text ready to go in an element, for Dash to read from
        document.getElementById("ketcher1_content").setAttribute("value", `{{ mol1_molfile }}`);
        document.getElementById("ketcher2_content").setAttribute("value", `{{ mol2_molfile }}`);

        // We aren't preventing an infinite loop here, but we take advantage of this variable to prevent another sketcher export/import round,
        // which would trigger a molchange event, which would erase the highlights
        ketcher1_molchange_infinite_loop_prevention = 'on';
        ketcher2_molchange_infinite_loop_prevention = 'on';

        ketcher1.setMolecule(`{{ mol1_molfile }}`)
            .then(function () {
                return ketcher2.setMolecule(`{{ mol2_molfile }}`);
            })
            .then(function () {

                generate_index_maps();
                ketcher1.editor.selection(null);
                ketcher2.editor.selection(null);

                highlights.mol1.variable.atoms = apply_index_map(highlights.mol1.inverse_maps.atoms, parseHighlights("{{ mol1_variable_atoms }}"));
                highlights.mol1.variable.bonds = apply_index_map(highlights.mol1.inverse_maps.bonds, parseHighlights("{{ mol1_variable_bonds }}"));
                highlights.mol1.environment.atoms = apply_index_map(highlights.mol1.inverse_maps.atoms, parseHighlights("{{ mol1_environment_atoms }}"));
                highlights.mol1.environment.bonds = apply_index_map(highlights.mol1.inverse_maps.bonds, parseHighlights("{{ mol1_environment_bonds }}"));
                highlights.mol2.variable.atoms = apply_index_map(highlights.mol2.inverse_maps.atoms, parseHighlights("{{ mol2_variable_atoms }}"));
                highlights.mol2.variable.bonds = apply_index_map(highlights.mol2.inverse_maps.bonds, parseHighlights("{{ mol2_variable_bonds }}"));
                highlights.mol2.environment.atoms = apply_index_map(highlights.mol2.inverse_maps.atoms, parseHighlights("{{ mol2_environment_atoms }}"));
                highlights.mol2.environment.bonds = apply_index_map(highlights.mol2.inverse_maps.bonds, parseHighlights("{{ mol2_environment_bonds }}"));

                createHighlight('ketcher1', {atoms: highlights.mol1.variable.atoms, bonds: highlights.mol1.variable.bonds, color: '#0FD730'});
                createHighlight('ketcher1', {atoms: highlights.mol1.environment.atoms, bonds: highlights.mol1.environment.bonds, color: '#fa6e6e'});
                createHighlight('ketcher2', {atoms: highlights.mol2.variable.atoms, bonds: highlights.mol2.variable.bonds, color: '#0FD730'});
                createHighlight('ketcher2', {atoms: highlights.mol2.environment.atoms, bonds: highlights.mol2.environment.bonds, color: '#fa6e6e'});

                // Initialize property selections
                if ("{{ REQUIRED_properties }}" !== "") {
                    setProps("req", "{{ REQUIRED_properties }}");
                }

                if ("{{ OPTIONAL_properties }}" !== "") {
                    setProps("opt", "{{ OPTIONAL_properties }}");
                }

                let imported_advanced_fields = {
                    variable_min_heavies : (("{{ variable_min_heavies }}" !== "0") ? "{{ variable_min_heavies }}" : ""),
                    variable_max_heavies : (("{{ variable_max_heavies }}" !== "0") ? "{{ variable_max_heavies }}" : ""),
                    compound_min_heavies : (("{{ compound_min_heavies }}" !== "0") ? "{{ compound_min_heavies }}" : ""),
                    compound_max_heavies : (("{{ compound_max_heavies }}" !== "0") ? "{{ compound_max_heavies }}" : ""),
                    aggregation_type : "{{ aggregation_type }}",
                };

                setAdvancedOptions(imported_advanced_fields);

                // Only load the Dash app after the entire snapshot content is initialized as above, because Dash app fires query on load
                document.getElementById("outputPlot").setAttribute("src", dash_path);
            })
    } else {
        window.setTimeout(initialize_content, 100);
    }
}

var mol1_entire_molecule_selected = false;
var mol2_entire_molecule_selected = false;

function clearVariableAction () {
    highlights.mol1.variable.atoms = [];
    highlights.mol1.variable.bonds = [];
    clearHighlights('ketcher1');
    highlights.mol2.variable.atoms = [];
    highlights.mol2.variable.bonds = [];
    clearHighlights('ketcher2');
    mol1_entire_molecule_selected = false;
    mol2_entire_molecule_selected = false;

    createHighlight('ketcher1', {atoms: highlights.mol1.environment.atoms, bonds: highlights.mol1.environment.bonds, color: '#fa6e6e'});
    createHighlight('ketcher2', {atoms: highlights.mol2.environment.atoms, bonds: highlights.mol2.environment.bonds, color: '#fa6e6e'});
}

function clearEnvironmentAction () {
    highlights.mol1.environment.atoms = [];
    highlights.mol1.environment.bonds = [];
    clearHighlights('ketcher1');
    highlights.mol2.environment.atoms = [];
    highlights.mol2.environment.bonds = [];
    clearHighlights('ketcher2');

    createHighlight('ketcher1', {atoms: highlights.mol1.variable.atoms, bonds: highlights.mol1.variable.bonds, color: '#0FD730'});
    createHighlight('ketcher2', {atoms: highlights.mol2.variable.atoms, bonds: highlights.mol2.variable.bonds, color: '#0FD730'});
}

var delayTime = 5;

function delayCall (func, args) {
    setTimeout(function () {
        func.apply(null, args);
    }, delayTime);
};

function handleSetVariableBtnClick () {
    delayCall(validateSelection, ['variable']);
};

function handleSetEnvironmentBtnClick  () {
    delayCall(validateSelection, ['environment']);
};

function handleClearVariableBtnClick () {
    delayCall(clearVariableAction, []);
};

function handleClearEnvironmentBtnClick () {
    delayCall(clearEnvironmentAction, []);
};

async function postRequest(url = '', data = {}) {

    const response = await fetch(url, {
        method: 'POST',
        mode: 'cors',
        cache: 'default',
        credentials: 'same-origin',
        headers: {
            'Content-Type': 'application/json'
        },
        redirect: 'follow',
        referrerPolicy: 'no-referrer',
        body: JSON.stringify(data)
    });
    return response.json();
}

async function getRequest(url = '') {

    const response = await fetch(url, {
        method: 'GET',
        mode: 'cors',
        cache: 'default',
        credentials: 'same-origin',
        headers: {
            'Content-Type': 'application/json'
        },
        redirect: 'follow',
        referrerPolicy: 'no-referrer',
    });
    return response.json();
}

function generate_index_maps () {
    ketcher1.editor.selection('all');
    ketcher1_selection = ketcher1.editor.selection();
    let mol1_atom_indexes = ketcher1_selection?.atoms;
    if (mol1_atom_indexes !== undefined) {
        mol1_atom_indexes.sort(function(x,y) {return x-y;});
        //let mol1_normalized_indexes = Array.from(Array(mol1_atom_indexes.length).keys());
        for (let i = 0; i < mol1_atom_indexes.length; i++) {
            highlights.mol1.maps.atoms[mol1_atom_indexes[i]] = i;
            highlights.mol1.inverse_maps.atoms[i] = mol1_atom_indexes[i];
        }
        let mol1_bond_indexes = ketcher1_selection?.bonds;
        mol1_bond_indexes.sort(function(x,y) {return x-y;});
        for (let i = 0; i < mol1_bond_indexes.length; i++) {
            highlights.mol1.maps.bonds[mol1_bond_indexes[i]] = i;
            highlights.mol1.inverse_maps.bonds[i] = mol1_bond_indexes[i];
        }
    }

    ketcher2.editor.selection('all');
    ketcher2_selection = ketcher2.editor.selection();
    let mol2_atom_indexes = ketcher2_selection?.atoms;
    if (mol2_atom_indexes !== undefined) {
        mol2_atom_indexes.sort(function(x,y) {return x-y;});
        //let mol2_normalized_indexes = Array.from(Array(mol2_atom_indexes.length).keys());
        for (let i = 0; i < mol2_atom_indexes.length; i++) {
            highlights.mol2.maps.atoms[mol2_atom_indexes[i]] = i;
            highlights.mol2.inverse_maps.atoms[i] = mol2_atom_indexes[i];
        }
        let mol2_bond_indexes = ketcher2_selection?.bonds;
        mol2_bond_indexes.sort(function(x,y) {return x-y;});
        for (let i = 0; i < mol2_bond_indexes.length; i++) {
            highlights.mol2.maps.bonds[mol2_bond_indexes[i]] = i;
            highlights.mol2.inverse_maps.bonds[i] = mol2_bond_indexes[i];
        }
    }
}

function apply_index_map(map_object, array) {
    // don't modify in place
    let newArray = JSON.parse(JSON.stringify(array));
    for (let i=0; i<newArray.length; i++) {
        newArray[i] = map_object[newArray[i]];
    }
    return newArray;
}

function getSketchedContent() {
    let ketcher1_content = document.getElementById('ketcher1_content').getAttribute('value');
    let ketcher2_content = document.getElementById('ketcher2_content').getAttribute('value');
    let sketched_content = {
        mol1_molfile: ketcher1_content,
        // Here we translate atom indices that ketcher generates, to atom indices that span from 0 ... n, where n is number of atoms or bonds
        mol1_variable_atoms: apply_index_map(highlights.mol1.maps.atoms, highlights.mol1.variable.atoms).join(","),
        mol1_variable_bonds: apply_index_map(highlights.mol1.maps.bonds, highlights.mol1.variable.bonds).join(","),
        mol1_environment_atoms: apply_index_map(highlights.mol1.maps.atoms, highlights.mol1.environment.atoms).join(","),
        mol1_environment_bonds: apply_index_map(highlights.mol1.maps.bonds, highlights.mol1.environment.bonds).join(","),
        mol2_molfile: ketcher2_content,
        mol2_variable_atoms: apply_index_map(highlights.mol2.maps.atoms, highlights.mol2.variable.atoms).join(","),
        mol2_variable_bonds: apply_index_map(highlights.mol2.maps.bonds, highlights.mol2.variable.bonds).join(","),
        mol2_environment_atoms: apply_index_map(highlights.mol2.maps.atoms, highlights.mol2.environment.atoms).join(","),
        mol2_environment_bonds: apply_index_map(highlights.mol2.maps.bonds, highlights.mol2.environment.bonds).join(","),
    };

    return sketched_content;
}

function getSelectionData(selection_type) {
    let ketcher1_selection = ketcher1.editor.selection();
    let mol1_selected_atoms = ketcher1_selection?.atoms;
    let mol1_selected_bonds = ketcher1_selection?.bonds;

    let ketcher2_selection = ketcher2.editor.selection();
    let mol2_selected_atoms = ketcher2_selection?.atoms;
    let mol2_selected_bonds = ketcher2_selection?.bonds;

    // ketcher does not always assign atom/bond indexes starting from 0 and going to n, where n is the number of atoms or bonds
    // instead, ketcher sometimes assigns seemingly arbitrary index sequences such as 41, 42, 43 (instead of assigning 0,1,2) to the atoms in a molecule with 3 atoms
    // Here we assign maps, so that for example in the above case, we can later map from 41 -> 0, 42 -> 1, 43 -> 2 (highlights.mol1.maps) and in the inverse direction (highlights.mol1.inverse_maps)
    generate_index_maps();

    let selection_data = {
        selection_type: selection_type,
        mol1_selected_atoms: (mol1_selected_atoms !== undefined) ? apply_index_map(highlights.mol1.maps.atoms, mol1_selected_atoms).join(",") : "",
        mol1_selected_bonds: (mol1_selected_bonds !== undefined) ? apply_index_map(highlights.mol1.maps.bonds, mol1_selected_bonds).join(",") : "",
        mol2_selected_atoms: (mol2_selected_atoms !== undefined) ? apply_index_map(highlights.mol2.maps.atoms, mol2_selected_atoms).join(",") : "",
        mol2_selected_bonds: (mol2_selected_bonds !== undefined) ? apply_index_map(highlights.mol2.maps.bonds, mol2_selected_bonds).join(",") : "",
    };

    return selection_data;
}

// 1. Check if variable selection will result in a legal query fragment - if not, automatically select the closest legal fragment
//    Pass both atom and bond indices. Return validated atom and bond indices.
// 2. Check if variable selection overlaps with environment selection - if so, delete environment selection

function validateSelection (selection_type) {

    if (selection_type !== 'variable' && selection_type !== 'environment') {
        throw new Error("argument passed to validateSelection must be either 'variable' or 'environment' to specify the selection type");
    }

    return Promise.resolve()
        .then(function () {

            let selection_data = getSelectionData(selection_type);
            ketcher1.editor.selection(null);
            ketcher2.editor.selection(null);
            if (!selection_data.mol1_selected_atoms && !selection_data.mol1_selected_bonds && !selection_data.mol2_selected_atoms && !selection_data.mol2_selected_bonds) {
                return Promise.resolve();
            }

            // Include current highlights below. Don't remove current highlights if they are on the other side of rxn arrow from selection
            let sketched_content = getSketchedContent();
            selection_data['sketched_content'] = sketched_content;

            return postRequest(
                backend_root + '/validateSelection/',
                selection_data
            );
        })
        .then(function (data) {

            if (data === undefined) {return;}

            highlights.mol1.variable.atoms = apply_index_map(highlights.mol1.inverse_maps.atoms, data["mol1_variable_atoms"]);
            highlights.mol1.variable.bonds = apply_index_map(highlights.mol1.inverse_maps.bonds, data["mol1_variable_bonds"]);
            highlights.mol1.environment.atoms = apply_index_map(highlights.mol1.inverse_maps.atoms, data["mol1_environment_atoms"]);
            highlights.mol1.environment.bonds = apply_index_map(highlights.mol1.inverse_maps.bonds, data["mol1_environment_bonds"]);

            highlights.mol2.variable.atoms = apply_index_map(highlights.mol2.inverse_maps.atoms, data["mol2_variable_atoms"]);
            highlights.mol2.variable.bonds = apply_index_map(highlights.mol2.inverse_maps.bonds, data["mol2_variable_bonds"]);
            highlights.mol2.environment.atoms = apply_index_map(highlights.mol2.inverse_maps.atoms, data["mol2_environment_atoms"]);
            highlights.mol2.environment.bonds = apply_index_map(highlights.mol2.inverse_maps.bonds, data["mol2_environment_bonds"]);

            if (data["mol1_entire_molecule_selected"] === "True") {
                mol1_entire_molecule_selected = true;
            }
            if (data["mol1_entire_molecule_selected"] === "False") {
                mol1_entire_molecule_selected = false;
            }
            if (data["mol2_entire_molecule_selected"] === "True") {
                mol2_entire_molecule_selected = true;
            }
            if (data["mol2_entire_molecule_selected"] === "False") {
                mol2_entire_molecule_selected = false;
            }

            clearHighlights('ketcher1');
            createHighlight('ketcher1', {atoms: highlights.mol1.variable.atoms, bonds: highlights.mol1.variable.bonds, color: '#0FD730'});
            createHighlight('ketcher1', {atoms: highlights.mol1.environment.atoms, bonds: highlights.mol1.environment.bonds, color: '#fa6e6e'});

            clearHighlights('ketcher2');
            createHighlight('ketcher2', {atoms: highlights.mol2.variable.atoms, bonds: highlights.mol2.variable.bonds, color: '#0FD730'});
            createHighlight('ketcher2', {atoms: highlights.mol2.environment.atoms, bonds: highlights.mol2.environment.bonds, color: '#fa6e6e'});
        })
}

var ketcher1_molchange_infinite_loop_prevention = 'off';
var ketcher2_molchange_infinite_loop_prevention = 'off';

function handleMolChange (ketcher_id) {
    clearHighlights(ketcher_id, true);
    let sketcher = (ketcher_id === 'ketcher2') ? ketcher2 : ketcher1;

    // Ensure validateRGroupsAction finishes before exiting this function
    return validateRGroupsAction(sketcher, ketcher_id)
        .then(function () {return;})
}

function clearHighlights (ketcher_id, clear_stored=false) {
    if (ketcher_id === 'ketcher1') {
        if (clear_stored===true) {highlights.mol1 = JSON.parse(JSON.stringify(empty_highlights));}
        ketcher1_legit_highlight_change = true;
        ketcher1.editor.highlights.clear();
        ketcher1_legit_highlight_change = false;
        mol1_entire_molecule_selected = false;
    } else if (ketcher_id === 'ketcher2') {
        if (clear_stored===true) {highlights.mol2 = JSON.parse(JSON.stringify(empty_highlights));}
        ketcher2_legit_highlight_change = true;
        ketcher2.editor.highlights.clear();
        ketcher2_legit_highlight_change = false;
        mol2_entire_molecule_selected = false;
    }
}

const num_patt = '([0-9][0-9][0-9]| [0-9][0-9]|  [0-9])';
const atom_line_patt = '([A-Z]|\\*)([a-z]| )' + num_patt.repeat(10);
//const bond_line_patt = '\\n' + num_patt.repeat(7) + '\\n';
const atom_line_regexp = new RegExp(atom_line_patt, 'g');
const null_molfile_signature = /  0  0  0  0  0  0            999 V2000/;

function storeMol(molfile, ketcher_id) {
    // Don't store null molecule
    let null_if_molecule = null_molfile_signature.exec(molfile);
    if (null_if_molecule !== null) {
        molfile = '';
    }
    document.getElementById(ketcher_id + '_content').setAttribute('value', molfile);
}

function validateRGroupsAction(sketcher, ketcher_id) {

    return sketcher.getMolfile()
        .then(function(sketched_molfile) {
            let M_conditions_start = /M  [A-Z]/.exec(sketched_molfile);
            let first_M_condition = sketched_molfile.slice(M_conditions_start.index + 3, M_conditions_start.index + 6);
            if (first_M_condition === 'END') {
                // This means there is nothing 'potentially wrong' with our drawn structure
                // Don't correct and reimport the structure, because if we do, the user can't undo changes using undo functions like Ctrl + z
                //document.getElementById(ketcher_id + '_content').setAttribute('value', sketched_molfile);
                storeMol(sketched_molfile, ketcher_id);
                return Promise.resolve();
            }
            sketched_molfile = sketched_molfile.slice(0, M_conditions_start.index) + "M  END\n";

            // We will get an infinite loop without using such a control measure
            if (ketcher_id === 'ketcher2') {
                ketcher2_molchange_infinite_loop_prevention = 'on';
            }
            else {ketcher1_molchange_infinite_loop_prevention = 'on';}

            //document.getElementById(ketcher_id + '_content').setAttribute('value', sketched_molfile);
            storeMol(sketched_molfile, ketcher_id);
            return sketcher.setMolecule(sketched_molfile);
        })
}

function attachMaps(rxnfile) {
    // Attach maps from 1 to n (where n is the total number of atoms) to the atoms, in order, across the two atom blocks in the rxnfile

    let beginning;
    let map_num_padding;
    let map_num = 0;
    let end;
    const matches = rxnfile.matchAll(atom_line_regexp);
    for (const match of matches) {

        beginning = rxnfile.slice(0, match.index + 29);
        map_num = map_num + 1;
        map_num_padding = ' '.repeat(3 - String(map_num).length);
        end = rxnfile.slice(match.index + 32);
        rxnfile = beginning + map_num_padding + String(map_num) + end;
    }

    return rxnfile;
}

function getProps (req_or_opt) {

    var selectedBoxes = document.getElementsByClassName(req_or_opt + "PropBox");
    let props = '';
    var j = 0;
    for (var i = 0; i < selectedBoxes.length; i++) {
        if (selectedBoxes[i].checked) {
            if (j > 0) {
                props += ',';
            }
            props += selectedBoxes[i].getAttribute("name");
            j += 1;
        }
    }
    return props
}

function setProps (req_or_opt, props_string) {

    var prop_boxes = document.getElementsByClassName(req_or_opt + "PropBox");
    props_array = props_string.split(",");

    // Create lookup table of property names vs. their index
    var prop_boxes_dict = {};
    for (var i = 0; i < prop_boxes.length; i++) {
        prop_boxes_dict[prop_boxes[i].getAttribute("name")] = i;
    }
    // Check property boxes that match the props_string inputs
    for (var j = 0; j < props_array.length; j++) {
        box_index = prop_boxes_dict[props_array[j]];
        prop_boxes[box_index].checked = true;
    }
}

function getAdvancedOptions () {
    let advanced_fields_boxes = document.getElementsByClassName("advancedFields");
    let advanced_fields = {};
    for (let i=0; i<advanced_fields_boxes.length; i++) {
        let field_name = advanced_fields_boxes[i].getAttribute("name");
        let field_number = parseInt(advanced_fields_boxes[i].value);
        if (!isNaN(field_number)) {
            advanced_fields[field_name] = field_number;
        }
    }
    let aggregation_type_radios = parent.document.getElementsByName("aggregationType");
    let aggregation_type;
    for (let i=0; i < aggregation_type_radios.length; i++) {
        if (aggregation_type_radios[i].checked) {
            aggregation_type = aggregation_type_radios[i].value;
        }
    }
    advanced_fields["aggregation_type"] = aggregation_type;
    return advanced_fields;
}

function setAdvancedOptions (advanced_fields) {
    let advanced_fields_boxes = document.getElementsByClassName("advancedFields");
    for (let i=0; i<advanced_fields_boxes.length; i++) {
        let field_name = advanced_fields_boxes[i].getAttribute("name");
        let field_number = advanced_fields[field_name]
        if (field_number !== undefined) {
            advanced_fields_boxes[i].setAttribute("value", String(field_number));
        }
    }

    document.getElementById(advanced_fields["aggregation_type"]).click();
    return;
}

function submitAction () {
    var sketchedPromise = inputController.sketcherInstance.exportStructure('mol', null);
    sketchedPromise.then(function(source) {
        inputController.submit(source);
    }, function(error) {
        alert(error);
    });
}

function push_original_url () {
    history.pushState("", document.title, window.location.pathname);
}

function jump_to_results () {
    // A delay is necessary here (when called from a Dash callback after the results are ready), otherwise the jump won't work
    setTimeout(function () {
        location.href = "#";
        location.href = "#outputPlot";
        push_original_url();
    }, 2000);
}

let root = document.documentElement;

function resize_iframe(iframe, desired_height="15000px") {

    setTimeout(function () {
        //iframe.style.height = iframe.contentWindow.document.documentElement.scrollHeight + "px";
        iframe.style.height = desired_height;
    }, 100);
}

function toggle_advanced_options_visibility() {
    let visibility = root.style.getPropertyValue('--advanced_options_div_visibility');
    if (visibility === '') {
        root.style.setProperty('--advanced_options_div_visibility', 'inline');
    } else {
        root.style.setProperty('--advanced_options_div_visibility', '');
    }
    return;
}

function bind_ketcher1() {
    ketcher1 = document.getElementById("ketcher_iframe1").contentWindow.ketcher;
    if (ketcher1 !== undefined) {
        ketcher1.editor.zoom(1.2);
        ketcher1.editor.subscribe('change', function (data) {
            // One possible outcome of this function is an automated import to the editor, which would trigger an infinite loop unless we use the below flag
            if (ketcher1_molchange_infinite_loop_prevention === 'on') {
                ketcher1_molchange_infinite_loop_prevention = 'off';
                return;
            }
            let execute_handle_mol_change = true;
            let op;
            let execute_undo = false;
            for (let i = 0; i < data.length; i++) {
                op = data[i]['operation'];
                // Don't count highlight operations as a change, because we remove highlights when change occurs
                if (op === 'Highlight' || op === 'Remove highlight' || op === 'Move atom' || op==='Move enhanced flag') {
                    execute_handle_mol_change = false;
                }
                if (op === 'Highlight' || op === 'Remove highlight') {
                    // See the definition of ketcher1_legit_highlight_change to understand why this is necessary
                    if (ketcher1_legit_highlight_change === false) {
                        execute_undo = true;
                    }
                }
            }
            if (execute_handle_mol_change === true) {handleMolChange('ketcher1');}
            if (execute_undo === true) {ketcher1.editor.undo();}
        });
        return;
    }
    window.setTimeout(bind_ketcher1, 100);
}

function bind_ketcher2() {
    ketcher2 = document.getElementById("ketcher_iframe2").contentWindow.ketcher;
    if (ketcher2 !== undefined) {
        ketcher2.editor.zoom(1.2);
        ketcher2.editor.subscribe('change', function (data) {
            // One possible outcome of this function is an automated import to the editor, which would trigger an infinite loop unless we use the below flag
            if (ketcher2_molchange_infinite_loop_prevention === 'on') {
                ketcher2_molchange_infinite_loop_prevention = 'off';
                return;
            }
            let execute_handle_mol_change = true;
            let op;
            let execute_undo = false;
            for (let i = 0; i < data.length; i++) {
                op = data[i]['operation'];
                // Don't count highlight operations as a change, because we remove highlights when change occurs
                if (op === 'Highlight' || op === 'Remove highlight' || op === 'Move atom' || op==='Move enhanced flag') {
                    execute_handle_mol_change = false;
                }
                if (op === 'Highlight' || op === 'Remove highlight') {
                    // See the definition of ketcher2_legit_highlight_change to understand why this is necessary
                    if (ketcher2_legit_highlight_change === false) {
                        execute_undo = true;
                    }
                }
            }
            if (execute_handle_mol_change === true) {handleMolChange('ketcher2');}
            if (execute_undo === true) {ketcher2.editor.undo();}
        });
        return;
    }
    window.setTimeout(bind_ketcher2, 100);
}

var property_metadata;

function get_prop_options() {

    //return getRequest(backend_root + '/propertyNames/')
    let query_schema = "";
    if (schema !== 'None') {
        query_schema = "?schema=" + schema;
    }
    return getRequest(backend_root + '/propertyMetadata' + query_schema)
    .then(function(data) {

        property_metadata = data;

        //let props = data['properties'];
        if (property_metadata.length !== 0) {
            let properties_dropdown = document.getElementById("optDropdown");
            // Remove the message saying that there are no properties to select
            properties_dropdown.removeChild(properties_dropdown.lastElementChild);
            // Populate the properties dropdown with all the props from the database
            for (const [key, value] of Object.entries(property_metadata)) {
                let checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.classList.add('optPropBox');
                checkbox.name = key;
                let checkbox_label = document.createElement('label');
                checkbox_label.textContent = ' ' + value['display_name'] + ' ';
                properties_dropdown.appendChild(checkbox);
                properties_dropdown.appendChild(checkbox_label);
                properties_dropdown.appendChild(document.createElement('br'));
            }
        }
    })
}

function bind_ketcher_editors() {
    bind_ketcher1();
    bind_ketcher2();
}

function bind_event_listeners() {
    document.getElementById("setVariableButton").addEventListener("click", handleSetVariableBtnClick);
    document.getElementById("setEnvironmentButton").addEventListener("click", handleSetEnvironmentBtnClick);
    document.getElementById("clearVariableButton").addEventListener("click", handleClearVariableBtnClick);
    document.getElementById("clearEnvironmentButton").addEventListener("click", handleClearEnvironmentBtnClick);
}

function starting_molecule_only_click () {
    root.style.setProperty('--ketcher_iframe2_div_visibility', 'none');
}

function starting_molecule_and_product_click() {
    root.style.setProperty('--ketcher_iframe2_div_visibility', 'inline-block');
}

/* When the user clicks on the button,
toggle between hiding and showing the dropdown content */
function showReqChoices() {
    document.getElementById("reqDropdown").classList.toggle("show");
}
function showOptChoices() {
    document.getElementById("optDropdown").classList.toggle("show");
}

// Close the dropdown if the user clicks outside of it
window.onclick = function(event) {
    if (!event.target.matches('.dropbtn') && !event.target.matches('.reqPropBox') && !event.target.matches('.optPropBox')) {
        var dropdowns = document.getElementsByClassName("dropdown-content");
        var i;
        for (i = 0; i < dropdowns.length; i++) {
            var openDropdown = dropdowns[i];
            if (openDropdown.classList.contains('show')) {
                openDropdown.classList.remove('show');
            }
        }
    }
}
