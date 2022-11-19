from pymol import cmd, stored
from time import time
import os
import re
from tkinter import *
from tkinter.ttk import Progressbar
from tkinter import filedialog
from tkinter.messagebox import showinfo



# Extract metadata for each DALI result
def extract_stats(dali_results_path, struc_range=[0]):
    dali_stats = {}
    with open(dali_results_path, 'rt') as f:
        data = f.readlines()
    for line in data:
        if line.__contains__('MOLECULE'):
            model = line.split()[1]
            dali_stats[model] = {
                'No':int(line.split()[0].replace(':', '')),
                'PDB':model.split('-')[0],
                'Chain':model.split('-')[1],
                'Z':float(line.split()[2]),
                'rmsd':float(line.split()[3]),
                'lali':int(line.split()[4]),
                'nres':int(line.split()[5]),
                '%%id':int(line.split()[6]),
                'Description':' '.join(line.split()[8:])
            }
    all_keys = list(dali_stats.keys())
    new_keys = []
    if struc_range != [0]:
        for struc_num in struc_range:
            new_keys.append(all_keys[struc_num-1])
        dali_stats = {key: dali_stats[key] for key in new_keys}
    return dali_stats


# Extract translation/rotation matrices for each DALI result
def extract_ttt(dali_results_path, struc_range=[0]):
    ttt_matrices = {}
    with open(dali_results_path, 'rt') as f:
        data = f.readlines()
    for line in data:
        
        if line.__contains__('U(1,.)'):
            model = line.split()[2]
            ttt_matrices[model] = [0,0,0,1]
            for position in [4,5,6,7]:
                ttt_value = line.split()[position]
                if '\"' in ttt_value:
                    ttt_value = ttt_value[:-1]
                ttt_value = float(ttt_value)
                ttt_matrices[model][position-4] = ttt_value
        
        if line.__contains__('U(2,.)'):
            model = line.split()[2]
            ttt_matrices[model] = ttt_matrices[model] + [0.0,0.0,0.0,1.0]
            for position in [4,5,6,7]:
                ttt_value = line.split()[position]
                if '\"' in ttt_value:
                    ttt_value = ttt_value[:-1]
                ttt_value = float(ttt_value)
                ttt_matrices[model][position] = ttt_value
        
        if line.__contains__('U(3,.)'):
            model = line.split()[2]
            ttt_matrices[model] = ttt_matrices[model] + [0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]
            for position in [4,5,6,7]:
                ttt_value = line.split()[position]
                if '\"' in ttt_value:
                    ttt_value = ttt_value[:-1]
                ttt_value = float(ttt_value)
                ttt_matrices[model][position+4] = ttt_value
    print('Extracted transformation matrices!')
    all_keys = list(ttt_matrices.keys())
    new_keys = []
    if struc_range != [0]:
        for struc_num in struc_range:
            new_keys.append(all_keys[struc_num-1])
        ttt_matrices = {key: ttt_matrices[key] for key in new_keys}
    return ttt_matrices


# Extract residue alignment ifo for each result
def extract_res_alignment(dali_results_path, struc_range=[0]):
    resi_mapping = {} #{ model : {query_res : model_res} }
    with open(dali_results_path, 'rt') as f:
        data = f.readlines()
    for line in data:
        if line.__contains__('<=>'):
            try:
                model = line.split()[2]
                query_start_res = int(line.split('(')[1].split()[1])
                query_end_res = int(line.split('(')[1].split()[4])
                query_range = range(query_start_res, query_end_res+1)
                pdb_start_res = int(re.sub(r'[^0-9]','', line.split('(')[1].split()[7]))
                pdb_end_res = int(re.sub(r'[^0-9]','', line.split('(')[1].split()[10]))
                pdb_range = range(pdb_start_res, pdb_end_res+1)
                if model not in resi_mapping:
                    resi_mapping[model] = {}
                for query_resi, pdb_resi in zip(query_range, pdb_range):
                    resi_mapping[model][pdb_resi] = query_resi
            except:
                print('Text processing error:', line)
    all_keys = list(resi_mapping.keys())
    new_keys = []
    if struc_range != [0]:
        for struc_num in struc_range:
            new_keys.append(all_keys[struc_num-1])
        resi_mapping = {key: resi_mapping[key] for key in new_keys}
    return resi_mapping


# Fetch structures from the PDB
# Save under ./pdb_downloads 
def load_pdbs(query_model_path, dali_stats, ttt_matrices, show_align=0):
    cmd.load(query_model_path)
    query_model = os.path.basename(query_model_path).split('.')[0]
    cmd.copy(query_model+'_sticks', query_model)
    if not os.path.exists('./pdb_downloads'):
        os.makedirs('./pdb_downloads')
    os.chdir('./pdb_downloads/')
    i=0
    window.update_idletasks()
    num_models = len(dali_stats)
    pdb_base_list = []
    for model in dali_stats.keys():
        i=i+1
        pdb_base = dali_stats[model]['PDB']
        pdb_chain = dali_stats[model]['Chain']
        output_text = 'Downloading PDB accession '+ pdb_base + ' ('+str(i)+'/'+str(num_models)+')' 
        print(output_text)
        current_task_label.configure(text=output_text)
        progress_bar['value'] = 25*(i/num_models)
        window.update_idletasks()
        cmd.fetch(pdb_base, model+'_full')
        cmd.transform_object(model+'_full', ttt_matrices[model])
        if show_align:
            cmd.distance(model+'_aln', 'none', 'none')
        cmd.create(model, model+'_full and polymer.protein and chain '+pdb_chain)
        if pdb_base not in pdb_base_list:
            pdb_base_list.append(pdb_base)
    for pdb_base in pdb_base_list:
        cmd.group(pdb_base, pdb_base+'*')
        cmd.disable(pdb_base)
    for model in dali_stats.keys():
        i=i+1
        pdb_base = model.split('-')[0]
        pdb_chain = model.split('-')[1]
        cmd.group(pdb_base+'.No'+str(dali_stats[model]['No'])+'-Chain_'+pdb_chain, model+'*')
    pdb_base_list = []
    for model in dali_stats.keys():
        i=i+1
        pdb_base = model.split('-')[0]
        pdb_chain = model.split('-')[1]
        if pdb_base in pdb_base_list:
            cmd.disable(pdb_base+'.No'+str(dali_stats[model]['No'])+'-Chain_'+pdb_chain)
            cmd.delete(model+'_full')
        else:
            pdb_base_list.append(pdb_base)
            cmd.copy(pdb_base+'_full', model+'_full')
            cmd.delete(model+'_full')
            cmd.group(pdb_base, pdb_base+'.* '+ pdb_base+'_full')
    os.chdir('../')
    return

sidechains = '(sc. + n. CA)'
sidechains_pro = '((sc. + n. CA) + (resn PRO and (sc. + n. CA + n. N)))'


# Store residue alignment info as atom properties
def apply_prop(query_model, resi_mapping, show_align = 0):
    print('Identifying equvalent residues')
    t_init = time()
    num_models = len(resi_mapping)
    i=0
    for model in resi_mapping:
        i=i+1
        output_text = 'Processing '+ model +' ('+ str(i) +'/'+ str(num_models)+ ')' 
        print('\t', output_text)
        current_task_label.configure(text=output_text)
        progress_bar['value'] = 25 + 75*(i/num_models)
        window.update_idletasks()
        t_init2 = time()
        pdb_chain = model.split('-')[1]
        selection = '/'+model+'*//'+pdb_chain+'/'
        resi_list = list(resi_mapping[model].keys())
        for resi in resi_list:
            selection = selection + str(resi)+'+'
            selection1 = '/'+model+'///'+str(resi)+'/CA'
            selection2 = '/'+query_model+'///'+str(resi_mapping[model][resi])+'/CA'
            if show_align:
                cmd.distance(model+'_aln', selection1, selection2, mode = 0, label=0)
        cmd.set_atom_property('aln', 1, selection)
        print('\t','done!', '%6.2f'% (time()-t_init2), 'sec' )
    print('done!', '%7.2f'% (time()-t_init), 'sec' )


# Set final output graphics
def set_graphics(query_model, fl_model=''):
    stored.trunc_start = ''
    stored.trunc_end = ''
    if fl_model != '':
        cmd.iterate('first '+query_model, "stored.trunc_start = resi")
        cmd.iterate('last '+query_model, "stored.trunc_end = resi")

    print('Coloring models...')
    t_init = time()
    cmd.color('grey')
    cmd.color('forest', query_model+'*')
    if fl_model != '':
        cmd.color('forest', fl_model+' and resi '+stored.trunc_start+'-'+stored.trunc_end)
    cmd.color('palegreen', 'elem C and p.aln = 1')
    cmd.color('atomic', 'not (bb. + elem C)')
    print('done!','%7.2f'% (time()-t_init), 'sec' )

    print('Setting atom representations...')
    t_init = time()
    cmd.hide('everything')
    cmd.show('cartoon', 'polymer and not *_sticks')
    cmd.show('spheres', 'metals')
    cmd.show('sticks', sidechains_pro +' and p.aln = 1')
    cmd.show('sticks', query_model +'* and '+sidechains_pro)
    cmd.show('sticks', query_model +'_sticks and bb. and not (elem O + elem H)')
    cmd.show('dashes')
    print('done!','%7.2f'% (time()-t_init), 'sec' )

    print('Setting stick radius...')
    t_init = time()
    cmd.set('stick_radius', 0.1, sidechains_pro)
    print('done!','%7.2f'% (time()-t_init), 'sec' )

    print('Finalizing...')
    t_init = time()
    #cmd.delete('*_aln')
    cmd.disable(query_model+'_sticks')
    if fl_model != '':
        cmd.disable(fl_model)
    cmd.center(query_model)
    print('done!','%7.2f'% (time()-t_init), 'sec' )


# Function that executes when "Run Daliview" button is clicked
def run_daliview(daliview_inputs):
    progress_bar['value'] = 0
    dali_results_path = os.path.normpath(daliview_inputs['dali_results_path'])
    query_model_path = os.path.normpath(daliview_inputs['query_model_path'])
    start_model = int(daliview_inputs['start_model'])
    end_model = int(daliview_inputs['end_model'])
    struc_range = range(start_model, end_model+1)
    query_model = os.path.basename(query_model_path).split('.')[0]
    show_align = daliview_inputs['show_align']

    dali_stats = extract_stats(dali_results_path, struc_range)
    ttt_matrices = extract_ttt(dali_results_path, struc_range) 
    resi_mapping = extract_res_alignment(dali_results_path, struc_range)#, struc_range)

    load_pdbs(query_model_path, dali_stats, ttt_matrices, show_align)
    apply_prop(query_model, resi_mapping, show_align)
    set_graphics(query_model)#, fl_model)#, resi_mapping)
    
    output_text = 'Done! You may now close this window.' 
    current_task_label.configure(text=output_text)
    window.update_idletasks()
    return



# Default input settings
daliview_inputs = {
    'dali_results_path':'example_data/dali_results/s001C-25.txt',
    'query_model_path':'example_data/query_model/query_model.pdb',
    'start_model':'1',
    'end_model':'10', 
    'show_align':0}


# Widget behaviour functions:
def browse_dali_results():
    filepath = filedialog.askopenfilename(initialdir = "/", title = "Select a File", filetypes = (("Text files", "*.txt*"), ("All files", "*.*")))
    daliview_inputs['dali_results_path'] = filepath
    filename = os.path.basename(filepath)
    # Change label contents
    if len(filename) > 32:
        filename = '[...]%s' % filename[-32:]
    file_label_1b.configure(text=filename)

def browse_query():
    filepath = filedialog.askopenfilename(initialdir = "/", title = "Select a File", filetypes = (("Atomic coordinates", ".pdb .cif"), ("All files", "*.*")))
    daliview_inputs['query_model_path'] = filepath
    filename = os.path.basename(filepath)
    # Change label contents
    if len(filename) > 32:
        filename = '[...]%s' % filename[-32:]
    file_label_2b.configure(text=filename)

def only_numbers(char):
    return char.isdigit()

def change_resi_align():
    if show_align.get() == 0:
        daliview_inputs['show_align'] = 0
    if show_align.get() == 1:
        daliview_inputs['show_align'] = 1


# Function that executes when "Run Daliview" button is clicked
# Perform checks, then if everything checks out, execute run_daliview()
def run_button_press():
    daliview_inputs['start_model'] = start_model.get()
    daliview_inputs['end_model'] = end_model.get()
    if (not os.path.exists(daliview_inputs['dali_results_path'])) or (daliview_inputs['dali_results_path'] == ''):
        showinfo(title='Error!', message='Incorrect DALI results filepath!\nPlease select a DALI text file. ')
        return
    if (not os.path.exists(daliview_inputs['query_model_path'])) or (daliview_inputs['query_model_path'] == ''):
        showinfo(title='Error!', message='Incorrect Query PDB filepath!\nPlease select a PDB file. ')
        return
    if not (daliview_inputs['start_model'].isnumeric() and daliview_inputs['end_model'].isnumeric()):
        showinfo(title='Error!', message='Incorrect model range!\nPlease examine the DALI output file and chose a range of structures to display in Daliview. ')
        return
    run_daliview(daliview_inputs)
    window.destroy()

# Create the root window
# tkinter-based GUI to get inputs from user
window = Tk()
window.title('Daliview')


# Create widgets and buttons
file_label_1a = Label(window, text = "DALI results: ", width = 20, height = 2, anchor ='e')
file_label_1b = Label(window, text = os.path.basename(daliview_inputs['dali_results_path']), width = 20, height = 2, anchor ='c')
explore_button_1 = Button(window, text = "Browse File", command = browse_dali_results)
file_label_2a = Label(window, text = "Query PDB model: ", width = 20, height = 2, anchor ='e')
file_label_2b = Label(window, text = os.path.basename(daliview_inputs['query_model_path']), width = 20, height = 2, anchor ='c')
explore_button_2 = Button(window, text = "Browse File", command = browse_query)
model_range = Label(window, text = 'Display results range: ', width = 20, height = 2, anchor ='e')
validation = window.register(only_numbers)
start_model = Entry(window,validate="key", validatecommand=(validation, '%S'), width=10)
end_model = Entry(window,validate="key", validatecommand=(validation, '%S'),width=10)
start_model.insert(0, daliview_inputs['start_model'])
end_model.insert(0, daliview_inputs['end_model'])
show_align = IntVar()
align_cb = Checkbutton(window, text="Show aligned residues (longer processing time)", variable=show_align, onvalue=1, offvalue=0, command=change_resi_align)
run_button = Button(window, text = "Run Daliview!", command = run_button_press)
current_task_label = Label(window, text = 'Idle...', anchor='w')
progress_bar = Progressbar(window, orient=HORIZONTAL)

# Poition widgets
file_label_1a.grid(column = 1, row = 1)
file_label_1b.grid(column = 2, row = 1)
explore_button_1.grid(column = 3, row = 1, padx=10)
file_label_2a.grid(column = 1, row = 2)
file_label_2b.grid(column = 2, row = 2)
explore_button_2.grid(column = 3, row = 2, padx=10)
model_range.grid(column = 1, row = 3)
start_model.grid(column = 2, row = 3)
end_model.grid(column = 3, row = 3)
align_cb.grid(column = 1, columnspan=3, row=4, sticky = W+E, padx=10 ,pady=10)
progress_bar.grid(column=1, columnspan=3, row=5, sticky = W+E, padx=10 ,pady=10)
current_task_label.grid(column=1, columnspan=3, row=6, sticky = W+E, padx=10)
run_button.grid(column = 2,row = 7, pady=10)

# Let the window wait for any events
window.mainloop()







