import pickle
from flask import Flask, redirect, render_template, request
import pandas as pd
import numpy as np
import os
import csv
import glob
from flask_wtf import FlaskForm
from wtforms import FileField, SubmitField, StringField, IntegerField, FloatField, validators
from padelpy import padeldescriptor
from werkzeug.utils import secure_filename
from wtforms.validators import InputRequired
from rdkit import Chem
from rdkit.Chem import Descriptors

from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp
from pubchempy import Compound, get_compounds
from prettytable import PrettyTable
import xgboost as xgb
import time as tm
from PIL import Image
from sklearn.preprocessing import StandardScaler
import base64
import io

#!/usr/bin/python

# Import modules for CGI handling 
import cgi, cgitb 


app = Flask(__name__)

app.config['SECRET_KEY'] = 'supersecretkey'
app.config['UPLOAD_FOLDER'] = 'static/assets/files'
app.config['EXTENSION'] = ['csv', 'txt']
app.config['STRUCTURE'] = 'static/assets/img/structure/'


# Create instance of FieldStorage 
form = cgi.FieldStorage() 

basedir = os.path.abspath(os.path.dirname(__file__))

input_file = 'input.csv'
input_file = os.path.join(basedir, app.config['UPLOAD_FOLDER'], input_file)

smi_file = 'input.smi'
smi_file = os.path.join(basedir, app.config['UPLOAD_FOLDER'], smi_file)

output_file = 'output.csv'
output_file = os.path.join(basedir, output_file)
output_file2 = 'output_2.csv'
output_file2 = os.path.join(basedir, output_file2)

result_file = os.path.join(basedir, 'PubChem.csv')
model_csv = os.path.join(basedir, 'model.csv')
xml = os.path.join(basedir, 'PubchemFingerprinter.xml')
logic_model = os.path.join(basedir, 'logic_model.pkl')


input_data = 'input_data.csv'
input_data = os.path.join(basedir, app.config['UPLOAD_FOLDER'], input_data)


default_img1=  "default1.png"
default_img1 = os.path.join(basedir, app.config['STRUCTURE'], default_img1)
default1 = Image.open(default_img1).convert('RGB')
data = io.BytesIO()
default1.save(data, 'JPEG')
def_img_1 = base64.b64encode(data.getvalue())
def_1 = def_img_1.decode("UTF-8")

default_img2=  "default2.png"
default_img2 = os.path.join(basedir, app.config['STRUCTURE'], default_img2)
default2 = Image.open(default_img2).convert('RGB')
data = io.BytesIO()
default2.save(data, 'JPEG')
def_img_2 = base64.b64encode(data.getvalue())
def_2 = def_img_2.decode("UTF-8")

structure_file = 'structure.png'
structure_file = os.path.join(basedir, app.config['UPLOAD_FOLDER'], structure_file)
struct1 = Image.open(structure_file)
data = io.BytesIO()
struct1.save(data, 'JPEG')
encode_img_1 = base64.b64encode(data.getvalue())
output_1 = encode_img_1.decode("UTF-8")


structure_prop_file = 'structure_prop.png'
structure_prop_file = os.path.join(basedir, app.config['UPLOAD_FOLDER'], structure_prop_file)
struct2 = Image.open(structure_prop_file)
struct2 = Image.open(structure_prop_file)
data = io.BytesIO()
struct2.save(data, 'JPEG')
encode_img_2 = base64.b64encode(data.getvalue())
output_2 = encode_img_2.decode("UTF-8")

def prediction_result(prediction):
    if prediction[0] == 0:
        return 'Active'
    else:
        return 'Inactive'
def prediction_result_all(prediction):
    for a in prediction:
        if a == 0:
           return 'Active'
        else:
            return 'Inactive'
            
def MOL_name(smiles):
    compounds = pcp.get_compounds(smiles, namespace='smiles')
    tm.sleep(0.1) 
    match = compounds[0]
    return match.iupac_name


def mol_W(smile_output):
    m = Chem.MolFromSmiles(smile_output)
    molwei = Descriptors.MolWt(m)
    return molwei
    

def cid_to_smiles(cid):
    c = pcp.Compound.from_cid(cid)
    gen_smiles = c.isomeric_smiles
    return gen_smiles


def name_to_smiles(name, identity):
    cs = get_compounds(name, identity )
    smiles = []
    for compound in cs:
        value = compound.isomeric_smiles
        smiles.append(value)
    if smiles == []:
        return 'invalid_input'
    else:
        return smiles[0]

def draw_structure(smile):
    if smile != 'invalid_input':
       mol_smiles = smile
       mol = Chem.MolFromSmiles(mol_smiles)
       Draw.MolToFile(mol, structure_file)

       for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
        Draw.MolToFile(mol, structure_prop_file)
    else:
        return 'invalid_input'

def allowed_file(name):
    if not '.' in name:
        return False
    ext = name.rsplit('.', 1)[1]
    if ext in app.config['EXTENSION']:
        return True
    else: 
        return False

def check_csv(name):
    ext = name.rsplit('.', 1)[1]
    if (ext in 'csv') & (ext in 'txt'):
        return True
    else: 
        return False


#pubchem calculation
def pubchem():

    padeldescriptor(mol_dir= smi_file,
                    d_file= result_file,  # 'Substructure.csv'
                    # descriptortypes='SubstructureFingerprint.xml',
                    descriptortypes= xml,
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=2,
                    removesalt=True,
                    log=True,
                    fingerprints=True)

    # Compare descriptor inputs
    load_desc = pd.read_csv(result_file)
    # read descriptors
    Xlist = list(pd.read_csv(model_csv).columns)

    subset = load_desc[np.intersect1d(load_desc.columns, Xlist)]

    sc = StandardScaler()
    desc_subset = sc.fit_transform(subset)
    print(desc_subset)

    return desc_subset

def disp_table(result):
    result_file = pd.read_csv(result_file)

@app.route('/tutorial')
def tutorial():
    return render_template('tutorial.html')


@app.route('/', methods=['GET', 'POST'])
def predict():
            
 # search file
    if request.method == 'POST' and "form2-submit" in request.form:
        default1 = def_1
        default2 = def_2 
        defauilt_text = 'N/A'
        smile = 'N/A'
        mol_name = 'N/A'

        uploaded_file = request.files['file']
        if uploaded_file.filename == '':
            return 'no file uploaded'  
             
        if allowed_file(uploaded_file.filename):
            print('file is allowed')
            uploaded_file.filename = 'input.csv'
            filename = secure_filename(uploaded_file.filename)
            uploaded_file.save(os.path.join(basedir, app.config['UPLOAD_FOLDER'], filename))
            txt_data = pd.read_csv(input_file)

            txt_data.to_csv(smi_file, sep='\t', header=False, index=False)
            smiles_file = pd.read_csv(smi_file)

            with open(smi_file) as f:
                list_smiles = [row.split()[0] for row in f]

            #mol names
            mol_list = []
            for i in list_smiles:
                i = i.upper()
                m_name = MOL_name(i)
                mol_list.append(m_name)
            
            # mol weight
            name_list = []
            for i in list_smiles:
                i = i.upper()
                names = mol_W(i)
                names = "%.2f" % names
                name_list.append(names)

            desc_subset = pubchem()
            load_model = pickle.load(open(logic_model , 'rb'))
            predict = load_model.predict(desc_subset)
            txt_data['Results'] = predict
            txt_data.loc[txt_data['Results'] == 1, 'Results'] = 'Inactive'
            txt_data.loc[txt_data['Results'] == 0, 'Results'] = 'Active'

            txt_data['molecule name'] = mol_list 
            txt_data['molecule weight'] = name_list
            
            txt_data.to_csv(output_file, sep=',', header=False)
            txt_data.to_csv(output_file2, sep=';', header=False)
            file_2 = open(output_file2) 
            
            #display first data
            with open(output_file, 'r') as f:
                reader = csv.reader(f)
                example_list = list(reader)
                print(example_list[0][1])
 
            ML_model = 'Logistic regression'
            mol_name = example_list[0][3]
            mw = example_list[0][4] 
            smile = list_smiles[0].upper()
            pred = example_list[0][2]
            draw_structure(smile) #draw structure

            structure1 = output_1
            structure2 = output_2
            
            return render_template('index.html', csv = file_2, pred = pred, ML_model = ML_model, mw =  mw, structure = structure1, structure_prop = structure2, smile = smile, mol_name=mol_name )
                  
        else: 
            return 'upload csv or txt file again'

    # single search
    elif request.method == 'POST' and "form1-submit" in request.form:
        input = request.form['selectinput']
        if(input == 'cid'):
            cid_value = request.form.get('value')
            output_value = int(cid_value)
            output_smile = cid_to_smiles(output_value)
            m_name = MOL_name(output_smile)
            identity = cid_value

        elif (input == 'smiles'):
            smile_value = request.form.get('value').upper()
            output_smile = smile_value
            identity = smile_value
            m_name = MOL_name(output_smile)

        elif (input == 'name'):
            identity = input
            mol_name = request.form.get('value')
            if mol_name != 'invalid_input':
                value = name_to_smiles(mol_name, identity)
                output_smile = value
                m_name = MOL_name(value)
            else:
                return render_template('message.html')
                
        else:
            return render_template('message.html')
            
        if request.form.get('value') == '':
            return 'no value entered'
    
       
        header = ['smile', 'identity']

        csvfile = open(input_data, 'w', newline='', encoding='utf-8')
        data = csv.writer(csvfile, delimiter=' ')
        data.writerow(header)
        data.writerow([output_smile, identity])
        csvfile.close()

        open_data = pd.read_csv(input_data)
        open_data.to_csv(smi_file, sep='\t', header=False, index=False)

        desc_subset = pubchem()
        mol = mol_W(output_smile) #molecular weight
        model = request.form['selectalgorithm']
        if model == 'LogisticRegression':
           load_model = pickle.load(open(logic_model , 'rb'))
           predict = load_model.predict(desc_subset)
           model = model
           molecular_weight = "%.2f" % mol
           result = prediction_result(predict)
           draw_structure(output_smile) #draw structure
           structure_prop_files = output_2
           structure_files = output_1
           return render_template('index.html',  pred= result, ML_model = model, mw =  molecular_weight, structure = structure_files, structure_prop = structure_prop_files, smile = output_smile, mol_name=m_name)
           

        elif model =='XGB':
            load_model = pickle.load(open(logic_model , 'rb'))
            predict = load_model.predict(desc_subset)
            print(prediction_result(predict))
            result = prediction_result(predict)
            model = model
            molecular_weight = "%.2f" % mol
            draw_structure(output_smile) #draw structure
            structure_prop_files = output_2
            structure_files = output_1
            return render_template('index.html',  pred= result, ML_model = model, mw =  molecular_weight, structure = structure_files, structure_prop = structure_prop_files, smile = output_smile, mol_name=m_name)


        elif model == 'RF':
            load_model = pickle.load(open(logic_model, 'rb'))
            predict = load_model.predict(desc_subset) 
            print(prediction_result(predict))
            result = prediction_result(predict)
            model= model
            molecular_weight = "%.2f" % mol
            draw_structure(output_smile) #draw structure
            structure_prop_files = output_2
            structure_files = output_1
            return render_template('index.html',  pred= result, ML_model = model, mw =  molecular_weight, structure = structure_files, structure_prop = structure_prop_files, smile = output_smile, mol_name=m_name)

        elif model == 'ANN':
            load_model = pickle.load(open(logic_model, 'rb'))
            predict = load_model.predict(desc_subset) 
            result = prediction_result(predict)
            molecular_weight = "%.2f" % mol
            model = model
            draw_structure(output_smile) #draw structure
            structure_prop_files = output_2
            structure_files = output_1
            return render_template('index.html',  pred = result, ML_model = model, mw =  molecular_weight, structure = structure_files, structure_prop = structure_prop_files, smile = output_smile, mol_name=m_name)

        else:
            return render_template('message.html')

    default1 = def_1
    default2 = def_2 
    defauilt_text = 'N/A'
    smile = 'N/A'
    mol_name = 'N/A'
    return render_template('index.html' , pred = defauilt_text, ML_model = defauilt_text, mw =  defauilt_text, structure = default1, structure_prop = default2, smile = smile, mol_name=mol_name)

if __name__ == "__main__":
    app.run(debug=True)       