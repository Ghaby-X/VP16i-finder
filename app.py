from flask import Flask, render_template, request, current_app

import numpy as np
import os
import pickle

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from myfunctions import compute_morganfps
from myfunctions import morgan_csv, rfc_csv_result
from myfunctions import format_smiles

from werkzeug.utils import secure_filename
from datetime import datetime

ALLOWED_EXTENSIONS = set(['csv'])

def allowed_file(filename):
    return '.' in filename and \
            filename.replit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def create_model(type):
    if type == "GBM":
        filepath = os.path.join(current_app.config['STATIC_FOLDER'], 'Models', 'gbm_pipeline.pkl')
    elif type == "SVM":
        filepath = os.path.join(current_app.config['STATIC_FOLDER'], 'Models', 'svm_pipeline.pkl')
    elif type == "RFC":
        filepath = os.path.join(current_app.config['STATIC_FOLDER'], 'Models', 'rfc_pipeline.pkl')
    else:
        return "Pick a valid model type"
    
    with open(filepath, 'rb') as file:
        model = pickle.load(file)
    
    return model


app = Flask(__name__)
app.config['STATIC_FOLDER'] = 'static'



@app.route('/')
def home():
    return render_template('index.html')

@app.route("/upload")
def upload():
    return render_template("upload.html")

@app.route('/results', methods = ['Post', 'Get'])
def results():
    if request.method == 'POST':
        smile = request.form.get('smile')
        model_type = request.form.get('model')

        mol = Chem.MolFromSmiles(smile)
        data = compute_morganfps(mol)


        #Saving the image to the directory --->> /static/Image/mol_image.png
        filepath = os.path.join(app.config['STATIC_FOLDER'], 'Image', "mol_image.png")
        try:
            Chem.Draw.MolToFile(mol,filepath, size=(300, 300))
        except:
            return "not a valid smile"
        
        model = create_model(model_type)
        
        activity = model.predict(data)
        if activity == 1:
            activity =  "Active" 
            confidence = model.predict_proba(data)[0][1]
        else:
            activity =  "Inactive"
            confidence = model.predict_proba(data)[0][0]


        try:
            return render_template("result.html", smile = smile, activity = activity, confidence = confidence)
        except:
            return "can't process image"
    else:
        return 'Method not allowed'


@app.route('/results_csv', methods = ['Post', 'Get'])
def results_csv():
    if request.method == "POST":
        file = request.files["smile_csv"]
        model_type = request.form.get('model')

        model = create_model(model_type)
        
        descriptors, smiles = morgan_csv(file)
        df_table = rfc_csv_result(descriptors, smiles, model)

        
        styled_df = df_table.style  \
                            .format({'Canonical Smiles' : format_smiles})

        data_table = styled_df.to_html(escape = False)

        return render_template("results_csv.html", df_html = data_table)
    else:
        return "Method not allowed"

if __name__ == '__main__':
    app.run()