from flask import Flask , redirect , url_for , render_template , flash , request
import os
import glob
import pandas as pd

app = Flask(__name__)
app.secret_key = "wololo"
param_path = r'Parameter_Files/'

res = []
param_select_list = []

for file in glob.glob(param_path+"*"):
    file_content_list = pd.read_csv(file).columns.values.tolist()
    for item in file_content_list: 
        if item not in param_select_list: 
            param_select_list.append(item) 

@app.route('/' , methods =["GET", "POST"])
def home():
    if request.method == "POST":
        # session.permanent = True  # To make session permanent
        param_selections = []
        for x in param_select_list:
            param_selections.append(request.form.getlist(x)) # Request the name from the form.
        param_selections = list(filter(None, param_selections))
        flash(str(param_selections))
        return render_template('form.html' , content_list = param_select_list)
    else:
        flash(f"Please check some values :")
        return render_template('form.html' , content_list = param_select_list)
    

@app.route("/")
def dnaprofiler():
    return render_template("DNAProfiler.html")

if __name__ == "__main__":
    app.run(debug=True)


