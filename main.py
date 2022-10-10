from flask import Flask , redirect , url_for , render_template , flash , request
import os
import glob
import pandas as pd
import sys
import itertools
import numpy as np
import plotly.express as px
from Bio import SeqIO

def VariableGenerator(length ,string):
    nucleotides = ["A" , "T", "G" , "C"]
    iterproduct = itertools.product(nucleotides, repeat = length)
    list = [''.join(iterproduct) for iterproduct in iterproduct]
    var_list = []
    for item in list:
        variable = str(item+string)
        var_list.append(variable)

    return var_list

def param_cleaner(param_input_df , param_selections):
    param_selections = set(param_selections).intersection(set(param_input_df.columns.values.tolist()))
    param_selections = list(param_selections)
    param_selections.insert(0,'Nucleotide')
    new_param_df = param_input_df[param_selections]
    return new_param_df
      
                 
# Parameter Check Function to create files based on parameters
def ParameterCheck(record_list,n):  
    global nString
    global param_input_df  
    global param_selections
    param_input_df = param_cleaner(param_input_df,param_selections)
    print(nString , " : " , param_input_df.columns)
    param_colnames = param_input_df.columns.values.tolist()
    Nucleic_acid_list = VariableGenerator(n,"")
    for record in record_list:
        seq = record.seq
        length = len(seq)            
        param_list = []
        na_list = []
        param_df = pd.DataFrame()
        print("\nSequence : " + record.name +"\nSequence Length : " + str(length) + "\nGenerating Parameter Outputs...")
        
        for param in param_colnames[1:] :
            param_list = []
            na_list = []
                   
            # Loop through entire sequence
            for i in range(length - n + 1):
                for na in Nucleic_acid_list:
                    if seq[i:i+n] == na:
                        #Append parameter and dinucleotide sequence in list.    
                        param_value = param_input_df.loc[param_input_df[param_colnames[0]] == na , param].iloc[0]
                        #print(param_value)
                        param_list.append(param_value)
                        na_list.append(na)

            #Write the parameter and dinucleotide sequences in a dataframe
            
            if "Nucleotide" not in param_df :
                param_df.insert(0, "Nucleotide" , na_list)
            param_df.insert(loc=param_input_df.columns.get_loc(param), column=param, value=param_list)
            param_df[param] = param_list
            
        # Export dataframe as csv          
        print(record.name + " writing to CSV ...")
        param_df.to_csv('Output/Parameters/Param_' + record.name + '_' + nString +"Nucleotide" + '.csv')
        
        spl = []
        for sp in range(len(param_df.index)):
            spl.append(int(sp + 1))
        #Plotting
        for col_val in param_df.columns.values[1:] :
            ## PLOTLY ##
            print("Parameters : Making Plotly Graph ...")
            fig = px.line(param_df, x=spl, y=col_val,)
            fig.update_xaxes(rangeslider_visible=True) 
            fig.update_layout(title=str(col_val)+" Concentration Plot",xaxis_title="Position",yaxis_title=str(col_val)+" Block Score")         
            fig.write_html('Output/Parameters/Plots/'+ col_val+ '_' + record.name + '_' + nString +"Nucleotide" + '.html')                                         
                                         
app = Flask(__name__)
app.secret_key = "wololo"
param_path = r'Parameter_Files/'
plot_path = r'static/Output/Parameters/Plots'
fasta_path = r'test.fasta'
record_list = list(SeqIO.parse(fasta_path, "fasta"))
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
        request.form.get("windowWidth")
        param_selections = list(filter(None, param_selections))
        flash(str(param_selections))
        return render_template('form.html' , content_list = param_select_list)
    else:
        flash(f"Please check some values :")
        return render_template('form.html' , content_list = param_select_list)
    
@app.route("/result")
def result():
    for n in range(2,4):
        ParameterCheck(record_list,n)
    plot_filenames = glob.glob(plot_path+"*")
    return render_template('form.html' , plots = plot_filenames)
        
    
@app.route("/")
def dnaprofiler():
    return render_template("DNAProfiler.html")

if __name__ == "__main__":
    app.run(debug=True)




