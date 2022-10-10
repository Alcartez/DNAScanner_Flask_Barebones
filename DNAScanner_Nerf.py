# Importing libraries
import os
import sys
import itertools
import numpy as np
import pandas as pd
import plotly.express as px
from Bio import SeqIO

# Setting up Environment and Import Files

print("\n Setting Up Envioronment ...")
print("\n Making Output Directories ...")

## Setting up environment ##

# Setting up directory tree
if not os.path.exists('static/Output/Parameters/Plots'):
    os.makedirs('Output/Parameters/Plots')

if not os.path.exists('static/Output/Nucleotide_Concentration/Plots'):
    os.makedirs('Output/Nucleotide_Concentration/Plots')

if not os.path.exists("Parameter_Files"):
    os.makedirs("Parameter_Files")

print("\nImporting and parsing input files ...")

#Import Fasta
fasta_path = sys.argv[1]
print(fasta_path)
record_list = list(SeqIO.parse(fasta_path, "fasta"))

# Import WindowWidth Variable
#windowWidth = 4
#param_selections = ['Arule', 'Trule', 'Grule', 'Bendability']


#Function to Generate Variables for every nucleotide
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

def Window_slidingBlock(record_list,windowWidth,n):
    global nString
    global repairIndex
    spl = [] # List for starting positions for window
    epl = [] # List for ending positions for window
    seq_list = [] # Declare empty list of sequences
    Block_size = int((windowWidth - n + 1) * len(record_list) )
    print("Block Size : " , Block_size )    
    Nucleic_acid_list = VariableGenerator(n,"")  

    for record in record_list:
        seq_list.append(record.seq) # We make a list of sequences as seq_list      

    #Looping through the list of sequences 
    seq = max(seq_list,key=len)      # Taking the longest sequence as the length of the DNA Scanner to cover most ground
    length = len(seq)
    Block_score_variables = VariableGenerator(n,"_Block_score")

    for var in Block_score_variables:
        globals()[var] = []        

    # Loop through entire sequence
    for i in range(length - windowWidth + 1):
        Count_variables = VariableGenerator(n,"Count")
        for var in Count_variables:
            globals()[var] = 0            

        # In each iteration, loop though window
        for j in range(i, windowWidth + i - repairIndex):
            for seq in seq_list :
                for na in Nucleic_acid_list:
                    if seq[j:j+n] == na:
                        globals()[na+"Count"] += 1                       

        # Update counts of each nucleotide in respective dictionary
        for z in Nucleic_acid_list:
            a = globals()[str(z+"_Block_score")]
            b = globals()[str(z+"Count")]
            a.append(b/Block_size)   
            start_pos = i+1
            end_pos = windowWidth+i                 
        spl.append(start_pos)
        epl.append(end_pos) 
        
    #Creating a dataframe 
    df = pd.DataFrame()
    df.insert(0, "Starting Position" ,spl)
    df.insert(1, "Ending Position" ,epl)
    print("Creating Dataframe")    

    for j in range(len(Nucleic_acid_list)):
        df.insert(j+2, str(Nucleic_acid_list[j]), globals()[str(Nucleic_acid_list[j]+"_Block_score")])

    # Exporting the output as csv 
    print("Writing to CSV ...")
    df.to_csv('Output/Nucleotide_Concentration/DNAScanner_'+str(nString)+'NucleotideRule__Output.csv')   

    # Plotting
    print("Nucleotides : Generating Plot(s) ...")

    for col_val in df.columns.values[2:] :
        ## PLOTLY ##
        print("Making Plotly Graph ...")
        fig = px.line(df, x=df.columns.values[0], y=col_val,)
        fig.update_xaxes(rangeslider_visible=True) 
        fig.update_layout(title=str(col_val)+" Concentration Plot",xaxis_title="Position",yaxis_title=str(col_val)+" Block Score")         
        fig.write_html("Output/Nucleotide_Concentration/Plots/"+str(col_val)+"_Nucleotide_Rule_Plot.html")               
    print("Done .... \n")                                
                   
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
                                         

for n in range(2,4):
    ## User Inputs ##
    print("\nImporting parameters ...")
    if n == 2:
        param_input_df = pd.read_csv(r'Parameter_Files\Parameter_Sheet_Dinucleotide - Sheet1.csv')
        #param_input_df = param_cleaner(param_input_df,param_selections)
        #print(param_input_df.columns)
        nString = "Di"
        repairIndex = 1
    elif n == 3:
        param_input_df = pd.read_csv(r'Parameter_Files\Parameter_Files_Trinucleotide - Sheet1.csv')
        #param_input_df = param_cleaner(param_input_df,param_selections)
        #print(param_input_df.columns)
        nString = "Tri"
        repairIndex = 1
    else :
        print("\n DON'T BREAK MY SCRIPT!!") 

    #Functions Called
    ParameterCheck(record_list,n)
    Window_slidingBlock(record_list,windowWidth,n)