from curses import flash
from tkinter import E
from unicodedata import name
from flask import Flask, redirect, url_for, render_template, request, redirect, Blueprint
import pandas as pd
from IPython.display import HTML

#import libraries needed to create and process forms
from flask_wtf import FlaskForm
from flask_wtf import Form
from wtforms import StringField, SubmitField
from wtforms.validators import Required
from wtforms import IntegerField

#create a flask application object
app = Flask(__name__)
#we need a scret key attribute for secure forms
app.config['SECRET_KEY'] = 'change this unsecure key'

#where to find the SNP information
snp_table_filename = 'snp_table.tsv'

#define the action for the top level route
#this is the front-page
@app.route("/")
def home():
    return render_template("home.html")

############################################
#SNP
#create a class to define the form
class   SNPForm(FlaskForm):
    snp_name = StringField('Enter a valid SNP name:', validators=[Required()])
    submit = SubmitField('Submit')

#Search for SNP
@app.route("/snpsearch", methods=['GET', 'POST'])
def index():
    #this route has been updated to use a template containing a form
    form = SNPForm() #create form to pass to template
    snp_name = None
    if form.validate_on_submit():
        snp_name = form.snp_name.data
        print('\n\n\n'+snp_name+'\n\n\n')
        return redirect(url_for('snp', snp_name=snp_name))
    return render_template('index_page.html', form=form, snp_name=snp_name)


#define a route called "snp" that accepts a snp name parameter
@app.route("/snp/<snp_name>")
def snp(snp_name):

    #load snp protein data from TSV file into pandas dataframe with snp name as index
    df = pd.read_csv(snp_table_filename, sep='\t', index_col=1)
    #snp_name = snp_name.upper() #ensure name is in captial letters

    try: #try to extract row for specific snp
        row = df.loc[snp_name]
        #if snp is found, return some info about it
        return render_template('snp_view.html', name=snp_name, geneid = row.GENEID, \
        genomic_coordinate= row.POS)
    except:
        #if protein not found a key error is thrown 
        return "We don't have any information about a snp called %s." % snp_name


        #return '<h1>' + snp_name + '</h1>' \
        #+ '<p>Full name: ' + row.ID + '</p>' \
        #+ '<p>Genomic coordinate: ' + row.POS + '</p' \
    e#xcept:
        #if protein is not found a key error is thrown and we end up here
        #return "We don't have any information about a snp called %s." % snp_name

#Gene name
#create a class to define the form
#@app.route("/genename", methods=['GET', 'POST'])
#def genename():
    geneform = gene()
    finding_gene = geneform.search.data 
    gene_aliases = get_all_aliases()

    if geneform.validate_on_submit():
        for x in range(len(gene_aliases)):  
            if finding_gene.upper() in gene_aliases[x]:  
                return redirect(url_for('generesult', finding_gene=finding_gene))
    else:
        "We don't have any information about a gene called %s." % genename
    return render_template("generesult.html", title= "list of genes", form=form)

##################################################################
##Gene

#create a class to define the form
class GeneForm(FlaskForm):
    gene_name = StringField('Enter a valid Gene ID:', validators=[Required()])
    submit = SubmitField('Submit')

#Search for gene
@app.route("/genesearch", methods=['GET', 'POST'])
def index_gene():
    #this route has been updated to use a template containing a form
    form = GeneForm() #create form to pass to template
    gene_name = None
    if form.validate_on_submit():
        gene_name = form.gene_name.data
        print('\n\n\n'+gene_name+'\n\n\n')
        return redirect(url_for('gene', gene_name=gene_name))
    return render_template('index_gene.html', form=form, gene_name=gene_name)


#define a route called "gene" that accepts a gene name parameter
@app.route("/gene/<gene_name>")
def gene(gene_name):

    #load snp protein data from TSV file into pandas dataframe with snp name as index
    df = pd.read_csv(snp_table_filename, sep='\t', index_col=4)
    #df2 = pd.read_csv(snp_table_filename, sep='\t', index_col=4)
    #alias = df[['ID', 'Gene Name']]
    #print(alias)
    #snp_name = snp_name.upper() #ensure name is in captial letters

    #row = df.loc[gene_name]
    #row2 = df2.loc[gene_name]
    
    #if gene_name in df["ID"]:
        #return render_template('gene_view.html', name=gene_name, rsvalue = row.ID, \
        #genomic_coordinate= row.POS)
    #if gene_name in row:
        #return render_template('gene_view.html', name=gene_name, rsvalue = row.ID, \
        #genomic_coordinate= row.POS)
    #else:
        #return "We don't have any information about a gene called %s." % gene_name

    try: #try to extract row for specific gene name
        row = df.loc[gene_name]
        #if snp is found, return some info about it
        return render_template('gene_view.html', name=gene_name, rsvalue = row.ID, \
        genomic_coordinate= row.POS)
    except:
        #if protein not found a key error is thrown 
        return "We don't have any information about a gene called %s." % gene_name




#############################################################
##genomic coordinate

class Coord(FlaskForm):
    startgenomic = IntegerField('Start Genomic coordinate:', validators=[Required()])
    endgenomic = IntegerField('End Genomic coordinate:', validators=[Required()])
    submit = SubmitField('Submit')

@app.route("/genomic", methods=['GET', 'POST'])
def index_genomic():
    form = Coord()
    startgenomic = None
    endgenomic = None

    if request.method == "POST":
        if form.validate_on_submit():
            startgenomic = form.startgenomic.data
            #print('\n\n\n'+str(startgenomic)+'\n\n\n')
            endgenomic = form.endgenomic.data
            #print('\n\n\n'+str(endgenomic)+'\n\n\n')
        #print('\n\n\n'+endgenomic+'\n\n\n')
            return redirect(url_for('genomic', startgenomic=startgenomic, endgenomic=endgenomic))
        else:
            #flash("We don't have info about this start genomic coordinate")
            return redirect(url_for('genomic'))

    return render_template('index_genomic.html', form=form, startgenomic=startgenomic, endgenomic=endgenomic)

@app.route("/genomic/<startgenomic>")
def genomic(startgenomic, endgenomic):

    #load snp protein data from TSV file into pandas dataframe with snp name as index
    df = pd.read_csv(snp_table_filename, sep='\t', index_col=2)
    #snp_name = snp_name.upper() #ensure name is in captial letters

    upper = df.query("POS > %s" % startgenomic)
    lower = upper.query("POS < %s" % endgenomic)
    result = upper.to_html()

    return render_template('genomic_view.html', name=startgenomic, rsvalue = lower.ID, geneid = lower.GENEID) 

    #coord_range = range(int(startgenomic), int(endgenomic))
    
    #try: #try to extract row for specific start gene 
        #start_row = df.loc[int(startgenomic)]
       # end_row = df.loc[int(endgenomic)]

        #if snp is found, return some info about it
        #return render_template('genomic_view.html', name=startgenomic, rsvalue = start_row.ID, geneid = start_row.GENEID,
        #endname=endgenomic, endrsvalue = end_row.ID, endgeneid = end_row.GENEID) 
   # except:
        #if protein not found a key error is thrown 
        #return "We don't have any information about a gene starting with %s." % startgenomic % endgenomic



#about page
@app.route("/about")
def about():
    #return "Hello! Welcome to SNPEuro! <h1> HELLO </h1>"
    return render_template("about.html")

#help page
@app.route("/help")
def help():
    return render_template("help.html")


@app.route("/<name>")
def user(name):
    return f"Hello {name}!"

@app.route("/admin")
def admin():
    return redirect(url_for("user", name="Admin!"))

#starts web sever
#runs the application
if __name__== "__main__":
    app.run(debug=True)