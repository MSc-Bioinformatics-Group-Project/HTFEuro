from curses import flash
from tkinter import E
from unicodedata import name
from flask import Flask, redirect, url_for, render_template, request, redirect
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
populations = 'chr21.five.populations'

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
    #df = pd.read_csv(snp_table_filename, sep='\t', index_col=1)
    pop = pd.read_csv(populations, sep='\t')
    #snp_name = snp_name.upper() #ensure name is in captial letters

    try: #try to extract row for specific snp
        row = pop.loc[snp_name]
        #if snp is found, return some info about it
        return render_template('snp_view.html', name=snp_name, geneid = row.GENEID, \
        genomic_coordinate= row.POS, gene_alias = row.ALIAS)
    except:
        #if protein not found a key error is thrown 
        return "We don't have any information about a snp called %s." % snp_name




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
#code for finding gene alias
@app.route("/gene/<gene_name>")
def gene(gene_name):

    #load snp protein data from TSV file into pandas dataframe with snp name as index
    df = pd.read_csv(snp_table_filename, sep='\t')

    #snp_name = snp_name.upper() #ensure name is in captial letters
    
    try:
        snp_name = df.query('GENEID == "%s"' % gene_name)
        snp_name2 = df.query('ALIAS == "%s"' % gene_name)

        snp_name3=[snp_name, snp_name2]

        snp_name4 = pd.concat(snp_name3)
        return render_template('gene_view2.html', tables=[snp_name4.to_html(classes='data', header="true")])
    
    except:
        #if gene or alias not found a key error is thrown 
        return "We don't have any information about a gene called %s." % gene_name

    #code for finding single gene name
    #try: #try to extract row for specific gene name
        #row = df.loc[gene_name]
        #if snp is found, return some info about it
        #return render_template('gene_view.html', name=gene_name, rsvalue = row.ID, \
        #genomic_coordinate= row.POS)
    #except:
        #if protein not found a key error is thrown 
       # return "We don't have any information about a gene called %s." % gene_name




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
            endgenomic = form.endgenomic.data
            return redirect(url_for('genomic', startgenomic=startgenomic, endgenomic=endgenomic))
        else:
            return "We don't have info about this these genomic coordinates"
            #return redirect(url_for('genomic'))

    return render_template('index_genomic.html', form=form, startgenomic=startgenomic, endgenomic=endgenomic)

@app.route("/genomic/<startgenomic>/<endgenomic>")
def genomic(startgenomic, endgenomic):

    #load snp protein data from TSV file into pandas dataframe with snp name as index
    df = pd.read_csv(snp_table_filename, sep='\t', index_col=2)
    #snp_name = snp_name.upper() #ensure name is in captial letters
    
    try:
        upper = df.query("POS >= %s" % startgenomic)
        lower = upper.query("POS <= %s" % endgenomic)

        return render_template('genomic_view.html', tables=[lower.to_html(classes='data', header="true")])
    
    except:
        return "We don't have any information about these coodinates."
    
#genomic coordinates results pahge
@app.route("/genomic/<startgenomic>/<endgenomic>/<coordinate_result>")
def coord_results(coordinate_result):
    return render_template('coordinate_result.html')

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