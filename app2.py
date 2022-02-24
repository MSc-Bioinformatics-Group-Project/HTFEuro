from crypt import methods
from curses import flash
import statistics
from tkinter import E
from unicodedata import name
from flask import Flask, redirect, url_for, render_template, request, redirect
import pandas as pd
from IPython.display import HTML

#import libraries needed to create and process forms
from flask_wtf import FlaskForm
from flask_wtf import Form
from wtforms import StringField, SubmitField, RadioField
from wtforms.validators import Required
from wtforms import IntegerField

#import R packages

import rpy2.robjects as robjects
r = robjects.r
r.source('big boy no variables.R')
positions = robjects.globalenv['getpositiondata']
SNPCount = robjects.globalenv['SNPCount']
sumdata = robjects.globalenv['sumdata']
neutrality = robjects.globalenv['neutrality']
slidingwindow = robjects.globalenv['SlidingWindow']
analysis = robjects.globalenv['analysis']
preplot = robjects.globalenv['pre_plot']
sum_pi = robjects.globalenv['sum_pi']
sum_fst = robjects.globalenv['sum_fst']
piPlot = robjects.globalenv['plot_pi']

population = robjects.StrVector(['British'])

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
        #R code goes here

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
    startgenomic = StringField('Start Genomic coordinate:', validators=[Required()])
    endgenomic = StringField('End Genomic coordinate:', validators=[Required()])
    submit = SubmitField('Submit')

@app.route("/genomic", methods=['GET', 'POST'])
def index_genomic():
    form = Coord()
    startgenomic = 0
    endgenomic = 0
    if request.method == "POST":
        if form.validate_on_submit():
            print('hello')
            startgenomic = int(form.startgenomic.data)
            endgenomic = int(form.endgenomic.data)
            return redirect(url_for('genomic', startgenomic=startgenomic, endgenomic=endgenomic))
        else:
            return "We don't have info about this these genomic coordinates"
            #return redirect(url_for('genomic'))
    return render_template('index_genomic.html', form=form, startgenomic=int(startgenomic), endgenomic=int(endgenomic))

class Stats(FlaskForm):
    statistics = StringField('Enter Stats', validators=[Required()])
    submit = SubmitField('Submit')

@app.route("/genomic/<startgenomic>/<endgenomic>", methods=['GET', 'POST'])
def genomic(startgenomic, endgenomic):
    #load snp protein data from TSV file into pandas dataframe with snp name as index
    df = pd.read_csv(snp_table_filename, sep='\t', index_col=2)
   
    positions(int(startgenomic),int(endgenomic))
    summary = list(sumdata())
   
    nsites=summary[0]
    bisites=summary[1]
    psites=summary[5]
   
    if bisites == 0 and psites == 0:
        #return no SNPS found
        return "We don't have info about this these genomic coordinates"
    if bisites == 1:    
        #return RSID search for that SNP.
        return "We don't have info about this these genomic coordinates"
    if bisites == 0 and psites == 1:
        #return RSID search for that SNP (polyallelic)
        return "We don't have info about this these genomic coordinates"
   
    slidingwindow(10000,2500,int(endgenomic)-int(startgenomic))
    analysis(10000)

    preplot()
    pre_pi = sum_pi(population)
    dict_pi = { key : pre_pi.rx2(key)[0] for key in pre_pi.names }
    pi = dict_pi['British_pi']


    try:
        #upper = df.query("POS >= %s" % startgenomic)
        #lower = upper.query("POS <= %s" % endgenomic)
        form = Stats()
        statistics = 0
        if request.method == 'POST':
            if form.validate_on_submit():
                print('how to train your brain')
                statistics = form.statistics.data
                return redirect(url_for("statspage",statistics = statistics))
        return render_template('genomic_view3.html', nsites=nsites, bisites=bisites ,psites=psites, pi = pi,form=form, statistics = statistics)
    except:
        return "We don't have any information about these coodinates."

#we need a page following this where user can select populations and then procedurally generate that form with results dxy, pi & fst means + graphs


@app.route("/statistic", methods=['GET'])
def statspage():
    statistics = request.args.get("statistics")
    print(statistics)
    if statistics == 'fst':
        #need to create stringvector in R to keep populations stored
        liam = robjects.StrVector(['British'])
        answer1 = sum_fst(liam)
        answer = answer1[0]

    return render_template("statistics.html",answer=answer, statistics=statistics)


#dummy page

# class GenerateForm(FlaskForm):
#     radio_fields = RadioField('', choices=[])
#     submit = SubmitField('submit')


# form = GenerateForm() # Instantiate it

# form.radio_fields.label = 'Label Example'
# form.radio_fields.choices = [('value_1', 'description'), ('value_2', 'description')]

# render_template('file.html', form=form)


# #Search for gene
# @app.route("/test", methods=['GET', 'POST'])
# def index_gene():  
#     #this route has been updated to use a template containing a form
#     form = GenerateForm() # Instantiate it

#     form.radio_fields.label = 'Label Example'
#     form.radio_fields.choices = [('value_1', 'description'), ('value_2', 'description')]

#     return render_template('file.html', form=form)


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