from crypt import methods
from curses import flash
import statistics
from tkinter import E
from unicodedata import name
from flask import Flask, redirect, url_for, render_template, request, redirect
import pandas as pd
from IPython.display import HTML
import json

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
sum_dxy = robjects.globalenv['sum_dxy']
piPlot = robjects.globalenv['plot_pi']

population = robjects.StrVector(['British'])

def get_rsid_data(population,rsid):
    with open(population+'.AFGeno.json','r') as f:
        data = json.load(f)
        result_dict = {}
        try:
            result = data[rsid]
            result_dict['Population'] = population
            result_dict['Position']=result[0][0]
            result_dict['Reference']=result[0][1]
            result_dict['Alternative']=result[0][2]
            result_dict['Genotypes']=result[0][3]
            result_dict['AF']=result[1][2]
            return result_dict
        except:
            KeyError
            return "We have no data on this rsid"

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
class SNPForm(FlaskForm):
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
    populations = ['British','CEPH','Iberian','Finnish','Toscani']
    items = []
    for pop in populations:
        rsid_data = get_rsid_data(pop,snp_name)
        items.append(rsid_data)
    #Take out the positions etc so they only appear once, they are consistent across all 5 dictionaries.
    Position = rsid_data['Position']
    Reference = rsid_data['Reference']
    Alternative = rsid_data['Alternative']

    return render_template('rsid_data.html', items=items,Position=Position,Reference=Reference,Alternative=Alternative)



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
        with open('genepositions.json','r') as f:
            genepositions = json.load(f)
        gene_name = form.gene_name.data
        gene_search = gene_name+'.positions'
        gene_range = genepositions[gene_search]
        startgen = gene_range[0]
        print(startgen)
        endgen = gene_range[1]
        print('\n\n\n'+gene_name+'\n\n\n')
        return redirect(url_for('gene', gene_name=gene_name, gene_search=gene_search, gene_range=gene_range, startgen=startgen, endgen=endgen))
        #R code goes here
    return render_template('index_gene.html', form=form, gene_name=gene_name)

#making stats form for gene search
class GeneStats(FlaskForm):
    statistics = StringField("Enter Statistical Tests (Fst, d_xy, Pi, Tajima's D)", validators=[Required()])
    enterpop = StringField('Enter Population (British, Finnish, Iberian, Toscani and CEPH)', validators=[Required()])
    slidwid = StringField('Enter Sliding Window Parameters', validators=[Required()])
    submit = SubmitField('Submit')

#define a route called "gene" that accepts a gene name parameter
#code for finding gene alias
@app.route("/gene/<startgen>/<endgen>/<gene_name>", methods=['GET', 'POST'])
def gene(startgen, endgen, gene_name):
    

    positions(int(startgen),int(endgen))
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

    try:
        form = GeneStats()
        statistics = 0
        enterpop = 0
        slidwid = 0
        if request.method == 'POST':
            print("hello")
            if form.validate_on_submit():
                print('how to train your brain')
                statistics = form.statistics.data
                enterpop = form.enterpop.data
                slidwid = form.slidwid.data

                return redirect(url_for("generesult",statistics = statistics, enterpop=enterpop, slidwid=slidwid, startgen=startgen, endgen=endgen))
        return render_template('gene_view2.html', nsites=nsites, bisites=bisites ,psites=psites,form=form, statistics = statistics, enterpop=enterpop, slidwid=slidwid)
    except:
        return "We don't know about this??."


@app.route("/generesult", methods=['GET']) #need GET to pass variable from one function to another
def generesult():
    statistics = request.args.get("statistics") # to pass variable from one function to another
    enterpop = request.args.get("enterpop")
    slidwid = request.args.get("slidwid")
    endgen = request.args.get("endgen")
    startgen = request.args.get("startgen")

    #convert pops to capitalise/CEPH
    populations = enterpop.split(", ")
    for i in range(len(populations)):
        populations[i] = populations[i].capitalize()
        if populations[i] == 'Ceph':
            populations[i] = populations[i].upper()

    #need to create stringvector in R to keep populations stored
    rpopulations = robjects.StrVector(populations)

    #convert stats to lower case
    statistics = statistics.split(', ')
    for i in range(len(statistics)):
        statistics[i] = statistics[i].lower()

    #sliding window
    slidwid = slidwid.split(', ')
    window_size=int(slidwid[0])
    window_jump=int(slidwid[1])
    vcf_size = int(endgen) - int(startgen)
    print(vcf_size)
    slidingwindow(window_size,window_jump,vcf_size)
    analysis(int(slidwid[0]))
    preplot()

    

    #items to store next results
    items = []

    #iterate through stats
    for stat in statistics:
        if stat == 'fst':
            #get result
            fst_result = sum_fst(rpopulations)
            fst_result_dict = { key : fst_result.rx2(key)[0] for key in fst_result.names }
           
            for key,value in fst_result_dict.items():
                a_pair = dict(result = key + ':' + str(value))
                items.append(a_pair)
               
        if stat == 'd_xy':
            #get result
            dxy_result = sum_dxy(rpopulations)
            dxy_result_dict = { key : dxy_result.rx2(key)[0] for key in dxy_result.names }
               
            for key,value in dxy_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)    

        if stat == 'pi':
            #get result
            pi_result = sum_pi(rpopulations)
            pi_result_dict = { key : pi_result.rx2(key)[0] for key in pi_result.names }
               
            for key,value in pi_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)              


    return render_template("generesult.html",items=items)

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
    slidwid = 0
    if request.method == "POST":
        if form.validate_on_submit():
            print('hello')
            startgenomic = int(form.startgenomic.data)
            endgenomic = int(form.endgenomic.data)
            return redirect(url_for('genomic', startgenomic=startgenomic, endgenomic=endgenomic, slidwid = slidwid))
        else:
            return "We don't have info about this these genomic coordinates"
            #return redirect(url_for('genomic'))
    return render_template('index_genomic.html', form=form, startgenomic=int(startgenomic), endgenomic=int(endgenomic))

class Stats(FlaskForm):
    statistics = StringField("Enter Statistical Tests (Fst, d_xy, Pi, Tajima's D)", validators=[Required()])
    enterpop = StringField('Enter Population (British, Finnish, Iberian, Toscani and CEPH)', validators=[Required()])
    slidwid = StringField('Enter Sliding Window Parameters', validators=[Required()])
    submit = SubmitField('Submit')

@app.route("/genomic/<startgenomic>/<endgenomic>", methods=['GET', 'POST'])
def genomic(startgenomic, endgenomic):
    #load snp protein data from TSV file into pandas dataframe with snp name as index
    #df = pd.read_csv(snp_table_filename, sep='\t', index_col=2)
   

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
   
  

    try:
        #upper = df.query("POS >= %s" % startgenomic)
        #lower = upper.query("POS <= %s" % endgenomic)
        form = Stats()
        statistics = 0
        enterpop = 0
        slidwid = 0
        if request.method == 'POST':
            if form.validate_on_submit():
                print('how to train your brain')
                statistics = form.statistics.data
                enterpop = form.enterpop.data
                slidwid = form.slidwid.data

                return redirect(url_for("statspage",statistics = statistics, enterpop=enterpop, slidwid=slidwid, startgenomic=startgenomic, endgenomic=endgenomic))
        return render_template('genomic_view3.html', nsites=nsites, bisites=bisites ,psites=psites,form=form, statistics = statistics, enterpop=enterpop, slidwid=slidwid)
    except:
        return "We don't have any information about these coodinates."

#we need a page following this where user can select populations and then procedurally generate that form with results dxy, pi & fst means + graphs


@app.route("/statistic", methods=['GET']) #need GET to pass variable from one function to another
def statspage():
    statistics = request.args.get("statistics") # to pass variable from one function to another
    enterpop = request.args.get("enterpop")
    slidwid = request.args.get("slidwid")
    endgenomic = request.args.get("endgenomic")
    startgenomic = request.args.get("startgenomic")

    #convert pops to capitalise/CEPH
    populations = enterpop.split(", ")
    for i in range(len(populations)):
        populations[i] = populations[i].capitalize()
        if populations[i] == 'Ceph':
            populations[i] = populations[i].upper()

    #need to create stringvector in R to keep populations stored
    rpopulations = robjects.StrVector(populations)

    #convert stats to lower case
    statistics = statistics.split(', ')
    for i in range(len(statistics)):
        statistics[i] = statistics[i].lower()

    #sliding window
    slidwid = slidwid.split(', ')
    window_size=int(slidwid[0])
    window_jump=int(slidwid[1])
    vcf_size = int(endgenomic) - int(startgenomic)
    print(vcf_size)
    slidingwindow(window_size,window_jump,vcf_size)
    analysis(int(slidwid[0]))
    preplot()

    

    #items to store next results
    items = []

    #iterate through stats
    for stat in statistics:
        if stat == 'fst':
            #get result
            fst_result = sum_fst(rpopulations)
            fst_result_dict = { key : fst_result.rx2(key)[0] for key in fst_result.names }
           
            for key,value in fst_result_dict.items():
                a_pair = dict(result = key + ':' + str(value))
                items.append(a_pair)
               
        if stat == 'd_xy':
            #get result
            dxy_result = sum_dxy(rpopulations)
            dxy_result_dict = { key : dxy_result.rx2(key)[0] for key in dxy_result.names }
               
            for key,value in dxy_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)    

        if stat == 'pi':
            #get result
            pi_result = sum_pi(rpopulations)
            pi_result_dict = { key : pi_result.rx2(key)[0] for key in pi_result.names }
               
            for key,value in pi_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)              


    return render_template("result.html",items=items)


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