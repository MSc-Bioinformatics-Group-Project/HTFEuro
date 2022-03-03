from crypt import methods
from curses import flash
from tkinter import E
from unicodedata import name
from flask import Flask, redirect, url_for, render_template, request, redirect, send_from_directory
from IPython.display import HTML
import json

#import libraries needed to create and process forms
from flask_wtf import FlaskForm, Form
from wtforms import StringField, SubmitField, IntegerField
from wtforms.validators import Required

#import R packages

import rpy2.robjects as robjects
r = robjects.r

r.source('SNPEuro.R')

#first analysis
positions = robjects.globalenv['getpositiondata']
SNPCount = robjects.globalenv['SNPCount']
sumdata = robjects.globalenv['sumdata']
neutrality = robjects.globalenv['neutrality']
slidingwindow = robjects.globalenv['SlidingWindow']
analysis = robjects.globalenv['analysis']
preplot = robjects.globalenv['pre_plot']

#stat summations
sum_pi = robjects.globalenv['sum_pi']
sum_fst = robjects.globalenv['sum_fst']
sum_dxy = robjects.globalenv['sum_dxy']
sum_hpw = robjects.globalenv['sum_hpw']
sum_tajd = robjects.globalenv['sum_tajd']

#stat plots
piPlot = robjects.globalenv['plot_pi']
dxyPlot = robjects.globalenv['plot_dxy']
fstPlot = robjects.globalenv['plot_fst']
hpwPlot = robjects.globalenv['plot_hpw']
tajdPlot = robjects.globalenv['plot_tajd']


#Python functions

#function for AF/GT presentation for ONE SNP search
def get_rsid_data(population,rsid):
    with open(population+'.AFGeno.json','r') as f:
        data = json.load(f)
        result_dict = {}
        #create dictionary with below results, pull through to HTML
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



#functions to present AF/GT for multiple SNPs returned from gene search


#step 1 - gene search, get all rsids

def get_gene_rsids(gene_name):
    with open('genes.json','r') as f:
        genedata = json.load(f)    
    rsids = genedata[gene_name]
    return rsids


#step 2 - put rsids in dict as keys, add ref, alt and then two lists of frequencies & GT in standard order

def create_gene_json_template(rsids):
    gene_json = {}
    for rsid in rsids:
        if rsid not in gene_json:
            gene_json[rsid] = None
    #add in the ref and alt, use one json to get this data
    with open('British.AFGeno.json','r') as f:
        AFGenodata = json.load(f)
    for rsid in gene_json.keys():
        position = AFGenodata[rsid][0][0]
        reference = AFGenodata[rsid][0][1]
        alternative = AFGenodata[rsid][0][2]
        info = [position,reference,alternative]
        gene_json[rsid] = info
    return gene_json

def add_gene_json_population_stats(gene_json, populations):
    #get the GT and AF for each population JSON
    for i in populations:
        with open(i+'.AFGeno.json','r') as f:
            data = json.load(f)
            #loop over the rsid keys to get the data out of the jsons
            for rsid in gene_json.keys():
                GT = data[rsid][0][3]
                AF = data[rsid][1][2]
                stats = [GT,AF]
                gene_json[rsid].append(stats)
    return gene_json  

#function to present RSID stats from coord search, will use with two functions above replacing get_gene_rsids

def get_coord_rsids(start,end):
    #will take the start and end coords, make a range and check if the position exists in our rsid/position json - then it will pull out the RSIDs
    with open('Position_rsid.json','r') as f:
        position_rsid = json.load(f)
    position_int = []
    for position in list(position_rsid.keys()):
        position_int.append(int(position))
    rsids = []
    for i in range(start,end):
        if i in position_int:
            #(str(i) is '124111' position key)
            rsids.append(position_rsid[str(i)])
    return rsids



#create a flask application object
app = Flask(__name__)
#we need a scret key attribute for secure forms
app.config['SECRET_KEY'] = 'change this unsecure key'

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

    #rsid_data only pulls the AF and GT out of the dictionaries in items

    return render_template('rsid_data.html', items=items,Position=Position,Reference=Reference,Alternative=Alternative)



##################################################################
##Gene

#create a class to define the form
class GeneForm(FlaskForm):
    gene_name = StringField('Enter a valid Gene ID or alias:', validators=[Required()])
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
        with open('gene_alias.json','r') as f:
            genealiases = json.load(f)
        gene_name = form.gene_name.data
        gene_search = gene_name+'.positions'

        #check for alias 
        try:
            gene_range = genepositions[gene_search]

        #if not in dict, try alias dict
        except: 
            KeyError
            try:
                gene_name = genealiases[gene_name]
                gene_search = gene_name+'.positions'
                gene_range = genepositions[gene_search]
            except:
                KeyError
                return "We do not have that alias or gene symbol stored"

        
        startgen = gene_range[0]
        print(startgen)
        endgen = gene_range[1]
        return redirect(url_for('gene', gene_name=gene_name, gene_search=gene_search, gene_range=gene_range, startgen=startgen, endgen=endgen))
        #R code goes here
    return render_template('index_gene.html', form=form, gene_name=gene_name)

#making stats form for gene search
class GeneStats(FlaskForm):
    statistics = StringField("Enter Statistical Tests of choice (Pi, Haplotype Diversity, Fst, d_xy, Tajima D); Pi = Nucleotide Diversity, Fst = Genetic Differentiation, d_xy = Divergence, Tajima D = neutrality. You may enter up to any 5 of the tests, separated by a comma.", validators=[Required()])
    enterpop = StringField('Enter Populations of choice (British, Finnish, Iberian, Toscani and CEPH). You may enter 1 or more of the populations, separated by a comma. (Downloaded file will output data in the order you have specified here)', validators=[Required()])
    slidwid = StringField('Enter Sliding Window Parameters. This must be 2 numbers (window size & window jump) separated by a comma.', validators=[Required()])
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
        return "We don't have info about this these genomic coordinates or gene"
    if bisites == 1:    
        #return RSID search for that SNP.
        return "We don't have info about this these genomic coordinates or gene"
    if bisites == 0 and psites == 1:
        #return RSID search for that SNP (polyallelic)
        return "We don't have info about this these genomic coordinates or gene"

    try:
        form = GeneStats()
        statistics = 0
        enterpop = 0
        slidwid = 0
        if request.method == 'POST':
            print("hello")
            if form.validate_on_submit():
                statistics = form.statistics.data
                enterpop = form.enterpop.data
                slidwid = form.slidwid.data

                return redirect(url_for("generesult",statistics = statistics, enterpop=enterpop, slidwid=slidwid, startgen=startgen, endgen=endgen,gene_name = gene_name))
        return render_template('gene_view.html', nsites=nsites, bisites=bisites ,psites=psites,form=form, statistics = statistics, enterpop=enterpop, slidwid=slidwid,gene_name=gene_name)
    except:
        return "We do not have information about this gene"


@app.route("/generesult", methods=['GET']) #need GET to pass variable from one function to another
def generesult():
    statistics = request.args.get("statistics") # to pass variable from one function to another
    enterpop = request.args.get("enterpop")
    slidwid = request.args.get("slidwid")
    endgen = request.args.get("endgen")
    startgen = request.args.get("startgen")

    gene_name = request.args.get("gene_name")

    #convert pops to capitalise/CEPH
    populations = enterpop.split(", ")
    for i in range(len(populations)):
        populations[i] = populations[i].capitalize()
        if populations[i] == 'Ceph':
            populations[i] = populations[i].upper()

    #need to create stringvector in R to keep populations stored
    rpopulations = robjects.StrVector(populations)
    print(rpopulations)

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

            fstPlot(rpopulations)
               
        if stat == 'd_xy':
            #get result
            dxy_result = sum_dxy(rpopulations)
            dxy_result_dict = { key : dxy_result.rx2(key)[0] for key in dxy_result.names }
               
            for key,value in dxy_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)

            dxyPlot(rpopulations)  

        if stat == 'pi':
            #get result
            pi_result = sum_pi(rpopulations)
            pi_result_dict = { key : pi_result.rx2(key)[0] for key in pi_result.names }
               
            for key,value in pi_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)

            piPlot(rpopulations)

        if stat == 'haplotype diversity':
            hpw_result = sum_hpw(rpopulations)
            hpw_result_dict = { key : hpw_result.rx2(key)[0] for key in hpw_result.names }
               
            for key,value in hpw_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)
            hpwPlot(rpopulations)

        if stat == 'tajima d':
            tajd_result = sum_tajd(rpopulations)
            tajd_result_dict = { key : tajd_result.rx2(key)[0] for key in tajd_result.names }
               
            for key,value in tajd_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)

            tajdPlot(rpopulations)        


    ## GET PLOTS

    plots = {}

    for stat in statistics:
        if stat == 'fst':
            #add plot
            plots['static/fstplot.jpg'] = ['Genetic Differentiation (FST) Plot','Population genetic variation/Fixation index (FST); proportion of the total genetic variance contained in a subpopulation relative to the total genetic variance. Population selection influences FST values, such that closely related groups are indistinguishable']


        if stat == 'd_xy':
            #add plot
            plots['static/dxyplot.jpg'] = ['Divergence (d_XY) Plot','Within genetic diversity, (d_XY), is the absolute nucleotide divergence between two populations.']

        if stat == 'pi':
            #add plot
            plots['static/piplot.jpg'] = ['Nucleotide Diversity (Pi) Plot','Nucleotide diversity (Pi); average number of nucleotide differences per site between two DNA sequences in all possible pairs within the same population, (this is denoted as Pi). Highly diverse libraries have approximately equal populations of all four nucleotides. Low diversity libraries have a high proportion of certain nucleotides.']

        if stat == 'haplotype diversity':
            plots['static/hpwplot.jpg'] = ['Haplotype Diversity Plot','Haplotype diversity represents the measure of the uniqueness of a particular haplotype in a given population.']

        if stat == 'tajima d':
            plots['static/tajdplot.jpg'] = ["Tajima's D Plot", 'The Tajima’s D test is a measure of neutrality, it is the difference between two measures of genetic diversity. A negative result signifies an excess of low frequency polymorphisms relative to expectation (population size expansion). A positive result signifies low levels of both low and high frequency polymorphisms, indicating a contraction in population size and/or balancing selection.']
    
    # Generate file for download
    
    gene_rsids = get_gene_rsids(gene_name)
    gene_json_first = create_gene_json_template(gene_rsids)
    gene_json_final = add_gene_json_population_stats(gene_json_first, rpopulations)
    
    with open('results.json','w') as f:
        json.dump(gene_json_final,f)
        f.close()
    return render_template("generesult.html",items=items,plots=plots)    
    

@app.route("/get-download")
def get_image():
    uploads = app.config["results"] = "/home/bt211033/flask/SNPEuro1"
    filename='results.json'
    try:
        return send_from_directory(uploads,filename, as_attachment=True)
    except FileNotFoundError:
        abort(404)




#############################################################
##genomic coordinate

class Coord(FlaskForm):
    startgenomic = StringField('Starting Genomic coordinate:', validators=[Required()])
    endgenomic = StringField('Ending Genomic coordinate:', validators=[Required()])
    submit = SubmitField('Submit')

@app.route("/genomic", methods=['GET', 'POST'])
def index_genomic():
    form = Coord()
    startgenomic = 0
    endgenomic = 0
    slidwid = 0
    if request.method == "POST":
        if form.validate_on_submit():
            startgenomic = int(form.startgenomic.data)
            endgenomic = int(form.endgenomic.data)
            return redirect(url_for('genomic', startgenomic=startgenomic, endgenomic=endgenomic, slidwid = slidwid))
        else:
            return "We don't have info about these genomic coordinates"
            #return redirect(url_for('genomic'))
    return render_template('index_genomic.html', form=form, startgenomic=int(startgenomic), endgenomic=int(endgenomic))

class Stats(FlaskForm):
    statistics = StringField("Enter Statistical Tests of choice (Pi, Haplotype Diversity, Fst, d_xy, Tajima D); Pi = Nucleotide Diversity, Fst = Genetic Differentiation, d_xy = Divergence, Tajima D = neutrality. You may enter up to any 5 of the tests, separated by a comma.", validators=[Required()])
    enterpop = StringField('Enter Populations of choice (British, Finnish, Iberian, Toscani and CEPH). You may enter 1 or more of the populations, separated by a comma. (Downloaded file will output data in the order you have specified here)', validators=[Required()])
    slidwid = StringField('Enter Sliding Window Parameters. This must be 2 numbers (window size & window jump) separated by a comma.', validators=[Required()])
    submit = SubmitField('Submit')


@app.route("/genomic/<startgenomic>/<endgenomic>", methods=['GET', 'POST'])
def genomic(startgenomic, endgenomic):
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
                statistics = form.statistics.data
                enterpop = form.enterpop.data
                slidwid = form.slidwid.data

                return redirect(url_for("statspage",statistics = statistics, enterpop=enterpop, slidwid=slidwid, startgenomic=startgenomic, endgenomic=endgenomic))
        return render_template('genomic_view.html', nsites=nsites, bisites=bisites ,psites=psites,form=form, statistics = statistics, enterpop=enterpop, slidwid=slidwid)
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
    print(rpopulations)

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

            fstPlot(rpopulations)
               
        if stat == 'd_xy':
            #get result
            dxy_result = sum_dxy(rpopulations)
            dxy_result_dict = { key : dxy_result.rx2(key)[0] for key in dxy_result.names }
               
            for key,value in dxy_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)

            dxyPlot(rpopulations)  

        if stat == 'pi':
            #get result
            pi_result = sum_pi(rpopulations)
            pi_result_dict = { key : pi_result.rx2(key)[0] for key in pi_result.names }
               
            for key,value in pi_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)

            piPlot(rpopulations)

        if stat == 'haplotype diversity':
            hpw_result = sum_hpw(rpopulations)
            hpw_result_dict = { key : hpw_result.rx2(key)[0] for key in hpw_result.names }
               
            for key,value in hpw_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)
            hpwPlot(rpopulations)

        if stat == 'tajima d':
            tajd_result = sum_tajd(rpopulations)
            tajd_result_dict = { key : tajd_result.rx2(key)[0] for key in tajd_result.names }
               
            for key,value in tajd_result_dict.items():
                    a_pair = dict(result = key + ':' + str(value))
                    items.append(a_pair)

            tajdPlot(rpopulations)        


    ## GET PLOTS

    plots = {}

    for stat in statistics:
        if stat == 'fst':
            #add plot
            plots['static/fstplot.jpg'] = ['Genetic Differentiation (FST) Plot','Population genetic variation/Fixation index (FST); proportion of the total genetic variance contained in a subpopulation relative to the total genetic variance. Population selection influences FST values, such that closely related groups are indistinguishable']


        if stat == 'd_xy':
            #add plot
            plots['static/dxyplot.jpg'] = ['Divergence (d_XY) Plot','Within genetic diversity, (d_XY), is the absolute nucleotide divergence between two populations.']

        if stat == 'pi':
            #add plot
            plots['static/piplot.jpg'] = ['Nucleotide Diversity (Pi) Plot','Nucleotide diversity (Pi); average number of nucleotide differences per site between two DNA sequences in all possible pairs within the same population, (this is denoted as Pi). Highly diverse libraries have approximately equal populations of all four nucleotides. Low diversity libraries have a high proportion of certain nucleotides.']

        if stat == 'haplotype diversity':
            plots['static/hpwplot.jpg'] = ['Haplotype Diversity Plot','Haplotype diversity represents the measure of the uniqueness of a particular haplotype in a given population.']

        if stat == 'tajima d':
            plots['static/tajdplot.jpg'] = ["Tajima's D Plot", 'The Tajima’s D test is a measure of neutrality, it is the difference between two measures of genetic diversity. A negative result signifies an excess of low frequency polymorphisms relative to expectation (population size expansion). A positive result signifies low levels of both low and high frequency polymorphisms, indicating a contraction in population size and/or balancing selection.']
    
    #Get data for download

    coord_rsids = get_coord_rsids(int(startgenomic),int(endgenomic))
    coord_json_first = create_gene_json_template(coord_rsids)
    coord_json_final = add_gene_json_population_stats(coord_json_first, rpopulations)
    
    with open('results.json','w') as f:
        json.dump(coord_json_final,f)
        f.close()

    return render_template("result.html",items=items, plots = plots)


#about page
@app.route("/about_statistics")
def about_statistics():
    return render_template("about_statistics.html")



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

#starts web sever
#runs the application
if __name__== "__main__":
    app.run(debug=True)
