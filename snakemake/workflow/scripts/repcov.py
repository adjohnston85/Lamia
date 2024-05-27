#!/usr/bin/python3

import ntpath
import pathlib
import os
import sys, getopt
import numpy as np
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Spacer, Paragraph
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from svglib.svglib import svg2rlg
import xlsxwriter 

def hotspot_report_generator(D_hotspots, D_hotspot_info, D_graphs, F_output):
    L_missing_hotspots = []
    D_sample_hotspots = {}
    L_output = []
    L_depths = []

    print(D_hotspots.keys())
    for chromo in D_hotspots:
            for exon in D_graphs[chromo]:
                                
                if D_graphs[chromo][exon][0] != chromo:
                    continue
                
                ex_start = D_graphs[chromo][exon][1]
                ex_stop = D_graphs[chromo][exon][2]
                L_ex_pos = D_graphs[chromo][exon][3]
                L_ex_depths = D_graphs[chromo][exon][4]
                
                for hotspot in D_hotspots[chromo]:
                    if ex_start <= hotspot <= ex_stop:
                        index = L_ex_pos.index(hotspot)
                        depth = L_ex_depths[index]
                        L_depths.append(depth)
                        
                        if chromo not in D_sample_hotspots:
                            D_sample_hotspots[chromo] = [(hotspot,depth,exon)]
                        else:
                            D_sample_hotspots[chromo].append((hotspot,depth,exon))  

    L = np.array([x for x in L_depths if x > 0])
    bins = int(np.max(L)/20)
    
    y, x, _ = plt.hist(L, density=True, bins=bins)
    plt.ylabel('Relative frequency')
    plt.xlabel('Read depth');
    plt.title("Targeted CpGs")

    mean_depth = 'mean depth = ' + str(int(np.mean(L)))
    max_depth = 'max depth = ' + str(int(np.max(L)))
    total_hotspots = 'CpGs targeted in panel = ' + str(len(L_depths))
    depth_0 = 'targeted CpGs with 0 reads = ' + str(L_depths.count(0))
    
    plt.annotate(mean_depth + "\n" + max_depth + "\n\n" 
                  + total_hotspots + "\n" + depth_0 + "\n\n" + 
                  "*CpGs with 0 read depth\nexcluded from graph",
                  xy=(np.max(L),y.max()+y.max()*0.01), color="black",
                  horizontalalignment='right', verticalalignment = "top",
                  fontsize = 10)
       
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.15, left=0.15)
    fig.set_size_inches(6.5, 4)
    fig.savefig("coverage_histogram.png", dpi=90)
    D_images["coverage_histogram.png"] =1
    
    plt.close()
    
    workbook = xlsxwriter.Workbook(F_output) 
    for chromo in D_sample_hotspots:
        for spot in D_sample_hotspots[chromo]:
            key = chromo + ":" + str(spot[0])
            info = D_hotspot_info[key]
            
            L= [chromo]
            L.extend(spot[0:2])
            L.extend(info)
            L_output.append(L)
    
    L_output = [["Chromosome","Position(hg38)",
                  "Read depth","CpG #","Amino_Acid_Position","Ref_Amino_Acid"]] + L_output
    
    worksheet = workbook.add_worksheet() 
    worksheet.set_column(0, 5, 12)
      
    row = 0
      
    for item in L_output:
        column = 0
        for x in item:
            worksheet.write(row, column, x) 
            column += 1

        row += 1
    
    worksheet.insert_image('I3', "coverage_histogram.png")
    workbook.close()
    return D_sample_hotspots


def coverage_histogram(L_depths, image_name):  
    x_axis_max = np.percentile(L_depths, 99.8)
    L_plt = [x for x in L_depths if x < x_axis_max]
    bins = int(np.max(L_plt))
    y, x, _ = plt.hist(L_plt, density=True, bins=bins)  # `density=False` would make counts

    plt.xlim(0, x_axis_max);
    plt.ylabel('Relative frequency')
    plt.xlabel('Read depth');
    plt.title("Histogram of all single base read depths")
    
    mean_depth = 'mean depth = ' + str(int(np.mean(L_depths)))
    max_depth = 'max depth = ' + str(int(np.max(L_depths)))
    total_bases = 'Number of bases targeted = ' + str(len(L_depths))
    depth_0 = 'Number of targed bases with 0 reads = ' + str(L_depths.count(0))
    
    plt.annotate(mean_depth + "\n" + max_depth + "\n\n" + total_bases + "\n" +
                 depth_0,
                 xy=(x_axis_max - x_axis_max*0.02,y.max()-y.max()*0.05),
                 color="black",
                 horizontalalignment='right', verticalalignment = "top",
                 fontsize = 8)

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.15)
    fig.set_size_inches(6.5, 4)   
    fig.savefig(image_name, dpi=90)
    
    image = svg2rlg(image_name)
    
    plt.close()
    
    return image


def gene_graph_generator(D_hotspots, min_threshold, mean_threshold,
                         F_bedgraph, F_outputPDF,D_graphs):  
              
    doc = SimpleDocTemplate(F_outputPDF,pagesize=letter,
                            rightMargin=15,leftMargin=15,
                            topMargin=30,bottomMargin=15)
    
    L_depths = []
    for gene in D_graphs:
        for exon in D_graphs[gene]:
            for depth in D_graphs[gene][exon][4]:
                L_depths.append(depth)
    print(L_depths)
    
    Story=[]
    
    image_name = F_outputPDF + '.svg'
    image = coverage_histogram(L_depths, image_name)

    D_images[image_name] = 1
        
    Story.append(image)
        
    count = 0
    for gene in sorted(D_graphs.keys()):
        
        L_y = []
        L_x = []
        for exon in D_graphs[gene]:
            chromo = D_graphs[gene][exon][0]
            L_y.extend(D_graphs[gene][exon][4])
            L_x.extend(D_graphs[gene][exon][3])
               
        y = np.array(L_y)
        x = np.array(L_x)
           
        mean_depth = int(np.mean(y))
        min_depth = np.min(y)
        
        count += 1
           
        plt.locator_params(axis='x', nbins=10)
        plt.yscale("log")
        
        plt.ylim(1, 1000);
        
        plt.scatter(x, y, s=2, color='blue')
        
        plt.xlabel("Nucleotide position on " + chromo + " (hg38)")
        plt.ylabel("Depth of coverage (log scale)")
        plt.title(gene)

        plt.axhline(mean_threshold, color='red', linestyle='dashed')
        plt.axhline(min_threshold, color='red', linestyle='dashed')
                
        spread = np.max(x) - np.min(x)
        offset_500 = spread * 0.11   
        offset_50 = spread * 0.09
        offset_stats = spread * 0.03

        plt.annotate(str(mean_threshold) + ' reads',
                     xy=(np.max(x)-offset_500,mean_threshold+mean_threshold*0.2),
                     color="red")
        plt.annotate(str(min_threshold) + ' reads',
                     xy=(np.max(x)-offset_50,min_threshold+min_threshold*0.2),
                     color="red")
        plt.annotate('mean = ' + str(mean_depth),
                     xy=(np.min(x)-offset_stats,1.2),
                     color="black")
        plt.annotate('min = ' + str(min_depth),
                     xy=(np.min(x)-offset_stats,1.9),
                     color="black")
        
        chromo = D_graphs[gene][exon][0]
        
        for exon in D_graphs[gene]:
            label_check = False
            hs_count = 0  
            for pos,cov in zip(D_graphs[gene][exon][3],D_graphs[gene][exon][4]):
        
                if pos in D_hotspots[chromo]:
                    hs_count+=1
                    label_check = True  
            
            L_exon_nos = []
            for ex in D_graphs[gene]:

                exon_no = float(exon.split("_exon_")[1].replace('a',''))
                L_exon_nos.append(float(ex.split("_exon_")[1].replace('a','')))

            if exon_no == min(L_exon_nos) or exon_no == max(L_exon_nos):
                
                plt.annotate("exon " + str(int(exon_no)) + "\n|", # this is the text
                          (np.mean(D_graphs[gene][exon][3]),4000), # this is the point to label
                          textcoords="offset points", # how to position the text
                          xytext=(-1,5), # distance from text to points (x,y)
                          ha='center', fontsize=7) # horizontal alignment can be left, right or center  
                    
            if label_check == True:
                
                label = "x" + str(hs_count)
                
                plt.annotate(label, # this is the text
                              (np.mean(D_graphs[gene][exon][3]),
                               np.max(D_graphs[gene][exon][4])), # this is the point to label
                              textcoords="offset points", # how to position the text
                              xytext=(-1,5), # distance from text to points (x,y)
                              ha='center', fontsize=7) # horizontal alignment can be left, right or center  

        plt.annotate('no. of mutation hotspots above each exon if >0 ',
                     xy=(np.max(x)+offset_stats,1.2),
                     color="black",
                     horizontalalignment='right')
        
        depth_plot = "gene_coverage_" + str(count) + ".svg"

        fig = plt.gcf()
        fig.subplots_adjust(bottom=0.15)
        fig.set_size_inches(6.5, 4)
        fig.savefig(depth_plot, dpi=90)

        image = svg2rlg("gene_coverage_" + str(count) + ".svg")
        D_images["gene_coverage_" + str(count) + ".svg"] = 1

        plt.close()
        Story.append(image)
        
    doc.build(Story)


def exon_graph_generator(D_hotspots, min_threshold, mean_threshold,
                         F_output,D_graphs):
   
    doc = SimpleDocTemplate(F_output,pagesize=letter,
                            rightMargin=15,leftMargin=15,
                            topMargin=30,bottomMargin=15)
    
    Story=[]
   
    count = 0
    for gene in D_graphs:
        for exon in D_graphs[gene]:
               
            y = np.array(D_graphs[gene][exon][4])
            x = np.array(D_graphs[gene][exon][3])
               
            mean_depth = round(np.mean(y),1)
            min_depth = np.min(y)
            


            if np.min(y) >= min_threshold and np.mean(y) >= mean_threshold:
                continue
               
            plt.locator_params(axis='x', nbins=10)
            plt.yscale("log")
            
            plt.ylim(1, 1000);
            
            plt.scatter(x, y, s=2, color='blue')
            
            plt.xlabel("Nucleotide position on " +
                       D_graphs[gene][exon][0] +
                       " (hg38)")
            plt.ylabel("Depth of coverage (log scale)")
            plt.title(exon)
      
            plt.axhline(mean_threshold, color='red', linestyle='dashed')
            plt.axhline(min_threshold, color='red', linestyle='dashed')
                    
            spread = np.max(x) - np.min(x)
            offset_500 = spread * 0.04
            offset_50 = spread * 0.04
            offset_stats = spread * 0.03
        
            plt.annotate(str(mean_threshold) + ' reads',
                         xy=(np.max(x)+offset_500,mean_threshold+mean_threshold*0.2),
                         color="red",
                         horizontalalignment='right')
            plt.annotate(str(min_threshold) + ' reads',
                         xy=(np.max(x)+offset_50,min_threshold+min_threshold*0.2),
                         color="red",
                         horizontalalignment='right')
            plt.annotate('mean = ' + str(mean_depth),
                         xy=(np.min(x)-offset_stats,1.2),
                         color="black")
            plt.annotate('min = ' + str(min_depth),
                         xy=(np.min(x)-offset_stats,1.9),
                         color="black")
            
            chromo = D_graphs[gene][exon][0]
            label_check = False
            for pos,cov in zip(D_graphs[gene][exon][3],D_graphs[gene][exon][4]):
         
                if pos in D_hotspots[chromo]:
                    
                    if label_check == False:
                        label_check = True
                    
                    label = "CpG"
                    
                    plt.annotate(label, # this is the text
                                 (pos,cov), # this is the point to label
                                 textcoords="offset points", # how to position the text
                                 xytext=(0,5), # distance from text to points (x,y)
                                 fontsize=2,
                                 ha='center') # horizontal alignment can be left, right or center       
      
            depth_plot = "exon_coverage_" + str(count) + ".svg"
            
            fig = plt.gcf()
            fig.subplots_adjust(bottom=0.15)
            fig.set_size_inches(6.5, 4)
            fig.savefig(depth_plot, dpi=90)
            
            image = svg2rlg("exon_coverage_" + str(count) + ".svg")
            D_images["exon_coverage_" + str(count) + ".svg"] = 1
            
            plt.close()
            Story.append(image)
                     
    doc.build(Story)


def main(argv):     

    global D_images
    
    D_bedgraph = {}
    D_hotspot_info = {}
    L_depths = []
    ref_directory = '/datasets/work/hb-meth-atlas/work/pipeline_data/Genomes/hg38/'
    D_images = {}
    D_sample_exons = {}
    
    # default full resolution per base coverage file name
    input_name = "XMP004020xrun54x13Ax102VGxT.qc-coverage-region-1_full_res"
    
    F_input = input_name + ".bed" # input is a BED/BEDgraph file format
    F_output = input_name # Output files will use the input file name as a prefix
    
    input_check = False # Check variable for whether or not input file name is specified by user
    output_check = False # Check variable for whether or not out file prefix is specified by user
    
    help_message =  " ".join([' RepCov v' + str(Version) + '\n',
                    'Executed using: ' + sys.executable + '\n\n',
                    'Usage:\tRepCov.py',
                    '-i <input file>\n\n',
                    'options:\n',
                    '  -i, --input_file\t\tinput file is a full resolution coverage BED file (i.e. begraph)\n',
                    '  -o, --output_prefix\t\toutput prefix. May include file path'
                    ', e.g. /home/coverage_reports/<output prefix>. Derived from input file if not specified\n',
                    '  -d, --depth_threshold\tminimum read depth threshold (50 by default)\n',
                    '  -m, --mean_thresold\t\tmean read depth threshold (500 by default)\n',
                    '  -s, --hotspots\t\thotspots bed file',
                    ])

    # The following takes input from the user when running script
    try:
       opts, args = getopt.getopt(argv,"hi:o:d:m:s:",
                                  ["input_file=", "output_prefix=",
                                   "depth_threshold=",
                                   "mean_threshold=",
                                   "hotspots="])
    except getopt.GetoptError:
       print (help_message)
       sys.exit(2)
    for opt, arg in opts:
        if opt == '-h': # -h = help, 
            print (help_message)
            sys.exit()
        elif opt in ("-i", "--input_file"): # input file name
            F_input = str(arg)
            input_check = True
        elif opt in ("-o", "--output_prefix"): # output file prefix
            F_output = str(arg)
            output_check = True
        elif opt in ("-d", "--depth_threshold"): # minimum coverage threshold
            min_depth_threshold = int(arg)
        elif opt in ("-m", "--mean_threshold"): # mean coverage threshold
            mean_depth_threshold = int(arg)
        elif opt in ("-s", "--hotspots"): # 
            hotspots = str(arg)

    # generate output prefix from input file name if no output specified
    prefix = ntpath.basename(F_input)
    if input_check == True:
        if output_check == False:
            directory = os.path.dirname(os.path.abspath(F_input))
            if "\\" in directory:
                F_output = directory + "\\" + prefix.split('.')[0]
            else:
                F_output = directory + "/" + prefix.split('.')[0]

    else:
        print('\nERROR: No input file provided\n')
        print (help_message)
        sys.exit()

    min_depth_threshold = 5 # default minimum Tumour coverage threshold
    mean_depth_threshold = 20 # default mean Tumour coverage threshold
        
    # print input parameters
    print ('Input per base coverage file (bedgraph):', F_input)
    print ('Output prefix:', F_output)
    print ('Minimum depth threshold: ', min_depth_threshold)
    print ('Mean depth threshold: ', mean_depth_threshold)

    # Create Dictinary of genes including chromosome position and per base coverage, split by exon
    L_input_bedgraph = []
    with open(F_input,"r") as F_input_bedgraph:
        for line in F_input_bedgraph:
            chunk = line.strip().split("\t")

            if int(chunk[2]) - int(chunk[1]) > 1:
                for x in range(int(chunk[1]),int(chunk[2])):
                    L_input_bedgraph.append(chunk[0] + "\t" + str(x) + "\t" +
                                            str(x+1) + "\t" + chunk[3] + "\n")
            else:
                L_input_bedgraph.append(line)
                
    for line in L_input_bedgraph:
        chunk = line.strip().split("\t")
        chromo = chunk[0]
        pos = int(chunk[1])
        depth = int(chunk[3])
        L_depths.append(int(chunk[3]))
        
        if chromo not in D_bedgraph:
            D_bedgraph[chromo] = [(pos,depth)]
        else:
            D_bedgraph[chromo].append((pos,depth))

    print(len(L_depths))
    print(D_bedgraph.keys())
    def data_structure(F_targets):
        D_graphs = {}
        D_exon_info={}
        L_exons = []
        with open(F_targets,"r") as F_input_exons:

            for line in F_input_exons:
                chunk = line.strip().split("\t")
                exon = chunk[0] + ':' + chunk[1] + '-' + chunk[2]
                info = (chunk[0],int(chunk[1]), int(chunk[2]))
                D_exon_info[exon] = info
                
                if exon not in L_exons:
                    L_exons.append(exon)        
        
        index = 0
        depthCount = 0
        lenBedgraph = 0
        
        for chromo in D_bedgraph:
            print(chromo)
            for pos in D_bedgraph[chromo]:
                if index >= len(L_exons):
                    break
                lenBedgraph += 1
                exon = L_exons[index]
                exon_info = D_exon_info[exon]
               
                if exon_info[0] != chromo:
                    if D_exon_info[L_exons[index-1]][0] != exon_info[0]:
                        break
                    else:
                        index += 1

                if exon_info[1] <= pos[0] < exon_info[2]:
                    depthCount += 1
                    if chromo not in D_graphs:            
                        D_graphs[chromo] = {exon:[chromo,exon_info[1],
                                                    exon_info[2], [pos[0]],
                                                    [pos[1]]]}
                        D_sample_exons[exon] = 1
                    elif exon not in D_graphs[chromo]:
                        D_graphs[chromo][exon] = [chromo,exon_info[1],
                                                    exon_info[2], [pos[0]],
                                                    [pos[1]]]
                        D_sample_exons[exon] = 1
                    else:
                        D_graphs[chromo][exon][3].append(pos[0])
                        D_graphs[chromo][exon][4].append(pos[1])
                        
                if pos[0] >= exon_info[2]-1:
                    index += 1

        return D_graphs, L_exons
    
    results = data_structure("/datasets/work/hb-meth-atlas/work/Users/andrew_temp/Brca_MRFF_probes_hg38.bed")
                             
    D_graphs = results[0]
    L_exons = results[1]
    
    coverage_histogram(L_depths, F_output + '._per_base_coverage_histogram.svg')
   
    with open(hotspots,"r") as F_input_hotspots:
        for line in F_input_hotspots:
            chunk = line.strip().split("\t")

            key = ":".join(chunk[0:2])

            hotspot =  [chunk[3]]

            if key not in D_hotspot_info:
               D_hotspot_info[key] = hotspot
            else:
                D_hotspot_info[key].append(chunk[3])

    D_hotspots={}
    for key in D_hotspot_info:
        new_key = key.split(":")
        if new_key[0] not in D_hotspots:
            D_hotspots[new_key[0]] = [int(new_key[1])]
        else:
            D_hotspots[new_key[0]].append(int(new_key[1]))


    #Hotspot Report Section
    print("hotspot_report_generator")
    hotspot_report_generator(D_hotspots, D_hotspot_info, D_graphs, F_output + ".CANCER_HOTSPOTS_REPORT.xlsx")

    #Exon Graph Section
    exon_graph_generator(D_hotspots, min_depth_threshold,
                          mean_depth_threshold, F_output +
                          ".EXON_COVERAGE_GRAPHS_FAILED_QC_72_gene_pancancer_panel.pdf",
                          D_graphs)

    #Hotspot Report Section 
    #print("hotspot_report_generator")
    #hotspot_report_generator(D_hotspots, D_hotspot_info, D_graphs, F_output + ".CANCER_HOTSPOTS_REPORT.xlsx")
    

    for x in D_images:
        os.remove(x)

if __name__ == "__main__":
   main(sys.argv[1:])


