import pysam
import numpy as np
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson

cd /home/rhythm/Desktop/NGSWORK13dec2014/4xdata/1.5&3.5chipbam

binsize = 200
genomic_space = range(0,6988209,binsize)
control_filename = '3.5Enewtrimse.sorted.bam'
chip_filename = '3.5Nnewtrimse.sorted.bam'

control_score = []
chip_score = []
for i in genomic_space:
    control_score.append(0)
    chip_score.append(0)

control_sam = pysam.AlignmentFile(control_filename, "rb")
for i in control_sam:
    if i.pos > 1 and i.flag != 4:
        bin_loc = i.pos/binsize
        control_score[bin_loc]+=1
chip_sam = pysam.AlignmentFile(chip_filename,"rb")
for i in chip_sam:
    if i.pos >1 and i.flag != 4:
        bin_loc = i.pos/binsize
        chip_score[bin_loc]+=1

control_score = np.asarray(control_score)
chip_score = np.asarray(chip_score)

df = pd.DataFrame(data={'control': control_score, 'chip': chip_score}, index=genomic_space)
df.plot(kind='line')

library_factor = float(chip_score.sum())/control_score.sum()
print "Library size factor over control ", library_factor
chip_score = chip_score/library_factor

df = pd.DataFrame(data={'control': control_score, 'chip': chip_score}, index=genomic_space)
df.plot(kind='line')

# to measure unmapped reads 0 forward read 16 reverse reads 4 unmapped reads)

#test_sam = pysam.AlignmentFile(chip_filename, "rb")
#    if i.flag in flag_counts:
        #flag_counts[i.flag]+=1
    #else:#test_sam = pysam.AlignmentFile(control_filename,"rb")
#flag_counts = {}
#for i in test_sam:
    #if i.flag in flag_counts:
        #flag_counts[i.flag]+=1
    #else:
        #flag_counts[i.flag]=1
#print flag_counts
        
        #flag_counts[i.flag]=1
#print flag_counts

#test_sam = pysam.AlignmentFile(control_filename,"rb")
#flag_counts = {}
#for i in test_sam:
    #if i.flag in flag_counts:
        #flag_counts[i.flag]+=1#test_sam = pysam.AlignmentFile(control_filename,"rb")
#flag_counts = {}
#for i in test_sam:
    #if i.flag in flag_counts:
        #flag_counts[i.flag]+=1
    #else:
        #flag_counts[i.flag]=1
#print flag_counts
        
    #else:
        #flag_counts[i.flag]=1
#print flag_counts

substracted = chip_score - control_score
positive_peaks = []
negative_peaks = []
for i in substracted:
    if i>0:
        positive_peaks.append(i)
        negative_peaks.append(0)
    elif i<0:
        positive_peaks.append(0)
        negative_peaks.append(i)
    else:
        positive_peaks.append(0)
        negative_peaks.append(0)


plt.plot(genomic_space, positive_peaks, lw=5)
plt.plot(genomic_space, negative_peaks, lw=5)

random_lambdas = []
for i in range(100000):
    positive_vals_1000_random = []
    while len(positive_vals_1000_random) < 400:
        random_val = positive_peaks[random.choice(xrange(len(positive_peaks)))]
        if random_val > 0:
            positive_vals_1000_random.append(random_val)        
    random_lambdas.append(np.mean(positive_vals_1000_random))

y,binEdges=np.histogram(random_lambdas,bins=400)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.plot(bincenters,y,'-')
plt.show()

total_positive_non_zero_vals = 0
for i in positive_peaks:
    if i > 0:
        total_positive_non_zero_vals+=1
print total_positive_non_zero_vals

alpha_val = 0.05
lambda_val = int(np.mean(random_lambdas))
#lambda_val = np.max(random_lambdas)
global_adjusted_alpha_val = alpha_val / total_positive_non_zero_vals
print "global adjusted cut off significance level set at ", global_adjusted_alpha_val
bins_of_interest = []
genomic_coordinates = []
max_reads_in_bin = int(max(positive_peaks))
print "Will test upto %i reads in bins for test of significance " %(max_reads_in_bin)
for i in range(lambda_val, max_reads_in_bin):
    num_bins = len([x for x in positive_peaks if x>=i])
    local_adjusted_alpha_val = alpha_val/num_bins
    status_noadjust = "Insignificant"
    status_global = "Insignificant"
    status_local = "Insignificant"
    poisson_probability = poisson.pmf(i,lambda_val)
    if poisson_probability <= alpha_val:
        status_noadjust = "SIGNIFICANT"
    if poisson_probability <= global_adjusted_alpha_val:
        status_global = "SIGNIFICANT"
    if poisson_probability <= local_adjusted_alpha_val:
        status_local = "SIGNIFICANT"
    print i, poisson_probability, status_noadjust , status_global, status_local, "-----  %s bins qualify" %(num_bins)
    if status_local == "SIGNIFICANT" and len(bins_of_interest) == 0:
        for j,m in enumerate(positive_peaks):
            if m >= i:
                bins_of_interest.append(j)
                genomic_coordinates.append("Chromosome\t" + str(j*binsize-binsize) + "\t" + str(j*binsize))
    if status_global == "SIGNIFICANT":
        break

OUT = open('promacs_significant_3.5mihf.bed', 'w')
for x,i in enumerate(bins_of_interest):
    #print control_score[i],  chip_score[i], substracted[i], genomic_coordinates[x]
    OUT.write(genomic_coordinates[x]+"\n")
OUT.close()

roi = range(6166447-500,6168759+500)6988$
bins_of_interest = [x/binsize for x in roi]
roi_control = []
roi_chip = []
roi_positive = []
for i,x in enumerate(bins_of_interest):
    roi_control.append(control_score[x])
    roi_chip.append(chip_score[x])
    roi_positive.append(positive_peaks[x])
roi_df = pd.DataFrame(data={'control': roi_control, 'chip' : roi_chip, 'positive': roi_positive}, index=roi)

roi_df.plot(kind='bar') 

for j,m in enumerate(positive_peaks):
    if m >= 300:
        print "Chromosome\t" + str(j*binsize-binsize) + "\t" + str(j*binsize)




