import ReconEvents as recon 
import numpy as np 
import glob,pickle,csv
label_path = './label/'
#filename = glob.glob(label_path+'labels_*.pkl')
runs = []
with open('phase2_runset.csv','r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[1] == 'wTh' and row[3] == 'test':
            runs.append(row[0])
filename = [label_path+'labels_%s-%03i.pkl'%(run,fileNumber) for run in runs for fileNumber in range(3)]
label = []
for file in filename:
    try:
        with open(file,'r') as f:
            l = pickle.load(f)
        for row in l:
            label.append(row)
    except:
        print 'No file:',file
label = np.array(label)
print len(label),label[0]
#recon.plot_corr([],[],label[:,3],label[:,4])
res,_,mu,_ = recon.fit_res(label[:,4],[],[4000,8000,50,100],[6000,10500],True)
label_files = glob.glob(label_path+'labels_*.pkl')
for file in label_files:
    with open(file,'r') as f:
        l = pickle.load(f)
    save_file = file.replace('labels','labels_cali')
    for i in range(len(l)):
        np.put(l[i],4,l[i][4] / mu * 2615)
    with open(save_file,'w') as f:
        pickle.dump(l,f)
cali_res,_,cali_mu,_ = recon.fit_res(label[:,4] / mu * 2615,[],[4000,2615,50,100],[2000,3000],True)