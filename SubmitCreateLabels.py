import csv,os

run_list = []
with open('phase2_runset.csv','r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[3] == 'test':
            run_list.append(int(row[0]))
body = 'bsub -R rhel60 -W 72:00 -o output/%i.out python create_labels.py --RunNumber %i --FileNumber %i'
for run in sorted(run_list):
    for i in range(3):
        cmd = body%(run,run,i)
        print cmd
        os.system(cmd)