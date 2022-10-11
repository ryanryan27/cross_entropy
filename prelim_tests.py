import os
import subprocess

ce = './cross_entropy.exe'

pth = './graphs/param_tests/'

dom_types = ['d', 't', 's', 'c', '2']
r_vals = range(5,10)
p_vals = ['0.0001', '0.001', '0.01', '0.1']
a_vals = ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1']
nm_vals = [['10','5'],['10','10'],['100','5'],['100','20'],['100','50'],['100','80'],['100','100']]

iters = '10'
n_def = '10'
m_def = '5'
r_def = '5'
p_def = '0.01'
a_def = '0.5'


for dt in dom_types:
    # vary r
    file_res = "./results/ce_"+dt+"_n"+n_def+"_m"+m_def+"_p"+p_def+"_a"+a_def+"_rvar.csv"
    for r in r_vals:
        for fn in os.listdir(pth):
            params = [ce, '-f', pth+fn, '1', '-n', n_def, '-m', m_def, '-R', p_def, '-a', a_def, '-i', iters, '-o', '-1', '-r', str(r)]
            result = subprocess.run(params, stdout=subprocess.PIPE).stdout.decode('utf-8')
            file_obj = open(file_res, 'a')
            file_obj.write(result[:-2] + ', ' + str(r) + '\n')
            file_obj.close
    
    #vary p
    file_res = "./results/ce_"+dt+"_n"+n_def+"_m"+m_def+"_pvar"+"_a"+a_def+"_r"+r_def+".csv"
    for p in p_vals:
        for fn in os.listdir(pth):
            params = [ce, '-f', pth+fn, '1', '-n', n_def, '-m', m_def, '-R', p, '-a', a_def, '-i', iters, '-o', '-1', '-r', r_def]
            result = subprocess.run(params, stdout=subprocess.PIPE).stdout.decode('utf-8')
            file_obj = open(file_res, 'a')
            file_obj.write(result[:-2] + ', ' + p + '\n')
            file_obj.close

    #vary a
    file_res = "./results/ce_"+dt+"_n"+n_def+"_m"+m_def+"_p"+p_def+"_avar"+"_r"+r_def+".csv"
    for a in a_vals:
        for fn in os.listdir(pth):
            params = [ce, '-f', pth+fn, '1', '-n', n_def, '-m', m_def, '-R', p_def, '-a', a, '-i', iters, '-o', '-1', '-r', r_def]
            result = subprocess.run(params, stdout=subprocess.PIPE).stdout.decode('utf-8')
            file_obj = open(file_res, 'a')
            file_obj.write(result[:-2] + ', ' + a + '\n')
            file_obj.close

    #vary mn
    file_res = "./results/ce_"+dt+"_nvar_mvar"+"_p"+p_def+"_a"+a_def+"_r"+r_def+".csv"
    for nm in nm_vals:
        for fn in os.listdir(pth):
            params = [ce, '-f', pth+fn, '1', '-n', nm[0], '-m', nm[1], '-R', p_def, '-a', a_def, '-i', iters, '-o', '-1', '-r', r_def]
            result = subprocess.run(params, stdout=subprocess.PIPE).stdout.decode('utf-8')
            file_obj = open(file_res, 'a')
            file_obj.write(result[:-2] + ', ' + nm[0] + ', ' + nm[1] + '\n')
            file_obj.close



#result = subprocess.run(['./cross_entropy.exe', '-f', file_name, '1', '-d', 'd', '-o', '-1'], stdout=subprocess.PIPE)
#print(result.stdout.decode('utf-8'))

