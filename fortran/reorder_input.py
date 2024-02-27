import numpy as np

files = ['abundances/K2_18b_1bar_conv_abund.csv',
          'abundances/K2_18b_1bar_eq_abund.csv',
           'abundances/K2_18b_2500bar_conv_abund.csv']

order = ['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'H2', 'He', 'HCN', 'C2H6', 'N2']
mmws  = [  18,    44,   28,   16,    17,     2,    4,    27,    30,    28]
for f in files:
    with open(f) as g:
        l = g.readline().split('\n')[0].split(',')
    data = np.genfromtxt(f, skip_header=1, delimiter=',')
    exit
    dic = {}
    for i,head in enumerate(l[1:]):
        dic[head] = data[:,i+1]

    i=0
    for spec in order:
        if spec in dic:
            i+=1

    reordered = np.zeros((data.shape[0], i+1))
    reordered[:,0] = np.flip(dic['Pressure (bars)'])*1.e5

    i=1
    print(dic.keys())
    for spec in order:
        if spec in dic:
            reordered[:,i]=dic[spec]
            print(spec)
            i+=1

    newf = f"{f.split('.')[0]}_ordered.csv"
    np.savetxt(newf, reordered)

    print(reordered.shape)
    with open(newf,'r') as g:
        d = g.read()
    with open(newf,'w') as g:
        g.write(f'{reordered.shape[0]}, {reordered.shape[1]-1}\n'+d)

