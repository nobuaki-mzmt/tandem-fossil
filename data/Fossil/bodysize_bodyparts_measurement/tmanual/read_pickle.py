with open( 'res.pickle', mode='rb') as f:
    tmanual_output = pickle.load(f)
x = []
y = []
for i in range(6):
    x.append( tmanual_output[2][0][4][i+2][0][0])
    y.append( tmanual_output[2][0][4][i+2][0][1])
