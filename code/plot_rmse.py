from os import listdir
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# data size
s1 = 3883
s2 = 6040
s3 = 10

outpath = './rmse'

rmse = []

files = listdir(outpath)
for i in range(len(files)):
    f = open(outpath + "/" + str(i) + ".txt", "r") 
    lines = f.readlines()
    for line in lines:
        rmse.append(float(line))

data = pd.DataFrame(list(zip(list(range(len(rmse))), rmse)))
data.columns = ["Iteration", "RMSE"]


ax = sns.lineplot(x="Iteration", y="RMSE", legend="full", data=data)
plt.show()


