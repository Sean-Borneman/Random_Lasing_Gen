import matplotlib.pyplot as plt
from pandas import read_csv
filename = "TEST-20231031T011812Z-002/TEST/Power_62[uJcm2]/flatLasngQ_0.062462[uJcm2].csv"
data = read_csv(filename,header=None)

plt.plot(data, color='magenta', marker='o',mfc='pink' ) #plot the data
plt.xticks(range(0,len(data)+1, 1)) #set the tick frequency on x-axis

plt.ylabel('data') #set the label for y axis
plt.xlabel('index') #set the label for x-axis
plt.title("Plotting a list") #set the title of the graph
plt.show() #display the graph
